
//#define PRINT_STUFF
#ifdef PRINT_STUFF
#define DBG(x) std::cout<<x
#else
#define DBG(x)
#endif

#define DO_TRACES 0

//-----------------------------------------------------------------------------
class ComputeAsWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn CoordsRZ,
                                ExecObject locatorRZ,
                                WholeCellSetIn<> cellSetRZ,
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,
                                WholeArrayInOut AsOut,
                                WholeArrayInOut dAsOut);

  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7);
  using InputDomain = _1;

  ~ComputeAsWorklet()
  {
  }

  ComputeAsWorklet(XGCParameters& xgcParams)
  {
    this->NumNodes = xgcParams.numNodes;
    this->NumPlanes = xgcParams.numPlanes;
  }

  template <typename CoordsType, typename LocatorType, typename CellSetType, typename AsPhiType, typename dAsPhiType, typename AsOutType, typename dAsOutType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const CoordsType& CoordRZ,
                            const LocatorType& locatorRZ,
                            const CellSetType& cellSetRZ,
                            const AsPhiType& As_phi_ff,
                            const dAsPhiType& dAs_phi_ff_RZP,
                            AsOutType& AsOut,
                            dAsOutType& dAsOut_RZP) const
  {

    vtkm::Id cellId;
    vtkm::Vec3f param;

    vtkm::ErrorCode status = locatorRZ.FindCell(CoordRZ, cellId, param);

    if (status != vtkm::ErrorCode::Success)
    {
      vtkm::Vec3f tmp(0,0,0);
      for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
      {
        vtkm::Id outIdx = n * this->Num2DPts + idx;
        AsOut.Set(outIdx, 0);
        dAsOut_RZP.Set(outIdx, tmp);
      }

      return;
    }

    auto tmp =  cellSetRZ.GetIndices(cellId);
    vtkm::Vec<vtkm::Id,3> vIds;
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];

    for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
    {
      vtkm::FloatDefault as = 0;
      vtkm::Vec3f das_rzp;

      vtkm::Id offset = n*this->NumNodes;
      vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
      vals.Append(As_phi_ff.Get(vIds[0]+offset));
      vals.Append(As_phi_ff.Get(vIds[1]+offset));
      vals.Append(As_phi_ff.Get(vIds[2]+offset));
      vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), as);

      vtkm::VecVariable<vtkm::Vec3f, 3> valv;
      valv.Append(dAs_phi_ff_RZP.Get(vIds[0]+offset));
      valv.Append(dAs_phi_ff_RZP.Get(vIds[1]+offset));
      valv.Append(dAs_phi_ff_RZP.Get(vIds[2]+offset));
      vtkm::exec::CellInterpolate(valv, param, vtkm::CellShapeTagTriangle(), das_rzp);

      vtkm::Id outIdx = n * this->Num2DPts + idx;
      AsOut.Set(outIdx, as);
      dAsOut_RZP.Set(outIdx, das_rzp);
    }
  }

  vtkm::Id Num2DPts;
  vtkm::Id NumNodes;
  vtkm::Id NumPlanes;
};
