
//#define PRINT_STUFF
#ifdef PRINT_STUFF
#define DBG(x) std::cout<<x
#else
#define DBG(x)
#endif

#define DO_TRACES 0

//-----------------------------------------------------------------------------
class ComputeAsCellWorklet : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  using ControlSignature = void(CellSetIn Cells,
                                WholeArrayIn Coords,
                                ExecObject locator,
                                WholeCellSetIn<> cellSet,
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,
                                WholeArrayInOut AsOut,
                                WholeArrayInOut dAsOut);

  using ExecutionSignature = void(InputIndex, PointCount, PointIndices, _2, _3, _4, _5, _6, _7, _8);
  using InputDomain = _1;

  ComputeAsCellWorklet(XGCParameters& xgcParams)
  {
    this->NumNodes = xgcParams.numNodes;
    this->NumPlanes = xgcParams.numPlanes;
  }

  template <typename CoordsType, typename VertexIndexType, typename LocatorType, typename CellSetType, typename AsPhiType, typename dAsPhiType, typename AsOutType, typename dAsOutType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const vtkm::IdComponent& numPoints,
                            const VertexIndexType& ptIndices,
                            const CoordsType& coordsRZN,
                            const LocatorType& locator,
                            const CellSetType& cellSet,
                            const AsPhiType& As_phi_ff,
                            const dAsPhiType& dAs_phi_ff_RZP,
                            AsOutType& AsOut,
                            dAsOutType& dAsOut_RZP) const
  {
    //Get the center point.
    vtkm::FloatDefault R=0, Z=0;
    for (vtkm::IdComponent i = 0; i < numPoints; i++)
    {
      auto coordRZ = coordsRZN.Get(ptIndices[i]);
      R += coordRZ[0];
      Z += coordRZ[1];
    }
    auto div = 1/static_cast<vtkm::FloatDefault>(numPoints);
    R *= div;
    Z *= div;

    //vtkm::Id N = static_cast<vtkm::Id>(coordsRZN.Get(ptIndices[0])[2]);

    /*
    if (idx == 0)
    {
      R = 3.0;
      Z = 0.1;
      N = 12;
    }
    */

    vtkm::Id cellId;
    vtkm::Vec3f ptRZ(R,Z,0), param;
    vtkm::Vec<vtkm::Id,3> vIds;
    vtkm::ErrorCode status = locator.FindCell(ptRZ, cellId, param);

    if (status != vtkm::ErrorCode::Success)
    {
      vtkm::Vec3f tmp(0,0,0);
      for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
      {
        vtkm::Id outIdx = n * this->Num2DCells + idx;
        AsOut.Set(outIdx, 0);
        dAsOut_RZP.Set(outIdx, tmp);
      }

      return;
    }

    for (vtkm::Id n = 0; n < this->NumPlanes*2; n++)
    {
      vtkm::FloatDefault as = 0;
      vtkm::Vec3f das_rzp;

      auto tmp =  cellSet.GetIndices(cellId);
      vIds[0] = tmp[0];
      vIds[1] = tmp[1];
      vIds[2] = tmp[2];
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

      vtkm::Id outIdx = n * this->Num2DCells + idx;
      AsOut.Set(outIdx, as);
      dAsOut_RZP.Set(outIdx, das_rzp);
      //std::cout<<idx<<" :: outIdx= "<<outIdx<<" out= "<<as<<" "<<das_rzp<<std::endl;
    }


#if 0
    AsOut = 0;
    dAsOut_RZP = {0,0,0};
    if (status != vtkm::ErrorCode::Success)
      return;

    auto tmp =  cellSet.GetIndices(cellId);
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];

    vtkm::Id offset = N*this->NumNodes;

    vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
    vals.Append(As_phi_ff.Get(vIds[0]+offset));
    vals.Append(As_phi_ff.Get(vIds[1]+offset));
    vals.Append(As_phi_ff.Get(vIds[2]+offset));
    vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), AsOut);

    vtkm::VecVariable<vtkm::Vec3f, 3> valv;
    valv.Append(dAs_phi_ff_RZP.Get(vIds[0]+offset));
    valv.Append(dAs_phi_ff_RZP.Get(vIds[1]+offset));
    valv.Append(dAs_phi_ff_RZP.Get(vIds[2]+offset));
    vtkm::exec::CellInterpolate(valv, param, vtkm::CellShapeTagTriangle(), dAsOut_RZP);

    if (idx == 0)
    {
      vtkm::Vec3f ptRZN(R,Z,N);
      std::cout<<" Pt: "<<ptRZN<<" N= "<<N<<" param= "<<param<<std::endl;
      std::cout<<"    cid: "<<cellId<<" vids= "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<std::endl;
      std::cout<<"    As:  "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<" --> "<<AsOut<<std::endl;
      std::cout<<"    dAs: "<<valv[0]<<" "<<valv[1]<<" "<<valv[2]<<" --> "<<dAsOut_RZP<<std::endl;
    }

    std::cout<<"idx: "<<idx<<" "<<AsOut<<std::endl;
#endif
  }

  vtkm::Id Num2DCells;
  vtkm::Id NumNodes;
  vtkm::Id NumPlanes;
};
