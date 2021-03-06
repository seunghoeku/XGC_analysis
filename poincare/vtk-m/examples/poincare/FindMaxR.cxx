#include "FindMaxR.h"

#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/CellLocatorTwoLevel.h>

class FindMaxRWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn Thetas,
                                ExecObject Locator,
                                FieldOut MaxR);
  using ExecutionSignature = void(_1, _2, _3);
  using InputDomain = _1;

  FindMaxRWorklet(vtkm::FloatDefault rmin,
                  vtkm::FloatDefault rmax,
                  vtkm::FloatDefault eq_r,
                  vtkm::FloatDefault eq_z)
    : RMin(rmin)
    , RMax(rmax)
    , EqR(eq_r)
    , EqZ(eq_z)
  {
  }

  template <typename LocatorType>
  VTKM_EXEC void operator()(const vtkm::FloatDefault& theta,
                            const LocatorType& locator,
                            vtkm::FloatDefault& val) const
  {
    const vtkm::FloatDefault cost = vtkm::Cos(theta), sint = vtkm::Sin(theta);

    vtkm::FloatDefault r0 = 0, r1 = this->RMax;// - this->EqR;
    //std::cout<<"Theta= "<<theta<<std::endl;
    //std::cout<<"R0/R1= "<<r0<<" "<<r1<<std::endl;

    vtkm::Vec3f pt(0,0,0), pcoords;
    vtkm::Id cellId;
    while (vtkm::Abs(r0-r1) > 1e-8)
    {
      vtkm::FloatDefault ri = (r0+r1) / 2.0;

      pt[0] = this->EqR + ri*cost;
      pt[1] = this->EqZ + ri*sint;
      //std::cout<<"  ri= "<<ri<<" ("<<r0<<" "<<r1<<")"<<std::endl;
      //std::cout<<"   pt= "<<pt<<std::endl;
      if (locator.FindCell(pt, cellId, pcoords) == vtkm::ErrorCode::Success)
      {
        // ri is INSIDE, set range to (ri, r1)
        r0 = ri;
      }
      else
      {
        // ri is OUTSIDE, set range to (r0, ri)
        r1 = ri;
      }
    }

    //Set the mid value.
    val = (r0+r1)/2.0;
    //std::cout<<"     MaxR= "<<val<<std::endl;
  }

  vtkm::FloatDefault RMin;
  vtkm::FloatDefault RMax;
  vtkm::FloatDefault EqR;
  vtkm::FloatDefault EqZ;
};


void
FindMaxR(const vtkm::cont::DataSet& ds,
         XGCParameters& xgcParams,
         std::vector<vtkm::FloatDefault>& thetas,
         vtkm::cont::ArrayHandle<vtkm::FloatDefault> &maxR)
{
  auto thetaArr = vtkm::cont::make_ArrayHandle(thetas, vtkm::CopyFlag::On);

  vtkm::cont::Invoker invoker;
  FindMaxRWorklet worklet(xgcParams.eq_min_r, xgcParams.eq_max_r, xgcParams.eq_axis_r, xgcParams.eq_axis_z);

  vtkm::cont::CellLocatorTwoLevel locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  invoker(worklet, thetaArr, locator, maxR);
}
