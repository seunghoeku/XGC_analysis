#include <vtkm/Particle.h>

#include "XGCParameters.h"
#include "RunPoincare2.h"
#include "Poincare2.h"


void
RunPoincare2(const vtkm::cont::DataSet& ds,
             vtkm::cont::ArrayHandle<vtkm::Particle>& seeds,
             XGCParameters& xgcParams,
             vtkm::FloatDefault& stepSize,
             vtkm::Id& numPuncs,
             bool useBOnly,
             bool useTraces,
             bool useLinearB,
             const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& As_ff,
             const vtkm::cont::ArrayHandle<vtkm::Vec3f>& dAs_ff_rzp,
             const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& coeff_1D,
             const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& coeff_2D,
             const vtkm::cont::ArrayHandle<vtkm::Vec3f>& B_RZP,
             const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& psi,
             vtkm::cont::ArrayHandle<vtkm::Vec3f>& tracesArr,
             vtkm::cont::ArrayHandle<vtkm::Vec2f>& outRZ,
             vtkm::cont::ArrayHandle<vtkm::Vec2f>& outTP,
             vtkm::cont::ArrayHandle<vtkm::Id>& outID)
{
  vtkm::cont::CellLocatorTwoLevel locator2L;
  locator2L.SetCellSet(ds.GetCellSet());
  locator2L.SetCoordinates(ds.GetCoordinateSystem());
  auto startL = std::chrono::steady_clock::now();
  locator2L.Update();
  std::chrono::duration<double> dt = std::chrono::steady_clock::now()-startL;
  std::cout<<"2L build= "<<dt.count()<<std::endl;

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> eq_psi_gridArr;
  ds.GetField("eq_psi_grid").GetData().AsArrayHandle(eq_psi_gridArr);
  auto dPsi = eq_psi_gridArr.ReadPortal().Get(1)-eq_psi_gridArr.ReadPortal().Get(0);

  vtkm::cont::Invoker invoker;

  PoincareWorklet2 worklet(numPuncs, 0.0f, stepSize, useTraces, xgcParams);
  worklet.UseBOnly = useBOnly;
  worklet.one_d_cub_dpsi_inv = 1.0/dPsi;
  worklet.UseLinearB = useLinearB;

  if (useTraces)
  {
    std::vector<vtkm::Vec3f> t;
    t.resize(seeds.GetNumberOfValues()*worklet.MaxIter, {-100, -100, -100});
    std::cout<<"Allocate TRACES: "<<t.size()<<std::endl;
    tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
  }

  invoker(worklet, seeds,
          locator2L,
          ds.GetCellSet(),
          ds.GetCoordinateSystem(),
          As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
          B_RZP, psi,
          tracesArr, outRZ, outTP, outID);

  outID.SyncControlArray();
}
