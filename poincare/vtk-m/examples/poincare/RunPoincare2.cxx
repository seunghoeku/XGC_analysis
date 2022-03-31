#include <vtkm/Particle.h>
#include <vtkm/cont/ArrayCopy.h>

#include "XGCParameters.h"
#include "RunPoincare2.h"
#include "Poincare2.h"

#include "perfstubs_api/timer.h"

void
RunPoincare2(const vtkm::cont::DataSet& ds,
             vtkm::cont::ArrayHandle<vtkm::Particle>& seeds,
             XGCParameters& xgcParams,
             std::map<std::string, std::vector<std::string>>& args,
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
  PERFSTUBS_SCOPED_TIMER_FUNC();
  //Get all the arguments...
  vtkm::FloatDefault stepSize = std::atof(args["--stepSize"][0].c_str());
  vtkm::Id numPunc = std::atoi(args["--numPunc"][0].c_str());

  bool useTraces = false;
  if (args.find("--traces") != args.end()) useTraces = std::atoi(args["--traces"][0].c_str());
  bool useBOnly = false;
  if (args.find("--useBOnly") != args.end()) useBOnly = true;
  bool useLinearB = false;
  if (args.find("--useLinearB") != args.end()) useLinearB = true;
  if (useLinearB)
  {
    useBOnly = true;
    std::cout<<"Warning: Using linear B, forcing UseBOnly = true."<<std::endl;
  }
  bool validateInterp = false;
  vtkm::Id validateInterpSkip = 1;
  if (args.find("--validateInterpolation") != args.end())
  {
    validateInterp = true;
    validateInterpSkip = static_cast<vtkm::Id>(std::stoi(args["--validateInterpolation"][0].c_str()));
  }

  auto cellSet = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  vtkm::cont::CellLocatorTwoLevel locator2L;
  locator2L.SetCellSet(cellSet);
  locator2L.SetCoordinates(ds.GetCoordinateSystem());
  auto startL = std::chrono::steady_clock::now();
  locator2L.Update();
  std::chrono::duration<double> dt = std::chrono::steady_clock::now()-startL;
  std::cout<<"2L build= "<<dt.count()<<std::endl;

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> eq_psi_gridArr;
  ds.GetField("eq_psi_grid").GetData().AsArrayHandle(eq_psi_gridArr);
  auto dPsi = eq_psi_gridArr.ReadPortal().Get(1)-eq_psi_gridArr.ReadPortal().Get(0);

  vtkm::cont::Invoker invoker;

  vtkm::Id maxItersPerPunc = 2500;
  if (args.find("--MaxItersPerPunc") != args.end())
    maxItersPerPunc = std::stoi(args["--MaxItersPerPunc"][0].c_str());

  PoincareWorklet2 worklet(numPunc, 0.0f, stepSize, maxItersPerPunc, useTraces, xgcParams);
  worklet.UseBOnly = useBOnly;

  worklet.one_d_cub_dpsi_inv = 1.0/dPsi;
  worklet.UseLinearB = useLinearB;
  worklet.ValidateInterpolation = validateInterp;
  worklet.ValidateInterpolationSkip = validateInterpSkip;

  if (useTraces)
  {
    std::vector<vtkm::Vec3f> t;
    t.resize(seeds.GetNumberOfValues()*worklet.MaxIter, {-100, -100, -100});
    std::cout<<"Allocate TRACES: "<<t.size()<<std::endl;
    tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
  }

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords;
  vtkm::cont::ArrayCopy(ds.GetCoordinateSystem().GetData(), coords);

  invoker(worklet, seeds,
          locator2L,
          cellSet, coords,
          As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
          B_RZP, psi,
          tracesArr, outRZ, outTP, outID);

  outID.SyncControlArray();
}
