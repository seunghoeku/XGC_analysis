#ifndef RunPoincare2_h
#define RunPoincare2_h

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ArrayHandle.h>

#include "XGCParameters.h"


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
             vtkm::cont::ArrayHandle<vtkm::Id>& outID);

#endif
