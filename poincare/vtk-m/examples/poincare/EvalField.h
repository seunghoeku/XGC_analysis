#ifndef EvalField_h
#define EvalField_h

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CellLocatorTwoLevel.h>
#include "XGCParameters.h"

vtkm::FloatDefault
EvalScalar(const vtkm::Vec3f& ptRZ,
           const vtkm::cont::CellLocatorTwoLevel& locator,
           const vtkm::cont::DataSet& ds,
           const std::string& fieldName);

vtkm::Vec3f
EvalVec(const vtkm::Vec3f& ptRZ,
        const vtkm::cont::CellLocatorTwoLevel& locator,
        const vtkm::cont::DataSet& ds,
        const std::string& fieldName);


#endif
