#ifndef EvalField_h
#define EvalField_h

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ArrayHandle.h>
#include "XGCParameters.h"

vtkm::FloatDefault
EvalField(const vtkm::Vec3f& ptRZ,
          const vtkm::cont::DataSet& ds,
          const std::string& fieldName);

#endif
