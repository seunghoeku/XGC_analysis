#ifndef FindMaxR_h
#define FindMaxR_h

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ArrayHandle.h>
#include "XGCParameters.h"

void
FindMaxR(const vtkm::cont::DataSet& ds,
         XGCParameters& xgcParams,
         const vtkm::Id& numThetas,
         vtkm::cont::ArrayHandle<vtkm::FloatDefault> &thetas,
         vtkm::cont::ArrayHandle<vtkm::FloatDefault> &maxR);



#endif
