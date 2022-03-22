#include "EvalField.h"

#include <vtkm/worklet/WorkletMapField.h>

//-----------------------------------------------------------------------------
class EvalFieldWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn points,
                                WholeCellSetIn<> cellSet,
                                ExecObject locator,
                                WholeArrayIn field,
                                FieldOut val);
  using ExecutionSignature = void(_1, _2, _3, _4, _5);
  using InputDomain = _1;

  template <typename LocatorType, typename CellSetType, typename FieldType, typename OutType>
  VTKM_EXEC void operator()(const vtkm::Vec3f& point,
                            const CellSetType& cellSet,
                            const LocatorType& locator,
                            const FieldType& field,
                            OutType& val) const
  {
    vtkm::Id cellId;
    vtkm::Vec3f pcoords;
    vtkm::ErrorCode status = locator.FindCell(point, cellId, pcoords);
    if (status != vtkm::ErrorCode::Success)
      this->RaiseError(vtkm::ErrorString(status));

    auto ptIndices = cellSet.GetIndices(cellId);
    vtkm::VecVariable<OutType, 3> vals;
    vals.Append(field.Get(ptIndices[0]));
    vals.Append(field.Get(ptIndices[1]));
    vals.Append(field.Get(ptIndices[2]));

    vtkm::exec::CellInterpolate(vals, pcoords, vtkm::CellShapeTagTriangle(), val);
  }
};

vtkm::FloatDefault
EvalScalar(const vtkm::Vec3f& ptRZ,
           const vtkm::cont::CellLocatorTwoLevel& locator,
           const vtkm::cont::DataSet& ds,
           const std::string& fieldName)
{
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> field, out;
  ds.GetField(fieldName).GetData().AsArrayHandle(field);

  auto ptArray = vtkm::cont::make_ArrayHandle(&ptRZ, 1, vtkm::CopyFlag::On);
  vtkm::cont::Invoker invoker;
  invoker(EvalFieldWorklet{}, ptArray, ds.GetCellSet(), locator, field, out);

  return  out.ReadPortal().Get(0);
}

vtkm::Vec3f
EvalVec(const vtkm::Vec3f& ptRZ,
        const vtkm::cont::CellLocatorTwoLevel& locator,
        const vtkm::cont::DataSet& ds,
        const std::string& fieldName)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> field, out;
  ds.GetField(fieldName).GetData().AsArrayHandle(field);

  auto ptArray = vtkm::cont::make_ArrayHandle(&ptRZ, 1, vtkm::CopyFlag::On);
  vtkm::cont::Invoker invoker;
  invoker(EvalFieldWorklet{}, ptArray, ds.GetCellSet(), locator, field, out);

  return  out.ReadPortal().Get(0);
}
