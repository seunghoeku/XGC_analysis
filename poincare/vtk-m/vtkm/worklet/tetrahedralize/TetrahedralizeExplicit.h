//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_worklet_TetrahedralizeExplicit_h
#define vtk_m_worklet_TetrahedralizeExplicit_h

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/cont/CellSetExplicit.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Field.h>

#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/ScatterCounting.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>

#include <vtkm/worklet/internal/TriangulateTables.h>

namespace vtkm
{
namespace worklet
{

/// \brief Compute the tetrahedralize cells for an explicit grid data set
class TetrahedralizeExplicit
{
public:
  TetrahedralizeExplicit() {}

  //
  // Worklet to count the number of tetrahedra generated per cell
  //
  class TetrahedraPerCell : public vtkm::worklet::WorkletVisitCellsWithPoints
  {
  public:
    using ControlSignature = void(CellSetIn cells, ExecObject tables, FieldOut tetrahedronCount);
    using ExecutionSignature = _3(CellShape, _2);
    using InputDomain = _1;

    VTKM_CONT
    TetrahedraPerCell() {}

    template <typename CellShapeTag>
    VTKM_EXEC vtkm::IdComponent operator()(
      CellShapeTag shape,
      const vtkm::worklet::internal::TetrahedralizeTablesExecutionObject& tables) const
    {
      return tables.GetCount(shape);
    }
  };

  //
  // Worklet to turn cells into tetrahedra
  // Vertices remain the same and each cell is processed with needing topology
  //
  class TetrahedralizeCell : public vtkm::worklet::WorkletVisitCellsWithPoints
  {
  public:
    using ControlSignature = void(CellSetIn cellset,
                                  ExecObject tables,
                                  FieldOutCell connectivityOut);
    using ExecutionSignature = void(CellShape, PointIndices, _2, _3, VisitIndex);
    using InputDomain = _1;

    using ScatterType = vtkm::worklet::ScatterCounting;

    template <typename CellArrayType>
    VTKM_CONT static ScatterType MakeScatter(const CellArrayType& cellArray)
    {
      return ScatterType(cellArray);
    }

    // Each cell produces tetrahedra and write result at the offset
    template <typename CellShapeTag, typename ConnectivityInVec, typename ConnectivityOutVec>
    VTKM_EXEC void operator()(
      CellShapeTag shape,
      const ConnectivityInVec& connectivityIn,
      const vtkm::worklet::internal::TetrahedralizeTablesExecutionObject& tables,
      ConnectivityOutVec& connectivityOut,
      vtkm::IdComponent visitIndex) const
    {
      vtkm::IdComponent4 tetIndices = tables.GetIndices(shape, visitIndex);
      connectivityOut[0] = connectivityIn[tetIndices[0]];
      connectivityOut[1] = connectivityIn[tetIndices[1]];
      connectivityOut[2] = connectivityIn[tetIndices[2]];
      connectivityOut[3] = connectivityIn[tetIndices[3]];
    }
  };

  template <typename CellSetType>
  vtkm::cont::CellSetSingleType<> Run(const CellSetType& cellSet,
                                      vtkm::cont::ArrayHandle<vtkm::IdComponent>& outCellsPerCell)
  {
    vtkm::cont::CellSetSingleType<> outCellSet;

    vtkm::cont::Invoker invoke;

    // Output topology
    vtkm::cont::ArrayHandle<vtkm::Id> outConnectivity;

    vtkm::worklet::internal::TetrahedralizeTables tables;

    // Determine the number of output cells each input cell will generate
    invoke(TetrahedraPerCell{}, cellSet, tables, outCellsPerCell);

    // Build new cells
    invoke(TetrahedralizeCell{},
           TetrahedralizeCell::MakeScatter(outCellsPerCell),
           cellSet,
           tables,
           vtkm::cont::make_ArrayHandleGroupVec<4>(outConnectivity));

    // Add cells to output cellset
    outCellSet.Fill(cellSet.GetNumberOfPoints(), vtkm::CellShapeTagTetra::Id, 4, outConnectivity);
    return outCellSet;
  }
};
}
} // namespace vtkm::worklet

#endif // vtk_m_worklet_TetrahedralizeExplicit_h
