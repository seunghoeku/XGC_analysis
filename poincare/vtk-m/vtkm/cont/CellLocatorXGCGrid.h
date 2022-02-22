//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtkm_cont_CellLocatorXGCGrid_h
#define vtkm_cont_CellLocatorXGCGrid_h

#include <vtkm/cont/internal/CellLocatorBase.h>
#include <vtkm/cont/CellLocatorTwoLevel.h>

#include <vtkm/exec/CellLocatorXGCGrid.h>

namespace vtkm
{
namespace cont
{

class VTKM_CONT_EXPORT CellLocatorXGCGrid
  : public vtkm::cont::internal::CellLocatorBase<CellLocatorXGCGrid>
{
  using Superclass = vtkm::cont::internal::CellLocatorBase<CellLocatorXGCGrid>;

public:
  VTKM_CONT vtkm::exec::CellLocatorXGCGrid PrepareForExecution(
    vtkm::cont::DeviceAdapterId device,
    vtkm::cont::Token& token) const;

private:

  vtkm::Id NumPlanes;
  vtkm::Id CellsPerPlane;
  bool IsCylindrical;
  vtkm::cont::CellLocatorTwoLevel TwoLevelLocator;
  vtkm::cont::CoordinateSystem PlaneCoords;
  vtkm::cont::CellSetSingleType<> PlaneCells;

  vtkm::Id3 CellDims;
  vtkm::Id3 PointDims;
  vtkm::Vec3f Origin;
  vtkm::Vec3f InvSpacing;
  vtkm::Vec3f MaxPoint;
  bool Is3D = true;

  friend Superclass;
  VTKM_CONT void Build();
};
}
} // vtkm::cont

#endif //vtkm_cont_CellLocatorXGCGrid_h
