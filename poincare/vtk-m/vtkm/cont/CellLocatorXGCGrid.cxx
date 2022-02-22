//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleXGCCoordinates.h>
#include <vtkm/cont/CellLocatorXGCGrid.h>
#include <vtkm/cont/CellSetExtrude.h>
#include <vtkm/cont/CellSetSingleType.h>

#include <vtkm/cont/CellLocatorTwoLevel.h>
#include <vtkm/exec/ConnectivityExtrude.h>


namespace vtkm
{
namespace cont
{

using XGCType = vtkm::cont::ArrayHandleXGCCoordinates<vtkm::FloatDefault>;
using ExtrudedCell = vtkm::cont::CellSetExtrude;

void CellLocatorXGCGrid::Build()
{
  vtkm::cont::CoordinateSystem coords = this->GetCoordinates();
  vtkm::cont::DynamicCellSet cellSet = this->GetCellSet();


  if (!coords.GetData().IsType<XGCType>())
    throw vtkm::cont::ErrorBadType("Coordinates are not XGC type.");
  if (!cellSet.IsSameType(ExtrudedCell()))
    throw vtkm::cont::ErrorBadType("Cells are not Extruded type.");

  //Create the coords/cellset for the plane.

  XGCType xgcCoordSys;
  //auto xgcCoordSys = coords.GetData().Cast<XGCType>();
  coords.GetData().AsArrayHandle(xgcCoordSys);

  auto xgcPts = xgcCoordSys.GetArray();
  this->IsCylindrical = xgcCoordSys.GetUseCylindrical();
  /*
  if (!this->IsCylindrical)
    throw vtkm::cont::ErrorBadType("XGC Coordinates are not cylindrical.");
  */

  auto xgcCellSet = cellSet.Cast<vtkm::cont::CellSetExtrude>();

  vtkm::Id ptsPerPlane = xgcCellSet.GetNumberOfPointsPerPlane();
  this->NumPlanes = xgcCellSet.GetNumberOfPlanes();

  vtkm::cont::ArrayHandle<vtkm::Vec3f> planePts;
  planePts.Allocate(ptsPerPlane);
  auto portal = xgcPts.ReadPortal();
  auto portal3d = planePts.WritePortal();
  for (vtkm::Id i = 0; i < ptsPerPlane; i++)
  {
    vtkm::FloatDefault R = portal.Get(i*2+0);
    vtkm::FloatDefault Z = portal.Get(i*2+1);
    vtkm::Vec3f pt(R, Z, 0);
    portal3d.Set(i, pt);
  }

  vtkm::cont::ArrayHandle<vtkm::Id> conn;
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandleCast<vtkm::Id>(xgcCellSet.GetConnectivityArray()),
                        conn);
  this->PlaneCells.Fill(ptsPerPlane, vtkm::CELL_SHAPE_TRIANGLE, 3, conn);
  this->CellsPerPlane = this->PlaneCells.GetNumberOfCells();

  this->PlaneCoords = vtkm::cont::CoordinateSystem("coords", planePts);
  this->TwoLevelLocator.SetCellSet(this->PlaneCells);
  this->TwoLevelLocator.SetCoordinates(this->PlaneCoords);
}

vtkm::exec::CellLocatorXGCGrid CellLocatorXGCGrid::PrepareForExecution(
  vtkm::cont::DeviceAdapterId device,
  vtkm::cont::Token& token) const
{
//  using CoordsPortalType =
//    typename vtkm::cont::CoordinateSystem::MultiplexerArrayType::ReadPortalType;

  this->TwoLevelLocator.Update();
  this->Update();

  auto locMux = this->TwoLevelLocator.PrepareForExecution(device, token);
  using CellLocatorType = vtkm::cont::CellSetSingleType<>::ExecConnectivityType<vtkm::TopologyElementTagCell, vtkm::TopologyElementTagPoint>;

  vtkm::cont::ArrayHandleXGCCoordinates<vtkm::FloatDefault> xgcCoords;
  this->GetCoordinates().GetData().AsArrayHandle(xgcCoords);

  auto coordsExec = xgcCoords.PrepareForInput(device, token);
  //CoordsPortalType coordsExec = xgcCoords.PrepareForInput(device, token);
//  coordsExec.meow();
  vtkm::internal::ArrayPortalXGCCoordinates<vtkm::internal::ArrayPortalBasicRead<vtkm::FloatDefault>> bum;
  bum = xgcCoords.PrepareForInput(device, token);
  vtkm::cont::DynamicCellSet cellSet = this->GetCellSet();
  auto xgcCellSet = cellSet.Cast<vtkm::cont::CellSetExtrude>();
  auto cellSetExec = xgcCellSet.PrepareForInput(device, vtkm::TopologyElementTagCell{}, vtkm::TopologyElementTagPoint{}, token);
  //cellSetExec.meow();

  auto locator = locMux.Locators.Get<vtkm::exec::CellLocatorTwoLevel<CellLocatorType>>();
  return vtkm::exec::CellLocatorXGCGrid(cellSetExec, coordsExec, locator, this->NumPlanes, this->CellsPerPlane, this->IsCylindrical);

  /*
  return vtkm::exec::CellLocatorXGCGrid(//xgcCellSet,
                                        this->GetCellSet(),
                                        locator,
                                        this->IsCylindrical,
                                        this->CellsPerPlane,
                                        this->NumPlanes);
//                                        device,
//                                        token);
//                                        cellSetExec);
*/
}

} //namespace cont
} //namespace vtkm
