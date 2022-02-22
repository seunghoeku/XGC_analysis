//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#include <vtkm/Math.h>
#include <vtkm/Matrix.h>
#include <vtkm/Transform3D.h>
#include <vtkm/cont/CellLocatorGeneral.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/testing/Testing.h>
#include <vtkm/exec/CellInterpolate.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/ScatterPermutation.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/cont/ArrayHandleXGCCoordinates.h>

#include <vtkm/io/VTKDataSetWriter.h>

#include <ctime>
#include <random>

namespace
{

std::default_random_engine RandomGenerator;

using PointType = vtkm::Vec3f;

//-----------------------------------------------------------------------------
vtkm::cont::DataSet MakeTestDataSetUniform()
{
  return vtkm::cont::DataSetBuilderUniform::Create(
    vtkm::Id3{ 32 }, PointType{ -32.0f }, PointType{ 1.0f / 64.0f });
}

vtkm::cont::DataSet MakeTestDataSetRectilinear()
{
  std::uniform_real_distribution<vtkm::FloatDefault> coordGen(1.0f / 128.0f, 1.0f / 32.0f);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> coords[3];
  for (int i = 0; i < 3; ++i)
  {
    coords[i].Allocate(16);
    auto portal = coords[i].WritePortal();

    vtkm::FloatDefault cur = 0.0f;
    for (vtkm::Id j = 0; j < portal.GetNumberOfValues(); ++j)
    {
      cur += coordGen(RandomGenerator);
      portal.Set(j, cur);
    }
  }

  return vtkm::cont::DataSetBuilderRectilinear::Create(coords[0], coords[1], coords[2]);
}

vtkm::cont::DataSet MakeTestDataSetCurvilinear()
{
  auto recti = MakeTestDataSetRectilinear();
  auto coords = recti.GetCoordinateSystem().GetDataAsMultiplexer();

  vtkm::cont::ArrayHandle<PointType> sheared;
  sheared.Allocate(coords.GetNumberOfValues());

  auto inPortal = coords.ReadPortal();
  auto outPortal = sheared.WritePortal();
  for (vtkm::Id i = 0; i < inPortal.GetNumberOfValues(); ++i)
  {
    auto val = inPortal.Get(i);
    outPortal.Set(i, val + vtkm::make_Vec(val[1], val[2], val[0]));
  }

  vtkm::cont::DataSet curvi;
  curvi.SetCellSet(recti.GetCellSet());
  curvi.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", sheared));

  return curvi;
}

vtkm::cont::DataSet MakeTestDataSetXGC()
{
  vtkm::Id numPlanes = 32;

  std::vector<vtkm::FloatDefault> rz = {1,0, 1,1, 2,0, 2,1};
  //std::vector<vtkm::FloatDefault> rz = {2,0, 1,1, 2,1};

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> pts;
  pts.Allocate(rz.size());
  auto portal = pts.WritePortal();
  for (vtkm::Id i = 0; i < static_cast<vtkm::Id>(rz.size()); i++)
    portal.Set(i, rz[i]);

  bool isCyl = false;
  auto coords = vtkm::cont::make_ArrayHandleXGCCoordinates(pts, numPlanes, isCyl);

  {
  std::cout<<"XGC Coords"<<std::endl;
  std::ofstream ptF; ptF.open("pts.txt");
  auto portal2 = coords.ReadPortal();
  for (int i = 0; i < coords.GetNumberOfValues(); i++)
    ptF<<portal2.Get(i)[0]<<", "<<portal2.Get(i)[1]<<", "<<portal2.Get(i)[2]<<std::endl;


  }

  std::vector<vtkm::Int32> nextNode = {0,1,2,3}; //,3,4,5,6};
  std::vector<vtkm::Int32> conn = {0,1,2,2,1,3};

  auto cellSet = vtkm::cont::make_CellSetExtrude(conn, coords, nextNode, isCyl);

  vtkm::cont::DataSet ds;
  ds.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", coords));
  ds.SetCellSet(cellSet);
  ds.PrintSummary(std::cout);


  vtkm::io::VTKDataSetWriter writer("xgc.vtk");
  writer.WriteDataSet(ds);

  return ds;
}

//-----------------------------------------------------------------------------
class ParametricToWorldCoordinates : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
  using ControlSignature = void(CellSetIn cellset,
                                FieldInPoint coords,
                                FieldInOutCell pcs,
                                FieldOutCell wcs);
  using ExecutionSignature = void(CellShape, _2, _3, _4);

  using ScatterType = vtkm::worklet::ScatterPermutation<>;

  VTKM_CONT
  static ScatterType MakeScatter(const vtkm::cont::ArrayHandle<vtkm::Id>& cellIds)
  {
    return ScatterType(cellIds);
  }

  template <typename CellShapeTagType, typename PointsVecType>
  VTKM_EXEC void operator()(CellShapeTagType cellShape,
                            PointsVecType points,
                            const PointType& pc,
                            PointType& wc) const
  {
    auto status = vtkm::exec::CellInterpolate(points, pc, cellShape, wc);
    if (status != vtkm::ErrorCode::Success)
    {
      this->RaiseError(vtkm::ErrorString(status));
    }
  }
};

void GenerateRandomInput(const vtkm::cont::DataSet& ds,
                         vtkm::Id count,
                         vtkm::cont::ArrayHandle<vtkm::Id>& cellIds,
                         vtkm::cont::ArrayHandle<PointType>& pcoords,
                         vtkm::cont::ArrayHandle<PointType>& wcoords)
{
  vtkm::Id numberOfCells = ds.GetNumberOfCells();

  std::uniform_int_distribution<vtkm::Id> cellIdGen(0, numberOfCells - 1);
  std::uniform_real_distribution<vtkm::FloatDefault> pcoordGen(0.0f, 1.0f);

  cellIds.Allocate(count);
  pcoords.Allocate(count);
  wcoords.Allocate(count);

  auto cellIdPortal = cellIds.WritePortal();
  auto pcoordsPortal = pcoords.WritePortal();
  for (vtkm::Id i = 0; i < count; ++i)
  {
    cellIdPortal.Set(i, cellIdGen(RandomGenerator));

    PointType pc{ pcoordGen(RandomGenerator),
                  pcoordGen(RandomGenerator),
                  pcoordGen(RandomGenerator) };
    pcoordsPortal.Set(i, pc);
  }

  vtkm::worklet::DispatcherMapTopology<ParametricToWorldCoordinates> dispatcher(
    ParametricToWorldCoordinates::MakeScatter(cellIds));
  dispatcher.Invoke(
    ds.GetCellSet(), ds.GetCoordinateSystem().GetDataAsMultiplexer(), pcoords, wcoords);
}

//-----------------------------------------------------------------------------
class FindCellWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn points,
                                ExecObject locator,
                                FieldOut cellIds,
                                FieldOut pcoords);
  using ExecutionSignature = void(_1, _2, _3, _4);

  template <typename LocatorType>
  VTKM_EXEC void operator()(const vtkm::Vec3f& point,
                            const LocatorType& locator,
                            vtkm::Id& cellId,
                            vtkm::Vec3f& pcoords) const
  {
    vtkm::ErrorCode status = locator.FindCell(point, cellId, pcoords);
    if (status != vtkm::ErrorCode::Success)
    {
      this->RaiseError(vtkm::ErrorString(status));
    }
  }
};

void TestWithDataSet(vtkm::cont::CellLocatorGeneral& locator, const vtkm::cont::DataSet& dataset)
{
  locator.SetCellSet(dataset.GetCellSet());
  locator.SetCoordinates(dataset.GetCoordinateSystem());
  locator.Update();

  vtkm::cont::ArrayHandle<vtkm::Id> expCellIds;
  vtkm::cont::ArrayHandle<PointType> expPCoords;
  vtkm::cont::ArrayHandle<PointType> points;
  GenerateRandomInput(dataset, 64, expCellIds, expPCoords, points);

  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<PointType> pcoords;

  vtkm::worklet::DispatcherMapField<FindCellWorklet> dispatcher;
  dispatcher.Invoke(points, locator, cellIds, pcoords);

  auto cellIdPortal = cellIds.ReadPortal();
  auto expCellIdsPortal = expCellIds.ReadPortal();
  auto pcoordsPortal = pcoords.ReadPortal();
  auto expPCoordsPortal = expPCoords.ReadPortal();
  for (vtkm::Id i = 0; i < 64; ++i)
  {
    VTKM_TEST_ASSERT(cellIdPortal.Get(i) == expCellIdsPortal.Get(i), "Incorrect cell ids");
    VTKM_TEST_ASSERT(test_equal(pcoordsPortal.Get(i), expPCoordsPortal.Get(i), 1e-3),
                     "Incorrect parameteric coordinates");
  }
}

PointType RotatePt(const PointType& pt, vtkm::FloatDefault deg)
{
  PointType pt2;
  auto mat = vtkm::Transform3DRotate(deg, {0,0,1});

  pt2 = vtkm::Transform3DVector(mat, pt);

  vtkm::FloatDefault r = vtkm::Sqrt(pt[0]*pt[0] + pt[2]*pt[2]);
//  vtkm::FloatDefault theta = vtkm::

//  pt2[0] =

  return pt2;
}

void TestCellLocatorGeneral()
{
  vtkm::cont::CellLocatorGeneral locator;

  /*
  TestWithDataSet(locator, MakeTestDataSetUniform());

  TestWithDataSet(locator, MakeTestDataSetRectilinear());

  TestWithDataSet(locator, MakeTestDataSetCurvilinear());
  */

  //TestWithDataSet(locator, MakeTestDataSetXGC());
  auto dataset = MakeTestDataSetXGC();
  locator.SetCellSet(dataset.GetCellSet());
  locator.SetCoordinates(dataset.GetCoordinateSystem());
  locator.Update();

  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<PointType> pcoords;
  vtkm::cont::ArrayHandle<PointType> points;
  std::vector<PointType> pts = {

    {1.2, 0, .1},
    {1.8, 0, .9},
    {1.5, 0, .1},
    {1.5, 0, .9},
    {1.5, 0, .3},
    {1.5, 0, .7},
    {1.5, 0, .501},

    {1.5, 0, .5},
    {1.5, 0, .49},

    {1.2, 0, .1},
    {1.75, 0, .9},
    {0,1.75, .9},
    {-1.75, 0, .9},
    {0,-1.01, .8}

  };
  int n = pts.size();

  for (int i = 0; i < n; i++)
    pts.push_back(RotatePt(pts[i], 45+1));
  for (int i = 0; i < n; i++)
    pts.push_back(RotatePt(pts[i], 90+1));
  for (int i = 0; i < n; i++)
    pts.push_back(RotatePt(pts[i], 135+1));
  for (int i = 0; i < n; i++)
    pts.push_back(RotatePt(pts[i], 180+1));
  for (int i = 0; i < n; i++)
    pts.push_back(RotatePt(pts[i], 225+1));
  for (int i = 0; i < n; i++)
    pts.push_back(RotatePt(pts[i], 270+1));

  for (int i = 0; i < n; i++)
    pts.push_back(RotatePt(pts[i], 315+1));

  std::ofstream ptF; ptF.open("PT.txt");
  for (int i = 0; i < pts.size(); i++)
    ptF<<pts[i][0]<<", "<<pts[i][1]<<", "<<pts[i][2]<<std::endl;



  points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

  vtkm::worklet::DispatcherMapField<FindCellWorklet> dispatcher;
  dispatcher.Invoke(points, locator, cellIds, pcoords);



}

} // anonymous namespace

int UnitTestCellLocatorGeneral(int argc, char* argv[])
{
  vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagSerial{});
  return vtkm::cont::testing::Testing::Run(TestCellLocatorGeneral, argc, argv);
}
