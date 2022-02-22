//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <typeinfo>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/testing/Testing.h>
#include <vtkm/filter/GhostCellClassify.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/EulerIntegrator.h>
#include <vtkm/worklet/particleadvection/Field.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Particles.h>
#include <vtkm/worklet/particleadvection/RK4Integrator.h>
#include <vtkm/worklet/particleadvection/Stepper.h>
#include <vtkm/worklet/testing/GenerateTestDataSets.h>

#include <vtkm/filter/Gradient.h>

#include <vtkm/io/VTKDataSetWriter.h>

#include <fides/DataSetReader.h>

#include <random>
#include <chrono>
#include <mpi.h>
#include "io.h"


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
    std::cout<<"Call Locator...."<<std::endl;
    vtkm::ErrorCode status = locator.FindCell(point, cellId, pcoords);
    std::cout<<"Loc: "<<point<<" cellId= "<<cellId<<std::endl;

    if (status != vtkm::ErrorCode::Success)
    {
      this->RaiseError(vtkm::ErrorString(status));
    }
  }
};

void RunPoinc(const vtkm::cont::DataSet& ds, vtkm::Id numSeeds, vtkm::Id maxPunctures)
{
#if 0
  std::cout<<"Locator test"<<std::endl;
  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  std::vector<vtkm::Vec3f> pts;
  /*
  pts.push_back({2.5, -3, 0});
  pts.push_back({3, 0, 0});
  pts.push_back({2.5, vtkm::Pi(), 0});
  pts.push_back({2.5, vtkm::Pi(), .5});
  pts.push_back({2.5, vtkm::Pi(), -.5});
  */
  pts.push_back({2.5, vtkm::Pi(), 0});
//  pts.push_back({2.5, 0, .5});
//  pts.push_back({2.5, 0.01, -.5});
  pts.push_back({2.5, vtkm::TwoPi()-0.001, 0});
  pts.push_back({2.5, vtkm::TwoPi(), 0});
  pts.push_back({2.5, 0, 0});
  pts.push_back({2.5, 0.001, 0});
  auto points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);
  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;

  vtkm::worklet::DispatcherMapField<FindCellWorklet> dispatcher;
  dispatcher.Invoke(points, locator, cellIds, pcoords);
  return;
#endif

  //poinc.
  using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec<double,3>>;
  using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
  using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;


  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  FieldHandle BField;
  ds.GetField("B").GetData().AsArrayHandle(BField);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  const vtkm::FloatDefault stepSize = 0.005;

  FieldType velocities(BField);
  GridEvalType eval(ds, velocities);
  Stepper rk4(eval, stepSize);
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  vtkm::Id maxSteps = maxPunctures * 10000;

  std::vector<vtkm::Particle> seeds;

  seeds.push_back({{2.5, 6, 0}, 0});
  vtkm::worklet::Poincare p;
  auto seedsArr = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);
  vtkm::Plane<> plane({0,0,0}, {0,1,0});
  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  auto res = p.Run(rk4, seedsArr, plane, maxSteps, maxPunctures, true);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
  std::cout<<"Timer= "<<dt.count()<<std::endl;
}

void
CalcBField(vtkm::cont::DataSet& ds)
{
    std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  //Calculate B_hat
  vtkm::cont::ArrayHandle<vtkm::Vec3f> b;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> A_s;
    std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  ds.GetField("B").GetData().AsArrayHandle(b);
  ds.GetField("apars").GetData().AsArrayHandle(A_s);
    std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  auto bPortal = b.ReadPortal();
  auto aPortal = A_s.ReadPortal();

  vtkm::Id n = b.GetNumberOfValues();
  std::vector<vtkm::Vec3f> As_bHat(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f b = vtkm::Normal(bPortal.Get(i));
    vtkm::FloatDefault As = aPortal.Get(i);
    As_bHat[i] = b * As;
  }
    std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  ds.AddField(vtkm::cont::make_FieldPoint("As_bHat", vtkm::cont::make_ArrayHandle(As_bHat, vtkm::CopyFlag::On)));

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetComputeVorticity(true);
  gradient.SetActiveField("As_bHat");
  gradient.SetOutputFieldName("grad_As_bHat");
  std::cout<<"Compute Grad"<<std::endl;
  ds = gradient.Execute(ds);
  std::cout<<"Compute Grad DONE"<<std::endl;
  ds.PrintSummary(std::cout);

  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> curl;
  ds.GetField("Vorticity").GetData().AsArrayHandle(curl);
  auto cPortal = curl.ReadPortal();
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  //This is in cartesian coordintes.....
  std::vector<vtkm::Vec3f> V(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f c = cPortal.Get(i);
    V[i] = bPortal.Get(i) + c;
  }
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  ds.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(V, vtkm::CopyFlag::On)));

  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec<double,3>>;
  using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
  using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;

  if (argc < 4)
  {
    std::cerr<<"Usage: "<<argv[0]<<" dataFile numSeeds maxPunctures [options]"<<std::endl;
    std::cerr<<config.Usage<<std::endl;
    return -1;
  }

  if (1)
  {
    bool fullGrid = true;
    std::string dataDir = std::string(argv[2]);
    vtkm::Id numSeeds = std::stoi(argv[3]);
    vtkm::Id maxPunctures = std::stoi(argv[4]);
    std::map<std::string, std::string> args;
    args["--dir"] = dataDir;

    adios = new adios2::ADIOS;
    adiosStuff["mesh"] = new adiosS(adios, "/xgc.mesh.bp", "mesh", args);
    adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", args);
    adiosStuff["diag"] = new adiosS(adios, "xgc.oneddiag.bp", "oneddiag", args);
    adiosStuff["units"] = new adiosS(adios, "xgc.units.bp", "units", args);
    adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "/node_data[0]/values", args);
    bool isXYZ = false, is2D = false, isExplicit = true;
    auto ds = ReadMesh(adiosStuff, fullGrid, extendToFull, isXYZ, is2D, isExplicit);
    std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
    ds.PrintSummary(std::cout);

    ReadVar("dpot", adiosStuff["data"], ds, is2D, isXYZ);
    ReadVar("apars", adiosStuff["data"], ds, is2D, isXYZ);
    ReadVar("/node_data[0]/values", adiosStuff["bfield"], ds, is2D, isXYZ, "B");

//    CalcBField(ds);

    /*
    vtkm::filter::Gradient gradient;
    gradient.SetComputePointGradient(true);
    gradient.SetActiveField("apars");
    gradient.SetOutputFieldName("grad_apars");
    std::cout<<"Compute Grad"<<std::endl;
    ds = gradient.Execute(ds);
    std::cout<<"Compute Grad DONE"<<std::endl;
    */

    std::cout<<"IO created dataset"<<std::endl;
    ds.PrintSummary(std::cout);

    std::cout<<"**** Dump file...."<<std::endl;
    vtkm::io::VTKDataSetWriter writer("xgc0.vtk");
    writer.WriteDataSet(ds);

    //RunPoinc(ds, numSeeds, maxPunctures);

    return 0;
  }

  //fides
  std::string dataModelFile = std::string(argv[1]);
  fides::io::DataSetReader reader(dataModelFile);

  std::unordered_map<std::string, std::string> paths;
  paths["mesh"] = std::string(argv[2]) + "/";
  paths["3d"] = std::string(argv[2]) + "/";
  paths["diag"] = std::string(argv[2]) + "/";
  paths["units"] = std::string(argv[2]) + "/";
  paths["bfield"] = std::string(argv[2]) + "/";

  auto metaData = reader.ReadMetaData(paths);
  int addPlanes = 0;

  fides::metadata::MetaData selections;

  /*
  fides::metadata::Index idx(2);
  selections.Set(fides::keys::STEP_SELECTION(), idx);
  */

  if (addPlanes > 0)
  {
    fides::metadata::Size numAddPlanes(2);
    selections.Set(fides::keys::fusion::PLANE_INSERTION(), numAddPlanes);
  }

  /*
  fides::metadata::Bool addR(true), addPhi(true), addPsi(true);
  selections.Set(fides::keys::fusion::ADD_R_FIELD(), addR);
  selections.Set(fides::keys::fusion::ADD_PHI_FIELD(), addPhi);
  selections.Set(fides::keys::fusion::ADD_PSI_FIELD(), addPsi);
  */

  vtkm::cont::PartitionedDataSet output = reader.ReadDataSet(paths, selections);
  if (output.GetNumberOfPartitions() == 0)
  {
    throw std::runtime_error("No output from XGC reader");
  }
  auto ds = output.GetPartition(0);



  vtkm::Id numPlanes = 48; //4 + 4*(addPlanes-1);
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> bfield;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> B3D;
  ds.GetField("bfield").GetData().AsArrayHandle(bfield);
  vtkm::Id n = bfield.GetNumberOfValues() / 3;
  B3D.Allocate(n*numPlanes);

  auto portalB3D = B3D.WritePortal();
  auto portalB = bfield.ReadPortal();
  vtkm::Id ii = 0;
  for (vtkm::Id p = 0; p < numPlanes; p++)
  {
    for (vtkm::Id i = 0; i < n; i++)
    {
      auto r = portalB.Get(i*3 + 0);
      auto z = portalB.Get(i*3 + 1);
      auto t = portalB.Get(i*3 + 2);

      r = 0.01;
      z = 0.01;
      t = 0.5;

      auto X = r*vtkm::Cos(t);
      auto Y = r*vtkm::Sin(t);
      auto Z = z;

      vtkm::Vec3f v(X,Z,Y); //v(X,Y,Z);
      portalB3D.Set(ii++, v);
    }
  }

  ds.AddPointField("B", B3D);

  ds.PrintSummary(std::cout);



  vtkm::io::VTKDataSetWriter writer("xgc.vtk");
  writer.WriteDataSet(ds);

  std::cout<<"Locator test"<<std::endl;
  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  std::vector<vtkm::Vec3f> pts = {{1.5, 1.5, 0}};
  auto points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);
  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;

  vtkm::worklet::DispatcherMapField<FindCellWorklet> dispatcher;
  dispatcher.Invoke(points, locator, cellIds, pcoords);

  /*
  for (const auto& p : pts)
  {
    vtkm::Id cellId = -1;
    vtkm::Vec3f pcoords;
    auto status = locator.FindCell(p, cellId, pcoords);

    std::cout<<p<<":   "<<cellId<<" "<<pcoords<<std::endl;
  }
  */

  MPI_Finalize();

#if 0

  std::string fname = argv[1];
  vtkm::Id numSeeds = std::stoi(argv[3]);
  vtkm::Id maxPunctures = std::stoi(argv[3]);

  /*
  std::string fname = "/media/dpn/disk2TB/proj/vtkm/poincare/xgc_bfield.vtk";
  fname = "/media/dpn/disk2TB/proj/vtkm/poincare/vtk-m/data/data/rectilinear/fusion.vtk";
  fname = "/home/dpn/proj/vtkm/poincare/fusion.vtk";
  fname = "/media/dpn/disk2TB/proj/vtkm/poincare/xgc_300.vtk";
  //fname = "/media/dpn/disk2TB/proj/vtkm/poincare/xgc_bfield.vtk";
  */

  vtkm::io::VTKDataSetReader reader(fname);

  auto ds = reader.ReadDataSet();

  FieldHandle BField;
  ds.GetField("B").GetData().AsArrayHandle(BField);

  const vtkm::FloatDefault stepSize = 0.005;

  FieldType velocities(BField);
  GridEvalType eval(ds, velocities);
  Stepper rk4(eval, stepSize);

  vtkm::FloatDefault x0 = 1.65, x1 = 2.25;
  x0 = 1.7; x1 = 2.0;
  vtkm::Id maxSteps = maxPunctures * 10000;

  vtkm::FloatDefault dx = (x1-x0) / (float)(numSeeds-1);
  vtkm::FloatDefault x = x0;

  std::vector<vtkm::Particle> pts;
  for (vtkm::Id id = 0; id < numSeeds; id++, x+=dx)
    pts.push_back({vtkm::Particle({x, 0, 0}, id)});

  //std::vector<vtkm::Particle> pts = {vtkm::Particle({0, 1.5, 0}, 0)};
  /*
  std::vector<vtkm::Particle> pts = {vtkm::Particle({1.4, 0, 0}, 0),
                                     vtkm::Particle({1.5, 0, 0}, 1),
                                     vtkm::Particle({1.6, 0, 0}, 2),
                                     vtkm::Particle({1.7, 0, 0}, 3)};
  */
  auto seedArray = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::Off);

  vtkm::worklet::Poincare p;
  auto seeds = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);
  vtkm::Plane<> plane({0,0,0}, {0,1,0});
  auto t1 = std::chrono::high_resolution_clock::now();
  auto res = p.Run(rk4, seeds, plane, maxSteps, maxPunctures, true);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
  std::cout<<"Timer= "<<dt.count()<<std::endl;


  vtkm::cont::DataSet outDS;
  outDS.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", res.Positions));
  outDS.SetCellSet(res.PolyLines);
  vtkm::io::VTKDataSetWriter writer("poincPts.vtk");
  writer.WriteDataSet(outDS);

  outDS.PrintSummary(std::cout);

  std::ofstream outPts;
  outPts.open("points.txt");
  int nPts = res.Positions.GetNumberOfValues();
  auto portal = res.Positions.ReadPortal();
  for (int i = 0; i < nPts; i++)
  {
    auto pt = portal.Get(i);
    outPts<<pt[0]<<", "<<pt[1]<<", "<<pt[2]<<std::endl;
  }
  outPts.close();
#endif
}

#if 0
int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr
      << "./xgc <name of the json file> <path of data source folder> <optional write data>\n";
    return 0;
  }

  bool writeVTK = false;
  if (argc >= 4)
    writeVTK = std::atoi(argv[3]) == 1;

  MPI_Init(&argc, &argv);

  std::string dataModelFile = std::string(argv[1]);
  fides::io::DataSetReader reader(dataModelFile);

  std::unordered_map<std::string, std::string> paths;
  paths["mesh"] = std::string(argv[2]) + "/";
  paths["3d"] = std::string(argv[2]) + "/";
  paths["diag"] = std::string(argv[2]) + "/";
  paths["units"] = std::string(argv[2]) + "/";
  paths["bfield"] = std::string(argv[2]) + "/";

  auto metaData = reader.ReadMetaData(paths);
  int addPlanes = 0;

  fides::metadata::MetaData selections;

  fides::metadata::Index idx(2);
  selections.Set(fides::keys::STEP_SELECTION(), idx);

  if (addPlanes > 0)
  {
    fides::metadata::Size numAddPlanes(2);
    selections.Set(fides::keys::fusion::PLANE_INSERTION(), numAddPlanes);
  }

  /*
  fides::metadata::Bool addR(true), addPhi(true), addPsi(true);
  selections.Set(fides::keys::fusion::ADD_R_FIELD(), addR);
  selections.Set(fides::keys::fusion::ADD_PHI_FIELD(), addPhi);
  selections.Set(fides::keys::fusion::ADD_PSI_FIELD(), addPsi);
  */

  vtkm::cont::PartitionedDataSet output = reader.ReadDataSet(paths, selections);
  if (output.GetNumberOfPartitions() == 0)
  {
    throw std::runtime_error("No output from XGC reader");
  }
  auto ds = output.GetPartition(0);

  vtkm::cont::ArrayHandle<int> tmp;
  ds.GetField("nphi").GetData().AsArrayHandle(tmp);
  int numPlanes = tmp.ReadPortal().Get(0);
  ds.GetField("n_n").GetData().AsArrayHandle(tmp);
  int ptsPerPlane = tmp.ReadPortal().Get(0);

  vtkm::cont::ArrayHandle<double> bfield;
  ds.GetField("bfield").GetData().AsArrayHandle(bfield);
  std::cout<<"BFIELD"<<std::endl;
  std::cout<<bfield.GetNumberOfValues()<<" nplanes= "<<numPlanes<<std::endl;
  int nPts = ds.GetNumberOfPoints();
  std::vector<vtkm::Vec3f> B(nPts);

  std::size_t bidx = 0;
  auto bfieldPortal = bfield.ReadPortal();
  for (int i = 0; i < numPlanes*(addPlanes+1); i++)
    for (int j = 0; j < ptsPerPlane; j++)
    {
      double r = bfieldPortal.Get(j*3 + 0);
      double z = bfieldPortal.Get(j*3 + 1);
      double phi = bfieldPortal.Get(j*3 + 2);
      double X = r*cos(phi);
      double Y = r*sin(phi);
      double Z = z;
      B[bidx][0] = r;
      B[bidx][1] = z;
      B[bidx][2] = phi;
      bidx++;
    }

  auto BArr = vtkm::cont::make_ArrayHandle(B, vtkm::CopyFlag::On);
  //ds.AddPointField("B", B);
  //ds.AddField(vtkm::cont::make_FieldPoint("B", B));
  ds.AddField(vtkm::cont::Field("B", vtkm::cont::Field::Association::POINTS, BArr));

  ds.PrintSummary(std::cout);

  vtkm::io::VTKDataSetWriter wrt("xgc-output.vtk");
  wrt.WriteDataSet(ds);

  MPI_Finalize();

#if 0
  MPI_Init(&argc, &argv);

  if (argc < 3)
  {
    std::cerr
      << "./xgc <name of the json file> <path of data source folder> <optional write data>\n";
    return 0;
  }

  std::string dataModelFile = std::string(argv[1]);
  fides::io::DataSetReader reader(dataModelFile);

  std::unordered_map<std::string, std::string> paths;
  paths["mesh"] = std::string(argv[2]) + "/";
  paths["3d"] = std::string(argv[2]) + "/";
//  paths["diag"] = std::string(argv[2]) + "/";
//  paths["units"] = std::string(argv[2]) + "/";
  paths["bfield"] = std::string(argv[2]) + "/";

  auto metaData = reader.ReadMetaData(paths);
  fides::metadata::MetaData selections;
  fides::metadata::Index idx(2);
  selections.Set(fides::keys::STEP_SELECTION(), idx);

  vtkm::cont::PartitionedDataSet output = reader.ReadDataSet(paths, selections);

  output.PrintSummary(std::cout);


  MPI_Finalize();
  return 0;
#endif
}
#endif

#if 0

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/filter/Streamline.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

// Example computing streamlines.
// An example vector field is available in the vtk-m data directory: magField.vtk
// Example usage:
//   this will advect 200 particles 50 steps using a step size of 0.01
//
// Particle_Advection <path-to-data-dir>/magField.vtk vec 200 50 0.01 output.vtk
//

int main(int argc, char** argv)
{
  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  if (argc < 8)
  {
    std::cerr << "Usage: " << argv[0]
              << "dataFile varName numSeeds numSteps stepSize outputFile [options]" << std::endl;
    std::cerr << "where options are: " << std::endl << config.Usage << std::endl;
    return -1;
  }

  std::string dataFile = argv[1];
  std::string varName = argv[2];
  vtkm::Id numSeeds = std::stoi(argv[3]);
  vtkm::Id numSteps = std::stoi(argv[4]);
  vtkm::FloatDefault stepSize = std::stof(argv[5]);
  std::string outputFile = argv[6];

  vtkm::cont::DataSet ds;

  if (dataFile.find(".vtk") != std::string::npos)
  {
    vtkm::io::VTKDataSetReader rdr(dataFile);
    ds = rdr.ReadDataSet();
  }
  else
  {
    std::cerr << "Unsupported data file: " << dataFile << std::endl;
    return -1;
  }

  //create seeds randomly placed withing the bounding box of the data.
  vtkm::Bounds bounds = ds.GetCoordinateSystem().GetBounds();
  std::vector<vtkm::Particle> seeds;

  for (vtkm::Id i = 0; i < numSeeds; i++)
  {
    vtkm::Particle p;
    vtkm::FloatDefault rx = (vtkm::FloatDefault)rand() / (vtkm::FloatDefault)RAND_MAX;
    vtkm::FloatDefault ry = (vtkm::FloatDefault)rand() / (vtkm::FloatDefault)RAND_MAX;
    vtkm::FloatDefault rz = (vtkm::FloatDefault)rand() / (vtkm::FloatDefault)RAND_MAX;
    p.Pos[0] = static_cast<vtkm::FloatDefault>(bounds.X.Min + rx * bounds.X.Length());
    p.Pos[1] = static_cast<vtkm::FloatDefault>(bounds.Y.Min + ry * bounds.Y.Length());
    p.Pos[2] = static_cast<vtkm::FloatDefault>(bounds.Z.Min + rz * bounds.Z.Length());
    p.ID = i;
    seeds.push_back(p);
  }
  auto seedArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::Off);

  //compute streamlines
  vtkm::filter::Streamline streamline;

  streamline.SetStepSize(stepSize);
  streamline.SetNumberOfSteps(numSteps);
  streamline.SetSeeds(seedArray);

  streamline.SetActiveField(varName);
  auto output = streamline.Execute(ds);

  vtkm::io::VTKDataSetWriter wrt(outputFile);
  wrt.WriteDataSet(output);

  return 0;
}

#endif
