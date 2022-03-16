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
#include <iomanip>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/Particle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/Algorithm.h>
#ifdef VTKM_CUDA
#include <vtkm/cont/cuda/internal/ScopedCudaStackSize.h>
#include <vtkm/cont/cuda/internal/DeviceAdapterAlgorithmCuda.h>
#endif

#include <adios2.h>
#include <random>
#include <chrono>
#include <variant>
#include <mpi.h>

#include "XGCParameters.h"
#include "FindMaxR.h"
#include "RunPoincare2.h"

adios2::ADIOS *adios = NULL;
class adiosS;
std::map<std::string, adiosS*> adiosStuff;

//#define VALGRIND
#ifdef VALGRIND
#include <valgrind/callgrind.h>
#endif

#define BUILD_POINC2

#ifdef BUILD_POINC1
#include "Poincare.h"
#endif
#ifdef BUILD_POINC2
#include "Poincare2.h"
#endif
#ifdef BUILD_POINC3
#include "Poincare3.h"
#include "ComputeB.h"
#include "ComputeBCell.h"
#endif

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v)
{
  out<<"[";
  for (const auto& x : v)
    out<<x<<" ";
  out<<"]";
  return out;
}

//using ArgumentType = std::variant<std::vector<int>, std::vector<float>, std::vector<std::string>>;
bool
ParseArgs(int argc, char **argv, std::map<std::string, std::vector<std::string>> &args)
{
  args.clear();
  /*
  args["--dir"] = {""};
  args["--numSeeds"] = {"10"};
  args["--numPunc"] = {"10"};
  args["--stepSize"] = {"0.01"};
  args["--varname"] = {"V"};
  args["--seed"] = {"0", "0", "0"};
  */

//  int i = 1;
  std::string a0;
  std::vector<std::string> a1;
  //while (i < argc)

  for (int i = 1; i < argc; i++)
  {
    std::string tmp(argv[i]);
    if (tmp.find("--") != std::string::npos)
    {
      if (!a0.empty())
      {
        args[a0] = a1;
        a1.clear();
      }

      a0 = tmp;
      continue;
    }
    else
      a1.push_back(tmp);
  }
  //last argument.
  if (!a0.empty())
    args[a0] = a1;

  std::cout<<"ARGS\n";
  for (const auto& it : args)
  {
    std::cout<<it.first<<" : {";
    for (const auto& jt : it.second)
      std::cout<<jt<<" ";
    std::cout<<"}\n";
  }

  return true;
}


class adiosS
{
public:
    adiosS() {}
    adiosS(adios2::ADIOS *adiosPtr,
           const std::string &fn,
           const std::string &ioNm,
           const std::map<std::string, std::string> &args) : ioName(ioNm)
    {
        std::string pathNm = ".";

        auto x = args.find("--dir")->second;
        if (x.size() > 0) pathNm = x;
        this->fileName = pathNm + "/" + fn;
        std::cout<<"Open: "<<this->fileName<<std::endl;
        this->io = adios2::IO(adiosPtr->DeclareIO(this->ioName));
        this->engine = io.Open(fileName, adios2::Mode::Read);
    }
    ~adiosS() { engine.Close(); }
    adiosS& operator=(const adiosS &a)
    {
        ioName = a.ioName;
        fileName = a.fileName;
        io = a.io;
        engine = a.engine;
        return *this;
    }

    std::string ioName, fileName;
    adios2::IO io;
    adios2::Engine engine;
};

void
ReadOther(adiosS* stuff,
          vtkm::cont::DataSet& ds,
          const std::string& vname,
          std::string fileName="")
{
  if (fileName.empty())
    fileName = vname;

  auto v = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

#ifdef VTKM_USE_DOUBLE_PRECISION
  ds.AddField(vtkm::cont::make_Field(vname,
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     tmp,
                                     vtkm::CopyFlag::On));
#else
  //convert from double to vtkm::FloatDefault
  std::vector<vtkm::FloatDefault> tmpF(tmp.size());
  for (std::size_t i = 0; i < tmp.size(); i++)
    tmpF[i] = static_cast<vtkm::FloatDefault>(tmp[i]);

  ds.AddField(vtkm::cont::make_Field(vname,
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     tmpF,
                                     vtkm::CopyFlag::On));
#endif
}

void
ReadScalar(adiosS* stuff,
           XGCParameters& xgcParams,
           vtkm::cont::DataSet& ds,
           const std::string& vname,
           std::string fileName="",
           bool add3D=false,
           bool addExtra=false)
{
  if (fileName.empty())
    fileName = vname;

  auto v = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

  vtkm::Id numPts = ds.GetNumberOfPoints();
  std::vector<double> tmpPlane(numPts);
  for (int i = 0; i < numPts; i++)
    tmpPlane[i] = tmp[i];

  //Normalize psi by eq_x_psi
  /*
  if (vname == "psi")
  {
    for (int i = 0; i < numPts; i++)
      tmpPlane[i] /= eq_x_psi;
  }
  */

  if (addExtra && add3D)
  {
    for (int i = 0; i < xgcParams.numNodes; i++)
      tmp.push_back(tmp[i]);
  }

  ds.AddField(vtkm::cont::make_FieldPoint(vname+"2D", vtkm::cont::make_ArrayHandle(tmpPlane, vtkm::CopyFlag::On)));
  if (add3D)
    ds.AddField(vtkm::cont::make_FieldPoint(vname, vtkm::cont::make_ArrayHandle(tmp, vtkm::CopyFlag::On)));
}

void
ReadB(adiosS* stuff,
      XGCParameters& xgcParams,
      vtkm::cont::DataSet& ds)
{
  std::string fileName = "/node_data[0]/values";

  auto Bvar = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(Bvar, tmp, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> b_rzp;
  for (int i = 0; i < xgcParams.numNodes; i++)
  {
    vtkm::Vec3f v(tmp[i*3+0], tmp[i*3+1], tmp[i*3+2]);
    b_rzp.push_back(v);
  }

  ds.AddField(vtkm::cont::make_FieldPoint("B_RZP",vtkm::cont::make_ArrayHandle(b_rzp, vtkm::CopyFlag::On)));
}

void
ReadPsiInterp(adiosS* eqStuff,
              adiosS* interpStuff,
              vtkm::cont::DataSet& ds,
              XGCParameters& xgcParams)
{
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mr"), &xgcParams.eq_mr, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mz"), &xgcParams.eq_mz, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_axis_r"), &xgcParams.eq_axis_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_axis_z"), &xgcParams.eq_axis_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_min_r"), &xgcParams.eq_min_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_max_r"), &xgcParams.eq_max_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_min_z"), &xgcParams.eq_min_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_max_z"), &xgcParams.eq_max_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_psi"), &xgcParams.eq_x_psi, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_r"), &xgcParams.eq_x_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_z"), &xgcParams.eq_x_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_z"), &xgcParams.eq_x_z, adios2::Mode::Sync);

  ReadOther(eqStuff, ds, "eq_I");
  ReadOther(eqStuff, ds, "eq_psi_grid");
  ReadOther(eqStuff, ds, "eq_psi_rz");
  ReadOther(eqStuff, ds, "eq_psi_grid");
  ReadOther(interpStuff, ds, "coeff_1D", "one_d_cub_psi_acoef");

  std::vector<double> tmp2D;
//  interpStuff->engine.Get(interpStuff->io.InquireVariable<double>("one_d_cub_psi_acoef"),
//                          tmp1D, adios2::Mode::Sync);
  interpStuff->engine.Get(interpStuff->io.InquireVariable<double>("psi_bicub_acoef"),
                          tmp2D, adios2::Mode::Sync);

  /*
  vtkm::Id n = ds.GetField("eq_psi_grid").GetData().GetNumberOfValues();
  std::vector<std::vector<double>> coef1D;
  coef1D.resize(n);
  for (int i = 0; i < n; i++)
    coef1D[i].resize(4);

  int idx = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < 4; j++)
    {
      coef1D[i][j] = tmp1D[idx];
      idx++;
    }
  */

  std::cout<<"****** READ: "<<tmp2D[0]<<" "<<tmp2D[16]<<" "<<tmp2D[32]<<std::endl;
  std::cout<<"          :: "<<tmp2D[1]<<" "<<tmp2D[17]<<std::endl;

  int idx = 0;
  int nr = xgcParams.eq_mr-1, nz = xgcParams.eq_mz-1;
  //int ni = 150, nj = 150;
  std::vector<std::vector<std::vector<std::vector<double>>>> coef2D;
  coef2D.resize(nz);
  for (int i = 0; i < nz; i++)
  {
    coef2D[i].resize(nr);
    for (int j = 0; j < nr; j++)
    {
      coef2D[i][j].resize(4);
      for (int k = 0; k < 4; k++)
      {
        coef2D[i][j][k].resize(4);
        for (int m = 0; m < 4; m++)
        {
          coef2D[i][j][k][m] = tmp2D[idx];
          idx++;
        }
      }
    }
  }

  idx = 0;
  std::vector<vtkm::FloatDefault> arr_coeff2D(nz*nr*4*4);
  for (int i = 0; i < nz; i++)
    for (int j = 0; j < nr; j++)
      for (int k = 0; k < 4; k++)
        for (int m = 0; m < 4; m++)
        {
          arr_coeff2D[idx] = coef2D[i][j][k][m];
          idx++;
        }

  ds.AddField(vtkm::cont::make_Field("coeff_2D",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr_coeff2D,
                                     vtkm::CopyFlag::On));


  //Reorder coeff2d to see if it's wrong...
  idx = 0;
  std::vector<std::vector<std::vector<std::vector<double>>>> coeff_2D;
  coeff_2D.resize(nz);
  for (int i = 0; i < nz; i++)
  {
    coeff_2D[i].resize(nr);
    for (int j = 0; j < nr; j++)
    {
      coeff_2D[i][j].resize(4);
      for (int k = 0; k < 4; k++)
      {
        coeff_2D[i][j][k].resize(4);
        for (int m = 0; m < 4; m++)
        {
          coeff_2D[i][j][k][m] = tmp2D[idx];
          idx++;
        }
      }
    }
  }


  //Put this on a 2D grid for debugging...
  vtkm::Vec2f origin2D(xgcParams.eq_min_r, xgcParams.eq_min_z);
  vtkm::Vec2f spacing2D((xgcParams.eq_max_r-xgcParams.eq_min_r)/double(xgcParams.eq_mr-1), (xgcParams.eq_max_z-xgcParams.eq_min_z)/double(xgcParams.eq_mz-1));
  auto ds2D = vtkm::cont::DataSetBuilderUniform::Create(vtkm::Id2(xgcParams.eq_mr, xgcParams.eq_mr),
                                                        origin2D, spacing2D);

  std::vector<vtkm::FloatDefault> c00;
  std::vector<std::vector<std::vector<vtkm::FloatDefault>>> cij(4);
  for (int k = 0; k < 4; k++) cij[k].resize(4);

  for (int i = 0; i < nz; i++)
    for (int j = 0; j < nr; j++)
    {
      c00.push_back(coeff_2D[i][j][0][0]);
      for (int k = 0; k < 4; k++)
        for (int m = 0; m < 4; m++)
          cij[k][m].push_back(coeff_2D[i][j][k][m]);
    }

  //ds2D.AddPointField("c00", c00);
  for (int k = 0; k < 4; k++)
  {
    for (int m = 0; m < 4; m++)
    {
      char nm[32];
      sprintf(nm, "c%d%d", k,m);
      //std::cout<<"Add cij: "<<nm<<" "<<cij[k][m].size()<<std::endl;
      ds2D.AddCellField(nm, cij[k][m]);
    }
  }

#if 0
  //ds2D.AddCellField("bum", cij[0][0]);
  //ds.PrintSummary(std::cout);
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> b3d;
  ds.GetField("eq_psi_rz").GetData().AsArrayHandle(arr);
  ds2D.AddPointField("eq_psi_rz", arr);
  vtkm::io::VTKDataSetWriter writer("psiGrid.vtk");
  writer.WriteDataSet(ds2D);
#endif


//  ds.GetField("B_RZP").GetData().AsArrayHandle(b3d);
//  std::vector<vtkm::Vec3f> B2D(

  /*
  vtkm::io::VTKDataSetWriter writer("psiGrid.vtk");
  writer.WriteDataSet(ds2D);
  vtkm::io::VTKDataSetWriter writer2("grid.vtk");
  writer2.WriteDataSet(ds);
  */


  /*
  auto v = stuff->io.InquireVariable<double>("eq_psi_rz");
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

  std::cout<<"eq_psi_rz: "<<tmp.size()<<std::endl;
  for (int i = eq_mr; i < 2*eq_mr; i++)
    std::cout<<"eq_psi_rz["<<i<<"] = "<<tmp[i]<<std::endl;
  */


  //Let's evaluate the b field.
  //vtkm::FloatDefault R = 2, Z = 0;
}

vtkm::cont::DataSet
ReadMesh(adiosS* meshStuff, XGCParameters& xgcParams)
{
  std::vector<double> rz;
  std::vector<int> conn, nextnode;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
  std::string connNm = "/cell_set[0]/node_connect_list";
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("nextnode"), nextnode, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> coords;
  std::vector<vtkm::Id> connIds;

  //points.
  for (int i = 0; i < xgcParams.numNodes; i++)
  {
    double R = rz[i*2 +0];
    double Z = rz[i*2 +1];

    vtkm::Vec3f ptRZ(R,Z,0);
    coords.push_back(ptRZ);
  }

  //cells
  for (int i = 0; i < xgcParams.numTri*3; i+=3)
  {
    int p0 = conn[i+0];
    int p1 = conn[i+1];
    int p2 = conn[i+2];
    connIds.push_back(p0);
    connIds.push_back(p1);
    connIds.push_back(p2);
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  return dsb.Create(coords, vtkm::CellShapeTagTriangle(), 3, connIds, "coords");
}

vtkm::cont::DataSet
ReadMesh3D(adiosS* meshStuff, bool addExtra, XGCParameters& xgcParams)
{
  std::vector<double> rz;
  std::vector<int> conn, nextnode;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
  std::string connNm = "/cell_set[0]/node_connect_list";
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("nextnode"), nextnode, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> coords;
  std::vector<vtkm::Id> connIds;

  double dPhi = vtkm::TwoPi()/static_cast<double>(xgcParams.numPlanes);

  //points.
  double phi = 0;
  int NP = xgcParams.numPlanes;
  if (addExtra) NP++;
  for (int p = 0; p < NP; p++)
  {
    std::cout<<"ReadMesh3D: phi= "<<phi<<std::endl;
    for (int i = 0; i < xgcParams.numNodes; i++)
    {
      double R = rz[i*2 +0];
      double Z = rz[i*2 +1];

      vtkm::Vec3f pt(R,Z,phi);
      coords.push_back(pt);
    }
    phi += dPhi;
  }

  //cells
  for (int p = 0; p < NP-1; p++)
  {
    for (int i = 0; i < xgcParams.numTri*3; i+=3)
    {
      int off = p*xgcParams.numNodes;
      int p0 = conn[i+0];
      int p1 = conn[i+1];
      int p2 = conn[i+2];
      connIds.push_back(p0+off);
      connIds.push_back(p1+off);
      connIds.push_back(p2+off);

      off = (p+1)*(xgcParams.numNodes);
      connIds.push_back(p0 + off);
      connIds.push_back(p1 + off);
      connIds.push_back(p2 + off);
    }
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  return dsb.Create(coords, vtkm::CellShapeTagWedge(), 6, connIds, "coords");
}

static bool Exists(const std::string& /*fname*/)
{
  return false;
  /*
  std::ifstream ifile;
  ifile.open(fname);
  if (ifile)
    return true;
  return false;
  */
}

void WriteHeader(const std::string& fname, const std::string& header)
{
  std::ofstream f(fname, std::ofstream::out);
  f<<header<<std::endl;
  f.close();
}

void
SaveOutput(const std::vector<std::vector<vtkm::Vec3f>>& traces,
           const vtkm::cont::ArrayHandle<vtkm::Vec2f>& outRZ,
           const vtkm::cont::ArrayHandle<vtkm::Vec2f>& outTP,
           const vtkm::cont::ArrayHandle<vtkm::Id>& outID,
           const std::string& outFileName = "")
{
  std::string tracesNm, puncNm, puncThetaPsiNm, adiosNm;
  if (outFileName.empty())
  {
    tracesNm = "./traces.txt";
    puncNm = "./punctures.vtk";
    puncThetaPsiNm = "./punctures.theta_psi.vtk";
    adiosNm = "./punctures.bp";
  }
  else
  {
    tracesNm = outFileName + ".traces.txt";
    puncNm = outFileName + ".punc.vtk";
    puncThetaPsiNm = outFileName + ".punc.theta_psi.vtk";
    adiosNm = outFileName + ".adios.bp";
  }

  bool tExists = Exists(tracesNm);
  bool pExists = Exists(puncNm);
  bool ptpExists = Exists(puncThetaPsiNm);
  std::cout<<"EXISTS: "<<tExists<<" "<<pExists<<" "<<ptpExists<<std::endl;

  //Write headers.
  if (!tExists)
    WriteHeader(tracesNm, "ID, R, Z, T");

  std::ofstream outTraces;
  outTraces.open(tracesNm, std::ofstream::app);
  //write traces
  for (int i = 0; i < (int)traces.size(); i++)
  {
    for (const auto& pt : traces[i])
    {
      auto R = pt[0];
      auto Z = pt[2];
      auto PHI_N = pt[1];
      while (PHI_N < 0)
        PHI_N += vtkm::TwoPi();

      outTraces<<i<<", "<<R<<", "<<Z<<", "<<PHI_N<<std::endl;
    }
  }
  /*
  vtkm::io::VTKDataSetWriter writer1(puncNm), writer2(puncThetaPsiNm);
  writer1.WriteDataSet(dsRZ);
  writer2.WriteDataSet(dsTP);
  */

  //save to adios
  std::vector<vtkm::FloatDefault> ptsRZ(100);
  adios2::ADIOS adiosW;
  adios2::IO io = adiosW.DeclareIO("io");

  std::size_t nPts = static_cast<std::size_t>(outRZ.GetNumberOfValues());

  //create a (R,Z,ID) and (Theta, Psi, ID) arrays.
  auto RZBuff = vtkm::cont::ArrayHandleBasic<vtkm::Vec2f>(outRZ).GetReadPointer();
  auto TPBuff = vtkm::cont::ArrayHandleBasic<vtkm::Vec2f>(outTP).GetReadPointer();
  auto IDBuff = vtkm::cont::ArrayHandleBasic<vtkm::Id>(outID).GetReadPointer();

  std::vector<std::size_t> shape = {nPts*2}, offset = {0}, size = {nPts*2};
  auto vRZ = io.DefineVariable<vtkm::FloatDefault>("RZ", shape, offset, size);
  auto vTP = io.DefineVariable<vtkm::FloatDefault>("ThetaPsi", shape, offset, size);
  std::vector<std::size_t> shape2 = {nPts}, size2 = {nPts};
  auto vID = io.DefineVariable<vtkm::Id>("ID", shape2, offset, size2);

  adios2::Engine bpWriter = io.Open(adiosNm, adios2::Mode::Write);
  bpWriter.Put<vtkm::FloatDefault>(vRZ, &(RZBuff[0][0]));
  bpWriter.Put<vtkm::FloatDefault>(vTP, &(TPBuff[0][0]));
  bpWriter.Put<vtkm::Id>(vID, IDBuff);
  bpWriter.Close();

  /*
  //write punctures
  for (int i = 0; i < (int)punctures.size(); i++)
  {
    for (const auto& p : punctures[i])
      outPunc<<std::setprecision(12)<<i<<", "<<p[0]<<", "<<p[2]<<std::endl;

    for (const auto& p : puncturesTP[i])
      outPuncThetaPsi<<std::setprecision(12)<<i<<", "<<p[0]<<", "<<p[1]<<std::endl;
  }
  */
}

/*
vtkm::cont::ArrayHandle<vtkm::Particle>
CreateParticles(const std::vector<vtkm::Vec3f>& pts,
                std::map<std::string, std::vector<std::string>>& args)
{
  vtkm::cont::ArrayHandle<vtkm::Particle> seeds;
  vtkm::Id n = static_cast<vtkm::Id>(pts.size());
  seeds.Allocate(n);
  auto sPortal = seeds.WritePortal();

  bool isPsiRange = false;
  int numPts = -1, numTheta = -1;
  if (args.find("--psiRange") != args.end())
  {
    auto vals = args["--psiRange"];
    numPts = std::atoi(vals[2].c_str());
    numTheta = std::atoi(vals[3].c_str());
    if (numTheta > 1)
      isPsiRange = true;
  }

  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Id pid = i;

    vtkm::Particle p(pts[i], pid);
    sPortal.Set(i, p);
  }

  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
    s.push_back(vtkm::Particle(pts[i], i));
  auto seeds = vtkm::cont::make_ArrayHandle(s, vtkm::CopyFlag::On);
}
*/

void
Poincare(const vtkm::cont::DataSet& ds,
         XGCParameters& xgcParams,
         vtkm::cont::ArrayHandle<vtkm::Particle>& seeds,
         std::map<std::string, std::vector<std::string>>& args)
{

  //Get the arguments.
  vtkm::FloatDefault stepSize = std::atof(args["--stepSize"][0].c_str());
  vtkm::Id numPunc = std::atoi(args["--numPunc"][0].c_str());

  bool useTraces = false;
  if (args.find("--traces") != args.end()) useTraces = std::atoi(args["--traces"][0].c_str());
  std::string outFileName = args["--output"][0];
  bool useBOnly = false;
  if (args.find("--useBOnly") != args.end()) useBOnly = true;
  bool useLinearB = false;
  if (args.find("--useLinearB") != args.end()) useLinearB = true;

  if (useLinearB)
  {
    useBOnly = true;
    std::cout<<"Warning: Using linear B, forcing UseBOnly = true."<<std::endl;
  }

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_ff, coeff_1D, coeff_2D, psi;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> B_rzp, B_Norm_rzp, dAs_ff_rzp;
  //ds.GetField("As_phi_ff").GetData().AsArrayHandle(As_ff);
  std::cout<<"DRP: Fix me?? "<<__FILE__<<" "<<__LINE__<<std::endl;
  ds.GetField("As_ff").GetData().AsArrayHandle(As_ff);
  ds.GetField("dAs_ff_rzp").GetData().AsArrayHandle(dAs_ff_rzp);
  ds.GetField("coeff_1D").GetData().AsArrayHandle(coeff_1D);
  ds.GetField("coeff_2D").GetData().AsArrayHandle(coeff_2D);
  ds.GetField("psi2D").GetData().AsArrayHandle(psi);
  ds.GetField("B_RZP").GetData().AsArrayHandle(B_rzp);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> tracesArr;
  vtkm::cont::ArrayHandle<vtkm::Vec2f> outRZ, outTP;
  vtkm::cont::ArrayHandle<vtkm::Id> outID;
  vtkm::Id nSeeds = seeds.GetNumberOfValues();

  outRZ.Allocate(numPunc*nSeeds);
  outTP.Allocate(numPunc*nSeeds);
  outID.Allocate(numPunc*nSeeds);

  auto start = std::chrono::steady_clock::now();
  RunPoincare2(ds, seeds, xgcParams,
               stepSize, numPunc, useBOnly, useTraces, useLinearB,
               As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
               B_rzp, psi,
               tracesArr, outRZ, outTP, outID);
  auto end = std::chrono::steady_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;

  std::cout<<"PoincareTime= "<<elapsed_seconds.count()<<std::endl;
  //std::cout<<"vtkm::cont::Timer= "<<timer.GetElapsedTime()<<std::endl;
  std::cout<<"outputRZ.size()= "<<outRZ.GetNumberOfValues()<<std::endl;
  std::cout<<"outputTP.size()= "<<outTP.GetNumberOfValues()<<std::endl;
  std::cout<<"punctureID.size()= "<<outID.GetNumberOfValues()<<std::endl;
  //vtkm::cont::printSummary_ArrayHandle(output, std::cout);
  //vtkm::cont::printSummary_ArrayHandle(punctureID, std::cout);

  //Save output
  std::vector<std::vector<vtkm::Vec3f>> traces;

  std::vector<std::vector<vtkm::Vec3f>> res;
  if (useTraces)
  {
    auto portal = tracesArr.ReadPortal();
    vtkm::Id n = portal.GetNumberOfValues();
    for (vtkm::Id i = 0; i < n; i++)
    {
      auto v = portal.Get(i);
      if (v[2] > -1)
        traces[0].push_back(v);
    }
  }

  SaveOutput(traces, outRZ, outTP, outID, outFileName);
}

vtkm::cont::ArrayHandle<vtkm::FloatDefault>
CalcAs(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& As, XGCParameters& xgcParams)
{
//  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As;

//  ds.GetField("As_phi_ff").GetData().AsArrayHandle(As);
  vtkm::Id numAs = As.GetNumberOfValues();
  auto asPortal = As.ReadPortal();

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::FloatDefault>>> As_arr;
  As_arr.resize(xgcParams.numPlanes);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    As_arr[p].resize(2);
    for (int i = 0; i < 2; i++)
      As_arr[p][i].resize(xgcParams.numNodes);
  }

  //Do some easy indexing.
  vtkm::Id idx = 0;
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    for (int n = 0; n < xgcParams.numNodes; n++)
    {
      //vtkm::Vec3f cBhat = cbPortal.Get(n);
      for (int i = 0; i < 2; i++)
      {
        auto as = asPortal.Get(idx);
        //std::cout<<"As: "<<as<<" cBhat= "<<cBhat<<" :: "<<vtkm::Magnitude(cBhat)<<"  val= "<<val<<std::endl;
        //As_curlBhat[p][i][n] = val;
        As_arr[p][i][n] = as;
        idx++;
      }
    }
  }

  //flatten to 1d index.
  idx = 0;
  std::vector<vtkm::FloatDefault> arrAs(numAs);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < xgcParams.numNodes; n++)
      {
        vtkm::FloatDefault asVal = As_arr[p][i][n];
        arrAs[idx] = asVal;
        idx++;
      }
    }
  }

  return vtkm::cont::make_ArrayHandle(arrAs, vtkm::CopyFlag::On);

  /*
  ds.AddField(vtkm::cont::make_Field("As_ff",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arrAs,
                                     vtkm::CopyFlag::On));
  */
}

vtkm::cont::ArrayHandle<vtkm::Vec3f>
Calc_dAs(const vtkm::cont::ArrayHandle<vtkm::FloatDefault>& dAsAH, XGCParameters& xgcParams)
{
  //Assumed that dAs_phi_ff is R,Z,Phi
  //vtkm::cont::ArrayHandle<vtkm::FloatDefault> dAsAH;
  //ds.GetField("dAs_phi_ff").GetData().AsArrayHandle(dAsAH);
  int nVals = dAsAH.GetNumberOfValues()/3;

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> dAs_ff;
  dAs_ff.resize(xgcParams.numPlanes);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    dAs_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      dAs_ff[p][i].resize(xgcParams.numNodes);
  }


  int idx = 0;
  auto dAsPortal = dAsAH.ReadPortal();

  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    for (int n = 0; n < xgcParams.numNodes; n++)
    {
      for (int c = 0; c < 3; c++) //vec coords
      {
        for (int i = 0; i < 2; i++) //idx
        {
          vtkm::FloatDefault val = dAsPortal.Get(idx);
          idx++;
          int cc = c;

          /*
          //Change to R,Phi,Z
          //swap phi and z
          if (c == 1) cc = 2;
          else if (c == 2) cc = 1;
          */

          dAs_ff[p][i][n][cc] = val;
        }
      }
    }
  }

  //
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> VEC_ff;
  VEC_ff.resize(xgcParams.numPlanes);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    VEC_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      VEC_ff[p][i].resize(xgcParams.numNodes);
  }

  idx = 0;
  std::vector<vtkm::Vec3f> arr_ff(nVals);
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < xgcParams.numNodes; n++)
      {
        vtkm::Vec3f val = VEC_ff[p][i][n];
        val = dAs_ff[p][i][n];
        arr_ff[idx] = val;
        idx++;
      }
    }
  }

  //do some testing..
  for (int p = 0; p < xgcParams.numPlanes; p++)
  {
    int off0 = p*xgcParams.numNodes*2;
    int off1 = p*xgcParams.numNodes*2 + xgcParams.numNodes;
    for (int n = 0; n < xgcParams.numNodes; n++)
    {
      auto V0 = dAs_ff[p][0][n];
      auto V1 = dAs_ff[p][1][n];

      auto X0 = arr_ff[off0 + n];
      auto X1 = arr_ff[off1 + n];
      auto diff = vtkm::Magnitude(V0-X0) + vtkm::Magnitude(V1-X1);
      if (diff > 0)
        std::cout<<"ERROR: "<<V0<<" : "<<V1<<"  diff0= "<<(V0-X0)<<" diff1= "<<(V1-X1)<<std::endl;

    }
  }

  return vtkm::cont::make_ArrayHandle(arr_ff, vtkm::CopyFlag::On);
  /*
  ds.AddField(vtkm::cont::make_Field("dAs_ff_rzp",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr_ff,
                                     vtkm::CopyFlag::On));
  */
}

template <typename T>
void
UpdateField(vtkm::cont::DataSet& ds,
            const std::string& vName,
            vtkm::cont::ArrayHandle<T> &var)
{
  if (ds.HasField(vName))
  {
    vtkm::cont::ArrayHandle<T> oldVar;
    ds.GetField(vName).GetData().AsArrayHandle(oldVar);

    vtkm::cont::Algorithm::Copy(var, oldVar);
  }
  else
    ds.AddField(vtkm::cont::Field(vName,
                                  vtkm::cont::Field::Association::WHOLE_MESH,
                                  var));
}

void
ReadTurbData(adiosS* turbStuff,
             std::map<std::string, std::vector<std::string>>& args,
             vtkm::cont::ArrayHandle<vtkm::FloatDefault>& As_arr,
             vtkm::cont::ArrayHandle<vtkm::FloatDefault>& dAs_arr)
{
  auto asV = turbStuff->io.InquireVariable<double>("As_phi_ff");
  auto dAsV = turbStuff->io.InquireVariable<double>("dAs_phi_ff");

  std::vector<double> arrAs, arrdAs;
  turbStuff->engine.Get(asV, arrAs, adios2::Mode::Sync);
  turbStuff->engine.Get(dAsV, arrdAs, adios2::Mode::Sync);

  bool useTurbulence = true;
  if (args.find("--turbulence") != args.end())
    useTurbulence = std::atoi(args["--turbulence"][0].c_str());
  if (useTurbulence == false)
  {
    for (auto& x : arrAs)
      x = 0.0;
    for (auto& x : arrdAs)
      x = 0.0;
  }

  As_arr = vtkm::cont::make_ArrayHandle(arrAs, vtkm::CopyFlag::On);
  dAs_arr = vtkm::cont::make_ArrayHandle(arrdAs, vtkm::CopyFlag::On);
  //auto AsPhi = vtkm::cont::make_ArrayHandle(arrAs, vtkm::CopyFlag::On);
  //auto dAsPhi = vtkm::cont::make_ArrayHandle(arrdAs, vtkm::CopyFlag::On);
  //UpdateField(ds, "As_phi_ff", AsPhi);
  //UpdateField(ds, "dAs_phi_ff", dAsPhi);
}

vtkm::cont::DataSet
ReadStaticData(std::map<std::string, std::vector<std::string>>& args,
               XGCParameters& xgcParams,
               const std::string& coeffFile)
{
  if (adios != nullptr)
  {
    std::cerr<<"Re-reading static data!!!"<<std::endl;
  }

  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  adios = new adios2::ADIOS;

  adiosStuff["mesh"] = new adiosS(adios, "xgc.mesh.bp", "mesh", adiosArgs);
  adiosStuff["equil"] = new adiosS(adios, "xgc.equil.bp", "equil", adiosArgs);
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "bfield", adiosArgs);
  adiosStuff["coeff"] = new adiosS(adios, coeffFile, "coeff", adiosArgs);

  auto meshStuff = adiosStuff["mesh"];
  auto equilStuff = adiosStuff["equil"];
  auto coeffStuff = adiosStuff["coeff"];
  auto bfieldStuff = adiosStuff["bfield"];

  meshStuff->engine.BeginStep();
  equilStuff->engine.BeginStep();
  coeffStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();

  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &xgcParams.numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &xgcParams.numTri, adios2::Mode::Sync);
  std::vector<double> psiVals;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("psi"), psiVals, adios2::Mode::Sync);
  xgcParams.psi_min = xgcParams.psi_max = psiVals[0];
  for (const auto& p : psiVals)
  {
    if (p < xgcParams.psi_min) xgcParams.psi_min = p;
    if (p > xgcParams.psi_max) xgcParams.psi_max = p;
  }
  std::cout<<"********                  PSI m/M = "<<xgcParams.psi_min<<" "<<xgcParams.psi_max<<std::endl;

  auto ds = ReadMesh(meshStuff, xgcParams);
  ReadPsiInterp(equilStuff, coeffStuff, ds, xgcParams);
  ReadScalar(meshStuff, xgcParams, ds, "psi");
  ReadB(bfieldStuff, xgcParams, ds);

  meshStuff->engine.EndStep();
  equilStuff->engine.EndStep();
  coeffStuff->engine.EndStep();
  bfieldStuff->engine.EndStep();

  return ds;
}

vtkm::cont::DataSet
ReadDataSet_ORIG(std::map<std::string, std::vector<std::string>>& args, XGCParameters& xgcParams)
{
  auto ds = ReadStaticData(args, xgcParams, "xgc.bfield-all.bp");

  //Get the data...
  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", adiosArgs);
  adiosStuff["bfield-all"] = new adiosS(adios, "xgc.bfield-all.bp", "Bs", adiosArgs);

  auto dataStuff = adiosStuff["data"];
  auto bfieldStuff = adiosStuff["bfield-all"];

  dataStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();

  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &xgcParams.numPlanes, adios2::Mode::Sync);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_arr, dAs_arr;
  ReadTurbData(bfieldStuff, args, As_arr, dAs_arr);

  dataStuff->engine.EndStep();
  bfieldStuff->engine.EndStep();

  auto As = CalcAs(As_arr, xgcParams);
  auto dAs = Calc_dAs(dAs_arr, xgcParams);
  UpdateField(ds, "As_ff", As);
  UpdateField(ds, "dAs_ff_rzp", dAs);

  return ds;


#if 0
  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
  adios = new adios2::ADIOS;
  adiosStuff["mesh"] = new adiosS(adios, "xgc.mesh.bp", "mesh", adiosArgs);
  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", adiosArgs);
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "/node_data[0]/values", adiosArgs);
  adiosStuff["bfield-all"] = new adiosS(adios, "xgc.bfield-all.bp", "Bs", adiosArgs);
  //adiosStuff["psi_bicub_acoef"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "psi_bicub_acoef", adiosArgs);
  //adiosStuff["one_d_cub_psi_acoef"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "one_d_cub_acoef", adiosArgs);
  adiosStuff["interp_coeff"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "interp_coeff", adiosArgs);
  adiosStuff["equil"] = new adiosS(adios, "xgc.equil.bp", "equil", adiosArgs);

  auto meshStuff = adiosStuff["mesh"];
  auto dataStuff = adiosStuff["data"];
  auto bfieldStuff = adiosStuff["bfield"];
  auto bfield_allStuff = adiosStuff["bfield-all"];
  auto interp_coeffStuff = adiosStuff["interp_coeff"];
  auto equilStuff = adiosStuff["equil"];
  //auto one_dcub_acoefStuff = adiosStuff["one_d_cub_psi_acoef"];

  meshStuff->engine.BeginStep();
  dataStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();
  bfield_allStuff->engine.BeginStep();
  interp_coeffStuff->engine.BeginStep();
  equilStuff->engine.BeginStep();
  //one_dcub_acoefStuff->engine.BeginStep();

  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &xgcParams.numPlanes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &xgcParams.numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &xgcParams.numTri, adios2::Mode::Sync);
  std::vector<double> psiVals;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("psi"), psiVals, adios2::Mode::Sync);
  xgcParams.psi_min = xgcParams.psi_max = psiVals[0];
  for (const auto& p : psiVals)
  {
    if (p < xgcParams.psi_min) xgcParams.psi_min = p;
    if (p > xgcParams.psi_max) xgcParams.psi_max = p;
  }
  std::cout<<"********                  PSI m/M = "<<xgcParams.psi_min<<" "<<xgcParams.psi_max<<std::endl;

  //Try and do everything in cylindrical coords and worklets.
  auto ds = ReadMesh(meshStuff, xgcParams);
  ReadScalar(dataStuff, xgcParams, ds, "dpot");
  ReadScalar(dataStuff, xgcParams, ds, "apars", "apars", true);
  ReadScalar(meshStuff, xgcParams, ds, "psi");
  ReadScalar(meshStuff, xgcParams, ds, "theta");

  ReadOther(bfieldStuff, ds, "As_phi_ff");
  ReadOther(bfieldStuff, ds, "dAs_phi_ff");
  ReadPsiInterp(equilStuff, interp_coeffStuff, ds, xgcParams);

  CalcAs(ds, xgcParams);
  Calc_dAs(ds, xgcParams);


//  vtkm::io::VTKDataSetWriter writer("grid.vtk");
//  writer.WriteDataSet(ds);

  return ds;
#endif
}

vtkm::cont::DataSet
ReadDataSet(std::map<std::string, std::vector<std::string>>& args, XGCParameters& xgcParams)
{
  if (args.find("--test") != args.end())
    return ReadDataSet_ORIG(args, xgcParams);

  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  auto ds = ReadStaticData(args, xgcParams, "xgc.poincare_init.bp");
  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", adiosArgs);

  auto dataStuff = adiosStuff["data"];

  //Go to 2nd step.
//  dataStuff->engine.BeginStep();
//  dataStuff->engine.EndStep();

  dataStuff->engine.BeginStep();
  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &xgcParams.numPlanes, adios2::Mode::Sync);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_arr, dAs_arr;
  ReadTurbData(dataStuff, args, As_arr, dAs_arr);
  dataStuff->engine.EndStep();

  auto As = CalcAs(As_arr, xgcParams);
  auto dAs = Calc_dAs(dAs_arr, xgcParams);

  //CalcAs(ds, xgcParams);
  //Calc_dAs(ds, xgcParams);
  UpdateField(ds, "As_ff", As);
  UpdateField(ds, "dAs_ff_rzp", dAs);

  return ds;

#if 0
  if (args.find("--Step") == args.end())
    return ReadDataSet_ORIG(args, xgcParams);

  int step = std::stoi(args["--Step"][0]);

  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  adios = new adios2::ADIOS;
  adiosStuff["init"] = new adiosS(adios, "xgc.poincare_init.bp", "init", adiosArgs);
  adiosStuff["mesh"] = new adiosS(adios, "xgc.mesh.bp", "mesh", adiosArgs);
  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", adiosArgs);
  adiosStuff["equil"] = new adiosS(adios, "xgc.equil.bp", "equil", adiosArgs);
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "bfield", adiosArgs);

  auto initStuff = adiosStuff["init"];
  auto meshStuff = adiosStuff["mesh"];
  auto dataStuff = adiosStuff["data"];
  auto equilStuff = adiosStuff["equil"];
  auto bfieldStuff = adiosStuff["bfield"];

  initStuff->engine.BeginStep();
  meshStuff->engine.BeginStep();
  equilStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();

  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &xgcParams.numPlanes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &xgcParams.numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &xgcParams.numTri, adios2::Mode::Sync);
  std::vector<double> psiVals;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("psi"), psiVals, adios2::Mode::Sync);
  xgcParams.psi_min = xgcParams.psi_max = psiVals[0];
  for (const auto& p : psiVals)
  {
    if (p < xgcParams.psi_min) xgcParams.psi_min = p;
    if (p > xgcParams.psi_max) xgcParams.psi_max = p;
  }
  std::cout<<"********                  PSI m/M = "<<xgcParams.psi_min<<" "<<xgcParams.psi_max<<std::endl;

  auto ds = ReadMesh(meshStuff, xgcParams);
  ReadPsiInterp(equilStuff, initStuff, ds, xgcParams);

  adios2::StepStatus status;
  for (int i = 0; i < step; i++)
  {
    status = dataStuff->engine.BeginStep();
    if (status != adios2::StepStatus::OK)
    {
      std::cout<<"********* Failed to get step."<<std::endl;
      return ds;
    }

    dataStuff->engine.EndStep();
  }

  ReadOther(dataStuff, ds, "As_phi_ff");
  ReadOther(dataStuff, ds, "dAs_phi_ff");
  ReadScalar(dataStuff, xgcParams, ds, "dpot");
  ReadScalar(meshStuff, xgcParams, ds, "psi");

  CalcAs(ds, xgcParams);
  Calc_dAs(ds, xgcParams);

  initStuff->engine.EndStep();
  meshStuff->engine.EndStep();
  equilStuff->engine.EndStep();
  bfieldStuff->engine.EndStep();

//  vtkm::io::VTKDataSetWriter writer("grid.vtk");
//  writer.WriteDataSet(ds);

  return ds;
#endif
}

int oneDBlocks = 16;
int threadsPerBlock = 16;
#ifdef VTKM_CUDA
vtkm::cont::cuda::ScheduleParameters
mySchedParams(char const* name,
              int major,
              int minor,
              int multiProcessorCount,
              int maxThreadsPerMultiProcessor,
              int maxThreadsPerBlock)
{
  vtkm::cont::cuda::ScheduleParameters p;
  p.one_d_blocks = oneDBlocks;
  p.one_d_threads_per_block = threadsPerBlock;

  return p;
}
#endif


template <typename Coeff_2DType>
VTKM_EXEC
void EvalBicub2(const vtkm::FloatDefault& x,
                const vtkm::FloatDefault& y,
                const vtkm::FloatDefault& xc,
                const vtkm::FloatDefault& yc,
                const vtkm::Id& offset,
                const Coeff_2DType& Coeff_2D,
                vtkm::FloatDefault &f00, vtkm::FloatDefault &f10, vtkm::FloatDefault &f01,
                vtkm::FloatDefault &f11, vtkm::FloatDefault &f20, vtkm::FloatDefault &f02)
{
  vtkm::FloatDefault dx = x - xc;
  vtkm::FloatDefault dy = y - yc;

  //fortran code.

  f00 = f01 = f10 = f11 = f20 = f02 = 0.0f;
  vtkm::FloatDefault xv[4] = {1, dx, dx*dx, dx*dx*dx};
  vtkm::FloatDefault yv[4] = {1, dy, dy*dy, dy*dy*dy};
  vtkm::FloatDefault fx[4] = {0,0,0,0};
  vtkm::FloatDefault dfx[4] = {0,0,0,0};
  vtkm::FloatDefault dfy[4] = {0,0,0,0};
  vtkm::FloatDefault dfx2[4] = {0,0,0,0};
  vtkm::FloatDefault dfy2[4] = {0,0,0,0};


  for (int j=0; j<4; j++)
  {
    for (int i=0; i<4; i++)
      fx[j] = fx[j] + xv[i]*Coeff_2D.Get(offset + i*4 + j); //acoeff[i][j];
    for (int i=1; i<4; i++)
      dfx[j] = dfx[j] + vtkm::FloatDefault(i)*xv[i-1]*Coeff_2D.Get(offset + i*4 + j); //acoeff[i][j];
    for (int i=2; i<4; i++)
      dfx2[j] = dfx2[j] + vtkm::FloatDefault(i*(i-1))*xv[i-2]*Coeff_2D.Get(offset + i*4 + j); //acoeff[i][j];
  }

  for (int j = 0; j < 4; j++)
  {
    f00 = f00 + fx[j]*yv[j];
    f10 = f10 + dfx[j]*yv[j];
    f20 = f20 + dfx2[j]*yv[j];
  }

  for (int j = 1; j < 4; j++)
  {
    dfy[j] = vtkm::FloatDefault(j)*yv[j-1];
    f01 = f01 + fx[j]*dfy[j];
    f11 = f11 + dfx[j]*dfy[j];
  }

  for (int j = 2; j < 4; j++)
  {
    dfy2[j] = vtkm::FloatDefault(j*(j-1))*yv[j-2];
    f02 = f02 + fx[j]*dfy2[j];
  }
}

int GetIndex(const vtkm::FloatDefault& x,
             const int& nx,
             const vtkm::FloatDefault& xmin,
             const vtkm::FloatDefault& dx_inv)
{
  int idx = std::max(1, std::min(nx  , 1 + int ((x-xmin)*dx_inv)) );
  return idx-1;
}


template <typename CoeffType>
vtkm::FloatDefault
InterpolatePsi(const vtkm::Vec2f& ptRZ,
               const CoeffType& coeff,
               const vtkm::Id& ncoeff,
               const vtkm::Id2& nrz,
               const vtkm::Vec2f& rzmin,
               const vtkm::Vec2f& drz)
{
  vtkm::FloatDefault psi = 0;

  int r_i = GetIndex(ptRZ[0], nrz[0], rzmin[0], 1.0/drz[0]);
  int z_i = GetIndex(ptRZ[1], nrz[1], rzmin[1], 1.0/drz[1]);
  vtkm::FloatDefault Rc = rzmin[0] + (vtkm::FloatDefault)(r_i)*drz[0];
  vtkm::FloatDefault Zc = rzmin[1] + (vtkm::FloatDefault)(z_i)*drz[1];
  auto Rc_1 = Rc + drz[0];
  auto Zc_1 = Zc + drz[1];
  Rc = (Rc + Rc_1) * 0.5;
  Zc = (Zc + Zc_1) * 0.5;

  vtkm::Id offset = (r_i * ncoeff + z_i) * 16;

  vtkm::FloatDefault dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
  EvalBicub2(ptRZ[0], ptRZ[1], Rc, Zc, offset, coeff, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);

  return psi;
}

std::vector<vtkm::Particle>
GeneratePsiRangeSeeds(vtkm::FloatDefault psiNorm0,
                      vtkm::FloatDefault psiNorm1,
                      int numPts,
                      int numTheta,
                      const vtkm::cont::DataSet& ds,
                      XGCParameters& xgcParams)

{
  std::vector<vtkm::Particle> seeds;

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxR, thetas;
  FindMaxR(ds, xgcParams, numTheta, thetas, maxR);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
  ds.GetField("coeff_2D").GetData().AsArrayHandle(arr);
  auto coeff = arr.ReadPortal();

  auto thetaPortal = thetas.ReadPortal();
  auto maxRPortal = maxR.ReadPortal();
  vtkm::Id ncoeff = xgcParams.eq_mr-1;
  vtkm::Id2 nrz(xgcParams.eq_mr, xgcParams.eq_mz);
  vtkm::Vec2f rzmin(xgcParams.eq_min_r, xgcParams.eq_min_z);
  vtkm::Vec2f drz((xgcParams.eq_max_r-xgcParams.eq_min_r)/static_cast<vtkm::FloatDefault>(xgcParams.eq_mr-1),
                  (xgcParams.eq_max_z-xgcParams.eq_min_z)/static_cast<vtkm::FloatDefault>(xgcParams.eq_mz-1));

  auto psiMin = psiNorm0 * xgcParams.eq_x_psi, psiMax = psiNorm1 * xgcParams.eq_x_psi;
  auto dPsi = (psiMax-psiMin) / static_cast<vtkm::FloatDefault>(numPts-1);
  auto psi = psiMin;

  for (int i = 0; i < numPts; i++)
  {
    //std::cout<<"Pt_"<<i<<" = psi= "<<psi<<" psiN= "<<(psi/xgcParams.eq_x_psi)<<std::endl;

    vtkm::Id pid = static_cast<vtkm::Id>(i);

    for (vtkm::Id j = 0; j < thetaPortal.GetNumberOfValues(); j++)
    {
      vtkm::FloatDefault rad0 = 1e-8, rad1 = maxRPortal.Get(j);
      auto theta = thetaPortal.Get(j);
      auto cost = vtkm::Cos(theta), sint = vtkm::Sin(theta);

      /*
      vtkm::Vec2f p0(eq_axis_r + rad0*cost, eq_axis_z + rad0*sint);
      vtkm::Vec2f p1(eq_axis_r + rad1*cost, eq_axis_z + rad1*sint);
      auto psi_0 = InterpolatePsi(p0,  coeff, ncoeff, nrz, rzmin, drz);
      auto psi_1 = InterpolatePsi(p1,  coeff, ncoeff, nrz, rzmin, drz);
      std::cout<<"  theta= "<<theta<<" rad= "<<rad0<<" "<<rad1<<" psi= "<<psi_0<<" "<<psi_1<<std::endl;
      std::cout<<"     p0= "<<p0<<" "<<p1<<std::endl;
      std::cout<<"      c/s = "<<cost<<" "<<sint<<std::endl;
      */

      do
      {
        vtkm::FloatDefault rmid = (rad0+rad1) / 2.0;

        vtkm::Vec2f pt(xgcParams.eq_axis_r + rmid*cost, xgcParams.eq_axis_z + rmid*sint);
        auto psimid = InterpolatePsi(pt,  coeff, ncoeff, nrz, rzmin, drz);
        //std::cout<<"     "<<rmid<<" ("<<rad0<<" "<<rad1<<") psimid= "<<psimid<<std::endl;

        if (psimid < psi) // mid is inside, set range to (radmid, rad1)
          rad0 = rmid;
        else
          rad1 = rmid;   // mid is outside, set range to (rad0, radmid)
      }
      while ((rad1-rad0) > 1e-8);

      vtkm::FloatDefault r = (rad0+rad1)/2.0;
      vtkm::Vec3f pt_rpz(xgcParams.eq_axis_r + r*cost, 0, xgcParams.eq_axis_z + r*sint);
      seeds.push_back({pt_rpz, pid});
      //std::cout<<"   PSI_pt: "<<psi<<" "<<theta<<" id= "<<pid<<std::endl;

      //auto psi_ = InterpolatePsi({pt[0], pt[1]},  coeff, ncoeff, nrz, rzmin, drz);
      //std::cout<<"   "<<pt<<"  psi= "<<psi_/eq_x_psi<<std::endl;
    }
    psi += dPsi;
  }

  return seeds;
}

vtkm::cont::ArrayHandle<vtkm::Particle>
GenerateSeeds(const vtkm::cont::DataSet& ds,
              XGCParameters& xgcParams,
              std::map<std::string, std::vector<std::string>>& args)
{
  std::vector<vtkm::Particle> seeds;

  if (args.find("--range") != args.end())
  {
    vtkm::Id numSeeds = std::stoi(args["--numSeeds"][0]);

    auto vals = args["--range"];
    vtkm::FloatDefault r0 = std::atof(vals[0].c_str());
    vtkm::FloatDefault r1 = std::atof(vals[1].c_str());

    vtkm::FloatDefault dr = (r1-r0) / (float)(numSeeds-1);
    vtkm::FloatDefault r = r0;

    for (vtkm::Id id = 0; id < numSeeds; id++, r+=dr)
      seeds.push_back({{r, .1, 0}, id});
    //std::cout<<"SEEDS= "<<seeds<<std::endl;
  }
  else if (args.find("--psiRange") != args.end())
  {
    auto vals = args["--psiRange"];
    vtkm::FloatDefault psi0 = std::atof(vals[0].c_str());
    vtkm::FloatDefault psi1 = std::atof(vals[1].c_str());
    int numPts = std::atoi(vals[2].c_str());
    int numTheta = std::atoi(vals[3].c_str());

    seeds = GeneratePsiRangeSeeds(psi0, psi1, numPts, numTheta, ds, xgcParams);
  }
  else if (args.find("--seed") != args.end())
  {
    auto vals = args["--seed"];
    vtkm::FloatDefault r = std::atof(vals[0].c_str());
    vtkm::FloatDefault z = std::atof(vals[1].c_str());
    vtkm::FloatDefault t = std::atof(vals[2].c_str());
    seeds.push_back({{r, t, z}, 0});
  }
  else if (args.find("--jong1") != args.end())
  {
    //first point in traces.v2
    //zone = 3132, incident nodes: 1521, 1612, 1613
    seeds = {{{3.029365, 6.183185, 0.020600}, 0}};
  }
  else if (args.find("--jong6") != args.end())
  {
    auto vals = args["--jong6"];
    std::cout<<"VALS= "<<vals<<std::endl;

    std::vector<vtkm::Vec3f> allSeeds;
    allSeeds = {
      {3.351443028564415449, 0.0, -0.451648806402756176}, //pt_0, ID= 10670   blue 3 islands
      {3.187329423521033878, 0.0, -0.665017624967372267}, //pt_1, ID= 12000
      {1.992020349316277139, 0.0, -0.126203396421661285}, //pt_2, ID= 13100
      {3.018666196722858963, 0.0, 0.073864239629065770}, //pt_3, ID= 0
      {3.176582679765305173, 0.0, -0.220557108925872658}, //pt_4, ID= 4000    stochastic problem area
      {2.129928300491922499, 0.0, -0.176508860570331577}, //pt_5, ID= 10153  semi stoch. xgc spread out, vtkm thin. good with 0.001
      {2.568671712782164995, 0.0, 0.050249128799423198}, //pt_6, ID= 100
      {2.934624677179501262, 0.0, 0.220686132855778427}, //pt_7, ID= 500
      {2.959288366738244580, 0.0, 0.448869653975662142}, //pt_8, ID= 5000
    };

    if (vals.size() == 0) //all the points.
    {
      std::size_t n = allSeeds.size();
      for (std::size_t i = 0; i < n; i++)
        seeds.push_back({allSeeds[i], static_cast<vtkm::Id>(i)});
    }
    else
    {
      int n = std::stoi(vals[0]);
      if (n >= (int)allSeeds.size())
        std::cout<<"Bad seed!!! #allseeds= "<<allSeeds.size()<<std::endl;

      seeds = {{allSeeds[n], 0}};
    }
  }
  else if (args.find("--jongrz") != args.end())
  {
    //Seed from the blue island.
    //seeds = {{3.351443,             0, -0.451649}};
    seeds = {{{3.351443028564415449, 0, -0.45164880640275617552}, 0}};
  }
  else if (args.find("--afterN") != args.end())
  {
    seeds = {
      {{3.321888620239255019, 0.0, 0.478933972623594384}, 0},
      {{2.568934684085571352, 0.0, 0.731290913908178353}, 1},
      {{3.493628202658771720, 0.0, 0.433951677589735296}, 2},
      {{2.862485694515508605, 0.0, 0.208737305948038576}, 3},
      {{2.905837753215041008, 0.0, -0.397811882628356817}, 4},
      {{3.391834939600261389, 0.0, -0.350011953142094434}, 5}
    };
  }
  else if (args.find("--parse") != args.end())
  {
    //Generaate the seed list by running jongAll.py
    std::cout<<"READING: "<<args["--parse"][0]<<std::endl;
    std::ifstream seedFile;
    seedFile.open(args["--parse"][0]);
    std::string line;
    vtkm::Id pid = 0;
    while (std::getline(seedFile, line))
    {
      vtkm::FloatDefault r, p, z;
#ifdef VTKM_USE_DOUBLE_PRECISION
      sscanf(line.c_str(), "%lf, %lf, %lf", &r, &p, &z);
#else
      sscanf(line.c_str(), "%f, %f, %f", &r, &p, &z);
#endif
      seeds.push_back({{r,p,z}, pid});
      pid++;
    }
  }

  auto seedsArray = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);

  return seedsArray;
}

void
StreamingPoincare(std::map<std::string, std::vector<std::string>>& args)
{
  XGCParameters xgcParams;

  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];
  auto ds = ReadStaticData(args, xgcParams, "xgc.poincare_init.bp");

  std::string fName = args["--streaming"][0];
  std::string outFileNameBase = args["--output"][0];

  adiosStuff["data"] = new adiosS(adios, fName, "3d", adiosArgs);
  auto dataStuff = adiosStuff["data"];

  adios2::StepStatus status;
  int step = 0;
  while (true)
  {
    status = dataStuff->engine.BeginStep(adios2::StepMode::Read, 1000.0);
    if (status != adios2::StepStatus::OK)
    {
      std::cerr<<"Failed to read "<<fName<<" step "<<step<<". Exiting now"<<std::endl;
      break;
    }

    //Initialize the num planes variable.
    if (step == 0)
      dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &xgcParams.numPlanes, adios2::Mode::Sync);

    int timeStep;
    dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("tindex"), &timeStep, adios2::Mode::Sync);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_arr, dAs_arr;
    ReadTurbData(dataStuff, args, As_arr, dAs_arr);
    auto As = CalcAs(As_arr, xgcParams);
    auto dAs = Calc_dAs(dAs_arr, xgcParams);
    UpdateField(ds, "As_ff", As);
    UpdateField(ds, "dAs_ff_rzp", dAs);

    std::cout<<step<<": Read data timeStep= "<<timeStep<<std::endl;
    dataStuff->engine.EndStep();

    auto seeds = GenerateSeeds(ds, xgcParams, args);

    std::string outputFile = outFileNameBase + "." + std::to_string(timeStep);
    std::cout<<"Dump to "<<outputFile<<std::endl;
    args["--output"][0] = outputFile;
    Poincare(ds, xgcParams, seeds, args);

    step++;
  }
}

int
main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank, numRanks;
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::map<std::string, std::vector<std::string>> args;
  ParseArgs(argc, argv, args);

  if (args.find("--gpuParams") != args.end())
  {
    oneDBlocks = std::atoi(args["--gpuParams"][0].c_str());
    threadsPerBlock = std::atoi(args["--gpuParams"][1].c_str());
  }
  std::cout<<"GPUParams: "<<oneDBlocks<<" "<<threadsPerBlock<<std::endl;

#ifdef VTKM_CUDA
  if (args.find("--gpu") != args.end())
    vtkm::cont::cuda::InitScheduleParameters(mySchedParams);
#endif

  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  if (argc < 7)
  {
    std::cerr<<"Usage: "<<argv[0]<<" dataFile numSeeds maxPunctures stepSize poincVar [options]"<<std::endl;
    std::cerr<<config.Usage<<std::endl;
    return -1;
  }

  if (args.find("--gpu") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagCuda{});
  if (args.find("--kokkos") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagKokkos{});
  else if (args.find("--openmp") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP{});
  else if (args.find("--serial") != args.end())
      vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagSerial{});
  else
    vtkm::cont::GetRuntimeDeviceTracker().ForceDevice(vtkm::cont::DeviceAdapterTagCuda{});

  if (args.find("--streaming") != args.end())
  {
    StreamingPoincare(args);
  }
  else
  {
    XGCParameters xgcParams;
    auto ds = ReadDataSet(args, xgcParams);
    //ds.PrintSummary(std::cout);
    //vtkm::io::VTKDataSetWriter writer2("grid.vtk");
    //writer2.WriteDataSet(ds);
    auto seeds = GenerateSeeds(ds, xgcParams, args);
    Poincare(ds, xgcParams, seeds, args);
  }

  return 0;
}

/*
//cyclone case

./examples/poincare/Poincare  --vField B --dir ../data/run_1920 --traces 0 --useHighOrder --turbulence 1 --jong6 --openmp  --output bumm --numPunc 303 --gpuParams 256 128 --stepSize 0.1 --test


//iter case
./examples/poincare/Poincare  --vField B --dir ../data/XGC_GB/test_GB_small_su455 --traces 0 --useHighOrder --turbulence 1 --psiRange .1 .9 4 2 --openmp  --output bumm --numPunc 303 --gpuParams 256 128 --stepSize 0.1

//streaming case
./examples/poincare/Poincare  --vField B --dir ../data/XGC_GB/test_GB_small_su455 --traces 0 --useHighOrder --turbulence 1 --psiRange .1 .9 4 2 --openmp  --output bumm --numPunc 500 --gpuParams 256 128 --stepSize 0.05 --useLinearB --streaming xgc.3d.panout.1.bp

*/
