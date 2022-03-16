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
#include <vtkm/Geometry.h>
#include <vtkm/Matrix.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleExtractComponent.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/Particle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/Contour.h>
#include <vtkm/filter/PointAverage.h>

#include <vtkm/cont/CellLocatorGeneral.h>
//#include <vtkm/cont/CellLocatorTwoLevelTriangle.h>
#include <vtkm/cont/CellLocatorBoundingIntervalHierarchy.h>
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

/*
radius values: max range: 2.5 3.7
for interesting psi: 3.5  3.7
*/

adios2::ADIOS *adios = NULL;
class adiosS;
std::map<std::string, adiosS*> adiosStuff;

bool useTurb=true;
int numPlanes = -1;
int numNodes = -1;
int numTri = -1;
float XScale = 1;
double eq_axis_r = -1, eq_axis_z = -1, eq_x_psi = -1, eq_x_r = -1, eq_x_z = -1;
double eq_min_r = -1, eq_max_r = -1;
double eq_min_z = -1, eq_max_z = 1;
double psi_min = -1, psi_max = -1;
int eq_mr = -1, eq_mz = -1;

using Ray3f = vtkm::Ray<vtkm::FloatDefault, 3, true>;

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
#include "Poincare2.1.h"
#include "ComputeAs.h"
#include "ComputeAsCell.h"
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

  if ((vname == "As_phi_ff" || vname == "dAs_phi_ff")&& useTurb == false)
  {
    std::cout<<"NO TURB. Zero everything out!"<<std::endl;
    for (auto& x : tmp)
      x = 0.0;
  }

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
    for (int i = 0; i < numNodes; i++)
      tmp.push_back(tmp[i]);
  }

  ds.AddField(vtkm::cont::make_FieldPoint(vname+"2D", vtkm::cont::make_ArrayHandle(tmpPlane, vtkm::CopyFlag::On)));
  if (add3D)
    ds.AddField(vtkm::cont::make_FieldPoint(vname, vtkm::cont::make_ArrayHandle(tmp, vtkm::CopyFlag::On)));
}

void
ReadPsiInterp(adiosS* eqStuff,
              adiosS* interpStuff,
              vtkm::cont::DataSet& ds)
{
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mr"), &eq_mr, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<int>("eq_mz"), &eq_mz, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_axis_r"), &eq_axis_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_axis_z"), &eq_axis_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_min_r"), &eq_min_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_max_r"), &eq_max_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_min_z"), &eq_min_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_max_z"), &eq_max_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_psi"), &eq_x_psi, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_r"), &eq_x_r, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_z"), &eq_x_z, adios2::Mode::Sync);
  eqStuff->engine.Get(eqStuff->io.InquireVariable<double>("eq_x_z"), &eq_x_z, adios2::Mode::Sync);

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
  int nr = eq_mr-1, nz = eq_mz-1;
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
  vtkm::Vec2f origin2D(eq_min_r, eq_min_z);
  vtkm::Vec2f spacing2D((eq_max_r-eq_min_r)/double(eq_mr-1), (eq_max_z-eq_min_z)/double(eq_mz-1));
  auto ds2D = vtkm::cont::DataSetBuilderUniform::Create(vtkm::Id2(eq_mr, eq_mr),
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
CreateRZPsiGrid(const vtkm::cont::DataSet& ds)
{
  vtkm::Vec2f origin2D(eq_min_r, eq_min_z);
  vtkm::Vec2f spacing2D((eq_max_r-eq_min_r)/double(eq_mr-1), (eq_max_z-eq_min_z)/double(eq_mz-1));
  auto ds2D = vtkm::cont::DataSetBuilderUniform::Create(vtkm::Id2(eq_mr, eq_mr),
                                                        origin2D, spacing2D);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> coeff_2D;
  ds.GetField("coeff_2D").GetData().AsArrayHandle(coeff_2D);

  std::vector<vtkm::FloatDefault> psi;
  vtkm::Id idx = 0;
  int nr = eq_mr-1, nz = eq_mz-1;
  for (int i = 0; i < nr*nz; i++)
  {
    psi.push_back(coeff_2D.ReadPortal().Get(idx));
    idx += 16;
  }

  ds2D.AddCellField("psiCell", psi);

  vtkm::filter::PointAverage avg;
  avg.SetOutputFieldName("psi");
  avg.SetActiveField("psiCell");
  return avg.Execute(ds2D);
}


//-----------------------------------------------------------------------------
class NormalizeWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn input,
                                FieldOut output);
  using ExecutionSignature = void(_1, _2);
  using InputDomain = _1;

  template <typename VectorType>
  VTKM_EXEC void operator()(const VectorType& input,
                            VectorType& output) const
  {
    output = vtkm::Normal(input);
  }
};

//Sets it as RZP
vtkm::cont::ArrayHandle<vtkm::Vec3f>
ComputeCurl(const vtkm::cont::DataSet& inDS)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B;
  inDS.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;
  vtkm::Id numPts = coords.GetNumberOfValues();
  curlBNorm.Allocate(numPts);
  return curlBNorm;

#if 0
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B;
  inDS.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  inDS.GetField("B_RZP").GetData().AsArrayHandle(B);

  //Get the gradients.
  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetActiveField("B_RZP");
  auto ds = gradient.Execute(inDS);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f, 3>> gradients;
  ds.GetField("Gradients").GetData().AsArrayHandle(gradients);

  vtkm::Id numPts = coords.GetNumberOfValues();
  curlBNorm.Allocate(numPts);

  auto gPortal = gradients.ReadPortal();
  auto cPortal = coords.ReadPortal();
  auto bPortal = B.ReadPortal();
  auto portal = curlBNorm.WritePortal();

  for (vtkm::Id i = 0; i < numPts; i++)
  {
    auto B0 = bPortal.Get(i);
    auto ptRPZ = cPortal.Get(i);
    auto R = ptRPZ[0];
    //auto Z = ptRPZ[2];

    auto inv_r = 1.0/R;
    auto Bmag = vtkm::Magnitude(B0);
    auto over_Bmag = 1.0/Bmag;
    auto over_Bmag2 = over_Bmag * over_Bmag;

    auto br = B0[0];
    auto bz = B0[1];
    auto bphi = B0[2];

    auto GRAD = gPortal.Get(i);

    auto dbrdr = GRAD[0][0];
    auto dbzdr = GRAD[0][1];
    auto dbpdr = GRAD[0][2];

    auto dbrdz = GRAD[1][0];
    auto dbzdz = GRAD[1][1];
    auto dbpdz = GRAD[1][2];

    auto dbrdp = GRAD[2][0];
    auto dbzdp = GRAD[2][1];
    //auto dbpdp = GRAD[2][2];

    //dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
    //dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
    //dbdphi=0D0  ! no B perturbation
    auto dbdr = (br*dbrdr + bphi*dbpdr + bz*dbzdr) * over_Bmag;
    auto dbdz = (br*dbrdz + bphi*dbpdz + bz*dbzdz) * over_Bmag;
    //auto dbdphi = 0;


    vtkm::Vec3f curl_B;
    //R curl_B(1)  = fld%dbzdp*inv_r - fld%dbpdz
    curl_B[0] =          dbzdp*inv_r - dbpdz;
    //Z curl_B(2)  = fld%bphi*inv_r + fld%dbpdr - fld%dbrdp*inv_r
    curl_B[1] =          bphi*inv_r +     dbpdr -     dbrdp*inv_r;
    //phi curl_B(3)  = fld%dbrdz - fld%dbzdr
    curl_B[2] =            dbrdz -     dbzdr;

    vtkm::Vec3f curl_nb;

    //R,Z,Phi
    //curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
    //curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
    //curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2

    curl_nb[0] = curl_B[0]*over_Bmag + (bphi * dbdz) * over_Bmag2;
    curl_nb[1] = curl_B[1]*over_Bmag + (-bphi * dbdr) * over_Bmag2;
    curl_nb[2] = curl_B[2]*over_Bmag + (bz * dbdr - br * dbdz) * over_Bmag2;

    portal.Set(i, curl_nb);
  }

  return curlBNorm;
#endif
}

void
ReadB(adiosS* stuff,
      vtkm::cont::DataSet& ds)
{
  std::string fileName = "/node_data[0]/values";

  auto Bvar = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(Bvar, tmp, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> b_rzp;
  for (int i = 0; i < numNodes; i++)
  {
    vtkm::Vec3f v(tmp[i*3+0], tmp[i*3+1], tmp[i*3+2]);
    b_rzp.push_back(v);
  }

  ds.AddField(vtkm::cont::make_FieldPoint("B_RZP",vtkm::cont::make_ArrayHandle(b_rzp, vtkm::CopyFlag::On)));

  vtkm::cont::ArrayHandle<vtkm::Vec3f> bhat_rzp, b;
  vtkm::cont::Invoker invoker;

  ds.GetField("B_RZP").GetData().AsArrayHandle(b);
  invoker(NormalizeWorklet{}, b, bhat_rzp);
  ds.AddField(vtkm::cont::make_FieldPoint("B_RZP_Norm",bhat_rzp));

  /*
  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm_rzp;
  curlBNorm_rzp = ComputeCurl(ds);
  ds.AddField(vtkm::cont::make_FieldPoint("curl_nb_rzp", curlBNorm_rzp));
  */
}

void
ReadVec(adiosS* stuff,
        vtkm::cont::DataSet& ds,
        const std::string& vname,
        std::string fileName="",
        bool add3D=false,
        bool addExtra=false)
{
  if (fileName.empty())
    fileName = vname;

  bool isB = (vname == "B");

  auto var = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(var, tmp, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> vec, vecRZP, vec2d;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < numNodes; i++)
    {
      //B is saved as: R,Z,Phi
      int vidx = (isB ? i : (p*numNodes+i));
      //Swap to R,Phi,Z
      vtkm::Vec3f v(tmp[vidx*3+0], tmp[vidx*3+2], tmp[vidx*3+1]);
      //vtkm::Vec3f v(tmp[vidx*3+0], tmp[vidx*3+1], tmp[vidx*3+2]);
      vec.push_back(v);
      if (p == 0)
      {
        vec2d.push_back(v);
        vecRZP.push_back(vtkm::Vec3f(tmp[vidx*3+0], tmp[vidx*3+1], tmp[vidx*3+2]));
      }
    }
  }

  if (addExtra && add3D)
  {
    for (int i = 0; i < numNodes; i++)
      vec.push_back(vec[i]);
  }

  if (add3D)
    ds.AddField(vtkm::cont::make_FieldPoint(vname,vtkm::cont::make_ArrayHandle(vec, vtkm::CopyFlag::On)));
  else
  {
    ds.AddField(vtkm::cont::make_FieldPoint(vname+"2D",vtkm::cont::make_ArrayHandle(vec2d, vtkm::CopyFlag::On)));
    ds.AddField(vtkm::cont::make_FieldPoint("B_RZP_2D",vtkm::cont::make_ArrayHandle(vecRZP, vtkm::CopyFlag::On)));

    if (isB)
    {
      vtkm::cont::ArrayHandle<vtkm::Vec3f> b, brzp, bhat, brzphat;
      vtkm::cont::Invoker invoker;
      ds.GetField("B2D").GetData().AsArrayHandle(b);
      invoker(NormalizeWorklet{}, b, bhat);
      ds.AddField(vtkm::cont::make_FieldPoint("B2D_Norm",bhat));

      ds.GetField("B_RZP_2D").GetData().AsArrayHandle(brzp);
      invoker(NormalizeWorklet{}, brzp, brzphat);
      ds.AddField(vtkm::cont::make_FieldPoint("BRZP_Norm", brzphat));


      vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;
      curlBNorm = ComputeCurl(ds);


      ds.AddField(vtkm::cont::make_FieldPoint("Curl_B2D_Norm", curlBNorm));
      ds.AddField(vtkm::cont::make_FieldPoint("curl_nb_rzp", curlBNorm));

      //Add B3D.
      std::vector<vtkm::Vec3f> B3D;
      auto portal = b.ReadPortal();
      for (vtkm::Id p = 0; p < numPlanes; p++)
        for (vtkm::Id n = 0; n < numNodes; n++)
          B3D.push_back(portal.Get(n));
      ds.AddField(vtkm::cont::make_Field("B3D",
                                         vtkm::cont::Field::Association::WHOLE_MESH,
                                         B3D,
                                         vtkm::CopyFlag::On));
    }
  }
}

vtkm::cont::DataSet
ReadMesh(adiosS* meshStuff)
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
  for (int i = 0; i < numNodes; i++)
  {
    double R = rz[i*2 +0];
    double Z = rz[i*2 +1];

    vtkm::Vec3f ptRZ(R,Z,0);
    coords.push_back(ptRZ);
  }

  //cells
  for (int i = 0; i < numTri*3; i+=3)
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
ReadMesh3D(adiosS* meshStuff, bool addExtra)
{
  std::vector<double> rz;
  std::vector<int> conn, nextnode;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
  std::string connNm = "/cell_set[0]/node_connect_list";
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("nextnode"), nextnode, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> coords;
  std::vector<vtkm::Id> connIds;

  double dPhi = vtkm::TwoPi()/static_cast<double>(numPlanes);

  //points.
  double phi = 0;
  int NP = numPlanes;
  if (addExtra) NP++;
  for (int p = 0; p < NP; p++)
  {
    std::cout<<"ReadMesh3D: phi= "<<phi<<std::endl;
    for (int i = 0; i < numNodes; i++)
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
    for (int i = 0; i < numTri*3; i+=3)
    {
      int off = p*numNodes;
      int p0 = conn[i+0];
      int p1 = conn[i+1];
      int p2 = conn[i+2];
      connIds.push_back(p0+off);
      connIds.push_back(p1+off);
      connIds.push_back(p2+off);

      off = (p+1)*(numNodes);
      connIds.push_back(p0 + off);
      connIds.push_back(p1 + off);
      connIds.push_back(p2 + off);
    }
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  return dsb.Create(coords, vtkm::CellShapeTagWedge(), 6, connIds, "coords");
}

//-----------------------------------------------------------------------------
class FindCellWorklet : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn points,
                                WholeCellSetIn<> cellSet,
                                ExecObject locator,
                                FieldOut cellIds,
                                FieldOut pcoords);
//                                FieldOut ptIndices);
  using ExecutionSignature = void(_1, _2, _3, _4, _5); //, _6);
  using InputDomain = _1;

  template <typename LocatorType, typename CellSetType> //, typename PtIndexType>
  VTKM_EXEC void operator()(const vtkm::Vec3f& point,
                            const CellSetType& vtkmNotUsed(cellSet),
                            const LocatorType& locator,
                            vtkm::Id& cellId,
                            vtkm::Vec3f& pcoords) const
  //                            PtIndexType& ptIndices ) const
  {
    vtkm::ErrorCode status = locator.FindCell(point, cellId, pcoords);
    if (status != vtkm::ErrorCode::Success)
    {
      //std::cout<<"Cell not found! "<<point<<std::endl;
//      std::cout<<"    ***** Try reducing the step size."<<std::endl;
      this->RaiseError(vtkm::ErrorString(status));
    }
    //ptIndices = cellSet.GetIndices(cellId);
    //auto x = cellSet.GetIndices(cellId);
    //ptIndices = x;
  }
};


std::vector<std::vector<vtkm::Vec3f>>
ConvertPuncturesToThetaPsi(const std::vector<std::vector<vtkm::Vec3f>>& puncturesRZ,
                           const vtkm::cont::DataSet& ds)
{
  std::cout<<"Convert Punctures: "<<puncturesRZ.size()<<std::endl;

  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  std::vector<std::vector<vtkm::Vec3f>> puncturesTP(puncturesRZ.size());


  std::vector<vtkm::Vec3f> tmp;
  auto allPuncRZ = vtkm::cont::make_ArrayHandle(tmp, vtkm::CopyFlag::On);


#if 1
  for (int p = 0; p < (int)puncturesRZ.size(); p++)
  {
    auto pts = puncturesRZ[p];
    for (auto& pt : pts)
    {
      auto R = pt[0];
      auto Z = pt[2];
      pt = vtkm::Vec3f(R,Z,0);
    }
    //std::cout<<"RZ_pts.size()= "<<pts.size()<<std::endl;
    auto RZpoints = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

    vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
    vtkm::cont::Invoker invoker;

    invoker(FindCellWorklet{}, RZpoints, ds.GetCellSet(), locator, cellIds, pcoords);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> psi, theta;
    ds.GetField("psi2D").GetData().AsArrayHandle(psi);
    ds.GetField("theta2D").GetData().AsArrayHandle(theta);

    auto cPortal = cellIds.ReadPortal();
    auto pPortal = pcoords.ReadPortal();
    auto psiPortal = psi.ReadPortal();
    auto thetaPortal = theta.ReadPortal();
    auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

    auto pRZpoints = RZpoints.ReadPortal();
    vtkm::Id numPts = RZpoints.GetNumberOfValues();
    puncturesTP[p].resize(numPts);
    for (vtkm::Id i = 0; i < numPts; i++)
    {
      auto R = pRZpoints.Get(i)[0];
      auto Z = pRZpoints.Get(i)[1];
      vtkm::Id vIds[3];
      vtkm::Id cid = cPortal.Get(i);
      cs.GetCellPointIds(cid, vIds);

      vtkm::VecVariable<vtkm::FloatDefault, 3> pVals, tVals;
      for (vtkm::Id j = 0; j < 3; j++)
      {
        pVals.Append(psiPortal.Get(vIds[j]));
        tVals.Append(thetaPortal.Get(vIds[j]));
      }

      vtkm::FloatDefault psiI, thetaI;
      vtkm::exec::CellInterpolate(pVals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), psiI);
      vtkm::exec::CellInterpolate(tVals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), thetaI);
      thetaI += vtkm::Pi();
      auto thetaVal = vtkm::ATan2(Z-eq_axis_z, R-eq_axis_r);
      if (thetaVal < 0)
        thetaVal += vtkm::TwoPi();

      //puncturesTP[p][i][0] = thetaI;
      puncturesTP[p][i][0] = thetaVal;
      puncturesTP[p][i][1] = psiI / eq_x_psi;
    }
  }
#endif

  return puncturesTP;
}

void
Poincare(const vtkm::cont::DataSet& ds,
         std::vector<vtkm::Vec3f>& pts,
         const std::map<std::string, std::vector<std::string>>& args,
         const std::string& /*vField*/,
         vtkm::FloatDefault h,
         int numPunc,
         int whichWorklet,
         bool useBOnly,
         bool useHighOrderB,
         vtkm::cont::ArrayHandle<vtkm::Vec2f>& outRZ,
         vtkm::cont::ArrayHandle<vtkm::Vec2f>& outTP,
         vtkm::cont::ArrayHandle<vtkm::Id>& outID,
         bool quickTest,
         vtkm::Id BGridSize,
         bool BGridCell,
         std::string& locName,
         std::string& locParam1,
         std::string& locParam2,
         std::vector<std::vector<vtkm::Vec3f>>* traces=nullptr)

{
  //vtkm::cont::CellLocatorGeneral locator;
  vtkm::cont::CellLocatorBoundingIntervalHierarchy locatorBIH;
  vtkm::cont::CellLocatorTwoLevel locator2L;

  bool locBIH = false;
  if (locName == "BIH")
  {
    locBIH = true;
    if (!locParam1.empty() && !locParam2.empty())
    {
      //defaults 4, 5
      locatorBIH.SetNumberOfSplittingPlanes(std::atoi(locParam1.c_str()));
      locatorBIH.SetMaxLeafSize(std::atoi(locParam2.c_str()));
    }
    std::cout<<"LocatorBIH: "<<locParam1<<" "<<locParam2<<std::endl;
    locatorBIH.SetCellSet(ds.GetCellSet());
    locatorBIH.SetCoordinates(ds.GetCoordinateSystem());
    auto startL = std::chrono::steady_clock::now();
    locatorBIH.Update();
    std::chrono::duration<double> dt = std::chrono::steady_clock::now()-startL;
    std::cout<<"BIH build= "<<dt.count()<<std::endl;
  }
  else
  {
    if (!locParam1.empty() && !locParam2.empty())
    {
      //defaults 32.0, 2.0
      locator2L.SetDensityL1(std::atof(locParam1.c_str()));
      locator2L.SetDensityL2(std::atof(locParam2.c_str()));
    }
    std::cout<<"Locator2Level: "<<locParam1<<" "<<locParam2<<std::endl;
    locator2L.SetCellSet(ds.GetCellSet());
    locator2L.SetCoordinates(ds.GetCoordinateSystem());
    auto startL = std::chrono::steady_clock::now();
    locator2L.Update();
    std::chrono::duration<double> dt = std::chrono::steady_clock::now()-startL;
    std::cout<<"2L build= "<<dt.count()<<std::endl;

  }

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> eq_psi_gridArr;
  ds.GetField("eq_psi_grid").GetData().AsArrayHandle(eq_psi_gridArr);
  auto dPsi = eq_psi_gridArr.ReadPortal().Get(1)-eq_psi_gridArr.ReadPortal().Get(0);

  vtkm::cont::Invoker invoker;
  std::vector<vtkm::Particle> s;
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
    s.push_back(vtkm::Particle(pts[i], i));
  auto seeds = vtkm::cont::make_ArrayHandle(s, vtkm::CopyFlag::On);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As_ff, coeff_1D, coeff_2D, psi;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> B_rzp, B_Norm_rzp, dAs_ff_rzp;//, AsCurlBHat_rzp, curl_nb_rzp;
  ds.GetField("As_ff").GetData().AsArrayHandle(As_ff);
  ds.GetField("B_RZP").GetData().AsArrayHandle(B_rzp);
  ds.GetField("B_RZP_Norm").GetData().AsArrayHandle(B_Norm_rzp);
  ds.GetField("dAs_ff_rzp").GetData().AsArrayHandle(dAs_ff_rzp);
  //ds.GetField("AsCurlBHat_RZP").GetData().AsArrayHandle(AsCurlBHat_rzp);
  //ds.GetField("curl_nb_rzp").GetData().AsArrayHandle(curl_nb_rzp);
  ds.GetField("coeff_1D").GetData().AsArrayHandle(coeff_1D);
  ds.GetField("coeff_2D").GetData().AsArrayHandle(coeff_2D);
  ds.GetField("psi2D").GetData().AsArrayHandle(psi);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> tracesArr;
  outRZ.Allocate(numPunc*pts.size());
  outTP.Allocate(numPunc*pts.size());
  outID.Allocate(numPunc*pts.size());

  vtkm::cont::Timer timer;
  auto start = std::chrono::steady_clock::now();
  std::cout<<"Using Worklet_"<<whichWorklet<<std::endl;
  if (whichWorklet == 1)
  {
#ifdef BUILD_POINC1
    PoincareWorklet worklet(numPunc, 0.0f, h, (traces!=nullptr), quickTest);
    worklet.UseBOnly = useBOnly;
    worklet.UseHighOrderB = useHighOrderB;

    if (traces != nullptr)
    {
      std::vector<vtkm::Vec3f> t;
      t.resize(pts.size()*worklet.MaxIter, {-100, -100, -100});
      std::cout<<"TRACES: "<<t.size()<<std::endl;
      tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
    }
#ifdef VTKM_CUDA
    // This worklet needs some extra space on CUDA.
    vtkm::cont::cuda::internal::ScopedCudaStackSize stack(16 * 1024);
    (void)stack;
#endif // VTKM_CUDA


    //timer.Start();
    start = std::chrono::steady_clock::now();
    invoker(worklet, seeds,
            locator2L,
            ds.GetCellSet(),
            B_rzp, B_Norm_rzp, /*curl_nb_rzp,*/ As_ff, dAs_ff_rzp,
            coeff_1D, coeff_2D,
            tracesArr, outRZ, outTP, outID);
#else
    std::cout<<"*************************************** Code NOT build with Worklet 1"<<std::endl;
#endif
  }
  else if (whichWorklet == 2)
  {
#ifdef BUILD_POINC2
    PoincareWorklet2 worklet(numPunc, 0.0f, h, (traces!=nullptr), quickTest);
    worklet.UseBOnly = useBOnly;
    worklet.UseHighOrderB = useHighOrderB;
    worklet.one_d_cub_dpsi_inv = 1.0/dPsi;

    if (traces != nullptr)
    {
      std::vector<vtkm::Vec3f> t;
      t.resize(pts.size()*worklet.MaxIter, {-100, -100, -100});
      std::cout<<"Allocate TRACES: "<<t.size()<<std::endl;
      tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
    }

    //timer.Start();
    start = std::chrono::steady_clock::now();
    if (locBIH)
      invoker(worklet, seeds,
              locatorBIH,
              ds.GetCellSet(),
              ds.GetCoordinateSystem(),
              As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
              tracesArr, outRZ, outTP, outID);
    else
      invoker(worklet, seeds,
              locator2L,
              ds.GetCellSet(),
              ds.GetCoordinateSystem(),
              As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
              tracesArr, outRZ, outTP, outID);
#else
    std::cout<<"*************************************** Code NOT build with Worklet 2"<<std::endl;
#endif
  }
  else if (whichWorklet == 21)
  {
#ifdef BUILD_POINC2

    //Compute the As values on a uniform grid.
    vtkm::Id3 dims(BGridSize, BGridSize, numPlanes*2);
    vtkm::FloatDefault dR = (eq_max_r-eq_min_r) / static_cast<vtkm::FloatDefault>(BGridSize-1);
    vtkm::FloatDefault dZ = (eq_max_z-eq_min_z) / static_cast<vtkm::FloatDefault>(BGridSize-1);
    vtkm::Vec3f origin(eq_min_r, eq_min_z, 0), maxPt(eq_max_r, eq_max_z, numPlanes*2);
    vtkm::Vec3f spacing(dR, dZ, 1);

    auto grid = vtkm::cont::DataSetBuilderUniform::Create(dims, origin, spacing);
    vtkm::cont::ArrayHandleUniformPointCoordinates uniform2DCoords({BGridSize, BGridSize, 1},
                                                                   origin, {dR, dZ, 0});
    vtkm::cont::CellSetStructured<2> uniform2DCells;
    uniform2DCells.SetPointDimensions(vtkm::Id2(BGridSize, BGridSize));

    std::cout<<"********** NumCells= "<<grid.GetCellSet().GetNumberOfCells()<<std::endl;
    grid.PrintSummary(std::cout);

    /*
    grid.AddField(vtkm::cont::make_FieldPoint("As", AsUniform));
    grid.AddField(vtkm::cont::make_FieldPoint("dAs", dAsUniform));
    vtkm::io::VTKDataSetWriter writer("asGrid.vtk");
    writer.WriteDataSet(grid);
    */

    PoincareWorklet2_1 worklet(numPunc, 0.0f, h, (traces!=nullptr), quickTest);
    worklet.SetAsGrid(origin, spacing, maxPt, dims);
    worklet.UseBOnly = useBOnly;
    worklet.UseHighOrderB = useHighOrderB;
    worklet.UseAsCell = BGridCell;

    /*
    auto cellCountP = locator2L.CellCount.ReadPortal();
    vtkm::Id maxCount = cellCountP.Get(0);
    for (vtkm::Id i = 0; i < cellCountP.GetNumberOfValues(); i++)
    {
      vtkm::Id val = cellCountP.Get(i);
      if (val > maxCount) maxCount = val;
    }
    worklet.MaxCellCount = maxCount;
    */

    if (traces != nullptr)
    {
      std::vector<vtkm::Vec3f> t;
      t.resize(pts.size()*worklet.MaxIter, {-100, -100, -100});
      std::cout<<"TRACES: "<<t.size()<<std::endl;
      tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
    }

    vtkm::cont::ArrayHandle<vtkm::Vec3f> dAsUniform;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> AsUniform;

    vtkm::cont::Timer timer2;
    timer2.Start();
    if (worklet.UseAsCell)
    {
      std::cout<<"********************* AS Cell"<<std::endl;
      vtkm::Id nCells = uniform2DCells.GetNumberOfCells() * numPlanes * 2;
      AsUniform.Allocate(nCells);
      dAsUniform.Allocate(nCells);
      ComputeAsCellWorklet computeAs;
      computeAs.Num2DCells = uniform2DCells.GetNumberOfCells();

      invoker(computeAs,
              uniform2DCells, uniform2DCoords,
              locator2L,
              ds.GetCellSet(),
              As_ff, dAs_ff_rzp,
              AsUniform, dAsUniform);
    }
    else
    {
      std::cout<<"********************* AS Point"<<std::endl;
      ComputeAsWorklet computeAs;
      vtkm::Id nPts = uniform2DCoords.GetNumberOfValues() * numPlanes * 2;
      AsUniform.Allocate(nPts);
      dAsUniform.Allocate(nPts);
      computeAs.Num2DPts = uniform2DCoords.GetNumberOfValues();

      invoker(computeAs,
              uniform2DCoords,
              locator2L,
              ds.GetCellSet(),
              As_ff, dAs_ff_rzp,
              AsUniform, dAsUniform);
    }
    timer2.Stop();
    std::cout<<"... BGrid compute done."<<std::endl;
    std::cout<<"Grid build timer= "<<timer2.GetElapsedTime()<<std::endl;

    start = std::chrono::steady_clock::now();
    invoker(worklet, seeds,
            locator2L,
            ds.GetCellSet(),
            ds.GetCoordinateSystem(),
            As_ff, dAs_ff_rzp, coeff_1D, coeff_2D,
            AsUniform, dAsUniform,
            tracesArr, outRZ, outTP, outID);

#else
    std::cout<<"*************************************** Code NOT build with Worklet 2"<<std::endl;
#endif
  }
  else if (whichWorklet == 3)
  {
#ifdef BUILD_POINC3
    /* N = 128
       B0: [-0.0203529,-0.00366229,-0.239883] [-0.0203388,-0.0036626,-0.239886] 1.42779e-05
       curlB0: [0,-1.38778e-17,-0.0750265] [0,-9.04278e-18,-0.0749931] 3.34226e-05
       curlNB0: [-0.0188931,-0.330013,-0.304966] [-0.0188503,-0.33007,-0.304813] 0.000168853
       gradPsi: [-0.0101055,0.0561605,0] [-0.010099,0.0561215,0] 3.96151e-05
    */

    //Precompute B
    vtkm::Id2 dims(BGridSize, BGridSize);
    vtkm::FloatDefault dR = (eq_max_r-eq_min_r) / static_cast<vtkm::FloatDefault>(BGridSize-1);
    vtkm::FloatDefault dZ = (eq_max_z-eq_min_z) / static_cast<vtkm::FloatDefault>(BGridSize-1);
    vtkm::Vec2f origin(eq_min_r, eq_min_z), maxPt(eq_max_r, eq_max_z);
    vtkm::Vec2f spacing(dR, dZ);
    auto grid = vtkm::cont::DataSetBuilderUniform::Create(dims, origin, spacing);
    std::cout<<"origin= "<<origin<<std::endl;
    std::cout<<"spacing= "<<spacing<<std::endl;
    std::cout<<"R: "<<eq_min_r<<" "<<eq_max_r<<" dR= "<<dR<<std::endl;
    std::cout<<"Z: "<<eq_min_z<<" "<<eq_max_z<<" dZ= "<<dZ<<std::endl;

    vtkm::cont::CellLocatorUniformGrid locatorB;
    locatorB.SetCellSet(grid.GetCellSet());
    locatorB.SetCoordinates(grid.GetCoordinateSystem());
    locatorB.Update();

    std::cout<<"Precomputing BGrid with resolution of "<<dims<<std::endl;
    ComputeBWorklet computeB;
    ComputeBCellWorklet computeBCell;

    vtkm::cont::Timer timer2;
    vtkm::cont::ArrayHandle<vtkm::Vec3f> B0, CurlB0, CurlNB0, GradPsi;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> Psi;

    if (!BGridCell)
    {
      timer2.Start();
      invoker(computeB,
              grid.GetCoordinateSystem(),
              coeff_1D,
              coeff_2D,
              Psi,
              B0,
              CurlB0,
              CurlNB0,
              GradPsi);
      timer2.Stop();
    }
    else
    {
      timer2.Start();
      invoker(computeBCell,
              grid.GetCellSet(),
              grid.GetCoordinateSystem(),
              coeff_1D,
              coeff_2D,
              Psi,
              B0,
              CurlB0,
              CurlNB0,
              GradPsi);
      timer2.Stop();
    }

    std::cout<<"... BGrid compute done."<<std::endl;
    std::cout<<"Grid build timer= "<<timer2.GetElapsedTime()<<std::endl;

    //vtkm::Id divGood = vtkm::cont::Algorithm::Reduce(DivergenceGood, vtkm::Id(0)); //vtkm::Sum());
    //std::cout<<"Div Good= "<<divGood<<std::endl;

    std::string fname;
    if (BGridCell)
    {
      grid.AddField(vtkm::cont::make_FieldCell("Psi", Psi));
      grid.AddField(vtkm::cont::make_FieldCell("B0", B0));
      grid.AddField(vtkm::cont::make_FieldCell("CurlB0", CurlB0));
      grid.AddField(vtkm::cont::make_FieldCell("CurlNB0", CurlNB0));
      grid.AddField(vtkm::cont::make_FieldCell("GradPsi", GradPsi));
      fname = "bgridCell.vtk";
    }
    else
    {
      grid.AddField(vtkm::cont::make_FieldPoint("Psi", Psi));
      grid.AddField(vtkm::cont::make_FieldPoint("B0", B0));
      grid.AddField(vtkm::cont::make_FieldPoint("CurlB0", CurlB0));
      grid.AddField(vtkm::cont::make_FieldPoint("CurlNB0", CurlNB0));
      grid.AddField(vtkm::cont::make_FieldPoint("GradPsi", GradPsi));
      fname = "bgridPt.vtk";
    }

    /*
    vtkm::io::VTKDataSetWriter writer(fname.c_str());
    writer.WriteDataSet(grid);
    grid.PrintSummary(std::cout);
    */

    PoincareWorklet3 worklet(numPunc, 0.0f, h, (traces!=nullptr), quickTest);
    worklet.SetBGrid(origin, spacing, maxPt, dims, BGridCell);
    worklet.UseBOnly = useBOnly;
    worklet.UseHighOrderB = useHighOrderB;

    if (traces != nullptr)
    {
      std::vector<vtkm::Vec3f> t;
      t.resize(pts.size()*worklet.MaxIter, {-100, -100, -100});
      std::cout<<"TRACES: "<<t.size()<<std::endl;
      tracesArr = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
    }

    //timer.Start();
    start = std::chrono::steady_clock::now();
    if (locBIH)
      invoker(worklet, seeds,
              locatorBIH, locatorB,
              ds.GetCellSet(), grid.GetCellSet(),
              B_rzp, /*B_Norm_rzp,*/ /*curl_nb_rzp,*/ As_ff, dAs_ff_rzp,
              Psi, B0, CurlB0, CurlNB0,  GradPsi,
              coeff_1D, coeff_2D,
              tracesArr, outRZ, outTP, outID);
    else
      invoker(worklet, seeds,
              locator2L, locatorB,
              ds.GetCellSet(), grid.GetCellSet(),
              B_rzp, /*B_Norm_rzp,*/ /*curl_nb_rzp,*/ As_ff, dAs_ff_rzp,
              Psi, B0, CurlB0, CurlNB0,  GradPsi,
              coeff_1D, coeff_2D,
              tracesArr, outRZ, outTP, outID);


    /* old worklet3
    invoker(worklet, seeds, locator, locatorB,
            ds.GetCellSet(), grid.GetCellSet(),
            Psi, B0, CurlB0, CurlNB0, GradPsi,
            As_ff, dAs_ff_rzp,
            coeff_1D, coeff_2D,
            tracesArr, outRZ, outTP, outID);
    */
#else
    std::cout<<"*************************************** Code NOT build with Worklet 3"<<std::endl;
#endif
  }
  outID.SyncControlArray();
  auto end = std::chrono::steady_clock::now();
  //timer.Stop();
  std::chrono::duration<double> elapsed_seconds = end-start;

  std::cout<<"PoincareTime= "<<elapsed_seconds.count()<<std::endl;
  //std::cout<<"vtkm::cont::Timer= "<<timer.GetElapsedTime()<<std::endl;
  std::cout<<"outputRZ.size()= "<<outRZ.GetNumberOfValues()<<std::endl;
  std::cout<<"outputTP.size()= "<<outTP.GetNumberOfValues()<<std::endl;
  std::cout<<"punctureID.size()= "<<outID.GetNumberOfValues()<<std::endl;
  //vtkm::cont::printSummary_ArrayHandle(output, std::cout);
  //vtkm::cont::printSummary_ArrayHandle(punctureID, std::cout);

  std::vector<std::vector<vtkm::Vec3f>> res;
  if (traces)
  {
    auto portal = tracesArr.ReadPortal();
    vtkm::Id n = portal.GetNumberOfValues();
    for (vtkm::Id i = 0; i < n; i++)
    {
      auto v = portal.Get(i);
      if (v[2] > -1)
        (*traces)[0].push_back(v);
    }
  }
}

void
CalcAs(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As;

  ds.GetField("As_phi_ff").GetData().AsArrayHandle(As);
  vtkm::Id numAs = As.GetNumberOfValues();
  auto asPortal = As.ReadPortal();

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::FloatDefault>>> As_arr;
  As_arr.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    As_arr[p].resize(2);
    for (int i = 0; i < 2; i++)
      As_arr[p][i].resize(numNodes);
  }

  //Do some easy indexing.
  vtkm::Id idx = 0;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int n = 0; n < numNodes; n++)
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
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::FloatDefault asVal = As_arr[p][i][n];
        arrAs[idx] = asVal;
        idx++;
      }
    }
  }

  ds.AddField(vtkm::cont::make_Field("As_ff",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arrAs,
                                     vtkm::CopyFlag::On));
}

void
Calc_dAs(vtkm::cont::DataSet& ds)
{
  //Assumed that dAs_phi_ff is R,Z,Phi
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> dAsAH;
  ds.GetField("dAs_phi_ff").GetData().AsArrayHandle(dAsAH);
  int nVals = dAsAH.GetNumberOfValues()/3;

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> dAs_ff;
  dAs_ff.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    dAs_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      dAs_ff[p][i].resize(numNodes);
  }


  int idx = 0;
  auto dAsPortal = dAsAH.ReadPortal();

  for (int p = 0; p < numPlanes; p++)
  {
    for (int n = 0; n < numNodes; n++)
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
  VEC_ff.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    VEC_ff[p].resize(2);
    for (int i = 0; i < 2; i++)
      VEC_ff[p][i].resize(numNodes);
  }

  idx = 0;
  std::vector<vtkm::Vec3f> arr_ff(nVals);
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::Vec3f val = VEC_ff[p][i][n];
        val = dAs_ff[p][i][n];
        arr_ff[idx] = val;
        idx++;
      }
    }
  }

  //do some testing..
  for (int p = 0; p < numPlanes; p++)
  {
    int off0 = p*numNodes*2;
    int off1 = p*numNodes*2 + numNodes;
    for (int n = 0; n < numNodes; n++)
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

  ds.AddField(vtkm::cont::make_Field("dAs_ff_rzp",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr_ff,
                                     vtkm::CopyFlag::On));
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

  auto RZBuff = vtkm::cont::ArrayHandleBasic<vtkm::Vec2f>(outRZ).GetReadPointer();
  auto TPBuff = vtkm::cont::ArrayHandleBasic<vtkm::Vec2f>(outTP).GetReadPointer();
  auto IDBuff = vtkm::cont::ArrayHandleBasic<vtkm::Id>(outID).GetReadPointer();

  std::vector<std::size_t> shape = {nPts*2}, offset = {0}, size = {nPts*2};
  auto vRZ = io.DefineVariable<vtkm::FloatDefault>("RZ", shape, offset, size);
  auto vTP = io.DefineVariable<vtkm::FloatDefault>("ThetaPsi", shape, offset, size);
  std::vector<std::size_t> shape2 = {nPts},size2 = {nPts};
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

vtkm::cont::DataSet
ReadData_ORIG(std::map<std::string, std::vector<std::string>>& args)
{
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

  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &numPlanes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &numTri, adios2::Mode::Sync);
  std::vector<double> psiVals;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("psi"), psiVals, adios2::Mode::Sync);
  psi_min = psi_max = psiVals[0];
  for (const auto& p : psiVals)
  {
    if (p < psi_min) psi_min = p;
    if (p > psi_max) psi_max = p;
  }
  std::cout<<"********                  PSI m/M = "<<psi_min<<" "<<psi_max<<std::endl;

  //Try and do everything in cylindrical coords and worklets.
  auto ds = ReadMesh(meshStuff);
  ReadScalar(dataStuff, ds, "dpot");
  ReadScalar(dataStuff, ds, "apars", "apars", true);
  ReadScalar(meshStuff, ds, "psi");
  ReadScalar(meshStuff, ds, "theta");
  ReadB(bfieldStuff, ds);

  /*
  ReadVec(bfieldStuff, ds, "B", "/node_data[0]/values");
  ReadVec(bfieldStuff, ds, "B3D", "/node_data[0]/values", true, false);
  ReadVec(bfield_allStuff, ds, "Bs", "Bs", true);
  */
  ReadOther(bfieldStuff, ds, "As_phi_ff");
  ReadOther(bfieldStuff, ds, "dAs_phi_ff");
  ReadPsiInterp(equilStuff, interp_coeffStuff, ds);

  CalcAs(ds);
  Calc_dAs(ds);

  if (0)
  {
    auto ds3d = ReadMesh3D(meshStuff, true);
    ReadScalar(dataStuff, ds3d, "apars", "apars", true, true);
    ReadVec(bfieldStuff, ds3d, "B", "/node_data[0]/values", true, true);
    ReadVec(bfield_allStuff, ds3d, "Bs", "Bs", true, true);
    //CalcX(ds3d);
    vtkm::io::VTKDataSetWriter writer("debug.vtk");
    writer.WriteDataSet(ds3d);
  }

  //ds.PrintSummary(std::cout);
  vtkm::io::VTKDataSetWriter writer("grid.vtk");
  writer.WriteDataSet(ds);

  return ds;
}

vtkm::cont::DataSet
ReadData(std::map<std::string, std::vector<std::string>>& args)
{
  if (args.find("--Step") == args.end())
    return ReadData_ORIG(args);

  int step = std::stoi(args["--Step"][0]);

  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  adios = new adios2::ADIOS;
  adiosStuff["init"] = new adiosS(adios, "xgc.poincare_init.bp", "init", adiosArgs);
  adiosStuff["mesh"] = new adiosS(adios, "xgc.mesh.bp", "mesh", adiosArgs);
  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", adiosArgs);
  adiosStuff["equil"] = new adiosS(adios, "xgc.equil.bp", "equil", adiosArgs);
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "bfield", adiosArgs);

  /*
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "/node_data[0]/values", adiosArgs);
  adiosStuff["bfield-all"] = new adiosS(adios, "xgc.bfield-all.bp", "Bs", adiosArgs);
  //adiosStuff["psi_bicub_acoef"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "psi_bicub_acoef", adiosArgs);
  //adiosStuff["one_d_cub_psi_acoef"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "one_d_cub_acoef", adiosArgs);
  adiosStuff["interp_coeff"] = new adiosS(adios, "xgc.bfield_with_coeff.bp", "interp_coeff", adiosArgs);
  adiosStuff["equil"] = new adiosS(adios, "xgc.equil.bp", "equil", adiosArgs);
  */

  auto initStuff = adiosStuff["init"];
  auto meshStuff = adiosStuff["mesh"];
  auto dataStuff = adiosStuff["data"];
  auto equilStuff = adiosStuff["equil"];
  auto bfieldStuff = adiosStuff["bfield"];

  /*
  auto bfield_allStuff = adiosStuff["bfield-all"];
  auto interp_coeffStuff = adiosStuff["interp_coeff"];
  //auto one_dcub_acoefStuff = adiosStuff["one_d_cub_psi_acoef"];
  */

  initStuff->engine.BeginStep();
  meshStuff->engine.BeginStep();
  equilStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();

  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &numPlanes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &numTri, adios2::Mode::Sync);
  std::vector<double> psiVals;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("psi"), psiVals, adios2::Mode::Sync);
  psi_min = psi_max = psiVals[0];
  for (const auto& p : psiVals)
  {
    if (p < psi_min) psi_min = p;
    if (p > psi_max) psi_max = p;
  }
  std::cout<<"********                  PSI m/M = "<<psi_min<<" "<<psi_max<<std::endl;

  auto ds = ReadMesh(meshStuff);
  ReadPsiInterp(equilStuff, initStuff, ds);
  ReadB(bfieldStuff, ds);

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
  ReadScalar(dataStuff, ds, "dpot");
  ReadScalar(meshStuff, ds, "psi");


/*
  ReadScalar(dataStuff, ds, "apars", "apars", true);
  ReadScalar(meshStuff, ds, "theta");
*/

  /*
  ReadVec(bfieldStuff, ds, "B", "/node_data[0]/values");
  ReadVec(bfieldStuff, ds, "B3D", "/node_data[0]/values", true, false);
  ReadVec(bfield_allStuff, ds, "Bs", "Bs", true);
  ReadOther(bfieldStuff, ds, "As_phi_ff");
  ReadOther(bfieldStuff, ds, "dAs_phi_ff");
  ReadPsiInterp(equilStuff, interp_coeffStuff, ds);
  */

  CalcAs(ds);
  Calc_dAs(ds);

  /*
  if (0)
  {
    auto ds3d = ReadMesh3D(meshStuff, true);
    ReadScalar(dataStuff, ds3d, "apars", "apars", true, true);
    ReadVec(bfieldStuff, ds3d, "B", "/node_data[0]/values", true, true);
    ReadVec(bfield_allStuff, ds3d, "Bs", "Bs", true, true);
    //CalcX(ds3d);
    vtkm::io::VTKDataSetWriter writer("debug.vtk");
    writer.WriteDataSet(ds3d);
  }
  */

  ds.PrintSummary(std::cout);
  vtkm::io::VTKDataSetWriter writer("grid.vtk");
  writer.WriteDataSet(ds);

  return ds;
}

void
SaveStuff(const vtkm::cont::DataSet& inDS)
{
  vtkm::cont::DataSet ds;
  ds.AddCoordinateSystem(inDS.GetCoordinateSystem());
  ds.SetCellSet(inDS.GetCellSet());

  //ds.AddField(inDS.GetField("B_RZP_Norm"));
  //ds.AddField(inDS.GetField("curl_nb_rzp"));
  ds.AddField(inDS.GetField("B_RZP"));

  //Add gradPsi
  vtkm::Id nPts = ds.GetNumberOfPoints();
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B0;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  inDS.GetField("B_RZP").GetData().AsArrayHandle(B0);
  auto cPortal = coords.ReadPortal();
  auto b0Portal = B0.ReadPortal();

  std::vector<vtkm::Vec3f> gradPsi(nPts);
  for (vtkm::Id i = 0; i < nPts; i++)
  {
    gradPsi[i][0] = b0Portal.Get(i)[2] * cPortal.Get(i)[0];
    gradPsi[i][1] = -b0Portal.Get(i)[0] * cPortal.Get(i)[0];
    gradPsi[i][2] = 0;
  }

  ds.AddField(vtkm::cont::make_FieldPoint("gradPsi", vtkm::cont::make_ArrayHandle(gradPsi, vtkm::CopyFlag::On)));

  vtkm::io::VTKDataSetWriter writer("stuff.vtk");
  writer.WriteDataSet(ds);
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


class FindMaxR : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn Thetas,
                                ExecObject Locator,
                                FieldOut MaxR);
  using ExecutionSignature = void(_1, _2, _3);
  using InputDomain = _1;

  FindMaxR(vtkm::FloatDefault rmin,
           vtkm::FloatDefault rmax,
           vtkm::FloatDefault eq_r,
           vtkm::FloatDefault eq_z)
    : RMin(rmin)
    , RMax(rmax)
    , EqR(eq_r)
    , EqZ(eq_z)
  {
  }

  template <typename LocatorType>
  VTKM_EXEC void operator()(const vtkm::FloatDefault& theta,
                            const LocatorType& locator,
                            vtkm::FloatDefault& val) const
  {
    const vtkm::FloatDefault cost = vtkm::Cos(theta), sint = vtkm::Sin(theta);

    vtkm::FloatDefault r0 = 0, r1 = this->RMax - this->EqR;
    //std::cout<<"Theta= "<<theta<<std::endl;

    vtkm::Vec3f pt(0,0,0), pcoords;
    vtkm::Id cellId;
    while (vtkm::Abs(r0-r1) > 1e-8)
    {
      vtkm::FloatDefault ri = (r0+r1) / 2.0;

      pt[0] = this->EqR + ri*cost;
      pt[1] = this->EqZ + ri*sint;
      //std::cout<<"  ri= "<<ri<<" ("<<r0<<" "<<r1<<")"<<std::endl;
      //std::cout<<"   pt= "<<pt<<std::endl;
      if (locator.FindCell(pt, cellId, pcoords) == vtkm::ErrorCode::Success)
      {
        // ri is INSIDE, set range to (ri, r1)
        r0 = ri;
      }
      else
      {
        // ri is OUTSIDE, set range to (r0, ri)
        r1 = ri;
      }
    }

    //Set the mid value.
    val = (r0+r1)/2.0;
  }

  vtkm::FloatDefault RMin;
  vtkm::FloatDefault RMax;
  vtkm::FloatDefault EqR;
  vtkm::FloatDefault EqZ;
};

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

void
GeneratePsiRangeSeeds(vtkm::FloatDefault psiNorm0,
                      vtkm::FloatDefault psiNorm1,
                      int numPts,
                      int numTheta,
                      const vtkm::cont::DataSet& ds,
                      std::vector<vtkm::Vec3f>& seeds)
{
  vtkm::cont::Invoker invoker;
  FindMaxR worklet(eq_min_r, eq_max_r, eq_axis_r, eq_axis_z);
  std::vector<vtkm::FloatDefault> t;
  vtkm::FloatDefault dTheta = vtkm::TwoPi() / static_cast<vtkm::FloatDefault>(numTheta);

  for (int j = 0; j < numTheta; j++)
    t.push_back(static_cast<vtkm::FloatDefault>(j) * dTheta);

  auto thetas = vtkm::cont::make_ArrayHandle(t, vtkm::CopyFlag::On);
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxR;

  vtkm::cont::CellLocatorTwoLevel locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();
  invoker(worklet, thetas, locator, maxR);
  std::cout<<"MaxR= ";
  vtkm::cont::printSummary_ArrayHandle(maxR, std::cout, true);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> arr;
  ds.GetField("coeff_2D").GetData().AsArrayHandle(arr);
  auto coeff = arr.ReadPortal();

  auto thetaPortal = thetas.ReadPortal();
  auto maxRPortal = maxR.ReadPortal();
  vtkm::Id ncoeff = eq_mr-1;
  vtkm::Id2 nrz(eq_mr, eq_mz);
  vtkm::Vec2f rzmin(eq_min_r, eq_min_z);
  vtkm::Vec2f drz((eq_max_r-eq_min_r)/static_cast<vtkm::FloatDefault>(eq_mr-1),
                  (eq_max_z-eq_min_z)/static_cast<vtkm::FloatDefault>(eq_mz-1));

  auto psiMin = psiNorm0 * eq_x_psi, psiMax = psiNorm1 * eq_x_psi;
  auto dPsi = (psiMax-psiMin) / static_cast<vtkm::FloatDefault>(numPts-1);
  auto psi = psiMin;

  for (int i = 0; i < numPts; i++)
  {
    std::cout<<"Pt_"<<i<<" = psi= "<<psi<<" psiN= "<<(psi/eq_x_psi)<<std::endl;

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

        vtkm::Vec2f pt(eq_axis_r + rmid*cost, eq_axis_z + rmid*sint);
        auto psimid = InterpolatePsi(pt,  coeff, ncoeff, nrz, rzmin, drz);
        //std::cout<<"     "<<rmid<<" ("<<rad0<<" "<<rad1<<") psimid= "<<psimid<<std::endl;

        if (psimid < psi) // mid is inside, set range to (radmid, rad1)
          rad0 = rmid;
        else
          rad1 = rmid;   // mid is outside, set range to (rad0, radmid)
      }
      while ((rad1-rad0) > 1e-8);

      vtkm::FloatDefault r = (rad0+rad1)/2.0;
      vtkm::Vec3f pt_rpz(eq_axis_r + r*cost, 0, eq_axis_z + r*sint);
      seeds.push_back(pt_rpz);

      //auto psi_ = InterpolatePsi({pt[0], pt[1]},  coeff, ncoeff, nrz, rzmin, drz);
      //std::cout<<"   "<<pt<<"  psi= "<<psi_/eq_x_psi<<std::endl;
    }
    psi += dPsi;
  }
}

int
main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
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

  vtkm::FloatDefault stepSize = std::atof(args["--stepSize"][0].c_str());
  int numPunc = std::atoi(args["--numPunc"][0].c_str());
  std::string vField = args["--vField"][0];

  std::vector<vtkm::Vec3f> seeds;
  int whichWorklet = std::atoi(args["--worklet"][0].c_str());
  bool useTraces = false;
  if (args.find("--traces") != args.end()) useTraces = std::atoi(args["--traces"][0].c_str());
  std::string outFileName = args["--output"][0];
  useTurb = true;
  if (args.find("--turbulence") != args.end())
    useTurb = std::atoi(args["--turbulence"][0].c_str());

  bool quickTest = false;
  if (args.find("--quickTest") != args.end())
    quickTest = true;
  vtkm::Id BGridSize = 512;
  bool BGridCell = true;
  if (args.find("--BGridSize") != args.end())
  {
    BGridSize = std::atoi(args["--BGridSize"][0].c_str());
    BGridCell = std::atoi(args["--BGridSize"][1].c_str());
  }
  std::cout<<"Grid: "<<BGridSize<<" "<<BGridCell<<std::endl;


  std::string locator, locParam1, locParam2;
  if (args.find("--Locator") != args.end())
  {
    const auto& la = args["--Locator"];
    locator = la[0];
    if (la.size() == 3)
    {
      locParam1 = la[1];
      locParam2 = la[2];
    }
  }

  auto ds = ReadData(args);

  //ds.PrintSummary(std::cout);
  //return 0;
  /*
  SaveStuff(ds);
  Debug(ds);
  return 0;
  */


  bool useBOnly = false, useHighOrderB = false;
  if (args.find("--useBOnly") != args.end()) useBOnly = true;
  if (args.find("--useHighOrderB") != args.end()) useHighOrderB = true;

  if (args.find("--range") != args.end())
  {
    vtkm::Id numSeeds = std::stoi(args["--numSeeds"][0]);

    auto vals = args["--range"];
    vtkm::FloatDefault r0 = std::atof(vals[0].c_str());
    vtkm::FloatDefault r1 = std::atof(vals[1].c_str());

    vtkm::FloatDefault dr = (r1-r0) / (float)(numSeeds-1);
    vtkm::FloatDefault r = r0;

    for (vtkm::Id id = 0; id < numSeeds; id++, r+=dr)
      seeds.push_back({r, .1, 0});
    //std::cout<<"SEEDS= "<<seeds<<std::endl;
  }
  else if (args.find("--psiRange") != args.end())
  {
    auto vals = args["--psiRange"];
    vtkm::FloatDefault psi0 = std::atof(vals[0].c_str());
    vtkm::FloatDefault psi1 = std::atof(vals[1].c_str());
    int numPts = std::atof(vals[2].c_str());
    int numTheta = std::atof(vals[3].c_str());

    GeneratePsiRangeSeeds(psi0, psi1, numPts, numTheta, ds, seeds);
  }
  else if (args.find("--seed") != args.end())
  {
    auto vals = args["--seed"];
    vtkm::FloatDefault r = std::atof(vals[0].c_str());
    vtkm::FloatDefault z = std::atof(vals[1].c_str());
    vtkm::FloatDefault t = std::atof(vals[2].c_str());
    seeds.push_back({r, t, z});
  }
  else if (args.find("--jong1") != args.end())
  {
    //first point in traces.v2
    //zone = 3132, incident nodes: 1521, 1612, 1613
    seeds = {{3.029365, 6.183185, 0.020600}};

    //take position at vid=1512
    //seeds = {{2.98872, 0, -0.113239}};

    //first point in log file w/ dpsi_dr
    //seeds = {{3.030292040947820, 0, 0}};

    //testing seed...
    //seeds = {{2,0,0}};
    //seeds = {{2,0,.5}};
    //seeds = {{2,0,-.5}};
    //seeds = {{2.8, 0, -.99}};

    /*
      Jong's changes in diagnosis.F90
      R,Z = 3.0, 0.0: psi,dpsi_dr,dpsi_dz: 1.0227020266024015E-002   4.6889284669612431E-002  -4.6889284669613417E-002
     */
//    seeds = {{3, 0, 0}};

    //first point in MEOW MEOW MEOW
    //Particle: RZP= 2.986310165629829 0.1353622587347799 6.183185307179587
    //seeds = {{2.986310165629829, 6.183185307179587, 0.1353622587347799}};
//    seeds = {{3.019020736951258, 6.183185307179587, 7.1164151319539654E-002}};

    //seed for runs on summit.
    //seeds = {{2.728835848680459, 6.183185307179587, 0.2190207369512611}};
/*
DRP:  derivs_sp MEOW MEOW MEOW ***********************************************
 Particle: RZP=     3.019020736951258        7.1164151319539654E-002     6.183185307179587
  B=   -5.8549499079761569E-003   1.8019644771892874E-002   -0.2192499017639860      mag=   0.2200670521901169   4.544069591735594
 ****************************************************************
 ****************************************************************
 Calc curl_B
   over_B    4.544069591735594
   inv_r   0.3312332332668388
   fld.br/bz/phi  -5.8549499079761569E-003   1.8019644771892874E-002   -0.2192499017639860
   fld.dbRdx:    7.5098966316331705E-003  -8.0453390321500493E-002     0.000000000000000
   fld.dbZdx:    5.9158099337405949E-002  -5.5705426429988473E-003     0.000000000000000
   fld.dbPdx:    7.2622853854721864E-002   -0.000000000000000     0.000000000000000
   dbdr =   -6.7708980323630388E-002
   dbdz =    1.6843565038781637E-003
   dbdp =     0.000000000000000
   As  =     2.0422619291435113E-006
   Bsa2 =   -1.5573090758542901E-008  -6.2601835972401922E-007   -1.3466549373597934E-006
   Bsa =    9.8643536496446089E-006  -7.1673502824034881E-005   -7.4497222472069995E-006
   dAs =    7.1300391181276443E-005   9.9533406652970541E-006   -4.4521462973275313E-007
fld stuff....
   fld.rzp=     3.019020736951258        7.1164151319539654E-002     6.183185307179587
   fld.psi=    6.9734500328945031E-003
   fld.B_rzp=   -5.8549499079761569E-003   1.8019644771892874E-002   -0.2192499017639860
   fld.As=    2.0422619291435113E-006
   fld.dAs_rzp=   -7.1300391181276443E-005  -9.9533406652970541E-006    4.4521462973275313E-007
   dpsi_dr=    5.4401681238839920E-002
   dpsi_dz=    1.7676215185990881E-002
   deltaB?=     0.000000000000000         0.000000000000000     0.000000000000000
 ****************************************************************
 ****************************************************************


 * bummy ***************************************************************
 * bummy ***************************************************************
 DRP: here in efidld_common.F90
   fld.rzp=     3.019020736951258        7.1164151319539654E-002
    6.183185307179587
   fld.psi=    6.9734500328945031E-003
   fld.B_rzp=   -5.8549499079761569E-003   1.8019644771892874E-002
  -0.2192499017639860
   fld.As=    2.0422619291435113E-006            ******** VTKm = As= 2.01509e-06
   fld.dAs_rzp=   -7.1300391181276443E-005  -9.9533406652970541E-006  4.4521462973275313E-007
vTKm=  gradAs_rpz= [-7.13004e-05,4.45215e-07,-9.95334e-06]

   dpsi_dr=    5.4401681238839920E-002
   dpsi_dz=    1.7676215185990881E-002
   gamma_psi=    17.48211272627623
   wp=    2.9502307908158221E-002
   iphi=            47
   node=          1617
   irhoX=            47           47
   irho0p1=            47
   wphi=    0.7639437268410916        0.2360562731589084
   rvec=    0.9510563239163462        0.3090175864554083
   zvec=   -0.3090175864554083        0.9510563239163462
   deltaB?=     0.000000000000000         0.000000000000000     0.000000000000000
   As calculation: wrho0=     1.000000000000000         0.000000000000000
                   wrho1=     1.000000000000000         0.000000000000000
                   irho0=            47
                   irho1=            47
                   irho0p1=            47
                   irho1p1=            47
                   node=          1617
                   VECPOTS(0,irho0, node)    2.4447478377585865E-006
                   VECPOTS(1,irho1, node)    2.8654162153181949E-006
                   VECPOTS(0,irho0p1,node)    2.4447478377585865E-006
                   VECPOTS(1,irho1p1,node)    2.8654162153181949E-006

 * zoommy ***************************************************************
 * zoommy ***************************************************************


 */





    /*
      B=  -1.1258751808180879E-002   1.5496320905027214E-002 -0.2216514572458676
      dAs   7.3166248873948126E-005   9.5118595155724718E-005   4.0429622421441507E-006
      PTL1(R,Z,P,rho)    2.986310165629829  0.1353622587347799         6.183185307179587       -4.5917006108385793E-004
      Bsa   9.5002789220787541E-005 -7.3646160159566971E-005  -1.1920069817798782E-005
      Bvec  -1.1258751808180879E-002 1.5496320905027214E-002  -0.2216514572458676
      As   3.0485639983994535E-006
      curl_nb  -1.4504793857055040E-002 -0.3136527384215291       -0.6593914097083449
      curl_B    0.000000000000000 -4.7905885599546799E-018  -0.1419851267560749

      Calc curl_B
        over_B    4.494835378368263
        inv_r   0.3348613990298943
        BVALUE: fld.br/bz/phi  -1.1258751808180879E-002   1.5496320905027214E-002 -0.2216514572458676
        fld.dbRdx:    1.2882152089190432E-002  -7.6545946509187848E-002  0.000000000000000
        fld.dbZdx:    6.5439180246887108E-002  -9.1120307073726311E-003  0.000000000000000
        fld.dbPdx:    7.4222517070366034E-002   -0.000000000000000  0.000000000000000

        dbdr =   -7.0040769970213385E-002
        dbdz =    3.2390182056756529E-003
        dbdp =     0.000000000000000

     //stuff from efield_common.F90
     DRP: here in efidld_common.F90
     fld.rzp=     3.019020736951258 7.1164151319539654E-002  6.183185307179587
     fld.psi=    6.9734500328945031E-003
     fld.B_rzp=   -5.8549499079761569E-003   1.8019644771892874E-002 -0.2192499017639860
     fld.As=    2.0422619291435113E-006
     fld.dAs_rzp=   -7.1300391181276443E-005  -9.9533406652970541E-006 4.4521462973275313E-007
     dpsi_dr=    5.4401681238839920E-002
     dpsi_dz=    1.7676215185990881E-002
     gamma_psi=    17.48211272627623

     wp=    2.9502307908158221E-002
     iphi=            47
     node=          1617
     irhoX=            47           47
     irho0p1=            47
     wphi=    0.7639437268410916        0.2360562731589084
     rvec=    0.9510563239163462        0.3090175864554083
     zvec=   -0.3090175864554083        0.9510563239163462
     deltaB?=     0.000000000000000         0.000000000000000    0.000000000000000
     */

    /*
      x_ff debugging:
      DRP: field_following_pos2()
      x=    3.019020736951258 7.1164151319539654E-002 -->
            3.021638837619697 6.2535188082394388E-002
      dphi=   1.7275076525105959E-002
      phi   6.183185307179587  -->   6.217735460229799

////use this one.
pIn     2.728835848680459        0.2190207369512611        6.183185307179587
 *************************************************
  i=             1
  x=    2.728835848680459        0.2190207369512611      phi=   6.183185307179587
     bvec_interpol:     2.728835848680459        0.2190207369512611       -->   -1.9935856993796582E-002  -6.4775663199158279E-003  -0.2425649752146412
 b/bphi=   8.2187698269940740E-002   2.6704458523675777E-002 x(1)=   2.728835848680459      dx=    0.2242767373595473      7.2872083739006915E-002
    bvec_interpol: dx=   0.2242767373595473        7.2872083739006915E-002
  dx1=   0.2242767373595473        7.2872083739006915E-002
 k1_rpz= 0.2242767373595448,     0.07287208373900525]                        GOOD

   x_tmp=    2.730773047580803      0.2196501723628287
    tmp= [2.730773047580803,6.181090142925015,0.2196501723628287]



     bvec_interpol:     2.730773047580803        0.2196501723628287       -->
  -1.9978785301185673E-002  -6.2967293967192660E-003  -0.2423929006426939
 b/bphi=   8.2423145431292824E-002   2.5977367241465283E-002 x(1)=
    2.730773047580803      dx=    0.2250789040406072
   7.0938294310101874E-002
    bvec_interpol: dx=   0.2250789040406072        7.0938294310101874E-002
  dx2=   0.2250789040406072        7.0938294310101874E-002
   x_tmp=    2.730779976326204        0.2196334691826448
     bvec_interpol:     2.730779976326204        0.2196334691826448       -->
  -1.9977575372099224E-002  -6.2961967126177248E-003  -0.2423922856247466
 b/bphi=   8.2418362946694565E-002   2.5975235541798633E-002 x(1)=
    2.730779976326204      dx=    0.2250664152164190
   7.0932653097900436E-002
    bvec_interpol: dx=   0.2250664152164190        7.0932653097900436E-002
  dx3=   0.2250664152164190        7.0932653097900436E-002
   x_tmp=    2.732723888226554        0.2202461039616561
     bvec_interpol:     2.732723888226554        0.2202461039616561       -->
  -2.0018867953906675E-002  -6.1149800853820989E-003  -0.2422198608691359
 b/bphi=   8.2647508268210379E-002   2.5245576739414608E-002 x(1)=
    2.732723888226554      dx=    0.2258528201469401
   6.8989190627854940E-002
    bvec_interpol: dx=   0.2258528201469401        6.8989190627854940E-002
  dx4=   0.2258528201469401        6.8989190627854940E-002
   ** x=    2.732723950718343        0.2202461248374213
 *************************************************


     */
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


      /*
      {3.351443028564415449, 0.0, -0.451648806402756176}, //blue 3 islands.
      {3.187329423521033878, 0.0, -0.665017624967372267},
      {1.992020349316277139, 0.0, -0.126203396421661285},
      {3.018666196722858963, 0.0, 0.073864239629065770},
      {3.176582679765305173, 0.0, -0.220557108925872658},  //stochastic region 4000
      {2.179226604128697176, 0.0, 0.291539359807166554},

      {2.552260904008052389, 0, -0.003112355795000767}, //stochastic: 100
      {2.904147293825020348, 0, 0.233676309266207888},  //stochastic: 500
      {2.834116132387578091, 0, 0.301602359825420996},  //stochastic: 1000
      {3.221313019897226848, 0, -0.169787551332608172}, //stochastic: 5000
      {2.1299283004919225 0, -0.17650886057033158},     //semi stoch: pt 10153: xgc spread out, vtkm thin.
      */

    if (vals.size() == 0) //all the points.
      seeds = allSeeds;
    else
    {
      int n = std::stoi(vals[0]);
      if (n >= (int)allSeeds.size())
        std::cout<<"Bad seed!!! #allseeds= "<<allSeeds.size()<<std::endl;

      seeds = {allSeeds[n]};
    }

    /*
allSeeds[0] 50 punctures
--stepSize 0.01   Total error=  0.0120200778585 maxErr= 0.000484904284968
--stepSize 0.005  Total error=  0.00764654588018 maxErr= 0.000306355820869
--stepSize 0.001  Total error=  0.00458809221615 maxErr= 0.000220362884491

-stepSize  0.00050 Total error=  0.00410665950475 maxErr= 0.000212761006711
-stepSize  0.00025 Total error=  0.00395742797169 maxErr= 0.000212790194159
-stepSize  0.00010 Total error=  0.00384771565615 maxErr= 0.000209736086882

allSeeds[1] 50 punctures
--stepSize 0.01   Total error=  0.0246758872498 maxErr= 0.00102710677529
--stepSize 0.005  Total error=  0.0204393203537 maxErr= 0.000897389164912
--stepSize 0.001  Total error=  0.0168112860277 maxErr= 0.00089728753628

allSeeds[2] 50 punctures
--stepSize 0.01  Total error=  0.020815736223 maxErr= 0.000694933702884
--stepSize 0.005 Total error=  0.0157798325015 maxErr= 0.000548627689596
--stepSize 0.001 Total error=  0.012512912312 maxErr= 0.000501483340821


allSeeds[3] 50 punctures
--stepSize 0.01  Total error=  0.00930387190718 maxErr= 0.000376516678999
--stepSize 0.005 Total error=  0.00660720169504 maxErr= 0.000267672744868
--stepSize 0.001 Total error=  0.00501206967247 maxErr= 0.000244691060404


allSeeds[4] 50 punctures ##Stochastic region
--stepSize 0.01    Total error=  3.11058237394 maxErr= 0.417674988212
--stepSize 0.005   Total error=  3.10614570903 maxErr= 0.422411249578
--stepSize 0.001   Total error=  3.0891389475 maxErr= 0.424464351612
--stepSize 0.0005  Total error=  3.10299880864 maxErr= 0.423057887001
--stepSize 0.00025 Total error=  3.10476342536 maxErr= 0.421858627004


10 punctures:
--stepSize 0.000001


allSeeds[5] 50 punctures
--stepSize 0.01  Total error=  0.0173742171359 maxErr= 0.000733842147223
--stepSize 0.005 Total error=  0.0130991456975 maxErr= 0.000589881209904
--stepSize 0.001 Total error=  0.00999171349076 maxErr= 0.000566749317839

allSeeds[5] PID=10153 50 punctures
--stepSize 0.001  Total error=  0.00630002849371 maxErr= 0.000329808548949
--stepSize 0.0005 Total error=  0.00625729353866 maxErr= 0.000329798904507



     */

    //traces.v2 pt near begining.
    //seeds = {{3.024768, 6.070249, 0.049700}};

    //seeds from data/sku_8000/jong.py
    //pts in: xgc_theta_psi.txt, xgc_punctures.txt
//    seeds = {
//      {3.351443028564415449, 0.0, -0.451648806402756176}, //blue 3 islands.
//      {3.187329423521033878, 0.0, -0.665017624967372267},
//      {1.992020349316277139, 0.0, -0.126203396421661285},
//      {3.018666196722858963, 0.0, 0.073864239629065770},
//      {3.176582679765305173, 0.0, -0.220557108925872658},  //stochastic region
//      {2.179226604128697176, 0.0, 0.291539359807166554},

      //stochastic region, i=4000
      /*
s         3.176582679765305173, -0.220557108925872658
p1        2.673415694283511446, -0.404035922651817592
p2        2.366230557551912916, -0.054420066925619182
p3        2.602919798501353910, 0.404861298966848748
p4        3.164266550057113658, 0.258207175353824703
p5        3.087268470563441003, -0.311153651480581217
p6        2.556488556701177028, -0.370963671520662341
p7        2.363954167247493743, 0.104848146584189700
p8        2.744892558803513793, 0.421606828485864227
p9        3.244749979529466088, 0.034975199218512609
p10       2.780834806912544810, -0.455027987386564192
p11       2.329460849125147615, -0.073678279004152566
       */
//    };
  }
  else if (args.find("--jongrz") != args.end())
  {
    //Seed from the blue island.
    //seeds = {{3.351443,             0, -0.451649}};
    seeds = {{3.351443028564415449, 0, -0.45164880640275617552}};
  }
  else if (args.find("--afterN") != args.end())
  {
    seeds = {
      {3.321888620239255019, 0.0, 0.478933972623594384},
      {2.568934684085571352, 0.0, 0.731290913908178353},
      {3.493628202658771720, 0.0, 0.433951677589735296},
      {2.862485694515508605, 0.0, 0.208737305948038576},
      {2.905837753215041008, 0.0, -0.397811882628356817},
      {3.391834939600261389, 0.0, -0.350011953142094434},
    };
  }
  else if (args.find("--parse") != args.end())
  {
//    ./examples/poincare/Simple2.3 --vField B --dir ../data/sku_8000/POINC --worklet 1 --traces 0 --useHighOrder --parse ../data/sku_8000/seeds.txt   --output bumm --numPunc 1000 --stepSize 0.001

    //Generaate the seed list by running jongAll.py
    std::cout<<"READING: "<<args["--parse"][0]<<std::endl;
    std::ifstream seedFile;
    seedFile.open(args["--parse"][0]);
    std::string line;
    while (std::getline(seedFile, line))
    {
      vtkm::FloatDefault r, p, z;
#ifdef VTKM_USE_DOUBLE_PRECISION
      sscanf(line.c_str(), "%lf, %lf, %lf", &r, &p, &z);
#else
      sscanf(line.c_str(), "%f, %f, %f", &r, &p, &z);
#endif
      seeds.push_back({r,p,z});
    }
    std::cout<<"NumSeeds= "<<seeds.size()<<std::endl;
  }

  std::vector<std::vector<vtkm::Vec3f>> traces(seeds.size());
  vtkm::cont::ArrayHandle<vtkm::Vec2f> outRZ, outTP;
  vtkm::cont::ArrayHandle<vtkm::Id> outID;
  std::cout<<"******************* NumSeeds= "<<seeds.size()<<std::endl;
  Poincare(ds, seeds, args, vField, stepSize, numPunc, whichWorklet, useBOnly, useHighOrderB, outRZ, outTP, outID, quickTest, BGridSize, BGridCell, locator, locParam1, locParam2, (useTraces ? &traces : nullptr));

  //std::cout<<"Convert to theta/psi"<<std::endl;
  //auto puncturesTP = ConvertPuncturesToThetaPsi(punctures, ds);

  std::cout<<"TRACES: "<<traces.size()<<std::endl;

  std::cout<<"SaveOutput"<<std::endl;
  SaveOutput(traces, outRZ, outTP, outID, outFileName);

#if 0
  std::ofstream outPts, outPtsPsiTheta;

  outPts.open("punctures.txt");
  outPts<<"ID,R,Z,T"<<std::endl;
  for (int i = 0; i < (int)punctures.size(); i++)
    for (const auto& p : punctures[i])
    {
      outPts<<i<<", "<<p[0]<<","<<p[2]<<","<<p[1]<<std::endl;
    }

  std::ofstream outTraces, RZ, thetaPsi;
  outTraces.open("traces.txt"), RZ.open("rz.txt"), thetaPsi.open("thetaPsi.txt");
  outTraces<<"ID,R,Z,T,THETA,PSI"<<std::endl;
  RZ<<"ID,R,Z,T"<<std::endl;
  thetaPsi<<"ID,theta,psi,Z"<<std::endl;
  for (int i = 0; i < (int)traces.size(); i++)
  {
    int idx = 0;
    for (const auto& p : traces[i])
    {
      auto R = p[0];
      auto Z = p[2];
      auto PHI = p[1]; //std::fmod(p[1], vtkm::TwoPi());
      int numRevs = 0;
      auto PHI_N = PHI;
      while (PHI_N < 0)
      {
        PHI_N += vtkm::TwoPi();
        numRevs++;
      }
      //auto PHI_N = PHI + (numRevs*vtkm::TwoPi());
//      PHI = p[1];
//      if (PHI < 0) PHI = -PHI;
      auto theta = vtkm::ATan2(Z-eq_axis_z, R-eq_axis_r);
      if (theta < 0) theta += vtkm::TwoPi();
      auto psi = vtkm::Sqrt(((R-eq_axis_r)*(R-eq_axis_r) + Z*Z));

      outTraces<<idx<<", "<<R<<", "<<Z<<", "<<PHI_N<<", "<<theta<<", "<<psi<<std::endl;
//      outTraces<<idx<<", "<<PHI<<" "<<PHI_N<<" nr= "<<numRevs<<std::endl;
      RZ<<idx<<", "<<p[0]<<", "<<p[2]<<", 0"<<std::endl;
      thetaPsi<<idx<<", "<<theta<<", "<<psi<<", 0"<<std::endl;
      idx++;
    }
  }
#endif

  return 0;
}


/*
XGC:
3.029365, 0.020600, 6.183185
3.029293, 0.021400, 6.180291



mine:
3.02936, 0.0206, 6.18318
3.02934, 0.0209196, 6.17947

delta= -0.000047, .0004804, .000821


 */
