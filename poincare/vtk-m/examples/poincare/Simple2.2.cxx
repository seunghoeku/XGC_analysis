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
#include <vtkm/Geometry.h>
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

#include <vtkm/cont/CellLocatorGeneral.h>

#include <vtkm/filter/Gradient.h>

#include <vtkm/io/VTKDataSetWriter.h>

#include <fides/DataSetReader.h>

#include <adios2.h>
#include <random>
#include <chrono>
#include <variant>
#include <filesystem>
#include <mpi.h>

/*
radius values: max range: 2.5 3.7
for interesting psi: 3.5  3.7
*/

/*
TODO:
Make sure that wrap around works.
Get the Bs in 3D working.

*/


class Rock
{
public:
  Rock() {}
  Rock(const vtkm::Id& id, const vtkm::Vec3f& p) : ID(id), Pos(p) {}

  vtkm::Id ID;
  vtkm::Vec3f Pos;
  vtkm::Id Punc=0;
  bool Valid=true;
};

adios2::ADIOS *adios = NULL;
class adiosS;
std::map<std::string, adiosS*> adiosStuff;

int numPlanes = -1;
int numNodes = -1;
int numTri = -1;
float XScale = 1;
vtkm::FloatDefault eq_axis_r = 2.8, eq_axis_z = 0.0;
//  vtkm::FloatDefault eq_x_psi = 0.0697345, eq_x_r = 2.8, eq_x_z = -0.99988;



using Ray3f = vtkm::Ray<vtkm::FloatDefault, 3, true>;

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v)
{
  out<<"[";
  for (const auto& x : v)
    out<<x<<" ";
  out<<"]";
  return out;
}

std::vector<vtkm::Vec3f>
ConvertToThetaPsi(const std::vector<vtkm::Vec3f>& pts)
{
  std::vector<vtkm::Vec3f> output;
  for (const auto& p : pts)
  {
    auto R = p[0];
    auto Z = p[2];
    auto psi = vtkm::Sqrt(((R-eq_axis_r)*(R-eq_axis_r) + Z*Z));
    auto theta = vtkm::ATan2(Z-eq_axis_z, R-eq_axis_r);

    output.push_back({theta, psi, 0});
  }

  return output;
}

void
GetPlaneIdx(const vtkm::FloatDefault& phi,
            const vtkm::Id& nPlanes,
            vtkm::FloatDefault& phiN,
            vtkm::Id& plane0,
            vtkm::Id& plane1,
            vtkm::FloatDefault& phi0,
            vtkm::FloatDefault& phi1,
            vtkm::Id& numRevs,
            vtkm::FloatDefault& T)
{
  vtkm::FloatDefault dPhi = vtkm::TwoPi()/static_cast<double>(nPlanes);
  std::vector<vtkm::FloatDefault> PHIs;

  /*
  vtkm::FloatDefault p = 0;
  for (int i = 0; i < nPlanes+1; i++, p+= dPhi)
    PHIs.push_back(p);
  std::cout<<"PHIs= "<<PHIs<<std::endl;
  */

  numRevs = vtkm::Floor(vtkm::Abs(phi / vtkm::TwoPi()));
  //rem = std::fmod(vtkm::Abs(phi), vtkm::TwoPi());
  phiN = phi;
  if (phi < 0)
  {
    //rem = -rem;
    phiN += ((1+numRevs) * vtkm::TwoPi());
  }

  plane0 = vtkm::Floor(phiN / dPhi);
  plane1 = plane0 + 1;
  phi0 = static_cast<vtkm::FloatDefault>(plane0)*dPhi;
  phi1 = static_cast<vtkm::FloatDefault>(plane1)*dPhi;
  if (plane1 == nPlanes)
    plane1 = 0;
  T = (phiN-phi0) / (phi1-phi0);

//  std::cout<<phi<<"  : ("<<phi0<<" "<<phi1<<") :: ["<<plane0<<" "<<plane1<<"]"<<std::endl;

  return;

  /*
  plane0 = -1;
  plane1 = -1;
  vtkm::FloatDefault phi0 = -1.0f, phi1 = -1.0f;
  for (int i = 0; i < PHIs.size()-1; i++)
  {
    if (phiN > PHIs[i] && phiN <= PHIs[i+1])
    {
      plane0 = i; //(i % nPlanes);
      plane1 = i+1; //((i+1) % nPlanes);
      break;
    }
  }
  auto S = phiN / dPhi;
  auto S0 = vtkm::Floor(S);
  auto S1 = S0 + 1;
  if (S1 == nPlanes)
    S1 = 0;
  std::cout<<" **** S= "<<S<<" ("<<S0<<" "<<S1<<")"<<std::endl;

  phi0 = (plane0*dPhi);
  phi1 = (plane1*dPhi);
  if (plane1 == nPlanes)
    plane1 = 0;
  T = (phiN-phi0) / (phi1-phi0);
  */


  std::cout<<"phiN= "<<phiN<<" ";
  std::cout<<"plane0,plane1= "<<plane0<<" "<<plane1<<std::endl;
  std::cout<<"Phi0, Phi1= "<<phi0<<" "<<phi1<<std::endl;
  std::cout<<"      T= "<<T<<std::endl;

/*
  vtkm::Id idx = vtkm::Floor(phiN / dPhi);
  std::cout<<"    idx= "<<idx<<" :: "<<rem<<" "<<dPhi<<" phiN= "<<phiN<<std::endl;
  if (idx < 0)
    idx += nPlanes;

  if (idx == nPlanes-1)
  {
    plane0 = 0;
    plane1 = idx;
  }
  else
  {
    plane0 = idx;
    plane1 = plane0+1;
  }
*/
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
ComputeV(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> b;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> A_s;
  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;

  ds.GetField("B").GetData().AsArrayHandle(b);
  ds.GetField("apars").GetData().AsArrayHandle(A_s);

  auto bPortal = b.ReadPortal();
  auto aPortal = A_s.ReadPortal();

  vtkm::Id n = b.GetNumberOfValues();
  std::vector<vtkm::Vec3f> As_bHat(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f bn = vtkm::Normal(bPortal.Get(i));
    vtkm::FloatDefault As = aPortal.Get(i);
    As_bHat[i] = bn * As;
  }

  ds.AddField(vtkm::cont::make_FieldPoint("As_bHat", vtkm::cont::make_ArrayHandle(As_bHat, vtkm::CopyFlag::On)));

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetComputeVorticity(true);
  gradient.SetActiveField("As_bHat");
  gradient.SetOutputFieldName("grad_As_bHat");
  std::cout<<"Compute Grad"<<std::endl;
  ds = gradient.Execute(ds);
  std::cout<<"Compute Grad DONE"<<std::endl;

#if 0
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
#endif

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f,3>> grad;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("grad_As_bHat").GetData().AsArrayHandle(grad);

  auto cPortal = coords.ReadPortal();
  auto gPortal = grad.ReadPortal();

  std::vector<vtkm::Vec3f> V(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f val = bPortal.Get(i);
    vtkm::FloatDefault R = cPortal.Get(i)[0];
    auto g = gPortal.Get(i);

    //From: https://www.therightgate.com/deriving-curl-in-cylindrical-and-spherical/
    //R: (1/R * dAz/dT  - dAT/dZ)
    //T: dAr/dZ - dAz/dr
    //Z: Az/R + dAt/dr - 1/R dAr/dT]
    vtkm::FloatDefault rv, tv, zv;
    rv = 1/R * g[2][1] - g[1][2];
    tv = g[0][2] - g[2][1];
    zv = As_bHat[i][1]/R + g[1][0] - 1/R*g[0][1];
    val[0] += rv;
    val[1] += tv;
    val[2] += zv;

    V[i] = val;
  }

  ds.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(V, vtkm::CopyFlag::On)));

  std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
}

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

  ds.AddField(vtkm::cont::make_Field(vname,
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     tmp,
                                     vtkm::CopyFlag::On));
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

  if (addExtra && add3D)
  {
    for (int i = 0; i < numNodes; i++)
      tmp.push_back(tmp[i]);
  }

  ds.AddField(vtkm::cont::make_FieldPoint(vname+"2D", vtkm::cont::make_ArrayHandle(tmpPlane, vtkm::CopyFlag::On)));
  if (add3D)
    ds.AddField(vtkm::cont::make_FieldPoint(vname, vtkm::cont::make_ArrayHandle(tmp, vtkm::CopyFlag::On)));
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

  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B;
  inDS.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  inDS.GetField("B_RZP_2D").GetData().AsArrayHandle(B);

  //Get the gradients.
  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  //gradient.SetActiveField("BRZP_Norm");
  gradient.SetActiveField("B_RZP_2D");
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
    auto Z = ptRPZ[2];

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
    auto dbpdp = GRAD[2][2];

    //dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
    //dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
    //dbdphi=0D0  ! no B perturbation
    auto dbdr = (br*dbrdr + bphi*dbpdr + bz*dbzdr) * over_Bmag;
    auto dbdz = (br*dbrdz + bphi*dbpdz + bz*dbzdz) * over_Bmag;
    auto dbdphi = 0;


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


  //old way...
#if 0
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("B_RZP_2D").GetData().AsArrayHandle(B);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;

  //wrong way...
  /*
    vtkm::filter::Gradient gradient;
    gradient.SetComputePointGradient(true);
    gradient.SetComputeVorticity(true);
    gradient.SetActiveField("B2D_Norm");
    auto tmpDS = gradient.Execute(ds);
    tmpDS.GetField("Vorticity").GetData().AsArrayHandle(curlBNorm);
   */

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetActiveField("BRZP_Norm");
  auto tmpDS = gradient.Execute(ds);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f, 3>> gradients;
  tmpDS.GetField("Gradients").GetData().AsArrayHandle(gradients);

  vtkm::Id numPts = coords.GetNumberOfValues();
  curlBNorm.Allocate(numPts);

  auto gPortal = gradients.ReadPortal();
  auto cPortal = coords.ReadPortal();
  auto bPortal = B.ReadPortal();
  auto portal = curlBNorm.WritePortal();

  for (vtkm::Id i = 0; i < numPts; i++)
  {
    vtkm::Vec3f ptRZ = cPortal.Get(i);
    auto R = ptRZ[0];
    auto invR = 1.0f / R;
    auto BPhi = bPortal.Get(i)[2]; //This needs to be the BPhi at the position ?
    //std::cout<<__LINE__<<" fix bPhi!"<<std::endl;

    vtkm::Vec<vtkm::Vec3f, 3> grad = gPortal.Get(i);
    // mesh is in R,Z space.
    auto _dBr = grad[0];
    auto _dBz = grad[1];
    auto _dBphi = grad[2];

//    if (R > 3.0 && R < 3.1)
//      std::cout<<"B= "<<bPortal.Get(i)<<" GRAD: "<<_dBr<<" "<<_dBz<<" "<<_dBphi<<std::endl;

    //((dBr/dR, dBz/dR, dBphi/dR) (dBr/dz, dBz/dz, dBphi/dz) (dBr/dPhi, dBz/dPhi, dBphi/dPhi))
    auto dBr_dR     = grad[0][0];
    auto dBz_dR     = grad[0][1];
    auto dBphi_dR   = grad[0][2];

    auto dBr_dZ     = grad[1][0];
    auto dBz_dZ     = grad[1][1];
    auto dBphi_dZ   = grad[1][2];

    auto dBr_dPhi   = grad[2][0];
    auto dBz_dPhi   = grad[2][1];
    auto dBphi_dPhi = grad[2][2];

    vtkm::Vec3f curl;
    //curl_R = 1/R * dBz/dPhi - dBphi/dZ
    curl[0] = invR* dBz_dPhi - dBphi_dZ;

    //curl_Phi = dBr/dZ - dBz/dR
    curl[1] = dBr_dZ - dBz_dR;

    //curl_Z = BPhi / R + dPhi_dR - dBr_dPhi *1/R
    //curl_B(2)  = fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r
    //curl_B.Z = BPhi * invR + dBphi_dR - dBr_dphi*invR
    curl[2] = BPhi * invR  + dBphi_dR  - dBr_dPhi * invR;

    //std::cout<<"dBr/dPhi= "<<dBr_dPhi<<std::endl;
    //curl[2] = invR*(BPhi + R*dBphi_dR - dBr_dPhi);
    //curl(2) = 1/r*(B_phi + dB_phi_dR * R - dB_r/dPhi);

    //std::cout<<"********CURL: "<<curl<<" "<<vtkm::Magnitude(curl)<<std::endl;
    portal.Set(i, curl);
  }

  return curlBNorm;
#endif
}

//This uses the curl_b
vtkm::cont::ArrayHandle<vtkm::Vec3f>
ComputeCurl2(const vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, Brzp;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("B_RZP_2D").GetData().AsArrayHandle(Brzp);

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetActiveField("B_RZP_2D");
  auto tmpDS = gradient.Execute(ds);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f, 3>> gradients;
  tmpDS.GetField("Gradients").GetData().AsArrayHandle(gradients);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBNorm;
  vtkm::Id numPts = coords.GetNumberOfValues();
  curlBNorm.Allocate(numPts);

  auto gPortal = gradients.ReadPortal();
  auto cPortal = coords.ReadPortal();
  auto bPortal = Brzp.ReadPortal();
  auto portal = curlBNorm.WritePortal();

  for (vtkm::Id i = 0; i < numPts; i++)
  {
    vtkm::Vec3f ptRZ = cPortal.Get(i);
    auto R = ptRZ[0];
    auto invR = 1.0f / R;
    auto B = bPortal.Get(i);
    auto BR = B[0];
    auto BZ = B[1];
    auto BPhi = B[2];
    auto BMag = vtkm::Magnitude(B);
    auto over_B = 1.0/BMag;
    auto over_B2 = over_B * over_B;

    vtkm::Vec<vtkm::Vec3f, 3> grad = gPortal.Get(i);
    if (R > 3.0 && R < 3.1)
      std::cout<<"B= "<<bPortal.Get(i)<<" GRAD: "<<grad[0]<<" "<<grad[1]<<" "<<grad[2]<<std::endl;

    //((dBr/dR, dBz/dR, dBphi/dR) (dBr/dz, dBz/dz, dBphi/dz) (dBr/dPhi, dBz/dPhi, dBphi/dPhi))
    auto dBr_dR   = grad[0][0];
    auto dBz_dR   = grad[0][1];
    auto dBphi_dR = grad[0][2];

    auto dBr_dZ   = grad[1][0];
    auto dBz_dZ   = grad[1][1];
    auto dBphi_dZ = grad[1][2];

    auto dBr_dPhi = grad[2][0];
    auto dBz_dPhi = grad[2][1];
    auto dBphi_dPhi = grad[2][2];


    //dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
    auto db_dr = (BR*dBr_dR + BPhi*dBphi_dR + BZ*dBz_dR) * over_B;

    //dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
    auto db_dz = (BR*dBr_dZ + BPhi*dBphi_dZ + BZ*dBz_dZ) * over_B;

    //dbdphi=0D0  ! no B perturbation
    auto db_dphi = 0;

    vtkm::Vec3f curl_B;

    //Setting curl in R,Phi,Z

    //curl_B.R = 1/R * dBz/dPhi - dBphi/dZ
    curl_B[0] = invR* dBz_dPhi - dBphi_dZ;

    //curl_B.Phi = dBr/dZ - dBz/dR
    curl_B[1] = dBr_dZ - dBz_dR;

    //curl_Z = BPhi / R + dPhi_dR - dBr_dPhi *1/R
    //curl_B(2)  = fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r
    //curl_B.Z = BPhi * invR + dBphi_dR - dBr_dphi*invR
    curl_B[2] = BPhi * invR  + dBphi_dR  - dBr_dPhi * invR;


    vtkm::Vec3f curl_nb;

    //curl_nb.R
    //curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
    curl_nb[0] = curl_B[0]*over_B + (BPhi * db_dz) * over_B2;

    //curl_nb.Phi
    //curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2
    curl_nb[1] = curl_B[1]*over_B + (BZ * db_dr - BR * db_dz) * over_B2;


    //curl_nb.Z
    //curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
    curl_nb[2] = curl_B[1]*over_B + (BPhi * db_dr) * over_B2;

    portal.Set(i, curl_nb);
  }

  return curlBNorm;
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

class CalculateABhat : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn B,
                                FieldIn apars,
                                FieldOut output);
  using ExecutionSignature = void(_1, _2, _3);

  VTKM_EXEC void operator()(const vtkm::Vec3f& bvec,
                            const vtkm::FloatDefault& apars,
                            vtkm::Vec3f& output) const
  {
    output = vtkm::Normal(bvec) * apars;
  }
};

class CalculateVecField : public vtkm::worklet::WorkletMapField
{
public:
  using ControlSignature = void(FieldIn coords,
                                FieldIn B,
                                FieldIn apars,
                                FieldIn gradient,
                                FieldOut output);
  using ExecutionSignature = void(_1, _2, _3, _4, _5);

  template <typename GradientType>
  VTKM_EXEC void operator()(const vtkm::Vec3f& point,
                            const vtkm::Vec3f& bvec,
                            const vtkm::Vec3f& a_bHat,
                            const GradientType& grad,
                            vtkm::Vec3f& output) const
  {
    const auto& R = point[0];
    output = bvec;

    //From: https://www.therightgate.com/deriving-curl-in-cylindrical-and-spherical/
    //R: (1/R * dAz/dT  - dAT/dZ)
    //T: dAr/dZ - dAz/dr
    //Z: Az/R + dAt/dr - 1/R dAr/dT]
    auto rv = 1/R * grad[2][1] - grad[1][2];
    auto tv = grad[0][2] - grad[2][1];
    auto zv = a_bHat[1]/R + grad[1][0] - 1/R*grad[0][1];

    output[0] += rv;
    output[1] += tv;
    output[2] += zv;
  }
};

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
      std::cout<<"Cell not found! "<<point<<std::endl;
      std::cout<<"    ***** Try reducing the step size."<<std::endl;
      this->RaiseError(vtkm::ErrorString(status));
    }
    //ptIndices = cellSet.GetIndices(cellId);
    //auto x = cellSet.GetIndices(cellId);
    //ptIndices = x;
  }
};

void
InterpVector(const vtkm::cont::DataSet& ds,
             const vtkm::cont::CellLocatorGeneral& locator,
             const std::vector<vtkm::Vec3f>& pts,
             const std::string& vName,
             std::vector<vtkm::Vec3f>& out,
             bool is3D = false)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> points;
  std::vector<vtkm::Id> offset(pts.size(), 0);
  if (is3D)
  {
    vtkm::FloatDefault phiSpacing = vtkm::TwoPi() / (vtkm::FloatDefault)(numPlanes);

    auto pts2 = pts;
    for (std::size_t i = 0; i < pts.size(); i++)
    {
      //std::cout<<"Wrap around check: line= "<<__LINE__<<std::endl;
      vtkm::Id off = static_cast<vtkm::Id>(vtkm::Floor(pts2[i][1] / phiSpacing));
      pts2[i][1] = pts2[i][2];
      pts2[i][2] = 0;
      offset[i] = off * numNodes;
//      std::cout<<"******* Offset:  "<<pts2[i]<<" "<<off<<" --> "<<offset[i]<<std::endl;
    }
    points = vtkm::cont::make_ArrayHandle(pts2, vtkm::CopyFlag::On);
  }
  else
    points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::Invoker invoker;

  //Find the cell on the RZ plane.
  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> V;
  ds.GetField(vName).GetData().AsArrayHandle(V);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto vPortal = V.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Id vIds[3];
    vtkm::Id cid = cPortal.Get(i);
    cs.GetCellPointIds(cid, vIds);
    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
//    std::cout<<"CID: "<<cid<<":: "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<std::endl;
//    std::cout<<"  vals: "<<vPortal.Get(vIds[0]+offset[i])<<" "<<vPortal.Get(vIds[1]+offset[i])<<" "<<vPortal.Get(vIds[2]+offset[i])<<std::endl;
    for (vtkm::Id j = 0; j < 3; j++)
      vals.Append(vPortal.Get(vIds[j]+offset[i]));

    vtkm::Vec3f v;
    vtkm::exec::CellInterpolate(vals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), v);
    out.push_back(v);
//    std::cout<<"    "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<std::endl;
//    std::cout<<" -----> "<<v<<std::endl;
  }
}

std::vector<std::vector<vtkm::Vec3f>>
EvalTensor(const vtkm::cont::DataSet& ds,
           const vtkm::cont::CellLocatorGeneral& locator,
           const std::vector<vtkm::Vec3f>& pts,
           const std::string& vName,
           const std::vector<int>& offset)
{
  for (std::size_t i = 0; i < pts.size(); i++)
    if (pts[i][2] != 0)
      std::cout<<"********************************************************** FIX ME: "<<__LINE__<<std::endl;

  auto points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::Off);
  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::Invoker invoker;

  //std::cout<<"EvalVector("<<vName<<"): "<<pts<<" offset= "<<offset<<std::endl;
  //Find the cell on the RZ plane.
  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f,3>> V;
  ds.GetField(vName).GetData().AsArrayHandle(V);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto vPortal = V.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  std::vector<std::vector<vtkm::Vec3f>> out;
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Id vIds[3];
    vtkm::Id cid = cPortal.Get(i);
    cs.GetCellPointIds(cid, vIds);

    auto v0 = vPortal.Get(vIds[0]+offset[i]);
    auto v1 = vPortal.Get(vIds[1]+offset[i]);
    auto v2 = vPortal.Get(vIds[2]+offset[i]);

    std::vector<vtkm::Vec3f> grads(3);
    for (int j = 0; j < 3; j++)
    {
      vtkm::VecVariable<vtkm::Vec3f, 3> vals;
      vals.Append(v0[j]);
      vals.Append(v1[j]);
      vals.Append(v2[j]);

      vtkm::Vec3f v;
      vtkm::exec::CellInterpolate(vals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), v);
      grads[j] = v;
    }
    out.push_back(grads);
  }

  return out;
}

std::vector<vtkm::Vec3f>
EvalVector(const vtkm::cont::DataSet& ds,
           const vtkm::cont::CellLocatorGeneral& locator,
           const std::vector<vtkm::Vec3f>& pts,
           const std::string& vName,
           const std::vector<int>& offset)
{
  for (std::size_t i = 0; i < pts.size(); i++)
    if (pts[i][2] != 0)
      std::cout<<"********************************************************** FIX ME: "<<__LINE__<<std::endl;

  auto points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::Off);
  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::Invoker invoker;

  //std::cout<<"EvalVector("<<vName<<"): "<<pts<<" offset= "<<offset<<std::endl;
  //Find the cell on the RZ plane.
  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> V;
  ds.GetField(vName).GetData().AsArrayHandle(V);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto vPortal = V.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  std::vector<vtkm::Vec3f> out;
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Id vIds[3];
    vtkm::Id cid = cPortal.Get(i);
    cs.GetCellPointIds(cid, vIds);
    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
    vals.Append(vPortal.Get(vIds[0]+offset[i]));
    vals.Append(vPortal.Get(vIds[1]+offset[i]));
    vals.Append(vPortal.Get(vIds[2]+offset[i]));

    //std::cout<<"CID: "<<cid<<":: "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<std::endl;
    //std::cout<<"     p= "<<pPortal.Get(i)<<std::endl;

    vtkm::Vec3f v;
    vtkm::exec::CellInterpolate(vals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), v);
    //v = vtkm::Vec3f(v[0], v[2], v[1]);
    out.push_back(v);
//    std::cout<<"    "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<" --- ("<<pPortal.Get(i)<<") ---> "<<v<<std::endl;
  }
  return out;
}

std::vector<vtkm::Vec3f>
EvalVector(const vtkm::cont::DataSet& ds,
           const vtkm::cont::CellLocatorGeneral& locator,
           const std::vector<vtkm::Vec3f>& pts,
           const std::string& vName)
{
  std::vector<int> offset(pts.size(), 0);
  return EvalVector(ds, locator, pts, vName, offset);
}


std::vector<vtkm::FloatDefault>
InterpScalar(const vtkm::cont::DataSet& ds,
             const vtkm::cont::CellLocatorGeneral& locator,
             const std::vector<vtkm::Vec3f>& pts,
             const std::string& vName,
             const std::vector<int>& offsets)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> points;

  points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::Invoker invoker;

  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> V;
  ds.GetField(vName).GetData().AsArrayHandle(V);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto vPortal = V.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();

  std::vector<vtkm::FloatDefault> out;
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Id vIds[3];
    vtkm::Id cid = cPortal.Get(i);
    cs.GetCellPointIds(cid, vIds);
    vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
    //std::cout<<cid<<":: "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<std::endl;
    for (vtkm::Id j = 0; j < 3; j++)
      vals.Append(vPortal.Get(vIds[0]+offsets[i]));

    vtkm::FloatDefault v;
    vtkm::exec::CellInterpolate(vals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), v);
    out.push_back(v);
    //std::cout<<"    "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<std::endl;
    //std::cout<<" -----> "<<interp<<std::endl;
  }
  return out;
}

void
Evaluate(const vtkm::cont::DataSet& ds,
         const vtkm::cont::CellLocatorGeneral& locator,
         const std::vector<vtkm::Vec3f>& pts,
         const std::string& vField,
         std::vector<vtkm::Vec3f>& output)
{
  /*
  for (std::size_t i = 0; i < pts.size(); i++)
    output.push_back({0, -.1, 0});
  return;
  */

  /*
  vtkm::FloatDefault phiSpacing = vtkm::TwoPi() / (vtkm::FloatDefault)(numPlanes);
  std::cout<<"\n\n********************************************************"<<std::endl;
  std::cout<<"phiSpacing= "<<phiSpacing<<std::endl;
  std::cout<<"NumPlanes= "<<numPlanes<<" spacing= "<<phiSpacing<<std::endl;
  std::cout<<"Plane, Phi"<<std::endl;
  for (int i = 0; i < numPlanes; i++)
    std::cout<<i<<", "<<((float)i * phiSpacing)<<"  deg= "<<(i*phiSpacing*57.92958)<<std::endl;
  std::cout<<"** pts= "<<pts<<std::endl;
  */

  bool isB = vField == "B3D";
  bool isV = vField == "V";
  bool isV2 = vField == "V2";
  bool isX = vField == "X";

  for (const auto& x : pts)
  {
    auto pt = x;
    vtkm::FloatDefault phi = pt[1];

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    GetPlaneIdx(phi, numPlanes, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);

    //std::cout<<pt<<" : "<<phiN<<" Pln: "<<planeIdx0<<" "<<planeIdx1<<" Phi: "<<Phi0<<" "<<Phi1<<std::endl;
    if (isX || isB)
    {
      vtkm::Vec3f ptRZ(pt[0], pt[2], 0);
      std::vector<vtkm::Vec3f> P = {ptRZ};
      std::vector<vtkm::Vec3f> B0 = EvalVector(ds, locator, P, "B2D");
      auto B = B0[0];
      B[1] = B[1] / pt[0];
      if (isB)
      {
        output.push_back(B);
        continue;
      }

      //Calculate X
      vtkm::Vec3f rayPt(pt[0], phiN, pt[2]);
      Ray3f ray0(rayPt, -B), ray1(rayPt, B);

      //std::cout<<"Ray: "<<rayPt<<" "<<B<<std::endl;
      vtkm::Plane<> Plane0({0,Phi0,0}, {0,-1,0}), Plane1({0,Phi1,0}, {0,-1,0});

      vtkm::Vec3f ptOnPlane0, ptOnPlane1;
      vtkm::FloatDefault T0, T1;
      bool tmp;
      Plane0.Intersect(ray0, T0, ptOnPlane0, tmp);
      Plane1.Intersect(ray1, T1, ptOnPlane1, tmp);

      auto dist01 = vtkm::Magnitude(ptOnPlane1-ptOnPlane0);
      auto dist0i = vtkm::Magnitude(pt-ptOnPlane0) / dist01;
      //auto disti1 = vtkm::Magnitude(pt-ptOnPlane1) / dist01;

      //Eval X(p0_rz, p1_rz)
      std::vector<vtkm::Vec3f> P2 = { {ptOnPlane0[0], ptOnPlane0[2], 0}, {ptOnPlane1[0], ptOnPlane1[2], 0} };
      std::vector<int> offsets(2);
      offsets[0] = (int)(planeIdx0 * numNodes);
      offsets[1] = (int)(planeIdx1 * numNodes);

      auto X = EvalVector(ds, locator, P2, vField, offsets);
      //std::cout<<"     Eval X @ "<<P2<<" ---> "<<X<<std::endl;

      auto res = vtkm::Lerp(X[0], X[1], dist0i);
      res[1] /= pt[0];
      res = res+B;
      output.push_back(res);
    }
    //For V
    else
    {
      //Wrap around case....
      if (planeIdx0 == numPlanes-1 && planeIdx1 == 0)
      {
        //std::cout<<"Wrap around: "<<planeIdx0<<" phi= "<<phiN<<std::endl;
      }

      vtkm::Vec3f particleRZ(pt[0], pt[2], 0);
      vtkm::Vec3f B0 = EvalVector(ds, locator, {particleRZ}, "B2D")[0];

      vtkm::FloatDefault PhiMid = Phi0 + (Phi1-Phi0)/2.0;
      vtkm::Plane<> midPlane({0, PhiMid, 0}, {0,1,0});
      Ray3f ray({pt[0], phiN, pt[2]}, B0);

      //Get point on mid plane.  Use the R,Z for this point for triangle finds.
      vtkm::FloatDefault RP_T;
      vtkm::Vec3f ptOnMidPlane;
      bool tmp;
      midPlane.Intersect(ray, RP_T, ptOnMidPlane, tmp);

      //Now, interpolate between Phi_i and Phi_i+1
      vtkm::FloatDefault T01 = (phiN - Phi0) / (Phi1-Phi0);
      vtkm::FloatDefault T10 = 1.0f - T01;

      //Get vec at Phi0 and Phi1.
      vtkm::Vec3f x_ff(ptOnMidPlane[0], ptOnMidPlane[2], 0);
      std::vector<int> offsets(2);
      offsets[0] = planeIdx0*numNodes*2;
      offsets[1] = planeIdx0*numNodes*2 + numNodes;

      vtkm::Vec3f vec_phi0, vec_phi1;
      vtkm::Vec3f vec;
      if (isV)
      {
        auto vecs = EvalVector(ds, locator, {x_ff, x_ff}, vField, offsets);
        vec_phi0 = vecs[0];
        vec_phi1 = vecs[1];
        vec = vec_phi0 * T01 + vec_phi1 * T10;
      }
      else if (isV2)
      {
        const vtkm::FloatDefault basis = 0.0f;

        //gradPsi: pt on mid plane?  (question)
        vtkm::Vec3f gradPsi(B0[2]*x_ff[0], -B0[0]*x_ff[0], 0);
        vtkm::FloatDefault gammaPsi = 1.0f/vtkm::Magnitude(gradPsi);

        vtkm::Vec2f rvec(0,0), zvec(0,0);
        rvec[0] = basis + (1.0-basis) * gammaPsi *   gradPsi[0];
        rvec[1] =         (1.0-basis) * gammaPsi *   gradPsi[1];
        zvec[0] =         (1.0-basis) * gammaPsi * (-gradPsi[1]);
        zvec[1] = basis + (1.0-basis) * gammaPsi *   gradPsi[0];

        //Get the vectors in the ff coordinates.
        auto dAs_ff = EvalVector(ds, locator, {x_ff, x_ff}, "dAs_ff", offsets);
        auto dAs_ff0 = dAs_ff[0];
        auto dAs_ff1 = dAs_ff[1];
        //std::cout<<x_ff<<":  dAs_ff= ::: "<<dAs_ff0<<" "<<dAs_ff1<<" mag: "<<vtkm::Magnitude(dAs_ff0)<<" "<<vtkm::Magnitude(dAs_ff1)<<std::endl;

        vtkm::FloatDefault wphi[2] = {T10, T01}; //{T01, T10};
        vtkm::Vec3f gradAs;
        //std::cout<<"T::: "<<phiN<<" ("<<Phi0<<" "<<Phi1<<") wphi: "<<wphi[0]<<" "<<wphi[1]<<std::endl;

        //vec.r = wphi[0]*( rvec[0]*V.r[0] + zvec[0]*V.z[0]) +
        //        wphi[1]*( rvec[0]*V.r[1] + zvec[0]*v.z[1]);
        //vec.p = wphi[0]*V.phi[0] +
        //        whpi[1]*V.phi[1];
        //vec.z = wphi[0]*( rvec[1]*V.r[0] + zvec[1]*V.z[0]) +
        //        wphi[1]*( rvec[1]*V.r[1] + zvec[1]*V.Z[1]);
        gradAs[0] = wphi[0]*(rvec[0]*dAs_ff0[0] + zvec[0]*dAs_ff0[2]) +
                    wphi[1]*(rvec[0]*dAs_ff1[0] + zvec[0]*dAs_ff1[2]);
        gradAs[2] = wphi[0]*(rvec[1]*dAs_ff0[0] + zvec[1]*dAs_ff0[2]) +
                    wphi[1]*(rvec[1]*dAs_ff1[0] + zvec[1]*dAs_ff1[2]);
        gradAs[1] = wphi[0] * dAs_ff0[1] +
                    wphi[1] * dAs_ff1[1];

        vtkm::FloatDefault BMag = vtkm::Magnitude(B0);
        //project using bfield.
        gradAs[1] = (gradAs[1]*BMag -gradAs[0]*B0[0] - gradAs[2]*B0[2]) / B0[1];

        //deltaB = AsCurl(bhat) + gradAs x bhat.
        std::vector<int> off = {planeIdx0*numNodes};
        vtkm::Vec3f AsCurl_bhat = EvalVector(ds, locator, {x_ff}, "AsCurlBHat", off)[0];

        vtkm::Vec3f curl_bhat = EvalVector(ds, locator, {particleRZ}, "Curl_B2D_Norm")[0];
        auto As_ff = InterpScalar(ds, locator, {x_ff, x_ff}, "As_ff", offsets);
        vtkm::FloatDefault As_ff0 = As_ff[0];
        vtkm::FloatDefault As_ff1 = As_ff[1];

        vtkm::FloatDefault As = wphi[0]*As_ff0 + wphi[1]*As_ff1;
        AsCurl_bhat = As * curl_bhat;

        //std::cout<<"As*Curl(Bhat):: "<<AsCurl_bhat<<" "<<vtkm::Magnitude(AsCurl_bhat)<<std::endl;
        vtkm::Vec3f bhat = EvalVector(ds, locator, {particleRZ}, "B2D_Norm")[0];
        //std::cout<<"bhat= "<<bhat<<" "<<vtkm::Magnitude(bhat)<<std::endl;

        vtkm::Vec3f deltaB = AsCurl_bhat + vtkm::Cross(gradAs, bhat);

        vec = deltaB;
        //std::cout<<"deltaB: "<<deltaB<<" :: "<<vtkm::Magnitude(deltaB)<<std::endl;

        std::cout<<"Meow: pt= "<<pt<<" "<<particleRZ<<" x= "<<x_ff<<" "<<phiN<<" "<<Phi0<<" "<<Phi1<<" wphi=("<<wphi[0]<<", "<<wphi[1]<<")"<<std::endl;
        std::cout<<"           B0= "<<B0<<"  /R= "<<vtkm::Vec3f(B0[0], B0[1]/particleRZ[0], B0[2])<<std::endl;
        std::cout<<"    curl_bhat= "<<curl_bhat<<std::endl;
        std::cout<<"           As= "<<As<<std::endl;
        std::cout<<"Bsa (As*curl_nb)= "<<AsCurl_bhat<<std::endl;
        std::cout<<"          dAs= "<<dAs_ff0<<" "<<dAs_ff1<<" --> "<<gradAs<<std::endl;
        std::cout<<"\n\n"<<std::endl;
      }

      B0[1] = B0[1] / particleRZ[0];
      vec[1] = vec[1] / particleRZ[0];
      auto result = vec + B0;
      output.push_back(result);
    }
  }
}

static int rk4_counter = 0;

std::vector<vtkm::Vec3f>
RK4(const vtkm::cont::DataSet& ds,
    const vtkm::cont::CellLocatorGeneral& locator,
    const std::vector<vtkm::Vec3f>& pts,
    const std::string& vField,
    const std::vector<bool>& /*pointMask*/,
    vtkm::FloatDefault h)
{
  std::vector<vtkm::Vec3f> k1;
  Evaluate(ds, locator, pts, vField, k1);

  std::vector<vtkm::Vec3f> tmp(pts.size());
  vtkm::FloatDefault h_2 = h/2.0;
  for (std::size_t i = 0; i < pts.size(); i++)
    tmp[i] = pts[i] + k1[i]*h_2;
  std::vector<vtkm::Vec3f> k2;
  Evaluate(ds, locator, tmp, vField, k2);

  std::vector<vtkm::Vec3f> k3;
  for (std::size_t i = 0; i < pts.size(); i++)
    tmp[i] = pts[i] + k2[i]*h_2;
  Evaluate(ds, locator, tmp, vField, k3);

  std::vector<vtkm::Vec3f> k4;
  for (std::size_t i = 0; i < pts.size(); i++)
    tmp[i] = pts[i] + k3[i]*h;
  Evaluate(ds, locator, tmp, vField, k4);

  vtkm::FloatDefault h_6 = h/6.0;
  std::vector<vtkm::Vec3f> newPts(pts.size());
  for (std::size_t i = 0; i < pts.size(); i++)
  {
    newPts[i] = pts[i] + h_6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    //auto rkv = h_6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    //std::cout<<"*******************  RK4_"<<rk4_counter<<" : "<<pts[i]<<" === "<<rkv<<" norm: "<<vtkm::Normal(rkv)<<" =======> "<<newPts[i]<<std::endl;

    /*
    //Wrap around.
    if (newPts[i][1] < 0)
      newPts[i][1] += vtkm::TwoPi();
    else if (newPts[i][1] > vtkm::TwoPi())
      newPts[i][1] -= vtkm::TwoPi();
    */
  }

  rk4_counter++;
  return newPts;
}

std::vector<std::vector<vtkm::Vec3f>>
Poincare(const vtkm::cont::DataSet& ds,
         std::vector<vtkm::Vec3f>& pts,
         const std::string& vField,
         vtkm::FloatDefault h,
         int numPunc,
         std::vector<std::vector<vtkm::Vec3f>>* traces=nullptr)

{
//  const vtkm::FloatDefault planeVal = 2.0f;
//  const vtkm::FloatDefault planeVal = vtkm::Pi();
  const vtkm::FloatDefault planeVal = 0; //vtkm::TwoPi();

  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  std::vector<std::vector<vtkm::Vec3f>> punctures(pts.size());
  std::vector<int> puncCount(pts.size(), 0);
  std::vector<bool> pointMask(pts.size(), true);

  if (traces)
    for (int i = 0; i < (int)pts.size(); i++)
      (*traces)[i].push_back(pts[i]);

  auto thetaPsi = ConvertToThetaPsi(pts);
  std::cout<<"Poincare: "<<pts<<"  theta,psi= "<<thetaPsi<<std::endl;

  vtkm::Vec3f p(pts[0][0], pts[0][2], 0);
  //B0: R,phi,Z
  auto B0 = (EvalVector(ds, locator, {p}, "B2D"))[0];
  std::cout<<"B0= "<<B0<<std::endl;
  auto B0_pol = vtkm::Sqrt(B0[0]*B0[0] + B0[2]*B0[2]);
  std::cout<<"B0_pol = "<<B0_pol<<std::endl;

  std::cout<<"B0_phi/R = qB0_pol/r_minor"<<std::endl;
  std::cout<<B0[1]<<" / "<<pts[0][0]<<" = q "<<B0_pol<<" / "<<thetaPsi[0][1]<<std::endl;
  auto x = B0[1] / pts[0][0];
  auto y = B0_pol / thetaPsi[0][1];
  std::cout<<"----> "<<x<<" "<<y<<std::endl;
  std::cout<<std::endl<<std::endl<<std::endl;


  int maxIter = numPunc*1000000;
  maxIter = 15000;
  for (int i = 0; i < maxIter; i++)
  {
    auto newPts = RK4(ds, locator, pts, vField, pointMask, h);

    for (std::size_t j = 0; j < pts.size(); j++)
    {
      //We puncture the plane if the points are on opposite sides of the plane.
      int nRevs0 = vtkm::Floor(vtkm::Abs(pts[j][1] / vtkm::TwoPi()));
      int nRevs1 = vtkm::Floor(vtkm::Abs(newPts[j][1] / vtkm::TwoPi()));
      //std::cout<<" PCHECK: "<<pts[j][1]<<" "<<newPts[j][1]<<" planeVal= "<<planeVal<<" nREVS0= "<<nRevs0<<" "<<nRevs1<<std::endl;

      if (i % 500 == 0) std::cout<<"Poinc iter: "<<i<<" "<<pts[j]<<" nRevs: "<<nRevs0<<" "<<nRevs1<<std::endl;
      //if (newPts[j][1] < 0) pointMask[j] = false;

      if (nRevs1 > nRevs0)
      {
        punctures[j].push_back(newPts[j]);
        puncCount[j]++;
        std::cout<<"PUNC: "<<pts[j]<<" --> "<<newPts[j]<<"  planeVal= "<<planeVal<<" punc= "<<puncCount[j]<<std::endl;
        std::cout<<j<<":    "<<newPts[j]<<" "<<puncCount[j]<<std::endl;
        if (puncCount[j] == numPunc)
          pointMask[j] = false;
      }

      if (traces)
        (*traces)[j].push_back(newPts[j]);
    }

    //All points are done.
    if (std::accumulate(pointMask.begin(), pointMask.end(), 0) == 0)
      break;

    pts = std::move(newPts);
  }

  return punctures;
}

void TestLocator(const vtkm::cont::DataSet& ds)
{
  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  std::vector<vtkm::Vec3f> pts = {{2, 0, 0}, {3,.1,0}, {2.4, -.5, 0}};

  auto points = vtkm::cont::make_ArrayHandle(pts, vtkm::CopyFlag::On);

  vtkm::cont::ArrayHandle<vtkm::Id> cellIds;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> pcoords;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Id, 3>> vertIds;
  vtkm::cont::Invoker invoker;

  invoker(FindCellWorklet{}, points, ds.GetCellSet(), locator, cellIds, pcoords);
  std::cout<<"Points:  "; vtkm::cont::printSummary_ArrayHandle(points, std::cout, true);
  std::cout<<"CellIds: "; vtkm::cont::printSummary_ArrayHandle(cellIds, std::cout, true);
  std::cout<<"Param:   "; vtkm::cont::printSummary_ArrayHandle(pcoords, std::cout, true);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> B;
  ds.GetField("B2D").GetData().AsArrayHandle(B);

  auto cPortal = cellIds.ReadPortal();
  auto pPortal = pcoords.ReadPortal();
  auto bPortal = B.ReadPortal();
  auto cs = ds.GetCellSet().Cast<vtkm::cont::CellSetSingleType<>>();
  for (vtkm::Id i = 0; i < (vtkm::Id)pts.size(); i++)
  {
    vtkm::Id vIds[3];
    vtkm::Id cid = cPortal.Get(i);
    cs.GetCellPointIds(cid, vIds);
    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
    std::cout<<cid<<":: "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<std::endl;
    for (vtkm::Id j = 0; j < 3; j++)
      vals.Append(bPortal.Get(vIds[0]));

    vtkm::Vec3f interp;
    vtkm::exec::CellInterpolate(vals, pPortal.Get(i), vtkm::CellShapeTagTriangle(), interp);
    std::cout<<"    "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<std::endl;
    std::cout<<" -----> "<<interp<<std::endl;
  }
}

void
CalcV(vtkm::cont::DataSet& ds)
{
  std::cout<<__FILE__<<" "<<__LINE__<<" fix me. only works for 1 index. (see below)"<<std::endl;
  //DeltaB = As * curl(Bhat) + grad_As x bhat
  //Advect: DeltaB + B0

  ds.PrintSummary(std::cout);
  vtkm::cont::ArrayHandle<vtkm::Vec3f> AsCurlBhat, gradAs, B0;
  ds.GetField("AsCurlBHat").GetData().AsArrayHandle(AsCurlBhat);
  ds.GetField("gradAs").GetData().AsArrayHandle(gradAs);
  ds.GetField("B2D").GetData().AsArrayHandle(B0);

  std::vector<vtkm::Vec3f> vField((numPlanes * numNodes * 2), vtkm::Vec3f(1,0,0));
  auto B = vField;
  auto Bn = vField;
  auto acb = vField;
  auto gas = vField;

  auto acbPortal = AsCurlBhat.ReadPortal();
  auto gasPortal = gradAs.ReadPortal();
  auto b0Portal = B0.ReadPortal();

  vtkm::Id idx = 0, idx2 = 0;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        auto b = b0Portal.Get(n);
        auto bn = vtkm::Normal(b);

        auto v1 = acbPortal.Get(idx);
        auto v2 = gasPortal.Get(idx);
        auto val = v1 + vtkm::Cross(v2, bn);
        vField[idx] = val;

#if 0
        if (i == 0) //<-------------------- Only do this for ONE plane.
        {
          auto v1 = acbPortal.Get(idx);
          auto v2 = gasPortal.Get(idx);
          //auto cross = vtkm::Cross(v2, bn);
          auto val  = v1 + vtkm::Cross(v2, bn); // + b;

          //R,Phi,Z.
          //vtkm::FloatDefault dPhi = (float)(p+1) * -0.01;
          //val = {0, dPhi, 0};
          //Need to store as R,Z,Phi
          //val = b;
          vField[idx2] = val;
          B[idx2] = b;
          Bn[idx2] = bn;
          acb[idx2] = v1;
          gas[idx2] = v2;

//          if (vtkm::Magnitude(v2) > 100)
//            std::cout<<"************************ "<<val<<" :: "<<v1<<" "<<v2<<" "<<bn<<" "<<b<<std::endl;

          idx2++;
        }
#endif
        idx++;
      }
    }
  }
  std::cout<<"Calc V: "<<idx2<<" "<<vField.size()<<std::endl;

#if 0
  //Now duplicate plane 0 to plane N.
  for (int n = 0; n < numNodes; n++)
  {
    vField.push_back(vField[n]);
    B.push_back(B[n]);
    Bn.push_back(Bn[n]);
    acb.push_back(acb[n]);
    gas.push_back(gas[n]);
  }


  ds3D.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(vField, vtkm::CopyFlag::On)));

  if (0)
  {
    ds3D.AddField(vtkm::cont::make_FieldPoint("B", vtkm::cont::make_ArrayHandle(B, vtkm::CopyFlag::On)));
    ds3D.AddField(vtkm::cont::make_FieldPoint("Bn", vtkm::cont::make_ArrayHandle(Bn, vtkm::CopyFlag::On)));
    ds3D.AddField(vtkm::cont::make_FieldPoint("ACB", vtkm::cont::make_ArrayHandle(acb, vtkm::CopyFlag::On)));
    ds3D.AddField(vtkm::cont::make_FieldPoint("GAS", vtkm::cont::make_ArrayHandle(gas, vtkm::CopyFlag::On)));

    std::cout<<"Dumping file"<<std::endl;
    vtkm::io::VTKDataSetWriter writer("dsWithV.vtk");
    writer.WriteDataSet(ds3D);
  }
#endif

  ds.AddField(vtkm::cont::make_Field("V",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     vField, vtkm::CopyFlag::On));
}

void
CalcX(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> B, Bs, X;

  bool isB2D = ds.HasField("B2D");

  if (isB2D)
    ds.GetField("B2D").GetData().AsArrayHandle(B);
  else
    ds.GetField("B").GetData().AsArrayHandle(B);

  ds.GetField("Bs").GetData().AsArrayHandle(Bs);

  X.Allocate(Bs.GetNumberOfValues());
  auto BP = B.ReadPortal();
  auto BsP = Bs.ReadPortal();
  auto XP = X.WritePortal();

  vtkm::Id idx = 0;
  for (vtkm::Id p = 0; p < numPlanes; p++)
    for (vtkm::Id i = 0; i < numNodes; i++, idx++)
    {
      vtkm::Id bidx = (isB2D ? i : idx);
      XP.Set(idx, BsP.Get(idx)-BP.Get(bidx));
    }

  ds.AddField(vtkm::cont::make_FieldPoint("X", X));
}

void
CalcAsCurlBHat(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> curlBhat;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> As;

  ds.GetField("Curl_B2D_Norm").GetData().AsArrayHandle(curlBhat);
  ds.GetField("As_phi_ff").GetData().AsArrayHandle(As);

  vtkm::Id numAs = As.GetNumberOfValues();

  auto cbPortal = curlBhat.ReadPortal();
  auto asPortal = As.ReadPortal();

  //Dim: (numPlanes, numNodes, idx)
  std::vector<std::vector<std::vector<vtkm::Vec3f>>> As_curlBhat;
  std::vector<std::vector<std::vector<vtkm::FloatDefault>>> As_arr;
  As_curlBhat.resize(numPlanes);
  As_arr.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    As_curlBhat[p].resize(2);
    As_arr[p].resize(2);
    for (int i = 0; i < 2; i++)
    {
      As_curlBhat[p][i].resize(numNodes);
      As_arr[p][i].resize(numNodes);
    }
  }

  //Do some easy indexing.
  vtkm::Id idx = 0;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int n = 0; n < numNodes; n++)
    {
      vtkm::Vec3f cBhat = cbPortal.Get(n);
      for (int i = 0; i < 2; i++)
      {
        auto as = asPortal.Get(idx);
        auto val = as * cBhat;
        //std::cout<<"As: "<<as<<" cBhat= "<<cBhat<<" :: "<<vtkm::Magnitude(cBhat)<<"  val= "<<val<<std::endl;
        As_curlBhat[p][i][n] = val;
        As_arr[p][i][n] = as;
        idx++;
      }
    }
  }

  //flatten to 1d index.
  idx = 0;
  std::vector<vtkm::Vec3f> arr(numAs);
  std::vector<vtkm::FloatDefault> arrAs(numAs);
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::FloatDefault asVal = As_arr[p][i][n];
        arrAs[idx] = asVal;
        vtkm::Vec3f val = As_curlBhat[p][i][n];
        arr[idx] = val;
        idx++;
      }
    }
  }

  ds.AddField(vtkm::cont::make_Field("AsCurlBHat",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr,
                                     vtkm::CopyFlag::On));
  ds.AddField(vtkm::cont::make_Field("As_ff",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arrAs,
                                     vtkm::CopyFlag::On));
}

void
CalcGradAs(vtkm::cont::DataSet& ds)
{
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

          //swap phi and z
          if (c == 1) cc = 2;
          else if (c == 2) cc = 1;

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

  ds.AddField(vtkm::cont::make_Field("dAs_ff",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr_ff,
                                     vtkm::CopyFlag::On));


  /*
  std::cout<<"Print Out dAs_phi_ff"<<std::endl;
  for (int p = 0; p < numPlanes; p++)
    for (int n = 0; n < numNodes; n++)
    {
      auto v0 = dAs_ff[p][0][n];
      auto v1 = dAs_ff[p][1][n];
      if (vtkm::Magnitude(v0) > 1e-7 && vtkm::Magnitude(v1) > 1e-7)
        std::cout<<p<<" :: dAs_phi_ff["<<n<<"]= "<<dAs_ff[p][0][n]<<" "<<dAs_ff[p][1][n]<<std::endl;
    }
  */


  //B is R,Z,Phi
  std::cout<<"****************************************************************   Get R,Z,Phi thing right....."<<std::endl;
  //calculate gradPsi;
  std::vector<vtkm::Vec3f> gradPsi;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B0;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("B2D").GetData().AsArrayHandle(B0);
  gradPsi.resize(numNodes);
  auto cPortal = coords.ReadPortal();
  auto b0Portal = B0.ReadPortal();
  for (int n = 0; n < numNodes; n++)
  {
    auto b = b0Portal.Get(n);
    auto pt = cPortal.Get(n);
    vtkm::Vec3f gv(b[2] * pt[0], -b[0] * pt[0], 0);
    gradPsi[n] = gv;
  }

  std::vector<std::vector<std::vector<vtkm::Vec3f>>> VEC;
  VEC.resize(numPlanes);
  for (int p = 0; p < numPlanes; p++)
  {
    VEC[p].resize(2);
    for (int i = 0; i < 2; i++)
      VEC[p][i].resize(numNodes);
  }

  vtkm::FloatDefault basis = 0.0f;
  std::vector<std::vector<vtkm::FloatDefault>> wphi = {{1.0f, 0.0f}, {0.0f, 1.0f}};

  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::Vec3f gPsi = gradPsi[n];
        vtkm::Vec3f dAs = dAs_ff[p][i][n];
        vtkm::FloatDefault gammaPsi = 1.0f/vtkm::Magnitude(gPsi);

        vtkm::Vec2f rvec(0,0), zvec(0,0);
        rvec[0] = basis + (1.0-basis) * gammaPsi *   gPsi[0];
        rvec[1] =         (1.0-basis) * gammaPsi *   gPsi[1];
        zvec[0] =         (1.0-basis) * gammaPsi * (-gPsi[1]);
        zvec[1] = basis + (1.0-basis) * gammaPsi *   gPsi[0];

        //R,Phi,Z
        //vec.r = vec[0] = rvec[0]*field.R[0] + zvec[0] * field.Z[0]  //plane 0
        //vec.r +=         rvec[0]*field.R[1] + zvec[0] * field.Z[1]  //plane 1
        //vec.z = vec[2] = rvec[1]*field.R[0] + zvec[1] * field.Z[0]  //plane 0
        //vec.z +=         rvec[1]*field.R[1] + zvec[1] * field.Z[1]  //plane 1
        //vec.phi = vec[1] = field.Phi[0]                             //plane 0
        //vec.phi += vec[1] = field.Phi[1]                            //plane 1
        vtkm::Vec3f vec;
        vec[0] = rvec[0] * dAs[0] + zvec[0] * dAs[2];
        vec[2] = rvec[1] * dAs[0] + zvec[1] * dAs[2];
        vec[1] = dAs[1];

        auto B = b0Portal.Get(n);
        auto Bmag = vtkm::Magnitude(B);

        //parallel derivative to phi derivative
        //vec.phi = (vec.phi*Bmag - vec.R*B.r - vec.z*b.z) / B.phi
        vec[1] = (vec[1]*Bmag - vec[0]*B[0] - vec[2]*B[2]) / B[1];
        VEC[p][i][n] = vec;
      }
    }
  }


  //flatten to 1d index
  idx = 0;
  std::vector<vtkm::Vec3f> arr(nVals);
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < 2; i++)
    {
      for (int n = 0; n < numNodes; n++)
      {
        vtkm::Vec3f val = VEC[p][i][n];
        arr[idx] = val;
        idx++;
      }
    }
  }

  ds.AddField(vtkm::cont::make_Field("gradAs",
                                     vtkm::cont::Field::Association::WHOLE_MESH,
                                     arr,
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
           const std::vector<std::vector<vtkm::Vec3f>>& punctures,
           const std::string& tracesNm="./traces.txt",
           const std::string& puncNm="./punctures.txt",
           const std::string& puncThetaPsiNm="./punctures.theta_psi.txt")
{
  bool tExists = Exists(tracesNm);
  bool pExists = Exists(puncNm);
  bool ptpExists = Exists(puncThetaPsiNm);
  std::cout<<"EXISTS: "<<tExists<<" "<<pExists<<" "<<ptpExists<<std::endl;

  //Write headers.
  if (!tExists)
    WriteHeader(tracesNm, "ID, R, Z, T");
  if (!pExists)
    WriteHeader(puncNm, "ID, R, Z, DUMMY");
  if (!ptpExists)
    WriteHeader(puncThetaPsiNm, "ID, THETA, PSI, DUMMY");

  std::ofstream outTraces, outPunc, outPuncPsiTheta;

  outTraces.open(tracesNm, std::ofstream::app);
  outPunc.open(puncNm, std::ofstream::app);
  outPuncPsiTheta.open(puncThetaPsiNm, std::ofstream::app);

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

  //write punctures
  for (int i = 0; i < (int)punctures.size(); i++)
  {
    for (const auto& p : punctures[i])
      outPunc<<i<<", "<<p[0]<<" "<<p[2]<<std::endl;

    auto thetaPsi = ConvertToThetaPsi(punctures[i]);
    for (const auto& p : thetaPsi)
      outPuncPsiTheta<<i<<", "<<p[0]<<", "<<p[1]<<", 0"<<std::endl;
  }
}

vtkm::cont::DataSet
ReadData(std::map<std::string, std::vector<std::string>>& args)
{
  std::map<std::string, std::string> adiosArgs;
  adiosArgs["--dir"] = args["--dir"][0];

  adios = new adios2::ADIOS;
  adiosStuff["mesh"] = new adiosS(adios, "xgc.mesh.bp", "mesh", adiosArgs);
  adiosStuff["data"] = new adiosS(adios, "xgc.3d.bp", "3d", adiosArgs);
  adiosStuff["bfield"] = new adiosS(adios, "xgc.bfield.bp", "/node_data[0]/values", adiosArgs);
  adiosStuff["bfield-all"] = new adiosS(adios, "xgc.bfield-all.bp", "Bs", adiosArgs);

  auto meshStuff = adiosStuff["mesh"];
  auto dataStuff = adiosStuff["data"];
  auto bfieldStuff = adiosStuff["bfield"];
  auto bfield_allStuff = adiosStuff["bfield-all"];
  meshStuff->engine.BeginStep();
  dataStuff->engine.BeginStep();
  bfieldStuff->engine.BeginStep();
  bfield_allStuff->engine.BeginStep();

  dataStuff->engine.Get(dataStuff->io.InquireVariable<int>("nphi"), &numPlanes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_n"), &numNodes, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("n_t"), &numTri, adios2::Mode::Sync);

  //Try and do everything in cylindrical coords and worklets.
  auto ds = ReadMesh(meshStuff);
  ReadScalar(dataStuff, ds, "dpot");
  ReadScalar(dataStuff, ds, "apars", "apars", true);
  ReadVec(bfieldStuff, ds, "B", "/node_data[0]/values");
  ReadVec(bfieldStuff, ds, "B3D", "/node_data[0]/values", true, false);
  ReadVec(bfield_allStuff, ds, "Bs", "Bs", true);
  ReadOther(bfieldStuff, ds, "As_phi_ff");
  ReadOther(bfieldStuff, ds, "dAs_phi_ff");
  CalcX(ds);

  CalcAsCurlBHat(ds);
  CalcGradAs(ds);
  CalcV(ds);

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

  ds.PrintSummary(std::cout);
//  vtkm::io::VTKDataSetWriter writer("debug.vtk");
//  writer.WriteDataSet(ds);

  return ds;
}

void
GradientTest()
{
  std::vector<vtkm::Vec3f> coords;

  coords.push_back(vtkm::Vec3f(0,0,0));
  coords.push_back(vtkm::Vec3f(1,0,0));
  coords.push_back(vtkm::Vec3f(0,1,0));
  coords.push_back(vtkm::Vec3f(1,1,0));

  std::vector<vtkm::Id> conn;

  //tri 0
  conn.push_back(0);
  conn.push_back(2);
  conn.push_back(1);

  //tri1
  conn.push_back(1);
  conn.push_back(2);
  conn.push_back(3);

  vtkm::cont::DataSetBuilderExplicit dsb;
  auto ds = dsb.Create(coords, vtkm::CellShapeTagTriangle(), 3, conn);

  std::vector<vtkm::Vec3f> vecs;
  vecs.push_back(vtkm::Vec3f(0,0,2));
  vecs.push_back(vtkm::Vec3f(1,0,1));
  vecs.push_back(vtkm::Vec3f(0,1,1));
  vecs.push_back(vtkm::Vec3f(1,1,1));

  ds.AddField(vtkm::cont::make_FieldPoint("V", vtkm::cont::make_ArrayHandle(vecs, vtkm::CopyFlag::On)));
  ds.PrintSummary(std::cout);

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetActiveField("V");
  gradient.SetOutputFieldName("gradV");
  auto out = gradient.Execute(ds);
  out.PrintSummary(std::cout);

}

void
Debug(const vtkm::cont::DataSet& inDS)
{

  auto ds = inDS;

  /*
  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  gradient.SetComputeVorticity(true);
  gradient.SetActiveField("B_RZP_2D");
  auto ds = gradient.Execute(inDS);
  */

  vtkm::cont::CellLocatorGeneral locator;
  locator.SetCellSet(ds.GetCellSet());
  locator.SetCoordinates(ds.GetCoordinateSystem());
  locator.Update();

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  auto cPortal = coords.ReadPortal();

  vtkm::Vec3f pt(3.029365, 0.020600, 0);

  //Use node position.
  //pt = cPortal.Get(1521);


  vtkm::Vec3f B0 = EvalVector(ds, locator, {pt}, "B_RZP_2D")[0];
  vtkm::Vec3f curlBh = EvalVector(ds, locator, {pt}, "Curl_B2D_Norm")[0];
  std::cout<<"B0_rzp("<<pt<<")= "<<B0<<std::endl;
  std::cout<<"curlBhat("<<pt<<")= "<<curlBh<<std::endl;
  std::cout<<"**************** DO Z,Phi need to swap????? ***************"<<std::endl;

  vtkm::filter::Gradient gradient;
  gradient.SetComputePointGradient(true);
  //gradient.SetActiveField("BRZP_Norm");
  gradient.SetActiveField("B_RZP_2D");
  auto tmpDS = gradient.Execute(ds);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f, 3>> gradients;
  tmpDS.GetField("Gradients").GetData().AsArrayHandle(gradients);
  auto GRAD = EvalTensor(tmpDS, locator, {pt}, "Gradients", {0})[0];
/*
  std::cout<<"dBr/dr= "<<GRAD[0][0]<<std::endl;
  std::cout<<"dBz/dr= "<<GRAD[0][1]<<std::endl;
  std::cout<<"dBp/dr= "<<GRAD[0][2]<<std::endl;

  std::cout<<"dBr/dz= "<<GRAD[1][0]<<std::endl;
  std::cout<<"dBz/dz= "<<GRAD[1][1]<<std::endl;
  std::cout<<"dBp/dz= "<<GRAD[1][2]<<std::endl;

  std::cout<<"dBr/dp= "<<GRAD[2][0]<<std::endl;
  std::cout<<"dBz/dp= "<<GRAD[2][1]<<std::endl;
  std::cout<<"dBp/dp= "<<GRAD[2][2]<<std::endl;
*/

  auto R = pt[0];
  auto Z = pt[1];
  auto inv_r = 1.0/R;
  auto Bmag = vtkm::Magnitude(B0);
  auto over_Bmag = 1.0/Bmag;
  auto over_Bmag2 = over_Bmag * over_Bmag;

  auto br = B0[0];
  auto bz = B0[1];
  auto bphi = B0[2];

  auto dbrdr = GRAD[0][0];
  auto dbzdr = GRAD[0][1];
  auto dbpdr = GRAD[0][2];

  auto dbrdz = GRAD[1][0];
  auto dbzdz = GRAD[1][1];
  auto dbpdz = GRAD[1][2];

  auto dbrdp = GRAD[2][0];
  auto dbzdp = GRAD[2][1];
  auto dbpdp = GRAD[2][2];

  std::cout<<"dbrdr= "<<dbrdr<<std::endl;
  std::cout<<"dbrdz= "<<dbrdz<<std::endl;
  std::cout<<"dbrdp= "<<dbrdp<<std::endl;

  std::cout<<"dbzdr= "<<dbzdr<<std::endl;
  std::cout<<"dbzdz= "<<dbzdz<<std::endl;
  std::cout<<"dbzdp= "<<dbzdp<<std::endl;

  std::cout<<"dbpdr= "<<dbpdr<<std::endl;
  std::cout<<"dbpdz= "<<dbpdz<<std::endl;
  std::cout<<"dbpdp= "<<dbpdp<<std::endl;


  //dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
  //dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
  //dbdphi=0D0  ! no B perturbation
  auto dbdr = (br*dbrdr + bphi*dbpdr + bz*dbzdr) * over_Bmag;
  auto dbdz = (br*dbrdz + bphi*dbpdz + bz*dbzdz) * over_Bmag;
  auto dbdphi = 0;

  auto div = dbrdr + br/R + dbzdz;
  std::cout<<"Check divervgence: "<<div<<std::endl;

  vtkm::Vec3f curl_B;
  //R curl_B(1)  = fld%dbzdp*inv_r - fld%dbpdz
  curl_B[0] =          dbzdp*inv_r - dbpdz;
  //Z curl_B(2)  = fld%bphi*inv_r + fld%dbpdr - fld%dbrdp*inv_r
  curl_B[1] =          bphi*inv_r +     dbpdr -     dbrdp*inv_r;
  std::cout<<"    curl_b.z: "<<(bphi*inv_r)<<" + "<<dbpdr<<" - "<<(dbrdp*inv_r)<<std::endl;
  //phi curl_B(3)  = fld%dbrdz - fld%dbzdr
  curl_B[2] =            dbrdz -     dbzdr;

  std::cout<<"curl_B= "<<curl_B<<std::endl;
  std::cout<<"  dbdr= "<<dbdr<<std::endl;
  std::cout<<"  dbdz= "<<dbdz<<std::endl;
  std::cout<<"  dbdphi= "<<dbdphi<<std::endl;


  vtkm::Vec3f curl_nb;

  //R,Z,Phi
  //curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
  //curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
  //curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2

  curl_nb[0] = curl_B[0]*over_Bmag + (bphi * dbdz) * over_Bmag2;
  curl_nb[1] = curl_B[1]*over_Bmag + (-bphi * dbdr) * over_Bmag2;
  curl_nb[2] = curl_B[2]*over_Bmag + (bz * dbdr - br * dbdz) * over_Bmag2;

  std::cout<<"curl_nb= "<<curl_nb<<std::endl;

#if 0
  vtkm::Id offset = 2673266;
  std::vector<vtkm::Id> vids = {1521, 1612, 1613};

  vtkm::cont::ArrayHandle<vtkm::Vec3f> b, coords, curl_bhat;
  ds.GetField("B_RZP_2D").GetData().AsArrayHandle(b);
  ds.GetField("Curl_B2D_Norm").GetData().AsArrayHandle(curl_bhat);
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  std::vector<vtkm::Vec3f> B_vals, curl_nb;
  for (const auto& id : vids)
  {
    B_vals.push_back(b.ReadPortal().Get(id));
    curl_nb.push_back(curl_bhat.ReadPortal().Get(id));
  }

  std::vector<vtkm::Vec<vtkm::Vec3f, 3>> grad;
  std::vector<vtkm::Vec3f> pts;
  for (const auto& id : vids)
  {
    grad.push_back(gradients.ReadPortal().Get(id));
    pts.push_back(coords.ReadPortal().Get(id));
  }

  for (int i = 0; i < 3; i++)
  {
    auto c = pts[i];
    auto g = grad[i];

    auto dBr_dR     = g[0][0];
    auto dBz_dR     = g[0][1];
    auto dBphi_dR   = g[0][2];

    auto dBr_dZ     = g[1][0];
    auto dBz_dZ     = g[1][1];
    auto dBphi_dZ   = g[1][2];

    auto dBr_dPhi   = g[2][0];
    auto dBz_dPhi   = g[2][1];
    auto dBphi_dPhi = g[2][2];

    std::cout<<i<<": "<<c<<std::endl;
    std::cout<<"    B= "<<B_vals[i]<<std::endl;
    std::cout<<"    Curl_bhat= "<<curl_nb[i]<<std::endl;
    std::cout<<"    dB*_dR= "<<g[0]<<std::endl;
    std::cout<<"    dB*_dZ= "<<g[1]<<std::endl;
    std::cout<<"    dB*_dP= "<<g[2]<<std::endl;

    std::cout<<"   dBp_dr= "<<dBphi_dR<<std::endl;
    std::cout<<"   dBr_dp= "<<dBr_dPhi<<std::endl;
    std::cout<<"   dBp_dz= "<<dBphi_dZ<<std::endl;
    std::cout<<"   dBz_dr= "<<dBz_dR<<std::endl;
    std::cout<<"   dBr_dz= "<<dBr_dZ<<std::endl;
    std::cout<<"\n\n";
  }
#endif
}

void
SaveStuff(const vtkm::cont::DataSet& inDS)
{
  vtkm::cont::DataSet ds;
  ds.AddCoordinateSystem(inDS.GetCoordinateSystem());
  ds.SetCellSet(inDS.GetCellSet());

  ds.AddField(inDS.GetField("B2D_Norm"));
  ds.AddField(inDS.GetField("Curl_B2D_Norm"));
  ds.AddField(inDS.GetField("B2D"));

  //Add gradPsi
  vtkm::Id nPts = ds.GetNumberOfPoints();
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, B0;
  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  inDS.GetField("B2D").GetData().AsArrayHandle(B0);
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

int
main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  auto config = vtkm::cont::Initialize(argc, argv, opts);

  std::map<std::string, std::vector<std::string>> args;
  ParseArgs(argc, argv, args);

  if (argc < 7)
  {
    std::cerr<<"Usage: "<<argv[0]<<" dataFile numSeeds maxPunctures stepSize poincVar [options]"<<std::endl;
    std::cerr<<config.Usage<<std::endl;
    return -1;
  }

  auto ds = ReadData(args);

  SaveStuff(ds);
  //Debug(ds);
  //return 0;

  vtkm::FloatDefault stepSize = std::atof(args["--stepSize"][0].c_str());
  int numPunc = std::atoi(args["--numPunc"][0].c_str());
  std::string vField = args["--vField"][0];

  std::vector<vtkm::Vec3f> seeds;
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
  }
  else if (args.find("--seed") != args.end())
  {
    auto vals = args["--seed"];
    vtkm::FloatDefault r = std::atof(vals[0].c_str());
    vtkm::FloatDefault z = std::atof(vals[1].c_str());
    vtkm::FloatDefault t = std::atof(vals[2].c_str());
    seeds.push_back({r, t, z});
  }
  else
  {
    //first point in traces.v2
    seeds = {{3.029365, 6.183185, 0.020600}};

    //traces.v2 pt near begining.
    //seeds = {{3.024768, 6.070249, 0.049700}};

    //seeds from data/sku_8000/jong.py
    /*
    seeds = {
      {3.351443, 0.0, -0.451649},
      {3.187329, 0.0, -0.665018},
      {1.992020, 0.0, -0.126203},
      {3.018666, 0.0, 0.073864},
      {3.176583, 0.0, -0.220557},
      {2.179227, 0.0, 0.291539},
    };
    */
  }



  std::vector<std::vector<vtkm::Vec3f>> traces(seeds.size());
  auto punctures = Poincare(ds, seeds, vField, stepSize, numPunc, &traces);

  SaveOutput(traces, punctures);

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


#if 0
  float phi = -0.01;
  int cnt = 0;
  numPlanes = 48;
  float dPhi = vtkm::TwoPi()/static_cast<double>(numPlanes);
  int idx = 0;
  std::cout<<"dPhi= "<<dPhi<<std::endl;
  int rev0 = 0;

  vtkm::FloatDefault phiSpacing = vtkm::TwoPi() / (vtkm::FloatDefault)(numPlanes);
  std::cout<<"\n\n********************************************************"<<std::endl;
  std::cout<<"phiSpacing= "<<phiSpacing<<std::endl;
  std::cout<<"NumPlanes= "<<numPlanes<<" spacing= "<<phiSpacing<<std::endl;
  std::cout<<"Plane, Phi"<<std::endl;
  for (int i = 0; i < numPlanes; i++)
    std::cout<<i<<", "<<((float)i * phiSpacing)<<"  deg= "<<(i*phiSpacing*57.2958)<<std::endl;
  std::cout<<std::endl<<std::endl;

  while (cnt < 10)
  {

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault T, phi0, phi1;
    GetPlaneIdx(phi, numPlanes, planeIdx0, planeIdx1, phi0, phi1, numRevs, T);
/*
    int  numRevs = vtkm::Floor(vtkm::Abs(phi / vtkm::TwoPi()));
    vtkm::FloatDefault rem = std::fmod(phi, vtkm::TwoPi());

    vtkm::Id planeIdx = static_cast<vtkm::Id>(vtkm::Floor(rem / dPhi));
    vtkm::Id planeIdx0, planeIdx1;
    if (planeIdx < 0)
      planeIdx += numPlanes;
    if (planeIdx == numPlanes-1)
    {
      planeIdx0 = 0;
      planeIdx1 = planeIdx;
    }
    else
    {
      planeIdx0 = planeIdx;
      planeIdx1 = planeIdx0+1; //no //B is going in the NEGATIVE phi direction.
    }
*/
    std::cout<<idx<<":  phi= "<<phi<<" #Rev= "<<numRevs<<" planes: ("<<planeIdx0<<" "<<planeIdx1<<") T= "<<T<<std::endl;
    std::cout<<"*******************************************\n\n\n"<<std::endl;

    if (numRevs > rev0)
    {
      rev0 = numRevs;
      cnt++;
    }

    phi -= 0.05;
    idx++;
  }
  return 0;
#endif
