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

//rendering
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/rendering/View2D.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/MapperVolume.h>
#include <vtkm/rendering/MapperPoint.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/cont/testing/MakeTestDataSet.h>
#include <vtkm/rendering/testing/RenderTest.h>
#include <vtkm/cont/ColorTable.h>

#include <fides/DataSetReader.h>

#include <adios2.h>
#include <random>
#include <chrono>
#include <variant>
#include <mpi.h>

adios2::ADIOS *adios = NULL;
class adiosS;
std::map<std::string, adiosS*> adiosStuff;

int numPlanes = -1;
int numNodes = -1;
int numTri = -1;

const double psi_x = 0.2661956235889000;
const double xpoint_r = 1.557545038402000;
const double xpoint_z = -1.177067412978000;
const int nWedge = 1;


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

  int i = 1;
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


class RenderScene
{
public:
    RenderScene(const std::string &imageName)
        : ImageName(imageName)
    {
        std::cout<<"Scene: "<<this->ImageName<<std::endl;
    }

    int Size() const { return this->DataSets.size(); }

    void Add(const RenderScene &rsc)
    {
        for (int i = 0; i < rsc.Size(); i++)
            Add(rsc.DataSets[i], rsc.ScalarRanges[i], rsc.VarNames[i], rsc.ColorTables[i], rsc.RenderTypes[i]);
    }

    void Add(const vtkm::cont::DataSet &ds,
             const vtkm::Range &range,
             const std::string &vname,
             const std::string &ct,
             const std::string &rt="RayTrace")

    {
        this->DataSets.push_back(ds);
        this->ScalarRanges.push_back(range);
        this->VarNames.push_back(vname);
        this->ColorTables.push_back(ct);
        this->RenderTypes.push_back(rt);
    }

    std::vector<vtkm::cont::DataSet> DataSets;
    std::vector<std::string> VarNames;
    std::vector<std::string> ColorTables;
    std::vector<vtkm::Range> ScalarRanges;
    std::vector<std::string> RenderTypes;
    std::string ImageName;
};

class Renderer
{
public:
    Renderer(int _sz) :
        imageSize(_sz),
        bg(0.0f, 0.0f, 0.0f, 1.0f),
        fg(1.0f, 1.0f, 1.0f, 1.0f),
        canvas(imageSize, imageSize)
    {
        std::vector<double> rgb;
        rgb.push_back(0); rgb.push_back(0.5); rgb.push_back(0); rgb.push_back(.5);
        rgb.push_back(.2); rgb.push_back(0); rgb.push_back(.2); rgb.push_back(1);
        rgb.push_back(.4); rgb.push_back(0); rgb.push_back(0.5); rgb.push_back(.7);
        rgb.push_back(.5); rgb.push_back(0.9); rgb.push_back(0.9); rgb.push_back(0.9);
        rgb.push_back(.8); rgb.push_back(1); rgb.push_back(0.2); rgb.push_back(0);
        rgb.push_back(1); rgb.push_back(1); rgb.push_back(1); rgb.push_back(0);

        std::vector<double> alpha = {.25, .25, .25, .25, .25, .25};

        auto xgcCT = vtkm::cont::ColorTable("xgc-ct",
                                            vtkm::cont::ColorSpace::RGB,
                                            {1,1,0}, rgb);
        this->colorTables["xgc-ct"] = xgcCT;

        auto xgcCTAlpha = vtkm::cont::ColorTable("xgc-ct",
                                            vtkm::cont::ColorSpace::RGB,
                                            {1,1,0}, rgb);

        xgcCTAlpha.AddPointAlpha(0.0, .0001);
        xgcCTAlpha.AddPointAlpha(1.0, .0001);
        this->colorTables["xgc-ct-alpha"] = xgcCTAlpha;

        auto xgcCTVolume = vtkm::cont::ColorTable("xgc-ct-volume",
                                                  vtkm::cont::ColorSpace::RGB,
                                                  {1,1,0}, rgb);

        this->colorTables["xgc-ct-volume"] = xgcCTVolume;
        xgcCTVolume.AddPointAlpha(0.0, .00001);
        xgcCTVolume.AddPointAlpha(1.0, .00001);

        /*
        xgcCTVolume.AddPointAlpha(0.0, 0); //.75);
        xgcCTVolume.AddPointAlpha(0.1, 0); //.050);
        xgcCTVolume.AddPointAlpha(0.15,0); // .0);
        xgcCTVolume.AddPointAlpha(0.5, 0);
        xgcCTVolume.AddPointAlpha(0.85, .0);
        xgcCTVolume.AddPointAlpha(0.9, .050);
        xgcCTVolume.AddPointAlpha(0.95, .50);
        xgcCTVolume.AddPointAlpha(1.0, 1.0);
        */

        vtkm::rendering::CanvasRayTracer canvas(imageSize, imageSize);
    }

    void RenderMultiMapper(const std::vector<vtkm::cont::DataSet> &dataSets,
                           const std::vector<std::string> &fieldNames,
                           const std::vector<vtkm::Range> &scalarRanges,
                           const std::vector<std::string> &ctTables,
                           const std::vector<std::string> &renderTypes,
                           const std::string &viewType,
                           const std::string &fname)
    {
        vtkm::Bounds bbox;
        for (std::size_t i = 0; i < dataSets.size(); i++)
        {
            auto ds = dataSets[i];
            bbox.Include(ds.GetCoordinateSystem().GetBounds());
        }
        this->SetCamera(bbox, viewType);


        for (std::size_t i = 0; i < dataSets.size(); i++)
        {
            auto ds = dataSets[i];

            if (renderTypes[i] == "RayTrace")
            {
                vtkm::rendering::MapperRayTracer rtMapper;
                rtMapper.SetCanvas(&this->canvas);
                rtMapper.SetActiveColorTable(this->colorTables[ctTables[i]]);
                if (i == 0) rtMapper.SetCompositeBackground(false);

                rtMapper.RenderCells(ds.GetCellSet(),
                                     ds.GetCoordinateSystem(),
                                     ds.GetField(fieldNames[i]),
                                     this->colorTables[ctTables[i]],
                                     this->camera,
                                     scalarRanges[i]);
            }
            else if (renderTypes[i] == "Volume")
            {
                vtkm::rendering::MapperVolume volMapper;
                volMapper.SetCanvas(&this->canvas);
                volMapper.SetActiveColorTable(this->colorTables[ctTables[i]]);
                if (i == 0) volMapper.SetCompositeBackground(false);
                //volMapper.SetSampleDistance(0.001);
                volMapper.RenderCells(ds.GetCellSet(),
                                      ds.GetCoordinateSystem(),
                                      ds.GetField(fieldNames[i]),
                                      this->colorTables[ctTables[i]],
                                      this->camera,
                                      scalarRanges[i]);
            }
            else if (renderTypes[i] == "Point")
            {
              vtkm::rendering::MapperPoint ptMapper;
              ptMapper.SetCanvas(&this->canvas);
              ptMapper.SetRadiusDelta(4.0f);
              ptMapper.UseVariableRadius(false);
              ptMapper.SetActiveColorTable(this->colorTables[ctTables[i]]);
              if (i == 0) ptMapper.SetCompositeBackground(false);
              ptMapper.RenderCells(ds.GetCellSet(),
                                    ds.GetCoordinateSystem(),
                                    ds.GetField(fieldNames[i]),
                                    this->colorTables[ctTables[i]],
                                    this->camera,
                                    scalarRanges[i]);
            }
        }
        this->canvas.SaveAs(fname);
    }

    void Render(const std::vector<vtkm::cont::DataSet> &dataSets,
                const std::string &fieldName,
                const vtkm::Range &scalarRange,
                const std::string &ctTable,
                const std::string &renderType,
                const std::string &viewType,
                const std::string &fname)
  {
    std::vector<std::string> fn, ct, rt;
    std::vector<vtkm::Range> sr;
    for (const auto& ds: dataSets)
    {
      fn.push_back(fieldName);
      ct.push_back(ctTable);
      rt.push_back(renderType);
      sr.push_back(scalarRange);
    }

    this->Render(dataSets, fn, sr, ct, rt, viewType, fname);
  }

    void Render(const std::vector<vtkm::cont::DataSet> &dataSets,
                const std::vector<std::string> &fieldNames,
                const std::vector<vtkm::Range> &scalarRanges,
                const std::vector<std::string> &ctTables,
                const std::vector<std::string> &renderTypes,
                const std::string &viewType,
                const std::string &fname)
    {
        if (dataSets.empty() || dataSets.size() != fieldNames.size() ||
            dataSets.size() != scalarRanges.size() || dataSets.size() != ctTables.size())
        {
            std::cout<<"vectors sizes for Init don't match"<<std::endl;
            return;
        }
        if (std::find(renderTypes.begin(), renderTypes.end(), "Volume") != renderTypes.end())
            return RenderMultiMapper(dataSets, fieldNames, scalarRanges,
                                     ctTables, renderTypes, viewType, fname);


        vtkm::rendering::Scene theScene;
        vtkm::Bounds bbox;
        for (std::size_t i = 0; i < dataSets.size(); i++)
        {
            auto ds = dataSets[i];
            auto actor = vtkm::rendering::Actor(ds.GetCellSet(),
                                                ds.GetCoordinateSystem(),
                                                ds.GetField(fieldNames[i]),
                                                this->colorTables[ctTables[i]]);
            actor.SetScalarRange(scalarRanges[i]);
            theScene.AddActor(actor);
            bbox.Include(ds.GetCoordinateSystem().GetBounds());
        }

        this->SetCamera(bbox, viewType);

        vtkm::rendering::View3D view(theScene, mapper, canvas, camera, fg, bg);
        //vtkm::rendering::View2D view(theScene, mapper, canvas, camera, fg, bg);
        view.SetWorldAnnotationsEnabled(false);
//        view.Initialize();
        view.Paint();
        view.SaveAs(fname);
    }

    void Render(const std::string &fname)
    {
        vtkm::rendering::View3D view(scene, mapper, canvas, camera, fg, bg);

        view.SetWorldAnnotationsEnabled(false);
//        view.Initialize();
        view.Paint();
        view.SaveAs(fname);
    }

    void SetCamera(const vtkm::Bounds &coordBounds, const std::string &viewType)
    {
        vtkm::Bounds b = coordBounds;
        b.Z.Min = 0;
        b.Z.Max = 4;
        camera.ResetToBounds(b);

        vtkm::Vec3f lookAt, pos, viewUp;
        vtkm::FloatDefault FOV, dolly;

        if (viewType == "" || viewType == "D3D" || viewType == "d3d")
        {
            lookAt = vtkm::Vec3f(0,.5,0); //orig
            lookAt = vtkm::Vec3f(0,.2,0);
            pos = vtkm::Vec3f(8,6, 2);
            viewUp = vtkm::Vec3f(0,0,1);
            FOV = static_cast<vtkm::FloatDefault>(30.0);
            dolly = static_cast<vtkm::FloatDefault>(1.2); //orig
            dolly = static_cast<vtkm::FloatDefault>(1.15);

            //view for the back chopping..
            pos = vtkm::Vec3f(-10, 3, 2.5);
            pos = vtkm::Vec3f(-10, 2.5, 2.5);
            lookAt = vtkm::Vec3f(0,-.1,0);
            lookAt = vtkm::Vec3f(0,0,0);
            dolly = static_cast<vtkm::FloatDefault>(1.17);
            dolly = static_cast<vtkm::FloatDefault>(2.0);
            dolly = static_cast<vtkm::FloatDefault>(2.5);
        }
        else if (viewType == "circular")
        {
            lookAt = vtkm::Vec3f(0,1.75,0);
            pos = vtkm::Vec3f(8.5,5.1, 1);
            pos = vtkm::Vec3f(5,3,0);
            viewUp = vtkm::Vec3f(0,0,1);
            FOV = static_cast<vtkm::FloatDefault>(30.0);
            dolly = static_cast<vtkm::FloatDefault>(1.8);

            lookAt = vtkm::Vec3f(0,0,0);
            pos = vtkm::Vec3f(-7.5, -6.0, 1.7);
            viewUp = vtkm::Vec3f(0,0,1);
            dolly = static_cast<vtkm::FloatDefault>(2.0);
        }
        else if (viewType == "poinc")
        {
          camera.SetModeTo2D();
          camera.SetViewport(0.f, 10.f, 0.f, 10.f);
          camera.Zoom(-2.0f);
          lookAt = vtkm::Vec3f(3,4,0);
          pos = vtkm::Vec3f(3,4,1);
          viewUp = vtkm::Vec3f(0,1,0);
          FOV = 30;
          dolly = 2;
        }
/*
        camera.SetLookAt(lookAt);
        camera.SetFieldOfView(FOV);
        camera.SetPosition(pos);
        camera.SetViewUp(viewUp);
        camera.Dolly(dolly);
*/
    }

    int imageSize;
    vtkm::rendering::Color bg, fg;
    vtkm::rendering::Scene scene;
    vtkm::rendering::Camera camera;
    vtkm::rendering::MapperRayTracer mapper;
    vtkm::rendering::CanvasRayTracer canvas;
    std::map<std::string, vtkm::cont::ColorTable> colorTables;
};



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

enum TransformType {CAR_TO_CYL, CYL_TO_CAR};


static vtkm::Vec3f
TransformCarToCyl(const vtkm::Vec3f& v)
{
  auto x = v[0];
  auto y = v[1];
  auto z = v[2];
  auto r = vtkm::Sqrt(x*x + y*y);
  //VTK way...
  //auto phi = vtkm::Pi() + vtkm::ATan2(-y, -x);
  auto phi = vtkm::ATan2(y,x);
  if (phi < 0)
    phi += vtkm::TwoPi();

  return vtkm::Vec3f(r,phi,z);
}

static vtkm::Vec3f
TransformCylToCar(const vtkm::Vec3f& v)
{
  auto r = v[0];
  auto t = v[1];
  auto z = v[2];

  return vtkm::Vec3f(r*vtkm::Cos(t),
                     r*vtkm::Sin(t),
                     z);
}

static vtkm::Vec3f
TransformVec(const vtkm::Vec3f& pt,
             const vtkm::Vec3f& vec,
             const TransformType& transType)
{
  const vtkm::FloatDefault eps = 1.e-5;
  const vtkm::FloatDefault epsInv = 1.e+5;

//  std::cout<<"XFORM: p= "<<pt<<" v= "<<vec<<" "<<transType<<std::endl;

  vtkm::Vec3f pt1;
  if (transType == CAR_TO_CYL)
    pt1 = TransformCarToCyl(pt);
  else
    pt1 = TransformCylToCar(pt);
//  std::cout<<"XPT: pt "<<pt<<" ----> "<<pt1<<std::endl;

  auto mag = vtkm::Sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  vtkm::Vec3f tmp(vec[0]/mag, vec[1]/mag, vec[2]/mag);
//  std::cout<<mag<<" tmp= |v|= "<<tmp<<std::endl;

  for (int k = 0; k < 3; k++)
    tmp[k] = pt[k] + (tmp[k]*eps);
//  std::cout<<"pt+|v|*eps: tmp= "<<tmp<<std::endl;

  vtkm::Vec3f v1;
  if (transType == CAR_TO_CYL)
  {
    v1 = TransformCarToCyl(tmp);
//    std::cout<<"  XF_car_cyl "<<tmp<<" --> "<<v1<<std::endl;
    auto dPhi = v1[1] - pt[1];
//    std::cout<<"******* dPhi= "<<dPhi<<std::endl;
    //This was the wee fix.
//    if (v1[1] > vtkm::Pi())
//    {
//      std::cout<<"     "<<v1<<std::endl;
//      v1[1] -= vtkm::TwoPi();
//    }
  }
  else
    v1 = TransformCylToCar(tmp);

//  std::cout<<"v1 "<<v1<<std::endl;
//  std::cout<<"pt1 "<<pt1<<std::endl;
//  std::cout<<"v1-pt1 "<<(v1-pt1)<<std::endl;

  auto _vt = v1-pt1;
//  std::cout<<"0 _vt "<<_vt<<std::endl;
  if (_vt[1] > 1)
    _vt[1] -= vtkm::TwoPi();
  else if (_vt[1] < -1)
    _vt[1] += vtkm::TwoPi();
//  std::cout<<"1 _vt "<<_vt<<std::endl;

  auto vt = _vt * (epsInv * mag);
/*
  vtkm::Vec3f vt((v1[0] - pt1[0]) * epsInv * mag,
                 (v1[1] - pt1[1]) * epsInv * mag,
                 (v1[2] - pt1[2]) * epsInv * mag);
*/
//  std::cout<<"vt "<<vt<<std::endl<<std::endl<<std::endl;

  return vt;
}


static vtkm::Vec3f Cyl2CarVec(const vtkm::Vec3f& pCyl,
                              const vtkm::Vec3f& vCyl)
{
  const vtkm::FloatDefault eps = 1.e-5;
  const vtkm::FloatDefault epsInv = 1.e+5;

  auto r = pCyl[0];
  auto t = pCyl[1];
  auto z = pCyl[2];

  vtkm::Vec3f pCar(r*vtkm::Cos(t),
                   r*vtkm::Sin(t),
                   z);
  auto mag = vtkm::Sqrt(r*r + t*t + z*z);
  vtkm::Vec3f tmp(vCyl[0]/mag, vCyl[1]/mag, vCyl[2]/mag);

  for (int k = 0; k < 3; k++) tmp[k] = pCyl[k] + tmp[k]*eps;

  vtkm::Vec3f vCar(tmp[0]*vtkm::Cos(tmp[1]),
                   tmp[0]*vtkm::Sin(tmp[1]),
                   tmp[2]);

  for (int k = 0; k < 3; k++) vCar[k] = (vCar[k] - pCar[k]) * mag * epsInv;

  return vCar;
}


void
CalcBFieldXYZ(vtkm::cont::DataSet& ds)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> b;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> A_s;

  ds.GetField("BXYZ").GetData().AsArrayHandle(b);
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
//  auto vec = grad.ReadPortal().Get(0);

  std::vector<vtkm::Vec3f> V(n);
  for (vtkm::Id i = 0; i < n; i++)
  {
    vtkm::Vec3f val = bPortal.Get(i);
    auto g = gPortal.Get(i);

    //X: dAz/dy - dAy/dz
    //Y: dAx/dz - dAz/dx
    //Z: dAy/dx - dAx/dy
    vtkm::FloatDefault xv, yv, zv;
    xv = g[2][1] - g[1][2];
    yv = g[0][2] - g[2][0];
    zv = g[1][0] - g[0][1];

    val[0] += xv;
    val[1] += yv;
    val[2] += zv;

    V[i] = val;
  }

  ds.AddField(vtkm::cont::make_FieldPoint("V_XYZ", vtkm::cont::make_ArrayHandle(V, vtkm::CopyFlag::On)));
}

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
ReadScalar(adiosS* stuff,
           vtkm::cont::DataSet& ds,
           const std::string& vname,
           bool /*isXYZ*/,
           bool cylAdd,
           std::string fileName="")
{
  if (fileName.empty())
    fileName = vname;

  auto v = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(v, tmp, adios2::Mode::Sync);

  if (cylAdd)
  {
    for (int i = 0; i < numNodes; i++)
      tmp.push_back(tmp[i]);
  }

  ds.AddField(vtkm::cont::make_FieldPoint(vname,
                                          vtkm::cont::make_ArrayHandle(tmp, vtkm::CopyFlag::On)));
}

void
ReadVec(adiosS* stuff,
        vtkm::cont::DataSet& ds,
        const std::string& vname,
        bool isXYZ,
        bool cylAdd,
        std::string fileName="")
{
  if (fileName.empty())
    fileName = vname;

  auto var = stuff->io.InquireVariable<double>(fileName);
  std::vector<double> tmp;
  stuff->engine.Get(var, tmp, adios2::Mode::Sync);

  std::vector<vtkm::Vec3f> vec;
  int NP = numPlanes;
  if (cylAdd)
    NP = NP+1;

  std::vector<vtkm::Vec3f> vecXYZ;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> coordsRTZ;
  ds.GetField("coordsRTZ").GetData().AsArrayHandle(coordsRTZ);
  auto coordsRTZPortal = coordsRTZ.ReadPortal();

  vtkm::cont::ArrayHandle<int> nextNode;
  ds.GetField("nextnode").GetData().AsArrayHandle(nextNode);
  auto nnPortal = nextNode.ReadPortal();

  bool isB = (vname == "B");
  vtkm::Id idx = 0;
  for (int p = 0; p < NP; p++)
  {
    //R,Z,T in file:
    for (int i = 0; i < numNodes; i++)
    {
      int vidx = (isB ? i : (p*numNodes+i));
      vtkm::Vec3f v(tmp[vidx*3+0], tmp[vidx*3+2], tmp[vidx*3+1]);
      vec.push_back(v);
      if (isXYZ)
      {
        auto ptRTZ = coordsRTZPortal.Get(idx);
        vecXYZ.push_back(Cyl2CarVec(ptRTZ, v));
//        vecXYZ.push_back(TransformVec(ptRTZ, v, CYL_TO_CAR));
        idx++;
      }
    }
  }

  //Map next node...
  /*
  if (!isB)
  {
    std::vector<vtkm::Vec3f> tmp(vec.size());
    for (int p = 0; p < NP; p++)
    {
      for (int i = 0; i < numNodes; i++)
      {
        int idx = p*numNodes + i;

        if (p == 0)
          tmp[idx] = vec[idx];
        else
        {
          tmp[idx] = vec[nnPortal.Get(i)];
        }
      }
      vec = tmp;
    }
  }
  */

  ds.AddField(vtkm::cont::make_FieldPoint(vname,
                                          vtkm::cont::make_ArrayHandle(vec, vtkm::CopyFlag::On)));
  if (isXYZ)
    ds.AddField(vtkm::cont::make_FieldPoint(vname+"XYZ",
                                            vtkm::cont::make_ArrayHandle(vecXYZ, vtkm::CopyFlag::On)));
}

vtkm::cont::DataSet
ReadMesh(adiosS* meshStuff, adiosS* /*dataStuff*/, bool isXYZ, bool cylAdd)
{
  std::vector<double> rz;
  std::vector<int> conn, nextnode;
  meshStuff->engine.Get(meshStuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
  std::string connNm = "/cell_set[0]/node_connect_list";
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);
  meshStuff->engine.Get(meshStuff->io.InquireVariable<int>("nextnode"), nextnode, adios2::Mode::Sync);

  double dPhi = vtkm::TwoPi()/static_cast<double>(numPlanes);

  std::vector<vtkm::Vec3f> coords, coordsXYZ, coordsRTZ;
  std::vector<vtkm::Id> connIds;
  double phi = 0.0;

  int NP = numPlanes;
  if (!isXYZ && cylAdd)
    NP++;

  //coords
  for (int p = 0; p < NP; p++)
  {
    if (p == NP-1)
      phi = vtkm::TwoPi();

    for (int i = 0; i < numNodes; i++)
    {
      double R = rz[i*2 +0];
      double Z = rz[i*2 +1];

      vtkm::Vec3f ptXYZ(R*cos(phi), R*sin(phi), Z);
      vtkm::Vec3f ptRTZ(R,phi,Z);

      if (isXYZ)
        coords.push_back(ptXYZ);
      else
        coords.push_back(ptRTZ);
      coordsXYZ.push_back(ptXYZ);
      coordsRTZ.push_back(ptRTZ);
    }
    phi += dPhi;
  }

  //cells
  for (int p = 0; p < numPlanes; p++)
  {
    if (!isXYZ && p == (numPlanes-1) && !cylAdd)
      break;

    for (int i = 0; i < numTri*3; i+=3)
    {
      int off = p*numNodes;
      int p0 = conn[i+0];
      int p1 = conn[i+1];
      int p2 = conn[i+2];
      connIds.push_back(p0 + off);
      connIds.push_back(p1 + off);
      connIds.push_back(p2 + off);

      off = (p+1)*(numNodes);
      if (isXYZ && p == (numPlanes-1))
        off = 0;
/*
      p0 = nextnode[p0];
      p1 = nextnode[p1];
      p2 = nextnode[p2];
*/

      connIds.push_back(p0 + off);
      connIds.push_back(p1 + off);
      connIds.push_back(p2 + off);
    }
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  auto grid = dsb.Create(coords, vtkm::CellShapeTagWedge(), 6, connIds, "coords");
  grid.AddField(vtkm::cont::make_FieldPoint("coordsXYZ",
                                          vtkm::cont::make_ArrayHandle(coordsXYZ, vtkm::CopyFlag::On)));
  grid.AddField(vtkm::cont::make_FieldPoint("coordsRTZ",
                                          vtkm::cont::make_ArrayHandle(coordsRTZ, vtkm::CopyFlag::On)));
  grid.AddField(vtkm::cont::make_Field("nextnode",
                                       vtkm::cont::Field::Association::WHOLE_MESH,
                                       nextnode,
                                       vtkm::CopyFlag::On));
  return grid;
}

void
CompareB(vtkm::cont::DataSet& ds)
{
  using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec<double,3>>;
  using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
  using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;

  FieldHandle VField, BsField;
  ds.GetField("V").GetData().AsArrayHandle(VField);
  ds.GetField("Bs").GetData().AsArrayHandle(BsField);

  vtkm::Id n = BsField.GetNumberOfValues();
  vtkm::cont::ArrayHandle<double> err;
  err.Allocate(n);

  auto VPortal = VField.ReadPortal();
  auto BsPortal = BsField.ReadPortal();
  auto portal = err.WritePortal();

  for (vtkm::Id i = 0; i < n; i++)
  {
    auto diff = VPortal.Get(i) - BsPortal.Get(i);
    auto err  = vtkm::Magnitude(diff);
    portal.Set(i, err);
  }
  ds.AddField(vtkm::cont::make_FieldPoint("ERR", err));

  ds.PrintSummary(std::cout);
  vtkm::io::VTKDataSetWriter writer("debug.vtk");
  writer.WriteDataSet(ds);
}

void
RunPoincare(const vtkm::cont::DataSet& ds,
            std::map<std::string, std::vector<std::string>> args)
{
  using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec<double,3>>;
  using FieldType = vtkm::worklet::particleadvection::VelocityField<FieldHandle>;
  using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldType>;
  using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;
  using Stepper = vtkm::worklet::particleadvection::Stepper<RK4Type, GridEvalType>;

  vtkm::Id maxPunctures = std::stoi(args["--numPunc"][0]);
  vtkm::FloatDefault stepSize = std::stof(args["--stepSize"][0]);
  std::string vname = args["--varname"][0];

  std::cout<<"RunPoinc: "<<vname<<std::endl;

  FieldHandle BField;
  ds.GetField(vname).GetData().AsArrayHandle(BField);
  FieldType velocities(BField);
  GridEvalType eval(ds, velocities);
  Stepper rk4(eval, stepSize);

  vtkm::worklet::Poincare p;
  vtkm::Plane<> plane({0,3,0}, {0,1,0});

  vtkm::Id maxSteps = maxPunctures*100000000;

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
      seeds.push_back({vtkm::Particle({r, 2.95, 0}, id)});
  }
  else if (args.find("--seed") != args.end())
  {
    auto vals = args["--seed"];
    vtkm::FloatDefault r = std::atof(vals[0].c_str());
    vtkm::FloatDefault z = std::atof(vals[1].c_str());
    vtkm::FloatDefault t = std::atof(vals[2].c_str());
    seeds.push_back({vtkm::Particle({r, t, z}, 0)});
  }


#if 0
  if (numSeeds > 0)
  {
    vtkm::FloatDefault x0 = 2.9, x1 = 3.5;
    x0 = 2.3; x1 = 3.75; //0,1
    x0 = 3.5; x1 = 3.75;
    x0 = 3.5; x1 = 3.55; //.5, .58
/*
    x0 = 3; x1 = 3.5; //.05, .5
    x0 = 3.5; x1 = 3.55; //.5, .58
    x0 = 3.7; x1 = 4; //.8, 1
    x0 = 3.63; x1 = 3.69; //
*/

    vtkm::FloatDefault dx = (x1-x0) / (float)(numSeeds-1);
    vtkm::FloatDefault x = x0;

    for (vtkm::Id id = 0; id < numSeeds; id++, x+=dx)
      seeds.push_back({vtkm::Particle({x, 2.95, 0}, id)});
  }
  else
  {
    plane = vtkm::Plane<>({0,vtkm::TwoPi(),0}, {0,1,0});
    seeds.push_back(vtkm::Particle({2.2225268334841566, vtkm::TwoPi()-.01f, 0.4918939452569903}, 0));
  }
#endif

  auto seedsArr = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::On);

  auto t1 = std::chrono::high_resolution_clock::now();
  auto res = p.Run(rk4, seedsArr, plane, maxSteps, maxPunctures, true);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
  std::cout<<"Timer= "<<dt.count()<<std::endl;

  vtkm::cont::DataSet out;
  out.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", res.Positions));

  vtkm::cont::CellSetSingleType<> cellSet;
  vtkm::cont::ArrayHandle<vtkm::Id> connIds;
  vtkm::Id nPts = out.GetNumberOfPoints();

#if 1
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandleIndex(nPts), connIds);
  cellSet.Fill(nPts, vtkm::CELL_SHAPE_VERTEX, 1, connIds);
  out.SetCellSet(cellSet);
  std::vector<vtkm::FloatDefault> ids(nPts);
  for (int i = 0; i < nPts; i++)
    ids[i] = i;
  out.AddField(vtkm::cont::make_FieldPoint("id", vtkm::cont::make_ArrayHandle(ids, vtkm::CopyFlag::On)));
  using M = vtkm::rendering::MapperPoint;
  using C = vtkm::rendering::CanvasRayTracer;
  using V3 = vtkm::rendering::View3D;
  using V2 = vtkm::rendering::View2D;

  vtkm::cont::ColorTable colorTable("inferno");

  M mapper;
  std::cout << "Testing uniform delta radius\n";
//  mapper.SetRadius(.01f);
//  mapper.UseVariableRadius(false);
  vtkm::rendering::testing::Render<M, C, V2>(
    mapper, out, "id", colorTable, "pts.png");

  /*
  Renderer ren(512);
  vtkm::Range scalarRange(0, 1);
  std::string colorTable = "inferno"; //"jet";
  std::string renderType = "Point";
  std::string viewType = "poinc"; //"circular"; //"d3d";
  ren.Render({out}, "id", scalarRange, colorTable, renderType, viewType, "out.png");
  out.PrintSummary(std::cout);
  */
#endif

  vtkm::FloatDefault eq_axis_r = 2.8, eq_axis_z = 0.0;
  vtkm::FloatDefault eq_x_psi = 0.0697345, eq_x_r = 2.8, eq_x_z = -0.99988;

  vtkm::Id numCells = res.PolyLines.GetNumberOfCells();

  std::vector<vtkm::Id> cpids, ptIds(nPts, -1);
  for (vtkm::Id i = 0; i < numCells; i++)
  {
    auto n = res.PolyLines.GetNumberOfPointsInCell(i);
    cpids.resize(n);
    res.PolyLines.GetCellPointIds(i, cpids.data());
    for (int j = 0; j < n; j++)
      ptIds[cpids[j]] = i;
  }

  if (1)
  {
    std::ofstream outPts, outPtsPsiTheta;
    outPts.open(vname+".points.txt");
    outPtsPsiTheta.open("pointsPsiTheta.txt");
    int nPts = res.Positions.GetNumberOfValues();
    outPtsPsiTheta<<"ID,THETA,PSI,zz,R,Z"<<std::endl;
    auto portal = res.Positions.ReadPortal();
    for (int i = 0; i < nPts; i++)
    {
      auto pid = ptIds[i];
      auto pt = portal.Get(i);
      auto R = pt[0];
      auto Z = pt[2];

      auto theta = vtkm::ATan2(Z-eq_axis_z, R-eq_axis_r);
      if (theta < 0)
        theta += vtkm::TwoPi();
      //auto psi = ((R-eq_x_r)*(R-eq_x_r) + Z*Z) / eq_x_psi;
      auto psi = ((R-eq_x_r)*(R-eq_x_r) + Z*Z);

      outPts<<pid<<", "<<R<<", "<<Z<<", 0"<<std::endl;
      outPtsPsiTheta<<pid<<", "<<theta<<", "<<psi<<", 0, "<<R<<", "<<Z<<std::endl;
    }
    outPts.close();
  }
}

void
RunDebug(const vtkm::cont::DataSet& ds,
         vtkm::cont::DataSet& dsCyl)
{
  vtkm::cont::ArrayHandle<vtkm::Vec3f> B_XYZ, B_RTZ, coordsXYZ;
  ds.GetField("B").GetData().AsArrayHandle(B_XYZ);
  ds.GetField("coordsXYZ").GetData().AsArrayHandle(coordsXYZ);

  auto bPortal = B_XYZ.ReadPortal();
  auto cPortal = coordsXYZ.ReadPortal();
  vtkm::Id n = bPortal.GetNumberOfValues();

  B_RTZ.Allocate(n);
  auto bRTZPortal = B_RTZ.WritePortal();
  for (vtkm::Id i = 0; i < n; i++)
  {
    auto pXYZ = cPortal.Get(i);
    auto bXYZ = bPortal.Get(i);

    vtkm::Vec3f v = TransformVec(pXYZ, bXYZ, CAR_TO_CYL);
    bRTZPortal.Set(i, v);
  }
  dsCyl.AddField(vtkm::cont::make_FieldPoint("B_RTZ", B_RTZ));

  vtkm::io::VTKDataSetWriter writer("debug.vtk");
  writer.WriteDataSet(dsCyl);
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


void testbumbum()
{
  using M = vtkm::rendering::MapperPoint;
  using C = vtkm::rendering::CanvasRayTracer;
  using V3 = vtkm::rendering::View3D;
  using V2 = vtkm::rendering::View2D;

  vtkm::cont::testing::MakeTestDataSet maker;
  vtkm::cont::ColorTable colorTable("inferno");

  vtkm::cont::DataSetBuilderExplicit dsb;
  std::vector<vtkm::Vec3f> pts = {{0,0,0},
                                  {0,1,0},
                                  {1,0,0},
                                  {.3,.75,0},
                                  {.9,.9,0},
                                  {.25,.25,0},
                                  {.5,.5,0},
                                  {12,2,0},
                                  {11,1,0},
                                  {2,2,0},
                                  {4,2,0},
                                  {1,1,0}};
  std::vector<vtkm::Id> connIds;
  std::vector<vtkm::FloatDefault> var;
  for (int i = 0; i < pts.size(); i++)
  {
    connIds.push_back(i);
    var.push_back(2.0f + i);
  }
  auto grid = dsb.Create(pts, vtkm::CellShapeTagVertex(), 1, connIds);//, "coords");
  grid.AddField(vtkm::cont::make_FieldPoint("var",
                                            vtkm::cont::make_ArrayHandle(var, vtkm::CopyFlag::On)));

  M mapper;
  std::cout << "Testing uniform delta raduis\n";
  //mapper.SetRadiusDelta(9.0f);
  mapper.SetRadius(.05f);
  mapper.UseCells();
  vtkm::rendering::testing::Render<M, C, V2>(mapper, grid, "var", colorTable, "pts.png");

  /*
  // restore defaults
  mapper.SetRadiusDelta(0.5f);
  mapper.UseVariableRadius(false);

  mapper.SetRadius(0.2f);
  vtkm::rendering::testing::Render<M, C, V3>(
    mapper, maker.Make3DUniformDataSet1(), "pointvar", colorTable, "points_reg3D.pnm");

  mapper.UseCells();
  mapper.SetRadius(1.f);
  vtkm::rendering::testing::Render<M, C, V3>(
    mapper, maker.Make3DExplicitDataSet7(), "cellvar", colorTable, "spheres.pnm");
  */

}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  std::cout<<std::endl<<std::endl;

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
  bool isXYZ = false;
  bool cylAdd = false;
  auto ds = ReadMesh(meshStuff, dataStuff, isXYZ, cylAdd);
  ReadScalar(dataStuff, ds, "dpot", isXYZ, cylAdd);
  ReadScalar(dataStuff, ds, "apars", isXYZ, cylAdd);
  ReadVec(bfieldStuff, ds, "B", isXYZ, cylAdd, "/node_data[0]/values");
  ReadVec(bfield_allStuff, ds, "Bs", isXYZ, cylAdd, "Bs");

  vtkm::cont::Invoker invoker;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> a_bhat, b;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> apars;
  ds.GetField("B").GetData().AsArrayHandle(b);
  ds.GetField("apars").GetData().AsArrayHandle(apars);
  invoker(CalculateABhat{}, b, apars, a_bhat);
  ds.AddField(vtkm::cont::make_FieldPoint("As_bHat", a_bhat));

  vtkm::filter::Gradient grad;
  grad.SetComputePointGradient(true);
  grad.SetActiveField("As_bHat");
  grad.SetOutputFieldName("gradient");
  ds = grad.Execute(ds);

  vtkm::cont::ArrayHandle<vtkm::Vec3f> coords, output;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Vec3f,3>> gradient;

  ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);
  ds.GetField("gradient").GetData().AsArrayHandle(gradient);
  invoker(CalculateVecField{}, coords, b, a_bhat, gradient, output);
  ds.AddField(vtkm::cont::make_FieldPoint("V", output));

  ds.PrintSummary(std::cout);
  vtkm::io::VTKDataSetWriter writer("debug.vtk");
  writer.WriteDataSet(ds);
  return 0;

  //CompareB(ds);

  RunPoincare(ds, args);

  return 0;
}
