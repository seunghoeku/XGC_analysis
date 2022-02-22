#include "io.h"
#include <fides/DataSetReader.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>

adios2::ADIOS *adios = NULL;
int totNumPlanes = -1;
int numNodes=-1;
int numTri=-1;
int numPlanesInFile=-1;

std::map<std::string, adiosS*> adiosStuff;


static vtkm::Vec3f Cyl2CarVec(const vtkm::Vec3f& pCyl,
                              const vtkm::Vec3f& vCyl,
                              const vtkm::Vec3f& pCar0)
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

static vtkm::cont::DataSet
CreateXGCGeom(const std::vector<double>& rz,
              const std::vector<int>& conn,
              int numPlanes,
              double totAngle,
              std::vector<vtkm::FloatDefault> &RArr,
              std::vector<vtkm::FloatDefault> &ZArr,
              std::vector<vtkm::FloatDefault> &PhiArr,
              std::vector<vtkm::FloatDefault> &PlaneArr)
{

  auto rzArr = vtkm::cont::make_ArrayHandle(rz, vtkm::CopyFlag::On);
  bool isCyl = true;
  auto coords = vtkm::cont::make_ArrayHandleXGCCoordinates(rzArr, numPlanes, isCyl);

  int numNodes = rz.size()/2;
  std::vector<vtkm::Int32> nextNode(numNodes);
  for (int i = 0; i <numNodes; i++)
    nextNode[i] = i;

  auto connArr = vtkm::cont::make_ArrayHandle(conn, vtkm::CopyFlag::On);
  auto nextNodeArr = vtkm::cont::make_ArrayHandle(nextNode, vtkm::CopyFlag::On);
  auto cellSet = vtkm::cont::make_CellSetExtrude(connArr, coords, nextNodeArr, isCyl);

  RArr.resize(numNodes * totNumPlanes);
  ZArr.resize(numNodes * totNumPlanes);
  PhiArr.resize(numNodes * totNumPlanes);
  PlaneArr.resize(numNodes * totNumPlanes);

  double phi = 0.0;
  double dPhi = totAngle/static_cast<double>(numPlanes);
  int idx = 0;
  for (int p = 0; p < numPlanes; p++)
  {
    for (int i = 0; i < numNodes; i++)
    {
      RArr[idx] = rz[i*2 + 0];
      ZArr[idx] = rz[i*2 + 1];
      PhiArr[idx] = phi;
      PlaneArr[idx] = p;
    }
    phi += dPhi;
  }

  vtkm::cont::DataSet ds;
  ds.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coords", coords));
  ds.SetCellSet(cellSet);

  return ds;
}

static void CreateGeom(bool fullTorus,
                       const std::vector<double>& rz,
                       const std::vector<int>& conn,
                       int numPlanes,
                       double totAngle,
                       std::vector<vtkm::Vec3f>& coords,
                       std::vector<vtkm::Id> &connIds,
                       std::vector<vtkm::FloatDefault> &RArr,
                       std::vector<vtkm::FloatDefault> &ZArr,
                       std::vector<vtkm::FloatDefault> &PhiArr,
                       std::vector<vtkm::FloatDefault> &PlaneArr,
                       bool isXYZ,
                       bool is2D)
{
    int nNodes = rz.size()/2;
    if (is2D)
      numPlanes = 1;

    coords.resize(nNodes * numPlanes);
    RArr.resize(coords.size());
    ZArr.resize(coords.size());
    PhiArr.resize(coords.size());
    PlaneArr.resize(coords.size());

    double phi = 0.0;
    double dPhi = totAngle/static_cast<double>(numPlanes);


    for (int p = 0; p < numPlanes; p++)
    {
      //std::cout<<"Plane= "<<p<<" phi= "<<phi<<" last = "<<(p==numPlanes-1)<<std::endl;
        vtkm::Vec3f pt;
        for (int i = 0; i < nNodes; i++)
        {
            double R = rz[i*2 +0];
            double Z = rz[i*2 +1];

            if (is2D)
            {
              pt[0] = R;
              pt[1] = Z;
              pt[2] = 1;
            }
            else
            {
              if (isXYZ)
              {
                pt[0] = R*cos(phi);
                pt[1] = R*sin(phi);
                pt[2] = Z;
              }
              else
              {
                pt[0] = R;
                pt[1] = phi;
                pt[2] = Z;
              }
            }
            int idx = i+p*nNodes;
            coords[idx] = pt;

            RArr[idx] = static_cast<vtkm::FloatDefault>(R);
            ZArr[idx] = static_cast<vtkm::FloatDefault>(Z);
            PhiArr[idx] = totAngle - -phi;
            PlaneArr[idx] = 1 - static_cast<vtkm::FloatDefault>(p)/static_cast<vtkm::FloatDefault>(totNumPlanes);
        }

        //phi -= dPhi;
        phi += dPhi;
    }
    std::cout<<"FullTorus= "<<fullTorus<<std::endl;
    std::cout<<"#coords= "<<coords.size()<<std::endl;
    std::cout<<"#conn= "<<conn.size()<<std::endl;
    std::cout<<"totNumPlanes= "<<totNumPlanes<<std::endl;

    int plane = 0;
    for (int p = 0; p < numPlanes; p++)
    {
      if (!isXYZ && p == numPlanes-1)
        break;

      for (int i = 0; i < numTri*3; i+=3)
      {
        int off = plane*nNodes;
        int p0 = conn[i+0];
        int p1 = conn[i+1];
        int p2 = conn[i+2];

        int badID = -1;
        connIds.push_back(p0 + off);
        connIds.push_back(p1 + off);
        connIds.push_back(p2 + off);
        if (is2D)
          continue;

        if (p0+off >= coords.size()) badID = 0;
        if (p1+off >= coords.size()) badID = 1;
        if (p2+off >= coords.size()) badID = 2;

        off = (plane+1)*(nNodes);
        //p0 = nnPtr[p0];
        //p1 = nnPtr[p1];
        //p2 = nnPtr[p2];

        //Connect back to the first plane.
        if (plane == (numPlanes-1) && fullTorus)
          off = 0;

        connIds.push_back(p0 + off);
        connIds.push_back(p1 + off);
        connIds.push_back(p2 + off);

        if (p0+off >= coords.size()) badID = 3;
        if (p1+off >= coords.size()) badID = 4;
        if (p2+off >= coords.size()) badID = 5;

        if (badID != -1)
        {
          std::cout<<"BAD ID: plane= "<<plane<<" offset= "<<badID<<" i= "<<i<<std::endl;
          std::cout<<"   "<<p0+off<<" "<<p1+off<<" "<<p2+off<<std::endl;
          std::cout<<"   off= "<<plane*nNodes<<" "<<off<<std::endl;
          exit(-1);
        }
      }
      plane++;
    }
}

vtkm::cont::DataSet
readMesh(adiosS *stuff, int nPlanesBetween, bool extendToFullTorus, bool isXYZ, bool is2D, bool isExplicit)
{
    stuff->engine.BeginStep();
    stuff->engine.Get(stuff->io.InquireVariable<int>("n_n"), &numNodes, adios2::Mode::Sync);
    stuff->engine.Get(stuff->io.InquireVariable<int>("n_t"), &numTri, adios2::Mode::Sync);

    std::vector<double> rz;
    std::vector<int> conn;
    stuff->engine.Get(stuff->io.InquireVariable<double>("/coordinates/values"), rz, adios2::Mode::Sync);
    std::string connNm = "/cell_set[0]/node_connect_list";
    connNm = "nd_connect_list";
    stuff->engine.Get(stuff->io.InquireVariable<int>(connNm), conn, adios2::Mode::Sync);

    std::vector<double> psi;
    stuff->engine.Get(stuff->io.InquireVariable<double>("psi"), psi, adios2::Mode::Sync);
    std::for_each(psi.begin(), psi.end(), [](double &v) {v /= psi_x;});

    vtkm::cont::DataSet grid;
    totNumPlanes = numPlanesInFile + (numPlanesInFile*nPlanesBetween);

    if (extendToFullTorus && nWedge != 1)
        totNumPlanes = (nWedge*numPlanesInFile) + ((nWedge*numPlanesInFile)*nPlanesBetween);

    double totAngle = (2*M_PI) / static_cast<double>(nWedge);
    if (extendToFullTorus && nWedge != 1)
        totAngle = 2*M_PI;

    double dPhi = totAngle/static_cast<double>(totNumPlanes);


    std::vector<vtkm::FloatDefault> RArr, ZArr, PhiArr, PlaneArr;

    if (isExplicit)
    {
      std::vector<vtkm::Vec3f> coords;
      std::vector<vtkm::Id> connIds;
      CreateGeom(extendToFullTorus, rz, conn, totNumPlanes, totAngle, coords, connIds, RArr, ZArr, PhiArr, PlaneArr, isXYZ, is2D);

      vtkm::cont::DataSetBuilderExplicit dsb;
      if (is2D)
        grid = dsb.Create(coords, vtkm::CellShapeTagTriangle(), 3, connIds, "coords");
      else
        grid = dsb.Create(coords, vtkm::CellShapeTagWedge(), 6, connIds, "coords");
    }
    else
    {
      grid = CreateXGCGeom(rz, conn, totNumPlanes, totAngle, RArr, ZArr, PhiArr, PlaneArr);
    }

    /*
    if (!is2D && isXYZ)
    {
      std::vector<vtkm::Vec3f> cylCoords(coords.size());
      vtkm::Id n = coords.size();
      for (vtkm::Id i = 0; i < n; i++)
      {
        auto pt = coords[i];
        vtkm::Vec3d ptCyl;
        ptCyl[0] = vtkm::Sqrt(pt[0]*pt[0] + pt[2]*pt[2]);
        ptCyl[1] = vtkm::ATan2(pt[2], pt[0]);
        if (ptCyl[1] < 0)
          ptCyl[1] += vtkm::TwoPi();
        ptCyl[2] = pt[2];
        cylCoords[i] = ptCyl;
      }
      grid.AddField(vtkm::cont::make_FieldPoint("coords
    }
    */

    if (!is2D)
    {
      int sz = psi.size();
      std::vector<vtkm::FloatDefault> psif(sz);
      for (int i = 0; i < sz; i++)
        psif[i] = psi[i];

      std::vector<vtkm::FloatDefault> psi3D(RArr.size(), -1.0f);
      int idx = 0;
      for (int i = 0; i < totNumPlanes; i++)
        for (int j = 0; j < sz; j++)
          psi3D[idx++] = static_cast<vtkm::FloatDefault>(psi[j]);


      grid.AddField(vtkm::cont::make_Field("psi", vtkm::cont::Field::Association::WHOLE_MESH,
                                           psif, vtkm::CopyFlag::On));
      grid.AddField(vtkm::cont::make_FieldPoint("psi3D", vtkm::cont::make_ArrayHandle(psi3D, vtkm::CopyFlag::On)));
      grid.AddField(vtkm::cont::make_FieldPoint("R", vtkm::cont::make_ArrayHandle(RArr, vtkm::CopyFlag::On)));
      grid.AddField(vtkm::cont::make_FieldPoint("Z", vtkm::cont::make_ArrayHandle(ZArr, vtkm::CopyFlag::On)));
      grid.AddField(vtkm::cont::make_FieldPoint("Phi", vtkm::cont::make_ArrayHandle(PhiArr, vtkm::CopyFlag::On)));
      grid.AddField(vtkm::cont::make_FieldPoint("Plane", vtkm::cont::make_ArrayHandle(PlaneArr, vtkm::CopyFlag::On)));
    }

    return grid;
}

vtkm::cont::DataSet
ReadMesh(std::map<std::string, adiosS*> &adiosStuff, bool fullGrid, bool extend, bool isXYZ, bool is2D, bool isExplicit)
{
  auto data = adiosStuff["data"];
  data->engine.Get(data->io.InquireVariable<int>("nphi"), &numPlanesInFile, adios2::Mode::Sync);

  return readMesh(adiosStuff["mesh"], planesBetween, extend, isXYZ, is2D, isExplicit);
}

std::vector<vtkm::FloatDefault>
InterpFieldToSubPlanes(const std::vector<vtkm::FloatDefault> &input, int numPhi, int nPlanesBetween, int numNodes)
{
  std::vector<vtkm::FloatDefault> result;
  if (nPlanesBetween == 0 && nWedge == 1)
    result = input;
  else
  {
    vtkm::FloatDefault dT = static_cast<vtkm::FloatDefault>(1)/static_cast<vtkm::FloatDefault>(nPlanesBetween+1);

    int cnt = 0;
    for (int i = 0; i < nWedge; i++)
    {
      //std::cout<<"N= "<<i<<std::endl;
      for (int j = 0; j < numPhi; j++)
      {
        //Add the first plane.
        //std::cout<<" J= "<<j<<" : cnt= "<<cnt<<std::endl;
        for (int k = 0; k < numNodes; k++)
          result.push_back(input[j*numNodes + k]);
        cnt++;

        int offsetP0 = j*numNodes;
        int offsetP1 = (j+1)*numNodes;
        if (j == numPhi-1)
          offsetP1 = 0;

        vtkm::FloatDefault t = dT;
        for (int p = 0; p < nPlanesBetween; p++)
        {
          //std::cout<<"  Interp: "<<t<<" ("<<offsetP0/numNodes<<" "<<offsetP1/numNodes<<") : cnt= "<<cnt<<std::endl;
          for (int k = 0; k < numNodes; k++)
          {
            auto v0 = input[offsetP0 + k];
            auto v1 = input[offsetP1 + k];
            auto v = v0 + t*(v1-v0);
            result.push_back(v);
          }
          t += dT;
          cnt++;
        }
      }
    }
    //std::cout<<"Interp: nW= "<<nWedge<<" "<<numNodes<<" "<<numPhi<<" x "<<nPlanesBetween<<" --> "<<result.size()<<std::endl;
    //std::cout<<"cnt= "<<cnt<<std::endl;
  }

  return result;
}
void
READVAR(const std::string &vname, std::vector<vtkm::FloatDefault> &arr, adiosS *data, bool doInterp=true, bool extendFull=true)
{
    auto v = data->io.InquireVariable<double>(vname);
    //auto x = v.AllStepsBlocksInfo();
    //std::cout<<" x: "<<x.size()<<" "<<x[0][0].Step<<std::endl;

    std::vector<double> tmp;
    data->engine.Get(v, tmp, adios2::Mode::Sync);

    arr.resize(tmp.size());
    std::transform(tmp.begin(), tmp.end(), arr.begin(), [](double x) {return static_cast<float>(x);});

    if (doInterp)
      arr = InterpFieldToSubPlanes(arr, numPlanesInFile, planesBetween, numNodes);
}

adios2::StepStatus
ReadVar(const std::string &vname, adiosS *data, vtkm::cont::DataSet &ds, bool is2D, bool isXYZ, std::string dataSetVarName)
{
    adios2::StepStatus status;

    status = data->engine.BeginStep();
    /*
    if (status != adios2::StepStatus::OK)
    {
      std::cout<<"Read failed for "<<vname<<std::endl;
      return status;
    }
    */

    std::cout<<"Reading Step= "<<data->engine.CurrentStep()<<std::endl;
    std::vector<vtkm::FloatDefault> arr;
    READVAR(vname, arr, data);

    if (dataSetVarName.empty())
      dataSetVarName = vname;

    std::cout<<"Reading: "<<vname<<std::endl;

    if (dataSetVarName == "B")
    {
      //For cylindrical coords: X=r, Y=theta Z=z (X=r=B0, Y=theta=B2, Z=Z=B1)
      //B is written out as: B0, B2, B1 (R,Theta,Z)
      // R=B[0]
      // Theta = B[2]
      // Z=B[1]
      vtkm::Id idx = 0;
      vtkm::cont::ArrayHandle<vtkm::Vec3f> bvec;

      if (is2D)
      {
        bvec.Allocate(numNodes);
        auto portal = bvec.WritePortal();
        for (int i = 0; i < numNodes; i++)
        {
          auto R = arr[i*3+0];
          auto T = arr[i*3+2];
          auto Z = arr[i*3+1];
          vtkm::Vec3f v(R,Z,0);
          portal.Set(i, v);
        }
      }
      else
      {
        using XGCCoordsType = vtkm::cont::ArrayHandleXGCCoordinates<vtkm::FloatDefault>;
        using CoordsType = vtkm::cont::ArrayHandle<vtkm::Vec3f>;

        bool isXGC = false;
        XGCCoordsType xgcCoords;
        CoordsType coords;
        if (ds.GetCoordinateSystem().GetData().IsType<XGCCoordsType>())
        {
          ds.GetCoordinateSystem().GetData().AsArrayHandle(xgcCoords);
          isXGC = true;
        }
        else
          ds.GetCoordinateSystem().GetData().AsArrayHandle(coords);

        vtkm::cont::ArrayHandle<vtkm::FloatDefault> R, Z, Phi;
        ds.GetField("R").GetData().AsArrayHandle(R);
        ds.GetField("Z").GetData().AsArrayHandle(Z);
        ds.GetField("Phi").GetData().AsArrayHandle(Phi);
        auto rPortal = R.ReadPortal();
        auto zPortal = Z.ReadPortal();
        auto tPortal = Phi.ReadPortal();

        bvec.Allocate(totNumPlanes * numNodes);
        auto portal = bvec.WritePortal();
        for (int i = 0; i < totNumPlanes; i++)
        {
          for (int j = 0; j < numNodes; j++)
          {
            vtkm::Vec3f pCar;
            if (isXGC) pCar = xgcCoords.ReadPortal().Get(idx);
            else pCar = coords.ReadPortal().Get(idx);
            //auto pCar = cPortal.Get(idx);
            auto R = arr[i*3+0];
            auto T = arr[i*3+2];
            auto Z = arr[i*3+1];

            vtkm::Vec3f vCyl(R, T, Z);

            if (isXYZ)
            {
              auto ptR = rPortal.Get(idx);
              auto ptZ = zPortal.Get(idx);
              auto ptT = tPortal.Get(idx);

              //R,T,Z
              vtkm::Vec3f pCyl(ptR,ptT,ptZ);
              //R,Z,T  <<<-------------
              //vtkm::Vec3f vCyl(arr[j*3+0], arr[j*3+2], arr[j*3+1]); //this works....
              auto vCar = Cyl2CarVec(pCyl, vCyl, pCar);
              portal.Set(idx, vCar);
            }
            else
            {
              portal.Set(idx, vCyl);
            }

            idx++;
          }
        }
      }
      std::vector<vtkm::FloatDefault> b0, b1, b2;
      for (int j = 0; j < totNumPlanes; j++)
        for (int i = 0; i < arr.size(); i+=3)
        {
          b0.push_back(arr[i+0]);
          b1.push_back(arr[i+1]);
          b2.push_back(arr[i+2]);
        }
      ds.AddField(vtkm::cont::make_FieldPoint("B0", vtkm::cont::make_ArrayHandle(b0, vtkm::CopyFlag::On)));
      ds.AddField(vtkm::cont::make_FieldPoint("B1", vtkm::cont::make_ArrayHandle(b1, vtkm::CopyFlag::On)));
      ds.AddField(vtkm::cont::make_FieldPoint("B2", vtkm::cont::make_ArrayHandle(b2, vtkm::CopyFlag::On)));

      ds.AddField(vtkm::cont::make_FieldPoint(dataSetVarName, bvec));
    }
    else if (!ds.HasField(dataSetVarName))
    {
        auto ah = vtkm::cont::make_ArrayHandle(arr, vtkm::CopyFlag::On);
        ds.AddField(vtkm::cont::make_FieldPoint(dataSetVarName, ah));
    }
    else
    {
        auto f = ds.GetField(dataSetVarName).GetData();
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> fieldAH;
        fieldAH = f.Cast<vtkm::cont::ArrayHandle<vtkm::FloatDefault>>();
        auto portal = fieldAH.WritePortal();
        for (std::size_t i = 0; i < arr.size(); i++)
          portal.Set((vtkm::Id)i, arr[i]);
    }

    return status;
}
