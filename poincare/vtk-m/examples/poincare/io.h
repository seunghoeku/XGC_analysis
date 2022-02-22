#ifndef vtk_m_examples_poincare_io_h
#define vtk_m_examples_poincare_io_h

#include <vtkm/cont/DataSet.h>
#include <string>
#include <vector>
#include <map>

#include <fides/DataSetReader.h>
#include <adios2.h>

extern adios2::ADIOS *adios;
extern int totNumPlanes;
extern int numNodes;
extern int numTri;
extern int numPlanesInFile;
const int planesBetween = 0;

const bool extendToFull = true;

class adiosS;
extern std::map<std::string, adiosS*> adiosStuff;
const double psi_x = 0.2661956235889000;
const double xpoint_r = 1.557545038402000;
const double xpoint_z = -1.177067412978000;
const int nWedge = 1;

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

vtkm::cont::DataSet
ReadMesh(std::map<std::string, adiosS*> &adiosStuff, bool fullGrid=true, bool extend=true, bool isXYZ=true, bool is2D=false, bool isExplicit=true);

adios2::StepStatus
ReadVar(const std::string &vname, adiosS *data, vtkm::cont::DataSet &ds, bool is2D=false, bool isXYZ=true, std::string dataSetVarName="");


#endif
