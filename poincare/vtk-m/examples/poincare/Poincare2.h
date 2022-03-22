#ifndef PoincareWorklet2_h
#define PoincareWorklet2_h

#include <vtkm/Geometry.h>
#include <vtkm/Matrix.h>
#include <vtkm/Particle.h>
#include <vtkm/cont/CellLocatorTwoLevel.h>
#include <vtkm/worklet/WorkletMapField.h>

#include "XGCParameters.h"
#include <iomanip>

using Ray3f = vtkm::Ray<vtkm::FloatDefault, 3, true>;



//
//./examples/poincare/Simple2.3 --vField B --dir ../data/sku_8000/POINC --worklet 1 --traces 1 --numPunc 2 --stepSize 0.01 --useHighOrder --jong1 --output bumm



//XGC code:  /gpfs/alpine/proj-shared/csc143/pugmire/XGC-Devel/XGC-Devel-poincare-debug
// make xgc-eem
// make xgc-eem-cpp
//run: bsub job-summit.sh
/*

Jong's branch w/ coeff and debug stuff.
/gpfs/alpine/proj-shared/csc143/jyc/summit/exp-xgc-poincare/exp-poincare-8000-dave

//run: bsub job-summit.sh


//old code (that runs...)
/gpfs/alpine/proj-shared/csc143/pugmire/XGC-Devel/XGC-Devel-poincare-debug
 launched job: 1714079
wrong number of procs.
fixed script and ran again. Looks like it works.

run with DDT

bsub -q debug -nnodes 6  -P PHY122 -W 0:30 -Is $SHELL -l

[ -z $JOBSIZE ] && JOBSIZE=$(((LSB_DJOB_NUMPROC-1)/42))
[ -z $ISTEP ] && ISTEP=3000
WDIR=exp-poincare-$ISTEP
SDIR=poincare_plot_for_su419.org
export OMP_NUM_THREADS=2
export COLUMNS=512
export OMPI_MCA_coll_ibm_collselect_mode_barrier=failsafe

date



ddt --connect jsrun -n $((JOBSIZE*32)) -a1 -c1 -g0 -r32 -brs /usr/bin/stdbuf -oL -eL ./xgc-eem 2>&1 | tee run-$JOBID.log
ddt --connect jsrun -n 192 -a1 -c1 -g0 -r32 -brs /usr/bin/stdbuf -oL -eL ./xgc-eem 2>&1 | tee run-$JOBID.log

jsrun -n 192 -a1 -c1 -g0 -r32 -brs /usr/bin/stdbuf -oL -eL ./xgc-eem-rel 2>&1 | tee run.log


//jongs code.
/gpfs/alpine/proj-shared/csc143/jyc/summit/exp-xgc-poincare/exp-poincare-8000-dave

*/

//#define PRINT_STUFF
#ifdef PRINT_STUFF
#define DBG(x) std::cout<<x
#else
#define DBG(x)
#endif

#define DO_TRACES 0

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
class PoincareWorklet2 : public vtkm::worklet::WorkletMapField
{
  using DimensionType = vtkm::Int16;
  using DimVec3 = vtkm::Vec<DimensionType, 3>;
  using FloatVec3 = vtkm::Vec3f;

  class ParticleInfo
  {
  public:
    ParticleInfo() = default;
    ParticleInfo(const ParticleInfo&) = default;
    ParticleInfo& operator=(const ParticleInfo&) = default;

    vtkm::Id3 PrevCell = {-1,-1,-1};
    vtkm::FloatDefault Psi = 0;
    vtkm::FloatDefault dpsi_dr=0, dpsi_dz=0, d2psi_drdz=0, d2psi_d2r=0, d2psi_d2z=0;
    vtkm::Vec3f gradPsi_rzp, B0_rzp, curlB_rzp, curl_nb_rzp;
  };

public:
  using ControlSignature = void(FieldInOut particles,
                                ExecObject locator,
                                WholeCellSetIn<> cellSet,
                                WholeArrayIn Coords,
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,
                                WholeArrayIn coeff_1D,
                                WholeArrayIn coeff_2D,
                                WholeArrayIn B_RZP,
                                WholeArrayIn Psi,
                                WholeArrayInOut traces,
                                WholeArrayInOut outputRZ,
                                WholeArrayInOut outputTP,
                                WholeArrayInOut punctureID);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14);
  using InputDomain = _1;

  PoincareWorklet2(vtkm::Id maxPunc,
                   vtkm::FloatDefault planeVal,
                   vtkm::FloatDefault stepSize,
                   bool saveTraces,
                   const XGCParameters& xgcParams)
    : MaxIter(maxPunc * 1000000)
    , MaxPunc(maxPunc)
    , PlaneVal(planeVal)
    , StepSize(stepSize)
    , SaveTraces(saveTraces)
  {
    this->NumPlanes = xgcParams.numPlanes;
    this->NumNodes = xgcParams.numNodes;
    this->dPhi = vtkm::TwoPi()/static_cast<vtkm::FloatDefault>(this->NumPlanes);
    this->StepSize_2 = this->StepSize / 2.0;
    this->StepSize_6 = this->StepSize / 6.0;


    this->nr = xgcParams.eq_mr-1;
    this->nz = xgcParams.eq_mz-1;
    this->rmin = xgcParams.eq_min_r;
    this->rmax = xgcParams.eq_max_r;
    this->zmin = xgcParams.eq_min_z;
    this->zmax = xgcParams.eq_max_z;
    this->EqAxisR = xgcParams.eq_axis_r;
    this->EqAxisZ = xgcParams.eq_axis_z;
    this->EqXPsi = xgcParams.eq_x_psi;
    this->dr = (xgcParams.eq_max_r - xgcParams.eq_min_r) / vtkm::FloatDefault(this->nr);
    this->dz = (xgcParams.eq_max_z - xgcParams.eq_min_z) / vtkm::FloatDefault(this->nz);
    this->dr_inv = 1.0/this->dr;
    this->dz_inv = 1.0/this->dz;

    this->ncoeff = xgcParams.eq_mr-1;
    this->min_psi = xgcParams.psi_min;
    this->max_psi = xgcParams.psi_max;
    this->one_d_cub_dpsi_inv = 1.0 / ((this->max_psi-this->min_psi)/vtkm::FloatDefault(this->ncoeff));
  }

  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  vtkm::Vec3f HighOrderBOnly(const vtkm::Vec3f& ptRPZ,
                             const Coeff_1DType& Coeff_1D,
                             const Coeff_2DType& Coeff_2D) const
  {
    vtkm::FloatDefault R = ptRPZ[0], Z = ptRPZ[2];
    vtkm::Vec3f ptRZ(R,Z,0);

    int r_i = this->GetIndex(R, this->nr, this->rmin, this->dr_inv);
    int z_i = this->GetIndex(Z, this->nz, this->zmin, this->dz_inv);

    // rc(i), zc(j)
    vtkm::FloatDefault Rc = this->rmin + (vtkm::FloatDefault)(r_i)*this->dr;
    vtkm::FloatDefault Zc = this->zmin + (vtkm::FloatDefault)(z_i)*this->dz;
    auto Rc_1 = Rc + this->dr;
    auto Zc_1 = Zc + this->dz;
    Rc = (Rc + Rc_1) * 0.5;
    Zc = (Zc + Zc_1) * 0.5;

    //Get the coeffcients (z,r,4,4)
    vtkm::Matrix<vtkm::FloatDefault, 4, 4> acoeff;
    //vtkm::Id offset = (r_i * this->ncoeff + z_i) * 16; //DRP
    vtkm::Id offset = (z_i * ncoeff + r_i) * 16;

    /*
    std::cout<<"InterpolatePsi: "<<vtkm::Vec2f(R,Z)<<std::endl;
    std::cout<<"  i/j= "<<r_i<<" "<<z_i<<std::endl;
    std::cout<<"  ncoeff= "<<ncoeff<<std::endl;
    std::cout<<"  offset= "<<offset<<std::endl;
    std::cout<<"  Rc/Zc= "<<Rc<<" "<<Zc<<std::endl;
    */

    vtkm::FloatDefault psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
    //this->eval_bicub_2(R, Z, Rc, Zc, acoeff, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);
    this->EvalBicub2(R, Z, Rc, Zc, offset, Coeff_2D, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);

    vtkm::FloatDefault fld_I = this->I_interpol(psi, 0, Coeff_1D);

    vtkm::FloatDefault over_r = 1/R;
    vtkm::FloatDefault Br = -dpsi_dz * over_r;
    vtkm::FloatDefault Bz = dpsi_dr * over_r;
    vtkm::FloatDefault Bp = fld_I * over_r;

    return vtkm::Vec3f(Br, Bz, Bp);
  }

  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool HighOrderB(const vtkm::Vec3f& ptRPZ,
                  ParticleInfo& pInfo,
                  const Coeff_1DType& Coeff_1D,
                  const Coeff_2DType& Coeff_2D) const
  {
    vtkm::FloatDefault R = ptRPZ[0], Z = ptRPZ[2];
    vtkm::Vec3f ptRZ(R,Z,0);

    int r_i = this->GetIndex(R, this->nr, this->rmin, this->dr_inv);
    int z_i = this->GetIndex(Z, this->nz, this->zmin, this->dz_inv);

    // rc(i), zc(j)
    vtkm::FloatDefault Rc = rmin + (vtkm::FloatDefault)(r_i)*this->dr;
    vtkm::FloatDefault Zc = zmin + (vtkm::FloatDefault)(z_i)*this->dz;
    auto Rc_1 = Rc + this->dr;
    auto Zc_1 = Zc + this->dz;
    Rc = (Rc + Rc_1) * 0.5;
    Zc = (Zc + Zc_1) * 0.5;

    //Get the coeffcients (z,r,4,4)
    vtkm::Matrix<vtkm::FloatDefault, 4, 4> acoeff;
    //offset = ri * nz + zi
    //vtkm::Id offset = (r_i * this->ncoeff + z_i) * 16; //DRP
    vtkm::Id offset = (z_i * ncoeff + r_i) * 16;

#if 0
    vtkm::Id idx = 0;
    //std::cout<<"Offset= "<<(offset/16)<<" 16: "<<offset<<std::endl;
    for (vtkm::Id ii = 0; ii < 4; ii++)
      for (vtkm::Id jj = 0; jj < 4; jj++)
      {
        acoeff[ii][jj] = Coeff_2D.Get(offset+idx);
        idx++;
      }
#endif

    /*
    std::cout<<"  i/j= "<<r_i<<" "<<z_i<<std::endl;
    std::cout<<"  ncoeff= "<<ncoeff<<std::endl;
    std::cout<<"  offset= "<<offset<<std::endl;
    std::cout<<"  Rc/Zc= "<<Rc<<" "<<Zc<<std::endl;
    */

    vtkm::FloatDefault psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
    //this->eval_bicub_2(R, Z, Rc, Zc, acoeff, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);
    this->EvalBicub2(R, Z, Rc, Zc, offset, Coeff_2D, psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_d2r,d2psi_d2z);
    pInfo.Psi = psi;
    //PSI = psi;
    pInfo.gradPsi_rzp[0] = dpsi_dr;
    pInfo.gradPsi_rzp[1] = dpsi_dz;
    pInfo.gradPsi_rzp[2] = 0;
    pInfo.dpsi_dr = dpsi_dr;
    pInfo.dpsi_dz = dpsi_dz;
    pInfo.d2psi_drdz = d2psi_drdz;
    pInfo.d2psi_d2r = d2psi_d2r;
    pInfo.d2psi_d2z = d2psi_d2z;

    //PSI = psi;
    //gradPsi_rzp[0] = dpsi_dr;
    //gradPsi_rzp[1] = dpsi_dz;
    //gradPsi_rzp[2] = 0;

    vtkm::FloatDefault fld_I = this->I_interpol(psi, 0, Coeff_1D);
    vtkm::FloatDefault fld_dIdpsi = this->I_interpol(psi, 1, Coeff_1D);

    vtkm::FloatDefault over_r = 1/R;
    vtkm::FloatDefault over_r2 = over_r*over_r;
    vtkm::FloatDefault Br = -dpsi_dz * over_r;
    vtkm::FloatDefault Bz = dpsi_dr * over_r;
    vtkm::FloatDefault Bp = fld_I * over_r;

    pInfo.B0_rzp = vtkm::Vec3f(Br, Bz, Bp);

    const vtkm::FloatDefault bp_sign = 1.0;
/*
    vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;
    //Set the jacobian.
    const int PIR = 0;
    const int PIZ = 1;
    const int PIP = 2;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        jacobian_rzp[i][j] = 0;

    // R-derivatives depend on the geometry
    jacobian_rzp[PIR][PIR] = ( dpsi_dz * over_r2 - d2psi_drdz * over_r) * bp_sign;
    jacobian_rzp[PIR][PIZ] = (-dpsi_dr * over_r2 + d2psi_d2r  * over_r) * bp_sign;
    jacobian_rzp[PIR][PIP] = dpsi_dr * fld_dIdpsi * over_r - fld_I * over_r2;

    // Z and phi derivatives do not change between toroidal and cylindrical geometry
    jacobian_rzp[PIZ][PIR] = -d2psi_d2z * over_r * bp_sign;
    jacobian_rzp[PIP][PIR] = 0.0 * bp_sign;

    jacobian_rzp[PIZ][PIZ] = d2psi_drdz * over_r * bp_sign;
    jacobian_rzp[PIP][PIZ] = 0.0 * bp_sign;

    jacobian_rzp[PIZ][PIP] = fld_dIdpsi * dpsi_dz * over_r;
    jacobian_rzp[PIP][PIP] = 0.0;
    auto dBr_dr = jacobian_rzp[0][0];
    auto dBr_dz = jacobian_rzp[1][0];
    auto dBr_dp = jacobian_rzp[2][0];

    auto dBz_dr = jacobian_rzp[0][1];
    auto dBz_dz = jacobian_rzp[1][1];
    auto dBz_dp = jacobian_rzp[2][1];

    auto dBp_dr = jacobian_rzp[0][2];
    auto dBp_dz = jacobian_rzp[1][2];
    //auto dBp_dp = jacobian_rzp[2][2];
    */

#if 1
    auto dBr_dr = ( dpsi_dz * over_r2 - d2psi_drdz * over_r) * bp_sign;
    auto dBr_dz = -d2psi_d2z * over_r * bp_sign;
    auto dBr_dp = 0.0 * bp_sign;

    auto dBz_dr = (-dpsi_dr * over_r2 + d2psi_d2r  * over_r) * bp_sign;
    auto dBz_dz = d2psi_drdz * over_r * bp_sign;
    auto dBz_dp = 0.0 * bp_sign;

    auto dBp_dr = dpsi_dr * fld_dIdpsi * over_r - fld_I * over_r2;
    auto dBp_dz = fld_dIdpsi * dpsi_dz * over_r;
    //auto dBp_dp = jacobian_rzp[2][2];
#endif

    //calculate curl_B
    /*
    ! R
    curl_B(1)  = fld%dbzdp*inv_r - fld%dbpdz
    ! Z
    if (sml_cylindrical) then
      curl_B(2)  = fld%dbpdr-fld%dbrdp*inv_r
    else
      curl_B(2)  = fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r
    endif
    ! phi
    curl_B(3)  = fld%dbrdz - fld%dbzdr
    */
    //vtkm::Vec3f curlB_rzp;
    pInfo.curlB_rzp[0] = dBz_dp * over_r - dBp_dz;
    pInfo.curlB_rzp[1] = Bp*over_r + dBp_dr - dBr_dp*over_r;
    pInfo.curlB_rzp[2] = dBr_dz - dBz_dr;
    //std::cout<<"curl_B_rzp= "<<curlB_rzp<<std::endl;

    //calculate curl_nb
    /*
    !curl of norm b
    curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
    curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
    curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2
    */

    vtkm::FloatDefault Bmag = vtkm::Magnitude(pInfo.B0_rzp);
    vtkm::FloatDefault over_B = 1./Bmag, over_B2 = over_B*over_B;

    /*
    dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
    dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
    dbdphi=0D0  ! no B perturbation
    */

    vtkm::FloatDefault dBdr, dBdz; /*, dBdp = 0;*/
    dBdr = (Br*dBr_dr + Bp*dBp_dr + Bz*dBz_dr) * over_B;
    dBdz = (Br*dBr_dz + Bp*dBp_dz + Bz*dBz_dz) * over_B;

#if 0
    //Check divergence.
    auto divergence = dBr_dr + Br/R + dBz_dz;
#ifdef VTKM_USE_DOUBLE_PRECISION
    static const vtkm::FloatDefault divEps = 1e-12;
#else
    static const vtkm::FloatDefault divEps = 1e-8;
#endif
    if (vtkm::Abs(divergence) > divEps)
    {
      std::cout<<std::endl;
      std::cout<<"****************************************************** DIVERGENCE= "<<divergence<<std::endl;
      std::cout<<std::endl;
    }
#endif

    //vtkm::Vec3f curl_nb_rzp;
    pInfo.curl_nb_rzp[0] = pInfo.curlB_rzp[0] * over_B + ( Bp * dBdz)*over_B2;
    pInfo.curl_nb_rzp[1] = pInfo.curlB_rzp[1] * over_B + (-Bp * dBdr)*over_B2;
    pInfo.curl_nb_rzp[2] = pInfo.curlB_rzp[2] * over_B + (Bz*dBdr - Br*dBdz)*over_B2;

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType, typename BFieldType, typename PsiType, typename OutputType, typename OutputType2D, typename IdType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            vtkm::Particle& particle,
                            const LocatorType& locator,
                            const CellSetType& cellSet,
                            const CoordsType& coords,
                            const AsFieldType& AsPhiFF,
                            const DAsFieldType& DAsPhiFF_RZP,
                            const Coeff_1DType& Coeff_1D,
                            const Coeff_2DType& Coeff_2D,
                            const BFieldType& B_RZP,
                            const PsiType& Psi,
                            OutputType& traces,
                            OutputType2D& outputRZ,
                            OutputType2D& outputTP,
                            IdType punctureID) const
  {
#ifdef VALGRIND
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_TOGGLE_COLLECT;
#endif

    DBG("Begin: "<<particle<<std::endl);

    ParticleInfo pInfo;

    if (this->ValidateInterpolation)
    {
      vtkm::Id numPts = coords.GetNumberOfValues();

      const vtkm::FloatDefault epsPsi = 1e-13, epsB0 = 1e-12;
      for (vtkm::Id i = 0; i < numPts; i += this->ValidateInterpolationSkip)
      {
        auto ptRZ = coords.Get(i);
        vtkm::Vec3f ptRPZ(ptRZ[0], 0, ptRZ[1]);
        this->HighOrderB(ptRPZ, pInfo, Coeff_1D, Coeff_2D);

        //Compare interpolated values at nodes to the array values.
        //The values in the array are computed by XGC.
        auto dPsi = vtkm::Abs(pInfo.Psi - Psi.Get(i));
        auto dB0 = vtkm::Magnitude(pInfo.B0_rzp - B_RZP.Get(i));
        if (dPsi > epsPsi)
        {
          printf("%lld: Psi difference detected. Error= %16.15lf\n", i, dPsi);
          printf("    Vals= %16.15lf %16.15lf\n", pInfo.Psi, Psi.Get(i));
        }
        if (dB0 > epsB0)
        {
          printf("%lld: B0 difference detected. Error= %16.15lf\n", i, dB0);
          printf("   Vals= (%16.15lf, %16.15lf, %16.15lf)\n", pInfo.B0_rzp[0], pInfo.B0_rzp[1], pInfo.B0_rzp[2]);
          printf("   Vals= (%16.15lf, %16.15lf, %16.15lf)\n", B_RZP.Get(i)[0], B_RZP.Get(i)[1], B_RZP.Get(i)[2]);
        }
      }

      return;
    }

    while (true)
    {
      vtkm::Vec3f newPos;
      DBG("\n\n\n*********************************************"<<std::endl);
      DBG("   "<<particle.Pos<<" #s= "<<particle.NumSteps<<std::endl);


      if (!this->TakeRK4Step(particle.Pos, pInfo, locator, cellSet, coords,
                             AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, B_RZP, newPos))
      {
        break;
      }

      DBG("     *** Step--> "<<newPos<<std::endl);
      vtkm::Id numRevs0 = vtkm::Floor(vtkm::Abs(particle.Pos[1] / vtkm::TwoPi()));
      vtkm::Id numRevs1 = vtkm::Floor(vtkm::Abs(newPos[1] / vtkm::TwoPi()));

      particle.Pos = newPos;
      particle.NumSteps++;

      if (this->SaveTraces)
        traces.Set(idx*this->MaxIter + particle.NumSteps, particle.Pos);

      //std::cout<<std::setprecision(12)<<" Step: "<<particle.Pos<<std::endl;
      if (numRevs1 > numRevs0)
      {
        auto R = particle.Pos[0], Z = particle.Pos[2];
        auto theta = vtkm::ATan2(Z-this->EqAxisZ, R-this->EqAxisR);
        if (theta < 0)
          theta += vtkm::TwoPi();

        //calcualte psi. need to stash psi on the particle somehow....
        vtkm::Vec3f ptRPZ = particle.Pos;
        //vtkm::FloatDefault psi;
        //vtkm::Vec3f B0_rzp, curlB_rzp, curl_nb_rzp, gradPsi_rzp;
        //vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;
        if (this->UseLinearB)
        {
          vtkm::Vec3f ptRZ(R, Z, 0), param;
          vtkm::Vec<vtkm::Id,3> vids;
          this->PtLoc(ptRZ, pInfo, locator, cellSet, coords, param, vids);
          pInfo.Psi = this->EvalS(Psi, 0, vids, param);
        }
        else
        {
          this->HighOrderB(ptRPZ, pInfo, Coeff_1D, Coeff_2D);
        }

        vtkm::Id i = (idx * this->MaxPunc) + particle.NumPunctures;
        outputRZ.Set(i, vtkm::Vec2f(R, Z));
        outputTP.Set(i, vtkm::Vec2f(theta, pInfo.Psi/this->EqXPsi));
        punctureID.Set(i, idx);
        particle.NumPunctures++;

#if !defined(VTKM_CUDA) && !defined(VTKM_HIP)
        if (idx == 0 && particle.NumPunctures%10 == 0 ) std::cout<<" ***** PUNCTURE n= "<<particle.NumPunctures<<std::endl;
#endif
        DBG("************* PUNCTURE n= "<<particle.NumPunctures<<std::endl);
      }

      if (particle.NumSteps >= this->MaxIter || particle.NumPunctures >= this->MaxPunc)
      {
#if !defined(VTKM_CUDA) && !defined(VTKM_HIP)
        std::cout<<"************************************* All done: id= "<<particle.ID<<" #Punc= "<<particle.NumPunctures<<std::endl;
#endif
        break;
      }
    }

#if !defined(VTKM_CUDA) && !defined(VTKM_HIP)
    //std::cout<<"Particle done: "<<idx<<" "<<particle<<std::endl;
#endif

#ifdef VALGRIND
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
#endif
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType, typename BFieldType>
  VTKM_EXEC
  bool TakeRK4Step(const vtkm::Vec3f& ptRPZ,
                   ParticleInfo& pInfo,
                   const LocatorType& locator,
                   const CellSetType& cellSet,
                   const CoordsType& coords,
                   const AsFieldType& AsPhiFF,
                   const DAsFieldType& DAsPhiFF_RZP,
                   const Coeff_1DType& Coeff_1D,
                   const Coeff_2DType& Coeff_2D,
                   const BFieldType& B_RZP,
                   vtkm::Vec3f& res) const
  {
    vtkm::Vec3f k1, k2, k3, k4;

    //k1 = F(p)
    //k2 = F(p+hk1/2)
    //k3 = F(p+hk2/2)
    //k4 = F(p+hk3)
    //Yn+1 = Yn + 1/6 h (k1+2k2+2k3+k4)
    vtkm::Vec3f p0 = ptRPZ, tmp = ptRPZ;

    DBG("    ****** K1"<<std::endl);
    bool v1, v2, v3, v4;
    v1 = this->Evaluate(tmp, pInfo, locator, cellSet, coords, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, B_RZP, k1);
    tmp = p0 + k1*this->StepSize_2;

    DBG("    ****** K2"<<std::endl);
    v2 = this->Evaluate(tmp, pInfo, locator, cellSet, coords, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, B_RZP, k2);
    tmp = p0 + k2*this->StepSize_2;

    DBG("    ****** K3"<<std::endl);
    v3 = this->Evaluate(tmp, pInfo, locator, cellSet, coords, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, B_RZP, k3);
    tmp = p0 + k3*this->StepSize;

    DBG("    ****** K4"<<std::endl);
    v4 = this->Evaluate(tmp, pInfo, locator, cellSet, coords, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, B_RZP, k4);

    vtkm::Vec3f vec = (k1 + 2*k2 + 2*k3 + k4)/6.0;
    res = p0 + this->StepSize * vec;

    if (!(v1&&v2&&v3&&v4))
    {
      //printf("RK4 step failed\n");
      return false;
    }

    return true;
  }

  template <typename PortalType>
  VTKM_EXEC
  vtkm::FloatDefault
  EvalS(const PortalType& sPortal,
        const vtkm::Id& offset,
        const vtkm::Vec<vtkm::Id, 3>& vId,
        const vtkm::Vec3f& param) const
  {
    vtkm::FloatDefault s;
#if 1
    //Hard code...
    const auto& v0 = sPortal.Get(vId[0]+offset);
    const auto& v1 = sPortal.Get(vId[1]+offset);
    const auto& v2 = sPortal.Get(vId[2]+offset);

    const auto& w1 = param[0];
    const auto& w2 = param[1];
    const auto w0 = 1 - (w1+w2);
    s = v0*w0 + v1*w1 + v2*w2;
#else
    vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
    vals.Append(sPortal.Get(vId[0]+offset));
    vals.Append(sPortal.Get(vId[1]+offset));
    vals.Append(sPortal.Get(vId[2]+offset));
    vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), s);
#endif
    return s;
  }

  template <typename PortalType>
  VTKM_EXEC
  vtkm::Vec3f
  EvalV(const PortalType& vPortal,
        const vtkm::Id& offset,
        const vtkm::Vec3f& param,
        const vtkm::Vec<vtkm::Id, 3>& vId) const
  {
    vtkm::Vec3f v;
#if 1
    //Hard code...
    const auto& v0 = vPortal.Get(vId[0]+offset);
    const auto& v1 = vPortal.Get(vId[1]+offset);
    const auto& v2 = vPortal.Get(vId[2]+offset);

    const auto& w1 = param[0];
    const auto& w2 = param[1];
    const auto w0 = 1 - (w1+w2);

    v[0] = v0[0]*w0 + v1[0]*w1 + v2[0]*w2;
    v[1] = v0[1]*w0 + v1[1]*w1 + v2[1]*w2;
    v[2] = v0[2]*w0 + v1[2]*w1 + v2[2]*w2;
    //std::cout<<std::setprecision(16)<<"EvalV: "<<v<<" "<<v2<<" "<<vtkm::Magnitude(v-v2)<<std::endl;
#else
    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
    vals.Append(vPortal.Get(vId[0]+offset));
    vals.Append(vPortal.Get(vId[1]+offset));
    vals.Append(vPortal.Get(vId[2]+offset));
    vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), v);
#endif

    return v;
  }

  VTKM_EXEC inline vtkm::Id ComputeFlatIndex(const DimVec3& idx, const DimVec3 dim) const
  {
    return idx[0] + (dim[0] * (idx[1] + (dim[1] * idx[2])));
  }
  VTKM_EXEC inline vtkm::exec::Grid ComputeLeafGrid(const DimVec3& idx, const DimVec3& dim, const vtkm::exec::Grid& l1Grid) const
  {
    return { dim,
             0,
             l1Grid.Origin + (static_cast<FloatVec3>(idx) * l1Grid.BinSize),
             l1Grid.BinSize / static_cast<FloatVec3>(dim) };
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType>
  VTKM_EXEC
  vtkm::ErrorCode FindCell2(const vtkm::Vec3f& point,
                            const LocatorType& locator,
                            const CellSetType& cellSet,
                            const CoordsType& coords,
                            vtkm::Id& cellId,
                            vtkm::Vec3f& parametric,
                            vtkm::Id3& prevCell) const
  {
    using namespace vtkm::internal::cl_uniform_bins;

    cellId = -1;

    auto TopLevel = locator.GetTopLevel();
    DimVec3 binId3 = static_cast<DimVec3>((point - TopLevel.Origin) / TopLevel.BinSize);
    if (binId3[0] >= 0 && binId3[0] < TopLevel.Dimensions[0] && binId3[1] >= 0 &&
        binId3[1] < TopLevel.Dimensions[1] && binId3[2] >= 0 &&
        binId3[2] < TopLevel.Dimensions[2])
    {
      vtkm::Id binId = this->ComputeFlatIndex(binId3, TopLevel.Dimensions);

      auto ldim = locator.GetLeafDimensions().Get(binId);
      if (!ldim[0] || !ldim[1] || !ldim[2])
      {
        printf("2L CellNotFound:: %d\n", __LINE__);
        return vtkm::ErrorCode::CellNotFound;
      }

      auto leafGrid = this->ComputeLeafGrid(binId3, ldim, TopLevel);

      DimVec3 leafId3 = static_cast<DimVec3>((point - leafGrid.Origin) / leafGrid.BinSize);
      // precision issues may cause leafId3 to be out of range so clamp it
      leafId3 = vtkm::Max(DimVec3(0), vtkm::Min(ldim - DimVec3(1), leafId3));

      vtkm::Id leafStart = locator.GetLeafStartIndex().Get(binId);
      vtkm::Id leafId = leafStart + this->ComputeFlatIndex(leafId3, leafGrid.Dimensions);

      vtkm::Id start = locator.GetCellStartIndex().Get(leafId);
      vtkm::Id end = start + locator.GetCellCount().Get(leafId);
      for (vtkm::Id i = start; i < end; ++i)
      {
        vtkm::Id cid = locator.GetCellIds().Get(i);
        auto indices = cellSet.GetIndices(cid);
        auto pts = vtkm::make_VecFromPortalPermute(&indices, coords);
        vtkm::Vec3f pc;
        bool inside;
        PointInsideCell(point, cellSet.GetCellShape(cid), pts, pc, inside);
        if (inside)
        {
          prevCell[0] = start;
          prevCell[1] = end;
          prevCell[2] = cid;
//          std::cout<<"FindCell: s/e= "<<start<<" "<<end<<" i= "<<i<<" prev= "<<prevCell<<std::endl;
          cellId = cid;
          parametric = pc;
          return vtkm::ErrorCode::Success;
        }
      }
    }

    printf("2L CellNotFound:: %d\n", __LINE__);
    return vtkm::ErrorCode::CellNotFound;
  }

  template <typename PointsVecType>
  VTKM_EXEC inline bool InCellBounds(const PointsVecType& points, const vtkm::Vec3f& pt) const
  {
    vtkm::Vec2f minP, maxP;
    minP[0] = maxP[0] = points[0][0];
    minP[1] = maxP[1] = points[0][1];

    //X,Y
    for (vtkm::IdComponent i = 0; i < 2; i++)
    {
      if (points[1][i] < minP[i]) minP[i] = points[1][i];
      if (points[1][i] > maxP[i]) maxP[i] = points[1][i];

      if (points[2][i] < minP[i]) minP[i] = points[2][i];
      if (points[2][i] > maxP[i]) maxP[i] = points[2][i];
    }

    return (pt[0] >= minP[0] && pt[0] <= maxP[0] &&
            pt[1] >= minP[1] && pt[1] <= maxP[1]);
  }

  vtkm::FloatDefault dot(const vtkm::Vec3f& a, const vtkm::Vec3f& b) const
  {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  template <typename CellShapeTag, typename CoordsType>
  VTKM_EXEC bool PointInsideCell(const vtkm::Vec3f& point,
                                 const CellShapeTag& /*cellShape*/,
                                 const CoordsType& cellPoints,
                                 vtkm::Vec3f& parametricCoordinates,
                                 bool& inside) const
  {
    inside = this->InCellBounds(cellPoints, point);
    if (inside)
    {
      parametricCoordinates[2] = 0;
      vtkm::Vec3f norm(0,0,1); //2D, so always in Z.
      const auto& p0 = cellPoints[0];
      for (vtkm::IdComponent i = 0; i < 2; i++)
      {
        vtkm::Vec3f p1 = cellPoints[i+1];
        vtkm::Vec3f p2 = cellPoints[2-i];
        auto plnN = vtkm::Cross(norm, p2-p0);
        parametricCoordinates[i] = vtkm::Dot(point - p0, plnN) /  vtkm::Dot(p1-p0, plnN);
      }

      inside = (parametricCoordinates[0] >= 0 &&
                parametricCoordinates[1] >= 0 &&
                (parametricCoordinates[0] + parametricCoordinates[1]) <= 1);

      //vtkm::exec::WorldCoordinatesToParametricCoordinates(cellPoints, point, cellShape, parametricCoordinates);
      //inside = vtkm::exec::CellInside(parametricCoordinates, cellShape);
    }
    return true;
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType>
  VTKM_EXEC
  bool FindCell(const vtkm::Vec3f& ptRZ,
                ParticleInfo& pInfo,
                const LocatorType& locator,
                const CellSetType& cellSet,
                const CoordsType& coords,
                vtkm::Vec3f& param,
                vtkm::Id& cellId) const
  {
    cellId = -1;

    //prevCell: (start, end, prevID)
    auto prevCell = pInfo.PrevCell;

    if (prevCell[0] != -1)
    {
      //Check last cellId first.
      auto indices = cellSet.GetIndices(prevCell[2]);
      auto pts = vtkm::make_VecFromPortalPermute(&indices, coords);
      vtkm::Vec3f pc;
      bool inside;
      this->PointInsideCell(ptRZ, cellSet.GetCellShape(prevCell[2]), pts, pc, inside);
      if (inside)
      {
        cellId = prevCell[2];
//        std::cout<<"  Quick:: FindCell: "<<prevCell<<"  cid= "<<cellId<<std::endl;
        param = pc;
        return true;
      }
//      std::cout<<"    QuickFail "<<prevCell<<std::endl;
      //Check the rest of the cells in the bin, excluding the cell we just checked.
      for (vtkm::Id i = prevCell[0]; i < prevCell[1]; i++)
      {
        vtkm::Id cid = locator.GetCellIds().Get(i);
        if (cid != prevCell[2])
        {
          auto indices2 = cellSet.GetIndices(cid);
          auto pts2 = vtkm::make_VecFromPortalPermute(&indices2, coords);
          this->PointInsideCell(ptRZ, cellSet.GetCellShape(prevCell[2]), pts2, pc, inside);
          if (inside)
          {
//            std::cout<<"  "<<(i-prevCell[0])<<" -FindCell: "<<prevCell<<"  cid= "<<cid<<" len= "<<prevCell[1]-prevCell[0]<<std::endl;
            prevCell[2] = cid;
            cellId = cid;
            param = pc;
            return true;
          }
        }
      }
    }

//    std::cout<<" ***** No find: "<<prevCell<<std::endl;

    return (this->FindCell2(ptRZ, locator, cellSet, coords, cellId, param, prevCell) == vtkm::ErrorCode::Success);
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType>
  VTKM_EXEC
  bool PtLoc2(const vtkm::Vec3f& ptRZ,
              ParticleInfo& pInfo,
              const LocatorType& locator,
              const CellSetType& cs,
              const CoordsType& coords,
              vtkm::Vec3f& param,
              vtkm::Vec<vtkm::Id, 3>& vIds) const
  {
    vtkm::Id cellId;
    if (!this->FindCell(ptRZ, pInfo, locator, cs, coords, param, cellId))
      return false;

    auto tmp =  cs.GetIndices(cellId);
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType>
  VTKM_EXEC
  bool PtLoc(const vtkm::Vec3f& ptRZ,
             ParticleInfo& pInfo,
             const LocatorType& locator,
             const CellSetType& cs,
             const CoordsType& coords,
             vtkm::Vec3f& param,
             vtkm::Vec<vtkm::Id, 3>& vIds) const
  {
    vtkm::Id cellId;
    vtkm::ErrorCode status = locator.FindCell(ptRZ, cellId, param, pInfo.PrevCell);
    if (status != vtkm::ErrorCode::Success)
    {
      //printf("Find Cell failed! pt= %lf %lf %lf\n", ptRZ[0], ptRZ[1], ptRZ[2]);
      return false;
    }

    if (0)
    {
      vtkm::Vec3f p2;
      vtkm::Id cid2;

      ParticleInfo p2Info = pInfo;
      p2Info.PrevCell = vtkm::Id3(-1,-1,-1);
      this->FindCell(ptRZ, p2Info, locator, cs, coords, p2, cid2);
      //if (cid2 != cellId) std::cout<<" **************************** WRONG"<<std::endl;
      //if (vtkm::Magnitude(p2-param) > 1e-10) std::cout<<" **************************** WRONG"<<std::endl;
    }

    //vtkm::VecVariable<vtkm::Id, 3> tmp;
    auto tmp =  cs.GetIndices(cellId);
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];

    return true;
  }

  VTKM_EXEC
  int GetIndex(const vtkm::FloatDefault& x,
               const int& nx,
               const vtkm::FloatDefault& xmin,
               const vtkm::FloatDefault& dx_inv) const
  {
    int idx = std::max(0, std::min(nx-1,
                                   int((x-xmin)*dx_inv)) );
    return idx;

  }

  template <typename Coeff_2DType>
  VTKM_EXEC
  void EvalBicub2(const vtkm::FloatDefault& x,
                  const vtkm::FloatDefault& y,
                  const vtkm::FloatDefault& xc,
                  const vtkm::FloatDefault& yc,
                  const vtkm::Id& offset,
                  const Coeff_2DType& Coeff_2D,
                  vtkm::FloatDefault &f00, vtkm::FloatDefault &f10, vtkm::FloatDefault &f01,
                  vtkm::FloatDefault &f11, vtkm::FloatDefault &f20, vtkm::FloatDefault &f02) const
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

    /*
    for (int j = 0; j < 4; j++)
    {
      std::cout<<"acoeff_"<<j<<": ";
      for (int i = 0; i < 4; i++)
        std::cout<<Coeff_2D.Get(offset + j*4 + i)<<" ";
      std::cout<<std::endl;
    }
    */

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

  VTKM_EXEC
  bool eval_bicub_2(const vtkm::FloatDefault& x,
                    const vtkm::FloatDefault& y,
                    const vtkm::FloatDefault& xc,
                    const vtkm::FloatDefault& yc,
                    const vtkm::Matrix<vtkm::FloatDefault, 4, 4>& acoeff,
                    vtkm::FloatDefault &f00, vtkm::FloatDefault &f10, vtkm::FloatDefault &f01,
                    vtkm::FloatDefault &f11, vtkm::FloatDefault &f20, vtkm::FloatDefault &f02) const
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
        fx[j] = fx[j] + xv[i]*acoeff[i][j];
      for (int i=1; i<4; i++)
        dfx[j] = dfx[j] + vtkm::FloatDefault(i)*xv[i-1]*acoeff[i][j];
      for (int i=2; i<4; i++)
        dfx2[j] = dfx2[j] + vtkm::FloatDefault(i*(i-1))*xv[i-2]*acoeff[i][j];
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






    //c++ code.
    /*
    double fx_i, dfx_i, dfx2_i;

    f00 = 0;
    f01 = 0;
    f02 = 0;

    fx_i = ((acoeff[0][3]*dx + acoeff[0][2])*dx + acoeff[0][1])*dx + acoeff[0][0];
    f00 = f00 + fx_i;

    fx_i = ((acoeff[1][3]*dx + acoeff[1][2])*dx + acoeff[1][1])*dx + acoeff[1][0];
    f00 = f00 + dy*fx_i;
    f01 = f01 +    fx_i;

    fx_i = ((acoeff[2][3]*dx + acoeff[2][2])*dx + acoeff[2][1])*dx + acoeff[2][0];
    f00 = f00 + dy*(dy*fx_i);
    f01 = f01 + 2.0*(dy*fx_i);
    f02 = f02 + 2.0*fx_i;

    fx_i = ((acoeff[3][3]*dx + acoeff[3][2])*dx + acoeff[3][1])*dx + acoeff[3][0];
    f00 = f00 + dy*(dy*(dy*fx_i));
    f01 = f01 + 3.0*(dy*(dy*fx_i));
    f02 = f02 + 6.0*(dy*fx_i);

    f10 = 0;
    f11 = 0;

    dfx_i = (acoeff[0][3]*3.0*dx + acoeff[0][2]*2.0)*dx + acoeff[0][1];
    f10 = f10 + dfx_i;

    dfx_i = (acoeff[1][3]*3.0*dx + acoeff[1][2]*2.0)*dx + acoeff[1][1];
    f10 = f10 + dy*dfx_i;
    f11 = f11 +    dfx_i;

    dfx_i = (acoeff[2][3]*3.0*dx + acoeff[2][2]*2.0)*dx + acoeff[2][1];
    f10 = f10 + dy*(dy*dfx_i);
    f11 = f11 + 2.0*(dy*dfx_i);

    dfx_i = (acoeff[3][3]*3.0*dx + acoeff[3][2]*2.0)*dx + acoeff[3][1];
    f10 = f10 + dy*(dy*dy*dfx_i);
    f11 = f11 + 3.0*(dy*dy*dfx_i);

    f20 = 0;

    dfx2_i = (3.*2.*acoeff[0][3]*dx + 2.*1.*acoeff[0][2]);
    f20 = f20 + dfx2_i;

    dfx2_i = (3.*2.*acoeff[1][3]*dx + 2.*1.*acoeff[1][2]);
    f20 = f20 + dy*dfx2_i;

    dfx2_i = (3.*2.*acoeff[2][3]*dx + 2.*1.*acoeff[2][2]);
    f20 = f20 + (dy*dy)*dfx2_i;

    dfx2_i = (3.*2.*acoeff[3][3]*dx + 2.*1.*acoeff[3][2]);
    f20 = f20 + (dy*dy)*dy*dfx2_i;
    */
    //printf("  worklet %d\n", __LINE__);
    return true;
  }

  template <typename Coeff_1DType>
  VTKM_EXEC
  vtkm::FloatDefault I_interpol(const vtkm::FloatDefault& psi,
                                const int& ideriv,
                                const Coeff_1DType& coeff_1D) const
  {
    vtkm::FloatDefault pn = psi * this->one_d_cub_dpsi_inv;
    int ip=floor(pn);
    ip=std::min(std::max(ip,0),this->ncoeff-1);
    vtkm::FloatDefault wp=pn-(vtkm::FloatDefault)(ip);

    int idx = ip*4;

    //vtkm::FloatDefault acoef[4];
    //acoef[0] = one_d_cub_acoef(ip).coeff[0];
    //acoef[1] = one_d_cub_acoef(ip).coeff[1];
    //acoef[2] = one_d_cub_acoef(ip).coeff[2];
    //acoef[3] = one_d_cub_acoef(ip).coeff[3];

    const vtkm::FloatDefault acoef[4] = {coeff_1D.Get(idx+0),
                                         coeff_1D.Get(idx+1),
                                         coeff_1D.Get(idx+2),
                                         coeff_1D.Get(idx+3)};

    vtkm::FloatDefault iVal = 0.0;
    if (ideriv==0)
      iVal = acoef[0]+(acoef[1]+(acoef[2]+acoef[3]*wp)*wp)*wp;
    else if (ideriv==1)
      iVal = (acoef[1]+(2.0*acoef[2]+3.0*acoef[3]*wp)*wp)*one_d_cub_dpsi_inv;

    return iVal * this->sml_bp_sign;
  }

  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  vtkm::Vec3f GetB(vtkm::Vec3f& pt_rpz,
                   ParticleInfo& /*pInfo*/,
                   const Coeff_1DType& coeff_1D,
                   const Coeff_2DType& coeff_2D) const
  {
#if 0
    this->HighOrderB(pt_rpz, pInfo, coeff_1D, coeff_2D);
    //This gives is the time derivative: Br = dR/dt, Bz= dZ/dt, B_phi/R = dphi/dt
    //We need with respect to phi:
    // dR/dphi = dR/dt / (dphi/dt) = Br / B_phi * R
    // same for z;
    vtkm::Vec3f B0_rpz(pInfo.B0_rzp[0], pInfo.B0_rzp[2], pInfo.B0_rzp[1]);
    B0_rpz[0] /= pInfo.B0_rzp[2];
    B0_rpz[2] /= pInfo.B0_rzp[2];
    B0_rpz[0] *= pt_rpz[0];
    B0_rpz[2] *= pt_rpz[0];
    return B0_rpz;
#else
    auto B0_rzp = this->HighOrderBOnly(pt_rpz, coeff_1D, coeff_2D);

    vtkm::Vec3f B0_rpz(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

    //This gives is the time derivative: Br = dR/dt, Bz= dZ/dt, B_phi/R = dphi/dt
    //We need with respect to phi:
    // dR/dphi = dR/dt / (dphi/dt) = Br / B_phi * R
    // same for z;

    B0_rpz[0] /= B0_rzp[2];
    B0_rpz[2] /= B0_rzp[2];
    B0_rpz[0] *= pt_rpz[0];
    B0_rpz[2] *= pt_rpz[0];
    return B0_rpz;
#endif
  }

  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool CalcFieldFollowingPt(const vtkm::Vec3f& pt_rpz,
                            const ParticleInfo& pInfo,
                            const vtkm::Vec3f& B0_rpz,
                            const vtkm::FloatDefault& Phi0,
                            const vtkm::FloatDefault& Phi1,
                            const Coeff_1DType& coeff_1D,
                            const Coeff_2DType& coeff_2D,
                            vtkm::Vec3f& x_ff_rpz) const
  {
    vtkm::FloatDefault R = pt_rpz[0];
    vtkm::FloatDefault Phi = pt_rpz[1];
    vtkm::FloatDefault Z = pt_rpz[2];

    vtkm::FloatDefault PhiMid = Phi0 + (Phi1-Phi0)/2.0;
    vtkm::Plane<> midPlane({0, PhiMid, 0}, {0,1,0});
    Ray3f ray_rpz({R, Phi, Z}, B0_rpz);

    //Get point on mid plane.  Use the R,Z for this point for triangle finds.
    vtkm::FloatDefault RP_T;
    vtkm::Vec3f ptOnMidPlane_rpz;
    bool b;
    midPlane.Intersect(ray_rpz, RP_T, ptOnMidPlane_rpz, b);

    //Now, do it using RK4 and two steps.
    vtkm::FloatDefault h = (PhiMid-Phi) / 2.0;
    //h = -h;
    vtkm::FloatDefault h_2 = h / 2.0;

    //k1 = F(p)
    //k2 = F(p+hk1/2)
    //k3 = F(p+hk2/2)
    //k4 = F(p+hk3)
    //Yn+1 = Yn + 1/6 h (k1+2k2+2k3+k4)
    vtkm::Vec3f p0 = {R,Phi,Z}; //pt_rpz;
    vtkm::Vec3f tmp, k1, k2, k3, k4;
    //std::cout<<"     p0 = "<<p0<<std::endl;
    ParticleInfo pInfo2 = pInfo;
    for (int i = 0; i < 2; i++)
    {
      k1 = this->GetB(p0, pInfo2, coeff_1D, coeff_2D);
      tmp = p0 + k1*h_2;

      k2 = this->GetB(tmp, pInfo2, coeff_1D, coeff_2D);
      tmp = p0 + k2*h_2;

      k3 = this->GetB(tmp, pInfo2, coeff_1D, coeff_2D);
      tmp = p0 + k3*h;

      k4 = this->GetB(tmp, pInfo2, coeff_1D, coeff_2D);

      vtkm::Vec3f vec = (k1 + 2*k2 + 2*k3 + k4) / 6.0;
      p0 = p0 + h * vec;
    }

    x_ff_rpz = p0;

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename CoordsType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType, typename BFieldType>
  VTKM_EXEC
  bool Evaluate(const vtkm::Vec3f& ptRPZ,
                ParticleInfo& pInfo,
                const LocatorType& locator,
                const CellSetType& cellSet,
                const CoordsType& coords,
                const AsFieldType& AsPhiFF,
                const DAsFieldType& DAsPhiFF_RZP,
                const Coeff_1DType& coeff_1D,
                const Coeff_2DType& coeff_2D,
                const BFieldType& B_RZP,
                vtkm::Vec3f& res) const
  {
    auto R = ptRPZ[0];
    auto Phi = ptRPZ[1];
    auto Z = ptRPZ[2];

    //res is R,P,Z
    if (this->UseLinearB)
    {
      vtkm::Vec3f ptRZ(ptRPZ[0], ptRPZ[2], 0), param;
      vtkm::Vec<vtkm::Id,3> vids;
      this->PtLoc(ptRZ, pInfo, locator, cellSet, coords, param, vids);
      pInfo.B0_rzp = this->EvalV(B_RZP, 0, param, vids);
    }
    else
    {
      if (!this->HighOrderB(ptRPZ, pInfo, coeff_1D, coeff_2D))
        return false;
    }
    if (this->UseBOnly)
    {
      pInfo.B0_rzp[2] /= R;
      //res is R,P,Z
      res = vtkm::Vec3f(pInfo.B0_rzp[0], pInfo.B0_rzp[2], pInfo.B0_rzp[1]);
      return true;
    }

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    this->GetPlaneIdx(Phi, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);
    if (planeIdx0 == 48)
    {
      printf("Problem.... plane wraparound\n");
    }
    vtkm::Vec3f B0_rpz(pInfo.B0_rzp[0], pInfo.B0_rzp[2], pInfo.B0_rzp[1]);
    vtkm::Vec3f ff_pt_rpz;
    this->CalcFieldFollowingPt({R,phiN,Z}, pInfo, B0_rpz, Phi0, Phi1, coeff_1D, coeff_2D, ff_pt_rpz);

    //Now, interpolate between Phi_i and Phi_i+1
    vtkm::FloatDefault T01 = (phiN - Phi0) / (Phi1-Phi0);
    vtkm::FloatDefault T10 = 1.0f - T01;

    //Get vec at Phi0 and Phi1.
    //x_ff is in rzp
    //vtkm::Vec3f x_ff_rzp(ptOnMidPlane_rpz[0], ptOnMidPlane_rpz[2], 0);
    vtkm::Vec3f x_ff_rzp(ff_pt_rpz[0], ff_pt_rpz[2], 0);

    int offsets[2];
    offsets[0] = planeIdx0*this->NumNodes*2;
    offsets[1] = planeIdx0*this->NumNodes*2 + this->NumNodes;

    const vtkm::FloatDefault basis = 0.0f;
    //auto B0_R = B0_rzp[0];
    //auto B0_Z = B0_rzp[1];
    //auto x_ff_R = x_ff_rzp[0];
    //auto x_ff_Z = x_ff_rzp[1];

    //gradPsi: pt on mid plane?  (question)
    //dPsi/dR = B0_Z * R
    //dPsi/dZ = -B0_R * R;
    //vtkm::Vec3f gradPsi_rzp(B0_Z * x_ff_R, -B0_R * x_ff_R, 0);
    //use high order...
    vtkm::Vec3f gradPsi_rzp = pInfo.gradPsi_rzp;
    vtkm::FloatDefault gammaPsi = 1.0f/vtkm::Magnitude(gradPsi_rzp);

    vtkm::Vec2f rvec(0,0), zvec(0,0);
    rvec[0] = basis + (1.0-basis) * gammaPsi *   gradPsi_rzp[0];
    rvec[1] =         (1.0-basis) * gammaPsi *   gradPsi_rzp[1];
    zvec[0] =         (1.0-basis) * gammaPsi * (-gradPsi_rzp[1]);
    zvec[1] = basis + (1.0-basis) * gammaPsi *   gradPsi_rzp[0];

    //Get the vectors in the ff coordinates.
    //auto dAs_ff_rzp = EvalVector(ds, locator, {x_ff_rzp, x_ff_rzp}, "dAs_ff_rzp", offsets);
    //auto dAs_ff0_rzp = dAs_ff_rzp[0];
    //auto dAs_ff1_rzp = dAs_ff_rzp[1];

    vtkm::Vec3f x_ff_param;
    vtkm::Vec<vtkm::Id,3> x_ff_vids;

    if (!this->PtLoc(x_ff_rzp, pInfo, locator, cellSet, coords, x_ff_param, x_ff_vids))
      return false;

    //this->PtLoc2(x_ff_rzp, pInfo, locator, cellSet, coords, x_ff_param, x_ff_vids);
    auto dAs_ff0_rzp = this->EvalV(DAsPhiFF_RZP, offsets[0], x_ff_param, x_ff_vids);
    auto dAs_ff1_rzp = this->EvalV(DAsPhiFF_RZP, offsets[1], x_ff_param, x_ff_vids);

    vtkm::FloatDefault wphi[2] = {T10, T01}; //{T01, T10};
    vtkm::Vec3f gradAs_rpz;

    //vec.r = wphi[0]*( rvec[0]*V.r[0] + zvec[0]*V.z[0]) +
    //        wphi[1]*( rvec[0]*V.r[1] + zvec[0]*v.z[1]);
    //vec.p = wphi[0]*V.phi[0] +
    //        whpi[1]*V.phi[1];
    //vec.z = wphi[0]*( rvec[1]*V.r[0] + zvec[1]*V.z[0]) +
    //        wphi[1]*( rvec[1]*V.r[1] + zvec[1]*V.Z[1]);
    gradAs_rpz[0] = wphi[0]*(rvec[0]*dAs_ff0_rzp[0] + zvec[0]*dAs_ff0_rzp[1]) +
                    wphi[1]*(rvec[0]*dAs_ff1_rzp[0] + zvec[0]*dAs_ff1_rzp[1]);
    gradAs_rpz[1] = wphi[0] * dAs_ff0_rzp[2] +
                    wphi[1] * dAs_ff1_rzp[2];
    gradAs_rpz[2] = wphi[0]*(rvec[1]*dAs_ff0_rzp[0] + zvec[1]*dAs_ff0_rzp[1]) +
                    wphi[1]*(rvec[1]*dAs_ff1_rzp[0] + zvec[1]*dAs_ff1_rzp[1]);

    vtkm::FloatDefault BMag = vtkm::Magnitude(pInfo.B0_rzp);
    //project using bfield.
    //gradAs.Phi = (gradAs.Phi * BMag - gradAs.R*B0_pos.R - gradAs.Z*B0_pos.Z) / B0_pos.Phi
    gradAs_rpz[1] = (gradAs_rpz[1]*BMag -gradAs_rpz[0]*pInfo.B0_rzp[0] - gradAs_rpz[2]*pInfo.B0_rzp[1]) / pInfo.B0_rzp[2];

    //deltaB = AsCurl(bhat) + gradAs x bhat.
    //std::vector<int> off = {planeIdx0*this->NumNodes};
    //vtkm::Vec3f AsCurl_bhat_rzp = EvalVector(ds, locator, {x_ff_rzp}, "AsCurlBHat_RZP", off)[0];
    //auto AsCurl_bhat_rzp = this->EvalV(AsCurlBHat_RZP, 0, x_ff_vids, x_ff_param);
    auto As_ff0 = this->EvalS(AsPhiFF, offsets[0], x_ff_vids, x_ff_param);
    auto As_ff1 = this->EvalS(AsPhiFF, offsets[1], x_ff_vids, x_ff_param);

    vtkm::FloatDefault As = wphi[0]*As_ff0 + wphi[1]*As_ff1;
    auto AsCurl_bhat_rzp = As * pInfo.curl_nb_rzp;

    auto bhat_rzp = vtkm::Normal(pInfo.B0_rzp);

    vtkm::Vec3f gradAs_rzp(gradAs_rpz[0], gradAs_rpz[2], gradAs_rpz[1]);
    vtkm::Vec3f deltaB_rzp = AsCurl_bhat_rzp + vtkm::Cross(gradAs_rzp, bhat_rzp);

    deltaB_rzp[2] /= R;
    pInfo.B0_rzp[2] /= R;

    //std::cout<<"Evaluate: "<<ptRPZ<<" : psi= "<<pInfo.Psi<<" B0= "<<pInfo.B0_rzp<<std::endl;

    vtkm::Vec3f vec_rzp = pInfo.B0_rzp + deltaB_rzp;
    vtkm::Vec3f vec_rpz(vec_rzp[0], vec_rzp[2], vec_rzp[1]);
    res = vec_rpz;
    return true;
  }

  VTKM_EXEC void
  GetPlaneIdx(const vtkm::FloatDefault& phi,
              vtkm::FloatDefault& phiN,
              vtkm::Id& plane0,
              vtkm::Id& plane1,
              vtkm::FloatDefault& phi0,
              vtkm::FloatDefault& phi1,
              vtkm::Id& numRevs,
              vtkm::FloatDefault& T) const
  {
    numRevs = vtkm::Floor(vtkm::Abs(phi / vtkm::TwoPi()));
    phiN = phi;
    if (phi < 0)
    {
      phiN += ((1+numRevs) * vtkm::TwoPi());
    }

    plane0 = vtkm::Floor(phiN / this->dPhi);
    if (plane0 == this->NumPlanes)
      plane0 = 0;

    plane1 = plane0 + 1;
    phi0 = static_cast<vtkm::FloatDefault>(plane0)*this->dPhi;
    phi1 = static_cast<vtkm::FloatDefault>(plane1)*this->dPhi;

    if (plane1 == this->NumPlanes)
      plane1 = 0;
    T = (phiN-phi0) / (phi1-phi0);
  }

  vtkm::Id MaxIter = 0;
  vtkm::Id MaxPunc = 0;
  vtkm::FloatDefault PlaneVal = 0.0f;
  vtkm::FloatDefault StepSize;
  vtkm::FloatDefault StepSize_2;
  vtkm::FloatDefault StepSize_6;

  vtkm::Id NumPlanes;
  vtkm::Id NumNodes;
  vtkm::FloatDefault dPhi;

  bool UseBOnly = false;
  bool SaveTraces = false;
  bool UseLinearB = false;

  int nr, nz;
  vtkm::FloatDefault rmin, zmin, rmax, zmax;
  vtkm::FloatDefault dr, dz, dr_inv, dz_inv;
  vtkm::FloatDefault EqAxisZ, EqAxisR, EqXPsi;

  int ncoeff;
  vtkm::FloatDefault min_psi, max_psi;
  vtkm::FloatDefault one_d_cub_dpsi_inv;
  vtkm::FloatDefault sml_bp_sign = -1.0f;

  bool ValidateInterpolation = false;
  vtkm::Id ValidateInterpolationSkip = 1;
};

#endif
