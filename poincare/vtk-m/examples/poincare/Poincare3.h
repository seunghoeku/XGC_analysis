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
class PoincareWorklet3 : public vtkm::worklet::WorkletMapField
{
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
                                ExecObject locatorRZ,
                                ExecObject locatorB,
                                WholeCellSetIn<> cellSetRZ,
                                WholeCellSetIn<> cellSet,
                                WholeArrayIn B_RZP,
                                WholeArrayIn As_phi_ff,
                                WholeArrayIn dAs_phi_ff_RZP,
                                WholeArrayIn Psi,
                                WholeArrayIn B0_RZP,
                                WholeArrayIn CurlB0_RZP,
                                WholeArrayIn Curl_NB0_RZP,
                                WholeArrayIn GradPsi_RZP,
                                WholeArrayIn coeff_1D,
                                WholeArrayIn coeff_2D,
                                WholeArrayInOut traces,
                                WholeArrayInOut outputRZ,
                                WholeArrayInOut outputTP,
                                WholeArrayInOut punctureID);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19);
  using InputDomain = _1;

  PoincareWorklet3(vtkm::Id maxPunc,
                   vtkm::FloatDefault planeVal,
                   vtkm::FloatDefault stepSize,
                   bool saveTraces,
                   bool quickTest)
    : MaxIter(maxPunc * 100000000)
    , MaxPunc(maxPunc)
    , PlaneVal(planeVal)
    , StepSize(stepSize)
    , SaveTraces(saveTraces)
    , QuickTest(quickTest)
  {
    this->NumPlanes = numPlanes;
    this->NumNodes = numNodes;
    this->dPhi = vtkm::TwoPi()/static_cast<vtkm::FloatDefault>(this->NumPlanes);
    this->StepSize_2 = this->StepSize / 2.0;
    this->StepSize_6 = this->StepSize / 6.0;


    this->nr = eq_mr-1;
    this->nz = eq_mz-1;
    this->rmin = eq_min_r;
    this->rmax = eq_max_r;
    this->zmin = eq_min_z;
    this->zmax = eq_max_z;
    this->EqAxisR = eq_axis_r;
    this->EqAxisZ = eq_axis_z;
    this->EqXPsi = eq_x_psi;
    this->dr = (eq_max_r - eq_min_r) / vtkm::FloatDefault(this->nr);
    this->dz = (eq_max_z - eq_min_z) / vtkm::FloatDefault(this->nz);
    this->dr_inv = 1.0/this->dr;
    this->dz_inv = 1.0/this->dz;

    this->ncoeff = eq_mr-1;
    this->min_psi = psi_min;
    this->max_psi = psi_max;
    this->one_d_cub_dpsi_inv = 1.0 / ((max_psi-min_psi)/vtkm::FloatDefault(this->ncoeff));
    #ifndef VTKM_CUDA
    std::cout<<"PSI min/max= "<<psi_min<<" "<<psi_max<<std::endl;
    #endif
  }

  void SetBGrid(const vtkm::Vec2f& origin,
                const vtkm::Vec2f& spacing,
                const vtkm::Vec2f& maxPt,
                const vtkm::Id2& dims,
                const bool& bcell)
  {
    this->Origin = origin;
    this->Spacing = spacing;
    this->MaxPt = maxPt;
    this->Dims = dims;
    this->CellDims = {this->Dims[0]-1, this->Dims[1]-1};
    this->InvSpacing = {1/this->Spacing[0], 1/this->Spacing[1]};
    this->BCell = bcell;
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
    vtkm::Id offset = (r_i * this->ncoeff + z_i) * 16;
    //offset = (z_i * this->ncoeff + r_i)*16;
    vtkm::Id idx = 0;
    //std::cout<<"Offset= "<<(offset/16)<<" 16: "<<offset<<std::endl;
    for (vtkm::Id ii = 0; ii < 4; ii++)
      for (vtkm::Id jj = 0; jj < 4; jj++)
      {
        acoeff[ii][jj] = Coeff_2D.Get(offset+idx);
        idx++;
      }

    vtkm::FloatDefault psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
    this->eval_bicub_2(R, Z, Rc, Zc, acoeff, psi, dpsi_dr, dpsi_dz, d2psi_drdz, d2psi_d2r, d2psi_d2z);
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

    vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;
    //Set the jacobian.
    const int PIR = 0;
    const int PIZ = 1;
    const int PIP = 2;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        jacobian_rzp[i][j] = 0;
    vtkm::FloatDefault bp_sign = 1.0;

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

    #ifndef VTKM_CUDA
    //Check divergence.
    auto divergence = dBr_dr + Br/R + dBz_dz;

#ifdef VTKM_USE_DOUBLE_PRECISION
    static const vtkm::FloatDefault divEps = 1e-16;
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

  template <typename LocatorRZType, typename LocatorBType, typename CellSetRZType, typename CellSetBType, typename VecFieldType, typename AsFieldType, typename DAsFieldType, typename ScalarFieldType, typename Coeff_1DType, typename Coeff_2DType, typename OutputType, typename OutputType2D, typename IdType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            vtkm::Particle& particle,
                            const LocatorRZType& locatorRZ,
                            const LocatorBType& locatorB,
                            const CellSetRZType& cellSetRZ,
                            const CellSetBType& cellSetB,
                            const VecFieldType& B_RZP,
                            const AsFieldType& AsPhiFF,
                            const DAsFieldType& DAsPhiFF_RZP,
                            const ScalarFieldType& Psi,
                            const VecFieldType& Bcell_RZP,
                            const VecFieldType& Bcell_Curl_RZP,
                            const VecFieldType& Curlcell_NB_RZP,
                            const VecFieldType& GradPsicell_RZP,
                            const Coeff_1DType& Coeff_1D,
                            const Coeff_2DType& Coeff_2D,
                            OutputType& traces,
                            OutputType2D& outputRZ,
                            OutputType2D& outputTP,
                            IdType punctureID) const
  {
#ifdef VALGRIND
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_TOGGLE_COLLECT;
#endif

    if (this->QuickTest)
    {
      for (vtkm::Id p = 0; p < this->MaxPunc; p++)
      {
        vtkm::Id i = (idx * this->MaxPunc) + p;
        outputRZ.Set(i, vtkm::Vec2f(0,0));
        outputTP.Set(i, vtkm::Vec2f(0,0));
        punctureID.Set(i, idx);
      }

      return;
    }

    DBG("Begin: "<<particle<<std::endl);

    ParticleInfo pInfo;
    while (true)
    {
      vtkm::Vec3f newPos;
      DBG("\n\n\n*********************************************"<<std::endl);
      DBG("   "<<particle.Pos<<" #s= "<<particle.NumSteps<<std::endl);
      if (!this->TakeRK4Step(particle.Pos, pInfo,
                             locatorRZ, locatorB, cellSetRZ, cellSetB,
                             B_RZP, AsPhiFF, DAsPhiFF_RZP,
                             Psi, Bcell_RZP, Bcell_Curl_RZP, Curlcell_NB_RZP, GradPsicell_RZP,
                             Coeff_1D, Coeff_2D, newPos))
      {
#ifndef VTKM_CUDA
        std::cout<<"*****************       All done: RK step failed."<<std::endl;
#endif
        break;
      }

      DBG("     *** Step--> "<<newPos<<std::endl);
      vtkm::Id numRevs0 = vtkm::Floor(vtkm::Abs(particle.Pos[1] / vtkm::TwoPi()));
      vtkm::Id numRevs1 = vtkm::Floor(vtkm::Abs(newPos[1] / vtkm::TwoPi()));

      particle.Pos = newPos;
      particle.NumSteps++;

      if (this->SaveTraces)
        traces.Set(idx*this->MaxIter + particle.NumSteps, particle.Pos);

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
        this->HighOrderB(ptRPZ, pInfo, Coeff_1D, Coeff_2D);
        vtkm::Id i = (idx * this->MaxPunc) + particle.NumPunctures;
        outputRZ.Set(i, vtkm::Vec2f(R, Z));
        outputTP.Set(i, vtkm::Vec2f(theta, pInfo.Psi/this->EqXPsi));
        punctureID.Set(i, idx);
        particle.NumPunctures++;

#ifndef VTKM_CUDA
        if (idx == 0 && particle.NumPunctures%10 == 0 ) std::cout<<" ***** PUNCTURE n= "<<particle.NumPunctures<<std::endl;
#endif
        DBG("************* PUNCTURE n= "<<particle.NumPunctures<<std::endl);
      }
      //printf("  worklet %d\n", __LINE__);
      if (particle.NumSteps >= this->MaxIter || particle.NumPunctures >= this->MaxPunc)
      {
#ifndef VTKM_CUDA
        std::cout<<"************************************* All done: "<<particle<<std::endl;
#endif
        break;
      }
    }
#ifndef VTKM_CUDA
    std::cout<<"Particle done: "<<idx<<std::endl;
#endif
#ifdef VALGRIND
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
#endif
  }

  template <typename LocatorRZType, typename CellSetRZType, typename LocatorBType, typename CellSetBType, typename VecFieldType, typename AsFieldType, typename DAsFieldType, typename ScalarFieldType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool TakeRK4Step(const vtkm::Vec3f& ptRPZ,
                   ParticleInfo& pInfo,
                   const LocatorRZType& locatorRZ,
                   const LocatorBType& locatorB,
                   const CellSetRZType& cellSetRZ,
                   const CellSetBType& cellSetB,
                   const VecFieldType& B_RZP,
                   const AsFieldType& AsPhiFF,
                   const DAsFieldType& DAsPhiFF_RZP,
                   const ScalarFieldType& Psi,
                   const VecFieldType& Bcell_RZP,
                   const VecFieldType& Bcell_Curl_RZP,
                   const VecFieldType& Curlcell_NB_RZP,
                   const VecFieldType& GradPsicell_RZP,
                   const Coeff_1DType& Coeff_1D,
                   const Coeff_2DType& Coeff_2D,
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
    if (!this->Evaluate(tmp, pInfo, locatorRZ, locatorB, cellSetRZ, cellSetB, B_RZP, AsPhiFF, DAsPhiFF_RZP, Psi, Bcell_RZP, Bcell_Curl_RZP, Curlcell_NB_RZP, GradPsicell_RZP, Coeff_1D, Coeff_2D, k1))
      return false;
    tmp = p0 + k1*this->StepSize_2;

    DBG("    ****** K2"<<std::endl);
    if (!this->Evaluate(tmp, pInfo, locatorRZ, locatorB, cellSetRZ, cellSetB, B_RZP, AsPhiFF, DAsPhiFF_RZP, Psi, Bcell_RZP, Bcell_Curl_RZP, Curlcell_NB_RZP, GradPsicell_RZP, Coeff_1D, Coeff_2D, k2))
      return false;
    tmp = p0 + k2*this->StepSize_2;

    DBG("    ****** K3"<<std::endl);
    if (!this->Evaluate(tmp, pInfo, locatorRZ, locatorB, cellSetRZ, cellSetB, B_RZP, AsPhiFF, DAsPhiFF_RZP, Psi, Bcell_RZP, Bcell_Curl_RZP, Curlcell_NB_RZP, GradPsicell_RZP, Coeff_1D, Coeff_2D, k3))
      return false;
    tmp = p0 + k3*this->StepSize;

    DBG("    ****** K4"<<std::endl);
    if (!this->Evaluate(tmp, pInfo, locatorRZ, locatorB, cellSetRZ, cellSetB, B_RZP, AsPhiFF, DAsPhiFF_RZP, Psi, Bcell_RZP, Bcell_Curl_RZP, Curlcell_NB_RZP, GradPsicell_RZP, Coeff_1D, Coeff_2D, k4))
      return false;

    vtkm::Vec3f vec = (k1 + 2*k2 + 2*k3 + k4)/6.0;
    res = p0 + this->StepSize * vec;

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
    vtkm::VecVariable<vtkm::FloatDefault, 3> vals;
    vals.Append(sPortal.Get(vId[0]+offset));
    vals.Append(sPortal.Get(vId[1]+offset));
    vals.Append(sPortal.Get(vId[2]+offset));

    /*
    std::cout<<"  EvalS idx: "<<vId[0]<<" "<<vId[1]<<" "<<vId[2]<<std::endl;
    std::cout<<"  EvalS params: "<<param<<std::endl;
    std::cout<<"  EvalS: "<<vals[0]<<" "<<vals[1]<<" "<<vals[2]<<std::endl;
    */

    vtkm::FloatDefault s;
    vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), s);
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
    //std::cout<<"    ******** vid= "<<vId[0]<<" "<<vId[1]<<" "<<vId[2]<<std::endl;
    //std::cout<<"    ******** vec= "<<vPortal.Get(vId[0])<<" "<<vPortal.Get(vId[1])<<" "<<vPortal.Get(vId[2])<<std::endl;
    vtkm::VecVariable<vtkm::Vec3f, 3> vals;
    vals.Append(vPortal.Get(vId[0]+offset));
    vals.Append(vPortal.Get(vId[1]+offset));
    vals.Append(vPortal.Get(vId[2]+offset));

    vtkm::Vec3f v;
    vtkm::exec::CellInterpolate(vals, param, vtkm::CellShapeTagTriangle(), v);
    return v;
  }

  template <typename LocatorType, typename CellSetType>
  VTKM_EXEC
  bool PtLoc(const vtkm::Vec3f& ptRZ,
             ParticleInfo& pInfo,
             const LocatorType& locator,
             const CellSetType& cs,
             vtkm::Vec3f& param,
             vtkm::Vec<vtkm::Id, 3>& vIds) const
  {
    vtkm::Id cellId;
    vtkm::ErrorCode status = locator.FindCell(ptRZ, cellId, param, pInfo.PrevCell);
    if (status != vtkm::ErrorCode::Success)
    {
      printf("Find Cell failed! pt= %lf %lf %lf\n", ptRZ[0], ptRZ[1], ptRZ[2]);
#ifndef VTKM_CUDA
      std::cout<<"PtLoc(): Point not found: "<<ptRZ<<std::endl;
#endif
      return false;
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
    //std::cout<<"GetIndex: "<<x<<" "<<nx<<" min= "<<xmin<<" dx_inv "<<dx_inv<<std::endl;
    //return std::max(0, std::min(nx-1, (int)((x-xmin)*dx_inv)) );
    //return std::max(0, std::min(nx-1,    (int)((x-xmin)*dx_inv)) );
    int idx = std::max(1, std::min(nx  , 1 + int ((x-xmin)*dx_inv)) );
    return idx-1;
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
    //printf("  worklet %d\n", __LINE__);
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

    vtkm::FloatDefault acoef[4] = {coeff_1D.Get(idx+0),
                                   coeff_1D.Get(idx+1),
                                   coeff_1D.Get(idx+2),
                                   coeff_1D.Get(idx+3)};

    vtkm::FloatDefault iVal = 0.0;
    if (ideriv==0)
      iVal = acoef[0]+(acoef[1]+(acoef[2]+acoef[3]*wp)*wp)*wp;
    else if (ideriv==1)
      iVal = (acoef[1]+(2.0*acoef[2]+3.0*acoef[3]*wp)*wp)*one_d_cub_dpsi_inv;

    //printf("  worklet %d\n", __LINE__);
    return iVal * this->sml_bp_sign;
  }

  template <typename PortalType>
  VTKM_EXEC
  vtkm::FloatDefault
  EvalS4(const PortalType& sPortal,
         const vtkm::Vec<vtkm::Id, 4>& vId,
         const vtkm::Vec3f& param) const
  {
    vtkm::FloatDefault vals[4] = {sPortal.Get(vId[0]),
                                  sPortal.Get(vId[1]),
                                  sPortal.Get(vId[2]),
                                  sPortal.Get(vId[3])};

    const vtkm::FloatDefault& u = param[0];
    const vtkm::FloatDefault& v = param[1];
    vtkm::FloatDefault a = vtkm::Lerp(vals[0], vals[1], u);
    vtkm::FloatDefault b = vtkm::Lerp(vals[3], vals[2], u);
    vtkm::FloatDefault c = vtkm::Lerp(a,b,v);
    return c;
  }

  template <typename PortalType>
  VTKM_EXEC
  vtkm::Vec3f
  EvalV4(const PortalType& vPortal,
         const vtkm::Vec<vtkm::Id, 4>& vId,
         const vtkm::Vec3f& param) const
  {
    vtkm::Vec3f vals[4] = {vPortal.Get(vId[0]),
                           vPortal.Get(vId[1]),
                           vPortal.Get(vId[2]),
                           vPortal.Get(vId[3])};

    const vtkm::FloatDefault& u = param[0];
    const vtkm::FloatDefault& v = param[1];
    vtkm::Vec3f a = vtkm::Lerp(vals[0], vals[1], u);
    vtkm::Vec3f b = vtkm::Lerp(vals[3], vals[2], u);
    vtkm::Vec3f c = vtkm::Lerp(a,b,v);
    return c;
  }

  template <typename LocatorBType, typename CellSetBType, typename VecFieldType>
  VTKM_EXEC
  vtkm::Vec3f GetB(const vtkm::Vec3f& pt_rpz,
                   ParticleInfo& pInfo,
                   const LocatorBType& locatorB,
                   const CellSetBType& cellSetB,
                   const VecFieldType& B_RZP) const
  {
    this->EvalB(pt_rpz, locatorB, cellSetB, B_RZP, pInfo);
    return pInfo.B0_rzp;
    /*
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
    */
  }

  template <typename LocatorBType, typename CellSetBType, typename VectorFieldType>
  VTKM_EXEC
  bool CalcFieldFollowingPt(const vtkm::Vec3f& pt_rpz,
                            const ParticleInfo& pInfo,
                            const vtkm::Vec3f& B0_rpz,
                            const vtkm::FloatDefault& Phi0,
                            const vtkm::FloatDefault& Phi1,
                            const LocatorBType& locatorB,
                            const CellSetBType& cellSetB,
                            const VectorFieldType& B_RZP,
                            vtkm::Vec3f& x_ff_rpz) const
  {
    vtkm::FloatDefault R = pt_rpz[0];
    vtkm::FloatDefault Phi = pt_rpz[1];
    vtkm::FloatDefault Z = pt_rpz[2];
    vtkm::FloatDefault PhiMid = Phi0 + (Phi1-Phi0)/2.0;
    vtkm::Plane<> midPlane({0, PhiMid, 0}, {0,1,0});
    Ray3f ray_rpz({R, Phi, Z}, B0_rpz);

    //Now, do it using RK4 and two steps.
    vtkm::FloatDefault h = (PhiMid-Phi) / 2.0;
    vtkm::FloatDefault h_2 = h / 2.0;

    //k1 = F(p)
    //k2 = F(p+hk1/2)
    //k3 = F(p+hk2/2)
    //k4 = F(p+hk3)
    //Yn+1 = Yn + 1/6 h (k1+2k2+2k3+k4)
    vtkm::Vec3f p0 = {R,Phi,Z}; //pt_rpz;
    vtkm::Vec3f tmp, k1, k2, k3, k4;
    ParticleInfo pInfo2 = pInfo;

    for (int i = 0; i < 2; i++)
    {
      k1 = this->GetB(p0, pInfo2, locatorB, cellSetB, B_RZP);
      tmp = p0 + k1*h_2;

      k2 = this->GetB(tmp, pInfo2, locatorB, cellSetB, B_RZP);
      tmp = p0 + k2*h_2;

      k3 = this->GetB(tmp, pInfo2, locatorB, cellSetB, B_RZP);
      tmp = p0 + k3*h;

      k4 = this->GetB(tmp, pInfo2, locatorB, cellSetB, B_RZP);

      vtkm::Vec3f vec = (k1 + 2*k2 + 2*k3 + k4) / 6.0;
      p0 = p0 + h * vec;
    }
    x_ff_rpz = p0;

    return true;
  }

  VTKM_EXEC
  vtkm::Id CellLocUniform(const vtkm::Vec3f& ptRPZ,
                          vtkm::Vec3f* param=nullptr) const
  {
    vtkm::Id cellId = -1;

    const vtkm::FloatDefault& R= ptRPZ[0], Z = ptRPZ[2];
    if (R < this->Origin[0] || R > this->MaxPt[0]) return cellId;
    if (Z < this->Origin[1] || Z > this->MaxPt[1]) return cellId;

    auto tmpR = (R-this->Origin[0])*this->InvSpacing[0];
    auto tmpZ = (Z-this->Origin[1])*this->InvSpacing[1];
    vtkm::Id i = vtkm::Id(tmpR);
    vtkm::Id j = vtkm::Id(tmpZ);

    cellId = j * this->CellDims[0] + i;

    if (param != nullptr)
    {
      (*param)[0] = tmpR - i;
      (*param)[1] = tmpZ - j;
      (*param)[2] = 0;
    }
    return cellId;
  }

  template <typename LocatorBType, typename CellSetBType>
  VTKM_EXEC
  bool PtLoc4(const vtkm::Vec3f& ptRPZ,
              ParticleInfo& /*pInfo*/,
              const LocatorBType& /*locator*/,
              const CellSetBType& cs,
              vtkm::Vec3f& param,
              vtkm::Vec<vtkm::Id, 4>& vIds) const
  {
#if 0
    vtkm::Id cellId;
    vtkm::ErrorCode status = locator.FindCell(ptRZ, cellId, param, pInfo.PrevCell);
    if (status != vtkm::ErrorCode::Success)
    {
      printf("Find Cell failed! pt= %lf %lf %lf\n", ptRZ[0], ptRZ[1], ptRZ[2]);
#ifndef VTKM_CUDA
      std::cout<<"PtLoc(): Point not found: "<<ptRZ<<std::endl;
#endif
      return false;
    }

    //vtkm::VecVariable<vtkm::Id, 3> tmp;
    auto tmp =  cs.GetIndices(cellId);
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];
    vIds[3] = tmp[3];

    return true;
#else
    vtkm::Id cellId = this->CellLocUniform(ptRPZ, &param);
    if (cellId < 0)
    {
#ifndef VTKM_CUDA
      std::cout<<"PtLoc4: Point not found: "<<ptRPZ<<std::endl;
#endif
      return false;
    }

    auto tmp =  cs.GetIndices(cellId);
    vIds[0] = tmp[0];
    vIds[1] = tmp[1];
    vIds[2] = tmp[2];
    vIds[3] = tmp[3];
    return true;
#endif

  }

  template <typename LocatorBType, typename CellSetBType, typename ScalarFieldType, typename VecFieldType>
  VTKM_EXEC
  bool EvalB(const vtkm::Vec3f& ptRPZ,
             const LocatorBType& locatorB,
             const CellSetBType& cellSetB,
             const ScalarFieldType& Psi,
             const VecFieldType& B_RZP,
             const VecFieldType& B_Curl_RZP,
             const VecFieldType& Curl_NB_RZP,
             const VecFieldType& GradPsi_RZP,
             ParticleInfo& pInfo) const
  {
    if (this->BCell)
    {
      vtkm::Id cellId = this->CellLocUniform(ptRPZ);

      if (cellId < 0)
      {
#ifndef VTKM_CUDA
        std::cout<<"CellId Not found!!"<<std::endl;
#endif
        return false;
      }

      pInfo.Psi = Psi.Get(cellId);
      pInfo.B0_rzp = B_RZP.Get(cellId);
      pInfo.curlB_rzp = B_Curl_RZP.Get(cellId);
      pInfo.curl_nb_rzp = Curl_NB_RZP.Get(cellId);
      pInfo.gradPsi_rzp = GradPsi_RZP.Get(cellId);
    }
    else
    {
      ParticleInfo pInfo2;
      vtkm::Vec3f param;
      vtkm::Vec<vtkm::Id, 4> vIds;
      this->PtLoc4(ptRPZ, pInfo, locatorB, cellSetB, param, vIds);

      pInfo.Psi = this->EvalS4(Psi, vIds, param);
      pInfo.B0_rzp = this->EvalV4(B_RZP, vIds, param);
      pInfo.curlB_rzp = this->EvalV4(B_Curl_RZP, vIds, param);
      pInfo.curl_nb_rzp = this->EvalV4(Curl_NB_RZP, vIds, param);
      pInfo.gradPsi_rzp = this->EvalV4(GradPsi_RZP, vIds, param);
/*
      //std::cout<<"  PTIDs= "<<vIds[0]<<" "<<vIds[1]<<" "<<vIds[2]<<" "<<vIds[3]<<std::endl;

      vtkm::VecVariable<vtkm::FloatDefault, 4> valsS;
      for (int i = 0; i < 4; i++) valsS.Append(Psi.Get(vIds[i]));
      vtkm::exec::CellInterpolate(valsS, param, vtkm::CellShapeTagQuad(), pInfo.Psi);

      vtkm::VecVariable<vtkm::Vec3f, 4> valsV;
      for (int i = 0; i < 4; i++) valsV.Append(B_RZP.Get(vIds[i]));
      vtkm::exec::CellInterpolate(valsV, param, vtkm::CellShapeTagQuad(), pInfo.B0_rzp);

      for (int i = 0; i < 4; i++) valsV[i] = B_Curl_RZP.Get(vIds[i]);
      vtkm::exec::CellInterpolate(valsV, param, vtkm::CellShapeTagQuad(), pInfo.curlB_rzp);
      //std::cout<<"   BCurl: "<<param<<std::endl;
      //for (int i = 0; i < 4; i++) std::cout<<"   bcurl_i= "<<valsV[i]<<std::endl;

      for (int i = 0; i < 4; i++) valsV[i] = Curl_NB_RZP.Get(vIds[i]);
      vtkm::exec::CellInterpolate(valsV, param, vtkm::CellShapeTagQuad(), pInfo.curl_nb_rzp);
      for (int i = 0; i < 4; i++) valsV[i] = GradPsi_RZP.Get(vIds[i]);
      vtkm::exec::CellInterpolate(valsV, param, vtkm::CellShapeTagQuad(), pInfo.gradPsi_rzp);
*/
    }
    return true;
  }

  template <typename LocatorBType, typename CellSetBType, typename VecFieldType>
  VTKM_EXEC
  bool EvalB(const vtkm::Vec3f& ptRPZ,
             const LocatorBType& locatorB,
             const CellSetBType& cellSetB,
             const VecFieldType& B_RZP,
             ParticleInfo& pInfo) const
  {
    if (this->BCell)
    {
      vtkm::Id cellId = this->CellLocUniform(ptRPZ);

      if (cellId < 0)
      {
#ifndef VTKM_CUDA
        std::cout<<"CellId Not found!!"<<std::endl;
#endif
        return false;
      }

      pInfo.B0_rzp = B_RZP.Get(cellId);
    }
    else
    {
      vtkm::Vec3f param;
      vtkm::Vec<vtkm::Id, 4> vIds;
      this->PtLoc4(ptRPZ, pInfo, locatorB, cellSetB, param, vIds);

      pInfo.B0_rzp = this->EvalV4(B_RZP, vIds, param);

      /*
      vtkm::VecVariable<vtkm::Vec3f, 4> valsV;
      for (int i = 0; i < 4; i++) valsV.Append(B_RZP.Get(vIds[i]));
      vtkm::exec::CellInterpolate(valsV, param, vtkm::CellShapeTagQuad(), pInfo.B0_rzp);
      */
    }

    return true;
  }

#ifndef VTKM_CUDA
  VTKM_EXEC
  void PrintDiffs(const vtkm::Vec3f& pt, const ParticleInfo& piA, const ParticleInfo& piB) const
  {
    std::cout<<std::setprecision(15)
    <<"Pt: "<<pt<<std::endl
    <<"  BGrid, RZ, Error "<<std::endl
    <<"  B-comp psi: "<<piA.Psi<<" "<<piB.Psi<<" : "<<vtkm::Abs(piA.Psi-piB.Psi)<<std::endl
    <<"  B-comp b0: "<<piA.B0_rzp<<" "<<piB.B0_rzp<<" "<<vtkm::Magnitude(piA.B0_rzp - piB.B0_rzp)<<std::endl
    <<"  B-comp c : "<<piA.curlB_rzp<<" "<<piB.curlB_rzp<<" "<<vtkm::Magnitude(piA.curlB_rzp - piB.curlB_rzp)<<std::endl
    <<"  B-comp cn: "<<piA.curl_nb_rzp<<" "<<piB.curl_nb_rzp<<" "<<vtkm::Magnitude(piA.curl_nb_rzp - piB.curl_nb_rzp)<<std::endl
    <<"  B-comp gs: "<<piA.gradPsi_rzp<<" "<<piB.gradPsi_rzp<<" "<<vtkm::Magnitude(piA.gradPsi_rzp - piB.gradPsi_rzp)<<std::endl;
  }
#endif

  template <typename LocatorRZType, typename LocatorBType, typename CellSetRZType, typename CellSetBType, typename VecFieldType, typename AsFieldType, typename DAsFieldType, typename ScalarFieldType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool Evaluate(const vtkm::Vec3f& ptRPZ,
                ParticleInfo& pInfo,
                const LocatorRZType& locatorRZ,
                const LocatorBType& locatorB,
                const CellSetRZType& cellSetRZ,
                const CellSetBType& cellSetB,
                const VecFieldType& /*B_RZP*/,
                const AsFieldType& AsPhiFF,
                const DAsFieldType& DAsPhiFF_RZP,
                const ScalarFieldType& Psi,
                const VecFieldType& Bcell_RZP,
                const VecFieldType& Bcell_Curl_RZP,
                const VecFieldType& Curlcell_NB_RZP,
                const VecFieldType& GradPsicell_RZP,
                const Coeff_1DType& /*coeff_1D*/,
                const Coeff_2DType& /*coeff_2D*/,
                vtkm::Vec3f& res) const
  {
    auto R = ptRPZ[0];
    auto Phi = ptRPZ[1];
    auto Z = ptRPZ[2];

    /*
    //expensive B...
    ParticleInfo pInfoOld;
    if (!this->HighOrderB(ptRPZ, pInfoOld, coeff_1D, coeff_2D))
      return false;
    */

    //uniform grid B
    if (!this->EvalB(ptRPZ, locatorB, cellSetB, Psi, Bcell_RZP, Bcell_Curl_RZP, Curlcell_NB_RZP, GradPsicell_RZP, pInfo))
      return false;

    /*
    std::cout<<"bummy"<<std::endl;
    this->PrintDiffs(ptRPZ, pInfo, pInfoOld);
    if (0)
    {
      std::cout<<"***************Time to compare"<<std::endl;
      this->HighOrderB(ptRPZ, pInfoOld, coeff_1D, coeff_2D);
      this->PrintDiffs(ptRPZ, pInfo, pInfoOld);
    }
    */

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
    vtkm::Vec3f B0_rpz(pInfo.B0_rzp[0], pInfo.B0_rzp[2], pInfo.B0_rzp[1]);

    vtkm::Vec3f ff_pt_rpz;
    this->CalcFieldFollowingPt({R,phiN,Z}, pInfo, B0_rpz, Phi0, Phi1, locatorB, cellSetB, Bcell_RZP, ff_pt_rpz);

    //Now, interpolate between Phi_i and Phi_i+1
    vtkm::FloatDefault T01 = (phiN - Phi0) / (Phi1-Phi0);
    vtkm::FloatDefault T10 = 1.0f - T01;

    //Get vec at Phi0 and Phi1.
    //x_ff is in rzp
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
    this->PtLoc(x_ff_rzp, pInfo, locatorRZ, cellSetRZ, x_ff_param, x_ff_vids);
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

    //std::cout<<" gradAs_rpz= "<<gradAs_rpz<<std::endl;

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
    plane1 = plane0 + 1;
    phi0 = static_cast<vtkm::FloatDefault>(plane0)*this->dPhi;
    phi1 = static_cast<vtkm::FloatDefault>(plane1)*this->dPhi;

    DBG("             ***** phi= "<<phi<<" phiN= "<<phiN<<" "<<plane0<<" "<<plane1<<" "<<phi0<<" "<<phi1<<std::endl);

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
  bool UseHighOrder = false;
  bool SaveTraces = false;
  bool QuickTest = false;

  int nr, nz;
  vtkm::FloatDefault rmin, zmin, rmax, zmax;
  vtkm::FloatDefault dr, dz, dr_inv, dz_inv;
  vtkm::FloatDefault EqAxisZ, EqAxisR, EqXPsi;

  int ncoeff;
  vtkm::FloatDefault min_psi, max_psi;
  vtkm::FloatDefault one_d_cub_dpsi_inv;
  vtkm::FloatDefault sml_bp_sign = -1.0f;

  //worklet3 stufff.
  vtkm::Vec2f Origin, Spacing, MaxPt;
  vtkm::Id2 Dims, CellDims;
  vtkm::Vec2f InvSpacing;
  bool BCell;
};
