
//#define PRINT_STUFF
#ifdef PRINT_STUFF
#define DBG(x) std::cout<<x
#else
#define DBG(x)
#endif

#define DO_TRACES 0

//-----------------------------------------------------------------------------
class ComputeBWorklet : public vtkm::worklet::WorkletMapField
{
  using FloatPortal = typename vtkm::cont::ArrayHandle<vtkm::FloatDefault>::ReadPortalType;
  using VecPortal = typename vtkm::cont::ArrayHandle<vtkm::Vec3f>::ReadPortalType;
public:
  using ControlSignature = void(FieldIn Coords,
                                WholeArrayIn coeff_1D,
                                WholeArrayIn coeff_2D,
                                FieldOut Psi,
                                FieldOut B0,
                                FieldOut CurlB0,
                                FieldOut CurlNB0,
                                FieldOut GradPsi);

  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6, _7, _8);
  using InputDomain = _1;

  ~ComputeBWorklet()
  {
  }

  ComputeBWorklet()
  {
    this->nr = eq_mr-1;
    this->nz = eq_mz-1;
    this->rmin = eq_min_r;
    this->rmax = eq_max_r;
    this->zmin = eq_min_z;
    this->zmax = eq_max_z;
    this->ncoeff = eq_mr-1;

    this->dr = (eq_max_r - eq_min_r) / vtkm::FloatDefault(this->nr);
    this->dz = (eq_max_z - eq_min_z) / vtkm::FloatDefault(this->nz);
    this->dr_inv = 1.0/this->dr;
    this->dz_inv = 1.0/this->dz;

    this->ncoeff = eq_mr-1;
    this->min_psi = psi_min;
    this->max_psi = psi_max;
    this->one_d_cub_dpsi_inv = 1.0 / ((this->max_psi-this->min_psi)/vtkm::FloatDefault(this->ncoeff));
  }

  template <typename CoordsType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            const CoordsType& CoordsRZ,
                            const Coeff_1DType& Coeff_1D,
                            const Coeff_2DType& Coeff_2D,
                            vtkm::FloatDefault& Psi,
                            vtkm::Vec3f& B0,
                            vtkm::Vec3f& CurlB0,
                            vtkm::Vec3f& CurlNB0,
                            vtkm::Vec3f& GradPsi) const
  {
    vtkm::FloatDefault R = CoordsRZ[0];
    vtkm::FloatDefault Z = CoordsRZ[1];

    //vtkm::Vec3f ptRPZ(3.029365, 6.183185, 0.020600);
    //ptRPZ[1] = 0;

    vtkm::Vec3f ptRPZ(R,0,Z);
    vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;

    bool printIt = false;
    /*
    if (idx == 19957837)
    {
      printIt = true;
      R = 3.3514430285644154;
      Z = -0.45164880640275618;
      ptRPZ[0] = R;
      ptRPZ[2] = Z;
    }
    */
    this->HighOrderB(ptRPZ, Coeff_1D, Coeff_2D,
                     B0, jacobian_rzp, CurlB0, CurlNB0, Psi, GradPsi, printIt);

#ifndef VTKM_CUDA
    if (printIt)
    {
      std::cout<<"IDX: "<<idx<<std::endl;
      std::cout<<std::setprecision(20)
               <<"RZ: "<<R<<" "<<Z<<std::endl
               <<"  psi: "<<Psi<<std::endl
               <<"  B0 : "<<B0<<std::endl
               <<"  CB : "<<CurlB0<<std::endl
               <<"  CNB : "<<CurlNB0<<std::endl
               <<"  GS  : "<<GradPsi<<std::endl;
    }
#endif

#ifndef VTKM_CUDA
    if (false)
    {
      ptRPZ[0] = 3.35144302856442;
      ptRPZ[2] = -0.451648806402756;
      this->HighOrderB(ptRPZ, Coeff_1D, Coeff_2D,
                       B0, jacobian_rzp, CurlB0, CurlNB0, Psi, GradPsi);
      std::cout<<"IDX========================: "<<idx<<std::endl;
      std::cout<<std::setprecision(20)
               <<"RZ: "<<R<<" "<<Z<<std::endl
               <<"  psi: "<<Psi<<std::endl
               <<"  B0 : "<<B0<<std::endl
               <<"  CB : "<<CurlB0<<std::endl
               <<"  CNB : "<<CurlNB0<<std::endl
               <<"  GS  : "<<GradPsi<<std::endl;
    }
#endif

    //std::cout<<"ComputeB: "<<ptRPZ<<" "<<Psi<<" "<<B0<<" "<<CurlB0<<" "<<CurlNB0<<" "<<GradPsi<<std::endl;
    //ComputeB: [3.02936,0,0.0206] 0.00697323 [-0.00168909,0.0188064,-0.218501] [0,0,-0.138882] [-0.00220759,-0.304296,-0.659427] [0.0569713,0.00511687,0]
    //exit(0);
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

    vtkm::FloatDefault acoef[4] = {coeff_1D.Get(idx+0),
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

  template <typename Coeff_2DType>
  VTKM_EXEC
  void GetCoeff2D(const vtkm::Id &offset,
                  const Coeff_2DType& Coeff_2D,
                  vtkm::Matrix<vtkm::FloatDefault, 4, 4>& coeff) const
  {
    vtkm::Id idx = offset;
    for (vtkm::Id ii = 0; ii < 4; ii++)
      for (vtkm::Id jj = 0; jj < 4; jj++)
        coeff[ii][jj] = Coeff_2D.Get(idx++);
  }

  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  void HighOrderB(const vtkm::Vec3f& ptRPZ,
                  const Coeff_1DType& Coeff_1D,
                  const Coeff_2DType& Coeff_2D,
                  vtkm::Vec3f& B0_rzp,
                  vtkm::Vec<vtkm::Vec3f, 3>& jacobian_rzp,
                  vtkm::Vec3f& curlB_rzp,
                  vtkm::Vec3f& curl_nb_rzp,
                  vtkm::FloatDefault& PSI,
                  vtkm::Vec3f& gradPsi_rzp,
                  bool printIt=false) const
  {
    vtkm::FloatDefault R = ptRPZ[0], Z = ptRPZ[2];
    //vtkm::Vec3f ptRZ(R,Z,0);

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

    this->GetCoeff2D(offset, Coeff_2D, acoeff);
/*
    vtkm::Id idx = 0;
    //std::cout<<"Offset= "<<(offset/16)<<" 16: "<<offset<<std::endl;
    for (vtkm::Id ii = 0; ii < 4; ii++)
      for (vtkm::Id jj = 0; jj < 4; jj++)
      {
        acoeff[ii][jj] = Coeff_2D.Get(offset+idx);
        //acoeff[jj][ii] = Coeff_2D.Get(offset+idx); //z,r
        //std::cout<<"c_"<<ii<<jj<<"= "<<Coeff_2D.Get(offset+idx)<<std::endl;
        idx++;
      }
*/

    vtkm::FloatDefault psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
    this->eval_bicub_2(R, Z, Rc, Zc, acoeff, psi, dpsi_dr, dpsi_dz, d2psi_drdz, d2psi_d2r, d2psi_d2z, printIt);

    /*
    pRPZ.Psi = psi;
    pRPZ.dpsi_dr = dpsi_dr;
    pRPZ.dpsi_dz = dpsi_dz;
    pRPZ.d2psi_drdz = d2psi_drdz;
    pRPZ.d2psi_d2r = d2psi_d2r;
    pRPZ.d2psi_d2z = d2psi_d2z;
    */
    PSI = psi;
    gradPsi_rzp[0] = dpsi_dr;
    gradPsi_rzp[1] = dpsi_dz;
    gradPsi_rzp[2] = 0;
    /*
    std::cout<<" psi= "<<psi<<std::endl;
    std::cout<<" dpsi_dr = "<<dpsi_dr<<std::endl;
    std::cout<<" dpsi_dz = "<<dpsi_dz<<std::endl;
    */

    vtkm::FloatDefault fld_I = this->I_interpol(PSI, 0, Coeff_1D);
    vtkm::FloatDefault fld_dIdpsi = this->I_interpol(PSI, 1, Coeff_1D);

    vtkm::FloatDefault over_r = 1/R;
    vtkm::FloatDefault over_r2 = over_r*over_r;
    vtkm::FloatDefault Br = -dpsi_dz * over_r;
    vtkm::FloatDefault Bz = dpsi_dr * over_r;
    vtkm::FloatDefault Bp = fld_I * over_r;

    B0_rzp = vtkm::Vec3f(Br, Bz, Bp);

    //Set the jacobian.
    const int PIR = 0;
    const int PIZ = 1;
    const int PIP = 2;
    //vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        jacobian_rzp[i][j] = 0;
    vtkm::FloatDefault bp_sign = 1.0;

    // R-derivatives depend on the geometry
    jacobian_rzp[PIR][PIR] = ( /*pRPZ.*/dpsi_dz * over_r2 - /*pRPZ.*/d2psi_drdz * over_r) * bp_sign;
    jacobian_rzp[PIR][PIZ] = (-/*pRPZ.*/dpsi_dr * over_r2 + /*pRPZ.*/d2psi_d2r  * over_r) * bp_sign;
    jacobian_rzp[PIR][PIP] = /*pRPZ.*/dpsi_dr * fld_dIdpsi * over_r - fld_I * over_r2;

    // Z and phi derivatives do not change between toroidal and cylindrical geometry
    jacobian_rzp[PIZ][PIR] = -/*pRPZ.*/d2psi_d2z * over_r * bp_sign;
    jacobian_rzp[PIP][PIR] = 0.0 * bp_sign;

    jacobian_rzp[PIZ][PIZ] = /*pRPZ.*/d2psi_drdz * over_r * bp_sign;
    jacobian_rzp[PIP][PIZ] = 0.0 * bp_sign;

    jacobian_rzp[PIZ][PIP] = fld_dIdpsi * /*pRPZ.*/dpsi_dz * over_r;
    jacobian_rzp[PIP][PIP] = 0.0;

    //std::cout<<"High order B: "<<std::endl;
    //vtkm::Vec3f ppp(R,Z,P);
    //std::cout<<"   p_rzp  = "<<ppp<<std::endl;
    //std::cout<<"   B0_rzp = "<<B0_rzp<<std::endl;
    //std::cout<<"   B0mag = "<<vtkm::Magnitude(B0_rzp)<<"  "<<(1.0/vtkm::Magnitude(B0_rzp))<<std::endl;
    //std::cout<<"   JAC_R = "<<jacobian_rzp[PIR]<<std::endl;
    //std::cout<<"   JAC_Z = "<<jacobian_rzp[PIZ]<<std::endl;
    //std::cout<<"   JAC_P = "<<jacobian_rzp[PIP]<<std::endl;

    auto dBr_dr = jacobian_rzp[0][0];
    auto dBr_dz = jacobian_rzp[1][0];
    auto dBr_dp = jacobian_rzp[2][0];

    auto dBz_dr = jacobian_rzp[0][1];
    auto dBz_dz = jacobian_rzp[1][1];
    auto dBz_dp = jacobian_rzp[2][1];

    auto dBp_dr = jacobian_rzp[0][2];
    auto dBp_dz = jacobian_rzp[1][2];
    //auto dBp_dp = jacobian_rzp[2][2];

    //these are correct.
    /*
    std::cout<<"dBr_dx: "<<dBr_dr<<" "<<dBr_dz<<" "<<dBr_dp<<std::endl;
    std::cout<<"dBz_dx: "<<dBz_dr<<" "<<dBz_dz<<" "<<dBz_dp<<std::endl;
    std::cout<<"dBp_dx: "<<dBp_dr<<" "<<dBp_dz<<" "<<dBp_dp<<std::endl;
    */

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
    /*pRPZ.*/curlB_rzp[0] = dBz_dp * over_r - dBp_dz;
    /*pRPZ.*/curlB_rzp[1] = Bp*over_r + dBp_dr - dBr_dp*over_r;
    /*pRPZ.*/curlB_rzp[2] = dBr_dz - dBz_dr;
    //std::cout<<"curl_B_rzp= "<<curlB_rzp<<std::endl;
#ifndef VTKM_CUDA
    if (printIt)
    {
      std::cout<<"  ComputeB"<<std::endl;
      std::cout<<"   acoeff= "<<acoeff<<std::endl;
      std::cout<<"   dBr_dz = "<<dBr_dz<<std::endl;
      std::cout<<"   dBz_dr = "<<dBz_dr<<std::endl;
      std::cout<<"    dpsi_dr= "<<dpsi_dr<<std::endl;
      std::cout<<"    d2psi_d2r= "<<d2psi_d2r<<std::endl;
      std::cout<<"    over_r/2 = "<<over_r<<" "<<over_r2<<std::endl<<std::endl;
    }
#endif

    //calculate curl_nb
    /*
    !curl of norm b
    curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
    curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
    curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2
    */

    vtkm::FloatDefault Bmag = vtkm::Magnitude(/*pRPZ.*/B0_rzp);
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
    /*pRPZ.*/curl_nb_rzp[0] = /*pRPZ.*/curlB_rzp[0] * over_B + ( Bp * dBdz)*over_B2;
    /*pRPZ.*/curl_nb_rzp[1] = /*pRPZ.*/curlB_rzp[1] * over_B + (-Bp * dBdr)*over_B2;
    /*pRPZ.*/curl_nb_rzp[2] = /*pRPZ.*/curlB_rzp[2] * over_B + (Bz*dBdr - Br*dBdz)*over_B2;

    /*
    std::cout<<"XXXXXX_B_rpz= "<<res<<std::endl;
    std::cout<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    */
  }


  VTKM_EXEC
  void eval_bicub_2(const vtkm::FloatDefault& x,
                    const vtkm::FloatDefault& y,
                    const vtkm::FloatDefault& xc,
                    const vtkm::FloatDefault& yc,
                    const vtkm::Matrix<vtkm::FloatDefault, 4, 4>& acoeff,
                    vtkm::FloatDefault &f00, vtkm::FloatDefault &f10, vtkm::FloatDefault &f01,
                    vtkm::FloatDefault &f11, vtkm::FloatDefault &f20, vtkm::FloatDefault &f02,
                    bool printIt = false) const
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

#ifndef VTKM_CUDA
    if (printIt) std::cout<<"RZ= "<<x<<" "<<y<<" RZc= "<<xc<<" "<<yc<<std::endl;
    if (printIt) std::cout<<"f11= "<<f11<<std::endl;
    if (printIt) std::cout<<"  fx= "<<fx[0]<<" "<<fx[1]<<" "<<fx[2]<<" "<<fx[3]<<std::endl;
    if (printIt) std::cout<<"  yv= "<<yv[0]<<" "<<yv[1]<<" "<<yv[2]<<" "<<yv[3]<<std::endl;
    if (printIt) std::cout<<" dfx= "<<dfx[0]<<" "<<dfx[1]<<" "<<dfx[2]<<" "<<dfx[3]<<std::endl;
    if (printIt) std::cout<<"  dx= "<<dx<<std::endl;
    if (printIt) std::cout<<"  dy= "<<dy<<std::endl;
#endif
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
  }


  int nr, nz;
  int ncoeff;
  vtkm::FloatDefault rmin, zmin, rmax, zmax;
  vtkm::FloatDefault dr, dz, dr_inv, dz_inv;
  vtkm::FloatDefault min_psi, max_psi;
  vtkm::FloatDefault one_d_cub_dpsi_inv;
  vtkm::FloatDefault sml_bp_sign = -1.0f;

#if 0
  template <typename Coeff_2DType>
  VTKM_EXEC
  void GetCoeff2D(const vtkm::Id &offset,
                  const Coeff_2DType& Coeff_2D,
                  vtkm::Matrix<vtkm::FloatDefault, 4, 4>& coeff) const
  {

    vtkm::Id idx = offset;
    for (vtkm::Id ii = 0; ii < 4; ii++)
      for (vtkm::Id jj = 0; jj < 4; jj++)
        coeff[ii][jj] = Coeff_2D.Get(idx++);
  }


  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool HighOrderB(vtkm::Particle& pRPZ,
                  const Coeff_1DType& Coeff_1D,
                  const Coeff_2DType& Coeff_2D) const
  /*
                  vtkm::Vec3f& B0_rzp,
                  vtkm::Vec<vtkm::Vec3f, 3>& jacobian_rzp,
                  vtkm::Vec3f& curlB_rzp,
                  vtkm::Vec3f& curl_nb_rzp,
                  vtkm::FloatDefault& PSI,
                  vtkm::Vec3f& gradPsi_rzp) const
  */
  {

    vtkm::FloatDefault R = pRPZ.Pos[0], Z = pRPZ.Pos[2], P = pRPZ.Pos[1];
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

    this->GetCoeff2D(offset, Coeff_2D, acoeff);
/*
    vtkm::Id idx = 0;
    //std::cout<<"Offset= "<<(offset/16)<<" 16: "<<offset<<std::endl;
    for (vtkm::Id ii = 0; ii < 4; ii++)
      for (vtkm::Id jj = 0; jj < 4; jj++)
      {
        acoeff[ii][jj] = Coeff_2D.Get(offset+idx);
        //acoeff[jj][ii] = Coeff_2D.Get(offset+idx); //z,r
        //std::cout<<"c_"<<ii<<jj<<"= "<<Coeff_2D.Get(offset+idx)<<std::endl;
        idx++;
      }
*/

    vtkm::FloatDefault psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z;
    this->eval_bicub_2(R, Z, Rc, Zc, acoeff, psi, dpsi_dr, dpsi_dz, d2psi_drdz, d2psi_d2r, d2psi_d2z, printIt);
    pRPZ.Psi = psi;
    pRPZ.dpsi_dr = dpsi_dr;
    pRPZ.dpsi_dz = dpsi_dz;
    pRPZ.d2psi_drdz = d2psi_drdz;
    pRPZ.d2psi_d2r = d2psi_d2r;
    pRPZ.d2psi_d2z = d2psi_d2z;
    //PSI = psi;
    pRPZ.gradPsi_rzp[0] = pRPZ.dpsi_dr;
    pRPZ.gradPsi_rzp[1] = pRPZ.dpsi_dz;
    pRPZ.gradPsi_rzp[2] = 0;
    /*
    std::cout<<" psi= "<<psi<<std::endl;
    std::cout<<" dpsi_dr = "<<dpsi_dr<<std::endl;
    std::cout<<" dpsi_dz = "<<dpsi_dz<<std::endl;
    */

    vtkm::FloatDefault fld_I = this->I_interpol(pRPZ.Psi, 0, Coeff_1D);
    vtkm::FloatDefault fld_dIdpsi = this->I_interpol(pRPZ.Psi, 1, Coeff_1D);

    vtkm::FloatDefault over_r = 1/R;
    vtkm::FloatDefault over_r2 = over_r*over_r;
    vtkm::FloatDefault Br = -pRPZ.dpsi_dz * over_r;
    vtkm::FloatDefault Bz = pRPZ.dpsi_dr * over_r;
    vtkm::FloatDefault Bp = fld_I * over_r;

    pRPZ.B0_rzp = vtkm::Vec3f(Br, Bz, Bp);

    //Set the jacobian.
    const int PIR = 0;
    const int PIZ = 1;
    const int PIP = 2;
    vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        jacobian_rzp[i][j] = 0;
    vtkm::FloatDefault bp_sign = 1.0;

    // R-derivatives depend on the geometry
    jacobian_rzp[PIR][PIR] = ( pRPZ.dpsi_dz * over_r2 - pRPZ.d2psi_drdz * over_r) * bp_sign;
    jacobian_rzp[PIR][PIZ] = (-pRPZ.dpsi_dr * over_r2 + pRPZ.d2psi_d2r  * over_r) * bp_sign;
    jacobian_rzp[PIR][PIP] = pRPZ.dpsi_dr * fld_dIdpsi * over_r - fld_I * over_r2;

    // Z and phi derivatives do not change between toroidal and cylindrical geometry
    jacobian_rzp[PIZ][PIR] = -pRPZ.d2psi_d2z * over_r * bp_sign;
    jacobian_rzp[PIP][PIR] = 0.0 * bp_sign;

    jacobian_rzp[PIZ][PIZ] = pRPZ.d2psi_drdz * over_r * bp_sign;
    jacobian_rzp[PIP][PIZ] = 0.0 * bp_sign;

    jacobian_rzp[PIZ][PIP] = fld_dIdpsi * pRPZ.dpsi_dz * over_r;
    jacobian_rzp[PIP][PIP] = 0.0;

    //std::cout<<"High order B: "<<std::endl;
    vtkm::Vec3f ppp(R,Z,P);
    //std::cout<<"   p_rzp  = "<<ppp<<std::endl;
    //std::cout<<"   B0_rzp = "<<B0_rzp<<std::endl;
    //std::cout<<"   B0mag = "<<vtkm::Magnitude(B0_rzp)<<"  "<<(1.0/vtkm::Magnitude(B0_rzp))<<std::endl;
    //std::cout<<"   JAC_R = "<<jacobian_rzp[PIR]<<std::endl;
    //std::cout<<"   JAC_Z = "<<jacobian_rzp[PIZ]<<std::endl;
    //std::cout<<"   JAC_P = "<<jacobian_rzp[PIP]<<std::endl;

    auto dBr_dr = jacobian_rzp[0][0];
    auto dBr_dz = jacobian_rzp[1][0];
    auto dBr_dp = jacobian_rzp[2][0];

    auto dBz_dr = jacobian_rzp[0][1];
    auto dBz_dz = jacobian_rzp[1][1];
    auto dBz_dp = jacobian_rzp[2][1];

    auto dBp_dr = jacobian_rzp[0][2];
    auto dBp_dz = jacobian_rzp[1][2];
    //auto dBp_dp = jacobian_rzp[2][2];

    //these are correct.
    /*
    std::cout<<"dBr_dx: "<<dBr_dr<<" "<<dBr_dz<<" "<<dBr_dp<<std::endl;
    std::cout<<"dBz_dx: "<<dBz_dr<<" "<<dBz_dz<<" "<<dBz_dp<<std::endl;
    std::cout<<"dBp_dx: "<<dBp_dr<<" "<<dBp_dz<<" "<<dBp_dp<<std::endl;
    */

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
    pRPZ.curlB_rzp[0] = dBz_dp * over_r - dBp_dz;
    pRPZ.curlB_rzp[1] = Bp*over_r + dBp_dr - dBr_dp*over_r;
    pRPZ.curlB_rzp[2] = dBr_dz - dBz_dr;
    //std::cout<<"curl_B_rzp= "<<curlB_rzp<<std::endl;

    //calculate curl_nb
    /*
    !curl of norm b
    curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
    curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
    curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2
    */

    vtkm::FloatDefault Bmag = vtkm::Magnitude(pRPZ.B0_rzp);
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
    if (vtkm::Abs(divergence) > 1e-16)
    {
      std::cout<<std::endl;
      std::cout<<"****************************************************** DIVERGENCE= "<<divergence<<std::endl;
      std::cout<<std::endl;
    }
    #endif


    //vtkm::Vec3f curl_nb_rzp;
    pRPZ.curl_nb_rzp[0] = pRPZ.curlB_rzp[0] * over_B + ( Bp * dBdz)*over_B2;
    pRPZ.curl_nb_rzp[1] = pRPZ.curlB_rzp[1] * over_B + (-Bp * dBdr)*over_B2;
    pRPZ.curl_nb_rzp[2] = pRPZ.curlB_rzp[2] * over_B + (Bz*dBdr - Br*dBdz)*over_B2;

    /*
    std::cout<<"XXXXXX_B_rpz= "<<res<<std::endl;
    std::cout<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    */

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType, typename OutputType, typename OutputType2D, typename IdType>
  VTKM_EXEC void operator()(const vtkm::Id& idx,
                            vtkm::Particle& particleRPZ,
                            const LocatorType& locator,
                            const CellSetType& cellSet,
                            const BFieldType& B_RZP,
                            const BFieldType& B_Norm_RZP,
                            //const BFieldType& Curl_NB_RZP,
                            const AsFieldType& AsPhiFF,
                            const DAsFieldType& DAsPhiFF_RZP,
                            const Coeff_1DType& Coeff_1D,
                            const Coeff_2DType& Coeff_2D,
                            OutputType& traces,
                            OutputType2D& outputRZ,
                            OutputType2D& outputTP,
                            IdType punctureID) const
  {

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

    while (true)
    {
      //printf("  worklet %d\n", __LINE__);
      vtkm::Vec3f newPos;
      DBG("\n\n\n*********************************************"<<std::endl);
      DBG("   "<<particleRPZ.Pos<<" #s= "<<particleRPZ.NumSteps<<std::endl);
      if (!this->TakeRK4Step(particleRPZ, locator, cellSet,
                             B_RZP, B_Norm_RZP,
                             AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, newPos))
      {
#ifndef VTKM_CUDA
        std::cout<<"*****************       All done: RK step failed."<<std::endl;
#endif
        break;
      }
      //printf("  worklet %d\n", __LINE__);

      DBG("     *** Step--> "<<newPos<<std::endl);
      vtkm::Id numRevs0 = vtkm::Floor(vtkm::Abs(particleRPZ.Pos[1] / vtkm::TwoPi()));
      vtkm::Id numRevs1 = vtkm::Floor(vtkm::Abs(newPos[1] / vtkm::TwoPi()));

      particleRPZ.Pos = newPos;
      particleRPZ.NumSteps++;

      //printf("  worklet %d\n", __LINE__);
      if (this->SaveTraces)
        traces.Set(idx*this->MaxIter + particleRPZ.NumSteps, particleRPZ.Pos);

      if (numRevs1 > numRevs0)
      {
        auto R = particleRPZ.Pos[0], Z = particleRPZ.Pos[2];
        auto theta = vtkm::ATan2(Z-this->EqAxisZ, R-this->EqAxisR);
        if (theta < 0)
          theta += vtkm::TwoPi();

        vtkm::FloatDefault psi = this->GetPsi(particleRPZ.Pos, Coeff_2D);

        //calcualte psi. need to stash psi on the particle somehow....
        //vtkm::Vec3f ptRPZ = particleRPZ.Pos;
        /*
        vtkm::FloatDefault psi;
        vtkm::Vec3f B0_rzp, curlB_rzp, curl_nb_rzp, gradPsi_rzp;
        vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;
        this->HighOrderB(particleRPZ, Coeff_1D, Coeff_2D, B0_rzp, jacobian_rzp, curlB_rzp, curl_nb_rzp, psi, gradPsi_rzp);
        */

        //std::cout<<" PUNC: "<<particleRPZ.Pos<<" "<<psi<<std::endl;
        vtkm::Id i = (idx * this->MaxPunc) + particleRPZ.NumPunctures;
        outputRZ.Set(i, vtkm::Vec2f(R, Z));
        outputTP.Set(i, vtkm::Vec2f(theta, psi/this->EqXPsi));
        punctureID.Set(i, idx);
        particleRPZ.NumPunctures++;

#ifndef VTKM_CUDA
        if (idx == 0 && particleRPZ.NumPunctures%10 == 0 ) std::cout<<" ***** PUNCTURE n= "<<particleRPZ.NumPunctures<<std::endl;
#endif
        DBG("************* PUNCTURE n= "<<particleRPZ.NumPunctures<<std::endl);
      }
      //printf("  worklet %d\n", __LINE__);
      if (particleRPZ.NumSteps >= this->MaxIter || particleRPZ.NumPunctures >= this->MaxPunc)
      {
#ifndef VTKM_CUDA
        std::cout<<"************************************* All done: "<<particleRPZ<<std::endl;
#endif
        break;
      }
    }
#ifndef VTKM_CUDA
    std::cout<<"Particle done: "<<idx<<std::endl;
#endif
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool TakeRK4Step(vtkm::Particle& particleRPZ,
                   const LocatorType& locator,
                   const CellSetType& cellSet,
                   const BFieldType& B_RZP,
                   const BFieldType& B_Norm_RZP,
                   const AsFieldType& AsPhiFF,
                   const DAsFieldType& DAsPhiFF_RZP,
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

    //std::cout<<"TakeRK4Step: "<<particleRPZ<<std::endl;

    vtkm::Particle p0, tmp;
    p0.Pos = particleRPZ.Pos;
    DBG("    ****** K1"<<std::endl);
    if (!this->Evaluate(p0, locator, cellSet, B_RZP, B_Norm_RZP, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, k1))
      return false;
    tmp.Pos = p0.Pos + k1*this->StepSize_2;
    //std::cout<<" K1: "<<p0.Pos<<" --> "<<tmp.Pos<<" k1= "<<k1<<std::endl;

    DBG("    ****** K2"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, k2))
      return false;
    tmp.Pos = p0.Pos + k2*this->StepSize_2;

    DBG("    ****** K3"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, k3))
      return false;
    tmp.Pos = p0.Pos + k3*this->StepSize;

    DBG("    ****** K4"<<std::endl);
    if (!this->Evaluate(tmp, locator, cellSet, B_RZP, B_Norm_RZP, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, k4))
      return false;

    vtkm::Vec3f vec = (k1 + 2*k2 + 2*k3 + k4)/6.0;
    res = p0.Pos + this->StepSize * vec;

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
             const LocatorType& locator,
             const CellSetType& cs,
             vtkm::Vec3f& param,
             vtkm::Vec<vtkm::Id, 3>& vIds) const
  {

    vtkm::Id cellId;
    vtkm::ErrorCode status = locator.FindCell(ptRZ, cellId, param);
    if (status != vtkm::ErrorCode::Success)
    {
      printf("Find Cell failed! pt= %lf %lf %lf\n", ptRZ[0], ptRZ[1], ptRZ[2]);
#ifndef VTKM_CUDA
      std::cout<<"Point not found: "<<ptRZ<<std::endl;
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


  template <typename Coeff_2DType>
  VTKM_EXEC
  vtkm::FloatDefault GetPsi(const vtkm::Vec3f& ptRPZ,
                            const Coeff_2DType& coeff_2D) const
  {
    vtkm::FloatDefault R = ptRPZ[0], Z = ptRPZ[2];

    int r_i = this->GetIndex(R, this->nr, this->rmin, this->dr_inv);
    int z_i = this->GetIndex(Z, this->nz, this->zmin, this->dz_inv);

    // rc(i), zc(j)
    vtkm::FloatDefault Rc = rmin + (vtkm::FloatDefault)(r_i)*this->dr;
    vtkm::FloatDefault Zc = zmin + (vtkm::FloatDefault)(z_i)*this->dz;
    auto Rc_1 = Rc + this->dr;
    auto Zc_1 = Zc + this->dz;
    Rc = (Rc + Rc_1) * 0.5;
    Zc = (Zc + Zc_1) * 0.5;

    vtkm::Matrix<vtkm::FloatDefault, 4, 4> acoeff;
    vtkm::Id offset = (r_i * this->ncoeff + z_i) * 16;

    this->GetCoeff2D(offset, coeff_2D, acoeff);

    //taken from eval_bicub_2;
    vtkm::FloatDefault dx = R - Rc, dy = Z - Zc;
    vtkm::FloatDefault fx[4] = {0,0,0,0};
    vtkm::FloatDefault xv[4] = {1, dx, dx*dx, dx*dx*dx};
    vtkm::FloatDefault yv[4] = {1, dy, dy*dy, dy*dy*dy};

    for (int j=0; j<4; j++)
      for (int i=0; i<4; i++)
        fx[j] = fx[j] + xv[i]*acoeff[i][j];

    vtkm::FloatDefault psi = 0;
    for (int j = 0; j < 4; j++)
      psi = psi + fx[j]*yv[j];

    return psi;
  }

  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  vtkm::Vec3f GetB(vtkm::Particle& pRPZ,
                   const Coeff_1DType& coeff_1D,
                   const Coeff_2DType& coeff_2D) const
  {
    /*
    vtkm::FloatDefault psi;
    vtkm::Vec3f B0_rzp, curlB_rzp, curl_nb_rzp, gradPsi_rzp;
    vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;
    */
    this->HighOrderB(pRPZ, coeff_1D, coeff_2D); //, B0_rzp, jacobian_rzp, curlB_rzp, curl_nb_rzp, psi, gradPsi_rzp);

    tt2 = new FuncTimer(&getBT0);
    //This gives is the time derivative: Br = dR/dt, Bz= dZ/dt, B_phi/R = dphi/dt
    //We need with respect to phi:
    // dR/dphi = dR/dt / (dphi/dt) = Br / B_phi * R
    // same for z;

    vtkm::Vec3f B0_rpz(pRPZ.B0_rzp[0], pRPZ.B0_rzp[2], pRPZ.B0_rzp[1]);
    B0_rpz[0]/=pRPZ.B0_rzp[2];
    B0_rpz[2]/=pRPZ.B0_rzp[2];
    B0_rpz[0] *= pRPZ.Pos[0];
    B0_rpz[2] *= pRPZ.Pos[0];

    delete tt2;

    return B0_rpz;
  }


#if 0
********************************************************
   field following code.
********************************************************
 pIn= 2.82104 0.1 -0.022034
CHANGE ME**********                               PhiMid= 0.15
  ***** plane method: [2.82237,0.15,-0.0261434]
(0.1 0.15) dP= -0.025
     p0 = [2.82104,0.1,-0.022034]
     p0_0 = [2.82098,0.105866,-0.0220856]
     p0_1 = [2.82093,0.111732,-0.0221371]
     p0_2 = [2.82087,0.117598,-0.0221884]
     p0_3 = [2.82082,0.123464,-0.0222396]
     p0_4 = [2.82077,0.12933,-0.0222907]
     p0_5 = [2.82071,0.135196,-0.0223417]
     p0_6 = [2.82066,0.141063,-0.0223925]
     p0_7 = [2.8206,0.146929,-0.0224432]
     p0_8 = [2.82055,0.152796,-0.0224937]
     p0_9 = [2.82049,0.158663,-0.0225442]
  **** rk4 method: [2.82049,0.158663,-0.0225442]
  **correct= [3.02164,0.15,0.0625352]
********************************************************


XGC
DRP: field_following_pos2() x=    2.821035000000000 -2.2034000000000002E-002
    phi_org/dest   0.1000000014901161  0.1500000059604645
 DRP: field_following_pos2() i=             1
 DRP: field_following_pos2() x=     2.821035000000000 -2.2034000000000002E-002 phi=   0.1000000014901161
 DRP: field_following_pos2() dx1=  -2.5891391172696388E-002 -2.4717688240964472E-002
 DRP: field_following_pos2() x_tmp=    2.820711357581406 -2.2342971130636227E-002
 DRP: field_following_pos2() dx2=  -2.6251378575081735E-002 -2.4334649913031497E-002
 DRP: field_following_pos2() x_tmp=    2.820706857738474 -2.2338183151108987E-002
 DRP: field_following_pos2() dx3=  -2.6245710166807108E-002 -2.4329322981324640E-002
 DRP: field_following_pos2() x_tmp=    2.820378857187166 -2.2642233128913392E-002
 DRP: field_following_pos2() dx4=  -2.6599791429478097E-002 -2.3941212662861113E-002
 DRP: field_following_pos2() x=    2.820378810940974 -2.2642278582269915E-002
DRP: field_following_pos2() i=             2
 DRP: field_following_pos2() x=     2.820378810940974 -2.2642278582269915E-002 phi=   0.1250000037252903
 DRP: field_following_pos2() dx1=  -2.6599844382876568E-002 -2.3941157948187672E-002
 DRP: field_following_pos2() x_tmp=    2.820046312856461 -2.2941543083378590E-002
 DRP: field_following_pos2() dx2=  -2.6948178306693129E-002 -2.3547812511435621E-002
 DRP: field_following_pos2() x_tmp=    2.820041958682024 -2.2936626264979592E-002
 DRP: field_following_pos2() dx3=  -2.6942360167079921E-002 -2.3542660313524748E-002
 DRP: field_following_pos2() x_tmp=    2.819705251876576 -2.3230845142729981E-002
 DRP: field_following_pos2() dx4=  -2.7284644996469994E-002 -2.3144425562066688E-002
 DRP: field_following_pos2() x=    2.819705204354387 -2.3230889173063193E-002
 DRP: field_following_pos2() x_dest=    2.819705204354387 -2.3230889173063193E-002



#endif

  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool CalcFieldFollowingPt(const vtkm::Vec3f& pt_rpz,
                            const vtkm::Vec3f& B0_rpz,
                            const vtkm::FloatDefault& Phi0,
                            const vtkm::FloatDefault& Phi1,
                            const Coeff_1DType& coeff_1D,
                            const Coeff_2DType& coeff_2D,
                            vtkm::Vec3f& x_ff_rpz) const
  {

    tt2 = new FuncTimer(&ffPtT0);
    vtkm::FloatDefault R = pt_rpz[0];
    vtkm::FloatDefault Phi = pt_rpz[1];
    vtkm::FloatDefault Z = pt_rpz[2];

    /*
    std::cout<<"********************************************************"<<std::endl;
    std::cout<<"   field following code."<<std::endl;
    std::cout<<"********************************************************"<<std::endl;
    std::cout<<std::setprecision(16);
    std::cout<<" pIn= "<<R<<" "<<Phi<<" "<<Z<<std::endl;
    */


    vtkm::FloatDefault PhiMid = Phi0 + (Phi1-Phi0)/2.0;
    //std::cout<<std::setprecision(16)<<"   PhiMid= "<<PhiMid<<std::endl;
    vtkm::Plane<> midPlane({0, PhiMid, 0}, {0,1,0});
    Ray3f ray_rpz({R, Phi, Z}, B0_rpz);

    //Get point on mid plane.  Use the R,Z for this point for triangle finds.
    vtkm::FloatDefault RP_T;
    vtkm::Vec3f ptOnMidPlane_rpz;
    bool b;
    midPlane.Intersect(ray_rpz, RP_T, ptOnMidPlane_rpz, b);
    delete tt2;

    tt2 = new FuncTimer(&ffPtT1);

    //Now, do it using RK4 and two steps.
    vtkm::FloatDefault h = (PhiMid-Phi) / 2.0;
    //h = -h;
    vtkm::FloatDefault h_2 = h / 2.0;
    //std::cout<<"("<<Phi<<" "<<PhiMid<<") dP= "<<h<<std::endl;

    //k1 = F(p)
    //k2 = F(p+hk1/2)
    //k3 = F(p+hk2/2)
    //k4 = F(p+hk3)
    //Yn+1 = Yn + 1/6 h (k1+2k2+2k3+k4)
    vtkm::Vec3f k1, k2, k3, k4;
    vtkm::Particle p0, tmp;

    p0.Pos = {R,Phi,Z}; //pt_rpz
    for (int i = 0; i < 2; i++)
    {
      k1 = this->GetB(p0, coeff_1D, coeff_2D);
      tmp.Pos = p0.Pos + k1*h_2;

      k2 = this->GetB(tmp, coeff_1D, coeff_2D);
      tmp.Pos = p0.Pos + k2*h_2;

      k3 = this->GetB(tmp, coeff_1D, coeff_2D);
      tmp.Pos = p0.Pos + k3*h;

      k4 = this->GetB(tmp, coeff_1D, coeff_2D);

      vtkm::Vec3f vec = (k1 + 2*k2 + 2*k3 + k4) / 6.0;
      p0.Pos = p0.Pos + h * vec;
    }

    //std::cout<<"  **** rk4 method: "<<p0<<std::endl;
    x_ff_rpz = p0.Pos;
    delete tt2;

    /*
    //correct
    //vtkm::Vec3f correct(2.736638733942683, PhiMid, 0.2214040972989031);
    vtkm::Vec3f correct(3.021638837619697, PhiMid, 6.2535188082394388E-002);
    std::cout<<"  **rk4 error=  "<<(correct-p0)<<std::endl;
    std::cout<<"********************************************************"<<std::endl;
    std::cout<<"********************************************************"<<std::endl;
    */

    return true;
  }

  template <typename LocatorType, typename CellSetType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  vtkm::Vec3f CalcDeltaB(vtkm::Particle& pRPZ,
                         vtkm::Vec3f& B0_rzp,
                         vtkm::Vec3f& GRADPSI_rzp,
                         vtkm::Vec3f& curl_nb_rzp,
                         const LocatorType& locator,
                         const CellSetType& cellSet,
                         const AsFieldType& AsPhiFF,
                         const DAsFieldType& DAsPhiFF_RZP,
                         const Coeff_1DType& coeff_1D,
                         const Coeff_2DType& coeff_2D) const
  {

    auto R = pRPZ.Pos[0];
    auto Phi = pRPZ.Pos[1];
    auto Z = pRPZ.Pos[2];

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    this->GetPlaneIdx(Phi, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);
    vtkm::Vec3f B0_rpz(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

    FuncTimer *ttt = new FuncTimer(&deltaBTaa);
    vtkm::Vec3f ff_pt_rpz;
    this->CalcFieldFollowingPt({R,phiN,Z}, B0_rpz, Phi0, Phi1, coeff_1D, coeff_2D, ff_pt_rpz);
    delete ttt;

    ttt = new FuncTimer(&deltaBTa);

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

    //gradPsi: pt on mid plane?  (question)
    //dPsi/dR = B0_Z * R
    //dPsi/dZ = -B0_R * R;
    //vtkm::Vec3f gradPsi_rzp(B0_Z * x_ff_R, -B0_R * x_ff_R, 0);
    //use high order...
    vtkm::Vec3f gradPsi_rzp = GRADPSI_rzp;
    vtkm::FloatDefault gammaPsi = 1.0f/vtkm::Magnitude(gradPsi_rzp);

    vtkm::Vec2f rvec(0,0), zvec(0,0);
    rvec[0] = basis + (1.0-basis) * gammaPsi *   gradPsi_rzp[0];
    rvec[1] =         (1.0-basis) * gammaPsi *   gradPsi_rzp[1];
    zvec[0] =         (1.0-basis) * gammaPsi * (-gradPsi_rzp[1]);
    zvec[1] = basis + (1.0-basis) * gammaPsi *   gradPsi_rzp[0];
    delete ttt;

    //Get the vectors in the ff coordinates.
    //auto dAs_ff_rzp = EvalVector(ds, locator, {x_ff_rzp, x_ff_rzp}, "dAs_ff_rzp", offsets);
    //auto dAs_ff0_rzp = dAs_ff_rzp[0];
    //auto dAs_ff1_rzp = dAs_ff_rzp[1];

    ttt = new FuncTimer(&deltaBTbb);
    vtkm::Vec3f x_ff_param;
    vtkm::Vec<vtkm::Id,3> x_ff_vids;
    this->PtLoc(x_ff_rzp, locator, cellSet, x_ff_param, x_ff_vids);
    auto dAs_ff0_rzp = this->EvalV(DAsPhiFF_RZP, offsets[0], x_ff_param, x_ff_vids);
    auto dAs_ff1_rzp = this->EvalV(DAsPhiFF_RZP, offsets[1], x_ff_param, x_ff_vids);
    delete ttt;

    ttt = new FuncTimer(&deltaBTb);

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

    vtkm::FloatDefault BMag = vtkm::Magnitude(B0_rzp);
    //project using bfield.
    //gradAs.Phi = (gradAs.Phi * BMag - gradAs.R*B0_pos.R - gradAs.Z*B0_pos.Z) / B0_pos.Phi
    gradAs_rpz[1] = (gradAs_rpz[1]*BMag -gradAs_rpz[0]*B0_rzp[0] - gradAs_rpz[2]*B0_rzp[1]) / B0_rzp[2];
    delete ttt;

    //auto As_ff = InterpScalar(ds, locator, {x_ff_rzp, x_ff_rzp}, "As_ff", offsets);
    //vtkm::FloatDefault As_ff0 = As_ff[0];
    //vtkm::FloatDefault As_ff1 = As_ff[1];
    //std::cout<<"****** Call EvalS on AsPhi_ff"<<std::endl;
    auto As_ff0 = this->EvalS(AsPhiFF, offsets[0], x_ff_vids, x_ff_param);
    auto As_ff1 = this->EvalS(AsPhiFF, offsets[1], x_ff_vids, x_ff_param);

    ttt = new FuncTimer(&deltaBTc);
    vtkm::FloatDefault As = wphi[0]*As_ff0 + wphi[1]*As_ff1;
    auto AsCurl_bhat_rzp = As * curl_nb_rzp;

    //auto bhat_rzp = this->EvalV(B_Norm_RZP, 0, particlePos_param, particlePos_vids);
    auto bhat_rzp = vtkm::Normal(B0_rzp);

    vtkm::Vec3f gradAs_rzp(gradAs_rpz[0], gradAs_rpz[2], gradAs_rpz[1]);
    vtkm::Vec3f deltaB_rzp = AsCurl_bhat_rzp + vtkm::Cross(gradAs_rzp, bhat_rzp);

    delete ttt;
    return deltaB_rzp;
  }


  template <typename LocatorType, typename CellSetType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool HighOrderEval(vtkm::Particle& pRPZ,
                     const LocatorType& locator,
                     const CellSetType& cellSet,
                     const AsFieldType& AsPhiFF,
                     const DAsFieldType& DAsPhiFF_RZP,
                     const Coeff_1DType& coeff_1D,
                     const Coeff_2DType& coeff_2D,
                     vtkm::Vec3f& res) const
  {

    //std::cout<<"    HighOrderEval: "<<pRPZ.Pos<<std::endl;
    auto R = pRPZ.Pos[0];

    //res is R,P,Z
    vtkm::Vec3f B0_rzp, curlB_rzp, curl_nb_rzp, GRADPSI_rzp;
    //vtkm::Vec<vtkm::Vec3f, 3> jacobian;
    if (!this->HighOrderB(pRPZ, coeff_1D, coeff_2D)) //, B0_rzp, jacobian, curlB_rzp, curl_nb_rzp, psi, GRADPSI_rzp))
      return false;
    //vtkm::FloatDefault psi = pRPZ.Psi;
    B0_rzp = pRPZ.B0_rzp;
    curlB_rzp = pRPZ.curlB_rzp;
    curl_nb_rzp = pRPZ.curl_nb_rzp;
    GRADPSI_rzp = pRPZ.gradPsi_rzp;


    //printf("  worklet %d\n", __LINE__);
    if (this->UseBOnly)
    {
      B0_rzp[2] /= R;
      //res is R,P,Z
      res = vtkm::Vec3f(B0_rzp[0], B0_rzp[2], B0_rzp[1]);
      return true;
    }

    vtkm::Vec3f deltaB_rzp;
    deltaB_rzp = this->CalcDeltaB(pRPZ, B0_rzp, GRADPSI_rzp, curl_nb_rzp,
                                  locator, cellSet,
                                  AsPhiFF, DAsPhiFF_RZP, coeff_1D, coeff_2D);
    deltaB_rzp[2] /= R;
    B0_rzp[2] /= R;

    //std::cout<<"****************** B0= "<<B0_rzp<<std::endl;
    //std::cout<<"****************** deltaB= "<<deltaB_rzp<<std::endl;
    vtkm::Vec3f vec_rzp = B0_rzp + deltaB_rzp;
    vtkm::Vec3f vec_rpz(vec_rzp[0], vec_rzp[2], vec_rzp[1]);
    res = vec_rpz;

    return true;

#if 0
    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    this->GetPlaneIdx(Phi, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);
    vtkm::Vec3f B0_rpz(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

    //printf("  worklet %d\n", __LINE__);
    vtkm::Vec3f ff_pt_rpz;
    this->CalcFieldFollowingPt({R,phiN,Z}, B0_rpz, Phi0, Phi1, coeff_1D, coeff_2D, ff_pt_rpz);

    //Now, interpolate between Phi_i and Phi_i+1
    vtkm::FloatDefault T01 = (phiN - Phi0) / (Phi1-Phi0);
    vtkm::FloatDefault T10 = 1.0f - T01;

    //Get vec at Phi0 and Phi1.
    //x_ff is in rzp
    //vtkm::Vec3f x_ff_rzp(ptOnMidPlane_rpz[0], ptOnMidPlane_rpz[2], 0);
    vtkm::Vec3f x_ff_rzp(ff_pt_rpz[0], ff_pt_rpz[2], 0);
    //std::cout<<"x_ff_rzp= "<<x_ff_rzp<<std::endl;
    //x_ff_rzp[0] = 3.021638837619697;
    //x_ff_rzp[1] = 6.2535188082394388E-002;
    //x_ff_rzp[2] = 6.217735460229799;
    //std::cout<<"****** XGC-fixed: x_ff_rzp= "<<x_ff_rzp<<std::endl;

    int offsets[2];
    offsets[0] = planeIdx0*this->NumNodes*2;
    offsets[1] = planeIdx0*this->NumNodes*2 + this->NumNodes;
    //std::cout<<"*** planeIdx = "<<planeIdx0<<" "<<planeIdx1<<std::endl;

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
    vtkm::Vec3f gradPsi_rzp = GRADPSI_rzp;
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
    this->PtLoc(x_ff_rzp, locator, cellSet, x_ff_param, x_ff_vids);
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

    vtkm::FloatDefault BMag = vtkm::Magnitude(B0_rzp);
    //project using bfield.
    //gradAs.Phi = (gradAs.Phi * BMag - gradAs.R*B0_pos.R - gradAs.Z*B0_pos.Z) / B0_pos.Phi
    gradAs_rpz[1] = (gradAs_rpz[1]*BMag -gradAs_rpz[0]*B0_rzp[0] - gradAs_rpz[2]*B0_rzp[1]) / B0_rzp[2];

    //std::cout<<" gradAs_rpz= "<<gradAs_rpz<<std::endl;

    //deltaB = AsCurl(bhat) + gradAs x bhat.
    //std::vector<int> off = {planeIdx0*this->NumNodes};
    //vtkm::Vec3f AsCurl_bhat_rzp = EvalVector(ds, locator, {x_ff_rzp}, "AsCurlBHat_RZP", off)[0];
    //auto AsCurl_bhat_rzp = this->EvalV(AsCurlBHat_RZP, 0, x_ff_vids, x_ff_param);

#if 0
    //vtkm::Vec3f curl_nb_rzp = EvalVector(ds, locator, {ptRZ}, "curl_nb_rzp")[0];
    //std::cout<<"    pos_ids= "<<particlePos_vids<<std::endl;
    //std::cout<<"    pos_parms= "<<particlePos_param<<std::endl;
    auto curl_nb_rzp = this->EvalV(Curl_NB_RZP, 0, particlePos_param, particlePos_vids);
#endif

    //auto As_ff = InterpScalar(ds, locator, {x_ff_rzp, x_ff_rzp}, "As_ff", offsets);
    //vtkm::FloatDefault As_ff0 = As_ff[0];
    //vtkm::FloatDefault As_ff1 = As_ff[1];
    //std::cout<<"****** Call EvalS on AsPhi_ff"<<std::endl;
    auto As_ff0 = this->EvalS(AsPhiFF, offsets[0], x_ff_vids, x_ff_param);
    auto As_ff1 = this->EvalS(AsPhiFF, offsets[1], x_ff_vids, x_ff_param);

    vtkm::FloatDefault As = wphi[0]*As_ff0 + wphi[1]*As_ff1;
    auto AsCurl_bhat_rzp = As * curl_nb_rzp;

    /*
    //As = 3.0485639983994535E-006;
    std::cout<<"    psi= "<<psi<<std::endl;
    std::cout<<"    Bsa= "<<AsCurl_bhat_rzp<<std::endl;
    std::cout<<"    As= "<<As<<std::endl;
    std::cout<<"       As0= "<<As_ff0<<std::endl;
    std::cout<<"       As1= "<<As_ff1<<std::endl;
    std::cout<<"  dAs0= "<<dAs_ff0_rzp<<std::endl;
    std::cout<<"  dAs1= "<<dAs_ff1_rzp<<std::endl;
    std::cout<<"    curl_nb_rzp= "<<curl_nb_rzp<<std::endl;
    std::cout<<"  over_B= "<<(1.0/BMag)<<std::endl;
    std::cout<<"  inv_r= "<<(1.0/R)<<std::endl;
    std::cout<<"  wphi: "<<wphi[0]<<" "<<wphi[1]<<std::endl;
    std::cout<<"  x_ff_vids= "<<x_ff_vids<<std::endl;
    std::cout<<"  rvec: "<<rvec<<std::endl;
    std::cout<<"  zvec: "<<zvec<<std::endl;
    std::cout<<"  gradPsi_rzp: "<<gradPsi_rzp<<std::endl;
    std::cout<<"  gammaPsi: "<<gammaPsi<<std::endl;
    std::cout<<"  GRADPSI_rzp: "<<GRADPSI_rzp<<std::endl;
    std::cout<<"  GAMMAPSI: "<<(1.0/vtkm::Magnitude(GRADPSI_rzp))<<std::endl;
    */


    //std::cout<<"    curl_nb_rzp.size()= "<<Curl_NB_RZP.GetNumberOfValues()<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v0)= "<<Curl_NB_RZP.Get(particlePos_vids[0])<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v1)= "<<Curl_NB_RZP.Get(particlePos_vids[1])<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v2)= "<<Curl_NB_RZP.Get(particlePos_vids[2])<<std::endl;
    //std::cout<<"    AsCurl_bhat_rzp= "<<AsCurl_bhat_rzp<<std::endl;

    //vtkm::Vec3f bhat_rzp = EvalVector(ds, locator, {ptRZ}, "B_RZP_Norm")[0];
    /////////////auto bhat_rzp = this->EvalV(B_Norm_RZP, 0, particlePos_param, particlePos_vids);
    auto bhat_rzp = vtkm::Normal(B0_rzp);
    //std::cout<<"    bhat_rzp= "<<bhat_rzp<<std::endl;


    vtkm::Vec3f gradAs_rzp(gradAs_rpz[0], gradAs_rpz[2], gradAs_rpz[1]);
    //std::cout<<"    gradAs_rzp= "<<gradAs_rzp<<std::endl;
    vtkm::Vec3f deltaB_rzp = AsCurl_bhat_rzp + vtkm::Cross(gradAs_rzp, bhat_rzp);
    //std::cout<<"    deltaB= "<<deltaB_rzp<<std::endl;

    deltaB_rzp[2] /= R;
    B0_rzp[2] /= R;

    //std::cout<<"****************** B0= "<<B0_rzp<<std::endl;
    //std::cout<<"****************** deltaB= "<<deltaB_rzp<<std::endl;
    vtkm::Vec3f vec_rzp = B0_rzp + deltaB_rzp;
    vtkm::Vec3f vec_rpz(vec_rzp[0], vec_rzp[2], vec_rzp[1]);
    res = vec_rpz;
    //std::cout<<"    vec_rpz= "<<vec_rpz<<std::endl<<std::endl;
    //printf("  worklet %d\n", __LINE__);
    return true;
#endif
  }

  template <typename LocatorType, typename CellSetType, typename BFieldType, typename AsFieldType, typename DAsFieldType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  bool Evaluate(vtkm::Particle& pRPZ,
                const LocatorType& locator,
                const CellSetType& cellSet,
                const BFieldType& /*B_RZP*/,
                const BFieldType& /*B_Norm_RZP*/,
                //const BFieldType& Curl_NB_RZP,
                const AsFieldType& AsPhiFF,
                const DAsFieldType& DAsPhiFF_RZP,
                const Coeff_1DType& Coeff_1D,
                const Coeff_2DType& Coeff_2D,
                vtkm::Vec3f& res) const
  {
    //printf("  worklet %d\n", __LINE__);
    //std::cout<<"  Evaluate: "<<pRPZ.Pos<<std::endl;
    if (this->UseHighOrder)
      return this->HighOrderEval(pRPZ, locator, cellSet, AsPhiFF, DAsPhiFF_RZP, Coeff_1D, Coeff_2D, res);

#if 0
    auto R = ptRPZ[0];
    auto Phi = ptRPZ[1];
    auto Z = ptRPZ[2];

    vtkm::Id planeIdx0, planeIdx1, numRevs;
    vtkm::FloatDefault phiN, Phi0, Phi1, T;
    this->GetPlaneIdx(Phi, phiN, planeIdx0, planeIdx1, Phi0, Phi1, numRevs, T);

    vtkm::Vec3f ptRZ(R, Z, 0);

    vtkm::Vec3f particlePos_param;
    vtkm::Vec<vtkm::Id,3> particlePos_vids;
    if (!this->PtLoc(ptRZ, locator, cellSet, particlePos_param, particlePos_vids))
      return false;

    auto B0_rzp = this->EvalV(B_RZP, 0, particlePos_param, particlePos_vids);
    //std::cout<<"Meow: "<<ptRPZ<<std::endl;
    //std::cout<<"    B0= "<<B0_rzp<<std::endl;

    if (this->UseBOnly)
    {
      B0_rzp[2] /= ptRZ[0];
      res = vtkm::Vec3f(B0_rzp[0], B0_rzp[2], B0_rzp[1]);
      return true;
    }

    vtkm::Vec3f B0_rpz(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

    vtkm::FloatDefault PhiMid = Phi0 + (Phi1-Phi0)/2.0;
    vtkm::Plane<> midPlane({0, PhiMid, 0}, {0,1,0});
    Ray3f ray_rpz({R, phiN, Z}, B0_rpz);

    //Get point on mid plane.  Use the R,Z for this point for triangle finds.
    vtkm::FloatDefault RP_T;
    vtkm::Vec3f ptOnMidPlane_rpz;
    bool tmp;
    midPlane.Intersect(ray_rpz, RP_T, ptOnMidPlane_rpz, tmp);

    //Now, interpolate between Phi_i and Phi_i+1
    vtkm::FloatDefault T01 = (phiN - Phi0) / (Phi1-Phi0);
    vtkm::FloatDefault T10 = 1.0f - T01;

    //Get vec at Phi0 and Phi1.
    //x_ff is in rzp
    vtkm::Vec3f x_ff_rzp(ptOnMidPlane_rpz[0], ptOnMidPlane_rpz[2], 0);
    int offsets[2];
    offsets[0] = planeIdx0*this->NumNodes*2;
    offsets[1] = planeIdx0*this->NumNodes*2 + this->NumNodes;


    const vtkm::FloatDefault basis = 0.0f;
    auto B0_R = B0_rzp[0];
    auto B0_Z = B0_rzp[1];
    auto x_ff_R = x_ff_rzp[0];
    //auto x_ff_Z = x_ff_rzp[1];

    //gradPsi: pt on mid plane?  (question)
    //dPsi/dR = B0_Z * R
    //dPsi/dZ = -B0_R * R;
    vtkm::Vec3f gradPsi_rzp(B0_Z * x_ff_R, -B0_R * x_ff_R, 0);
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
    this->PtLoc(x_ff_rzp, locator, cellSet, x_ff_param, x_ff_vids);
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

    vtkm::FloatDefault BMag = vtkm::Magnitude(B0_rzp);
    //project using bfield.
    //gradAs.Phi = (gradAs.Phi * BMag - gradAs.R*B0_pos.R - gradAs.Z*B0_pos.Z) / B0_pos.Phi
    gradAs_rpz[1] = (gradAs_rpz[1]*BMag -gradAs_rpz[0]*B0_rzp[0] - gradAs_rpz[2]*B0_rzp[1]) / B0_rzp[2];

    //deltaB = AsCurl(bhat) + gradAs x bhat.
    //std::vector<int> off = {planeIdx0*this->NumNodes};
    //vtkm::Vec3f AsCurl_bhat_rzp = EvalVector(ds, locator, {x_ff_rzp}, "AsCurlBHat_RZP", off)[0];
    //auto AsCurl_bhat_rzp = this->EvalV(AsCurlBHat_RZP, 0, x_ff_vids, x_ff_param);

    //vtkm::Vec3f curl_nb_rzp = EvalVector(ds, locator, {ptRZ}, "curl_nb_rzp")[0];
    //std::cout<<"    pos_ids= "<<particlePos_vids<<std::endl;
    //std::cout<<"    pos_parms= "<<particlePos_param<<std::endl;
    auto curl_nb_rzp = this->EvalV(Curl_NB_RZP, 0, particlePos_param, particlePos_vids);

    //auto As_ff = InterpScalar(ds, locator, {x_ff_rzp, x_ff_rzp}, "As_ff", offsets);
    //vtkm::FloatDefault As_ff0 = As_ff[0];
    //vtkm::FloatDefault As_ff1 = As_ff[1];
    auto As_ff0 = this->EvalS(AsPhiFF, offsets[0], x_ff_vids, x_ff_param);
    auto As_ff1 = this->EvalS(AsPhiFF, offsets[1], x_ff_vids, x_ff_param);

    vtkm::FloatDefault As = wphi[0]*As_ff0 + wphi[1]*As_ff1;
    auto AsCurl_bhat_rzp = As * curl_nb_rzp;
    //std::cout<<"    As= "<<As<<std::endl;
    //std::cout<<"    curl_nb_rzp= "<<curl_nb_rzp<<std::endl;
    //std::cout<<"    curl_nb_rzp.size()= "<<Curl_NB_RZP.GetNumberOfValues()<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v0)= "<<Curl_NB_RZP.Get(particlePos_vids[0])<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v1)= "<<Curl_NB_RZP.Get(particlePos_vids[1])<<std::endl;
    //std::cout<<"    curl_nb_rzp.Get(v2)= "<<Curl_NB_RZP.Get(particlePos_vids[2])<<std::endl;
    //std::cout<<"    AsCurl_bhat_rzp= "<<AsCurl_bhat_rzp<<std::endl;

    //vtkm::Vec3f bhat_rzp = EvalVector(ds, locator, {ptRZ}, "B_RZP_Norm")[0];
    auto bhat_rzp = this->EvalV(B_Norm_RZP, 0, particlePos_param, particlePos_vids);
    //std::cout<<"    bhat_rzp= "<<bhat_rzp<<std::endl;

    vtkm::Vec3f gradAs_rzp(gradAs_rpz[0], gradAs_rpz[2], gradAs_rpz[1]);
    //std::cout<<"    gradAs_rzp= "<<gradAs_rzp<<std::endl;
    vtkm::Vec3f deltaB_rzp = AsCurl_bhat_rzp + vtkm::Cross(gradAs_rzp, bhat_rzp);

    //std::cout<<"    deltaB= "<<deltaB_rzp<<std::endl;

    deltaB_rzp[2] /= R;
    B0_rzp[2] /= R;

    vtkm::Vec3f vec_rzp = B0_rzp + deltaB_rzp;
    vtkm::Vec3f vec_rpz(vec_rzp[0], vec_rzp[2], vec_rzp[1]);
    res = vec_rpz;
    //std::cout<<"    vec_rpz= "<<vec_rpz<<std::endl<<std::endl;
#endif
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
    //rem = std::fmod(vtkm::Abs(phi), vtkm::TwoPi());
    phiN = phi;
    if (phi < 0)
    {
      //rem = -rem;
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
#endif
};
