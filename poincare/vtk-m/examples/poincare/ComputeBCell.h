#include <vtkm/worklet/WorkletMapTopology.h>

//#define PRINT_STUFF
#ifdef PRINT_STUFF
#define DBG(x) std::cout<<x
#else
#define DBG(x)
#endif

#define DO_TRACES 0

//-----------------------------------------------------------------------------
class ComputeBCellWorklet : public vtkm::worklet::WorkletVisitCellsWithPoints
{
  using FloatPortal = typename vtkm::cont::ArrayHandle<vtkm::FloatDefault>::ReadPortalType;
  using VecPortal = typename vtkm::cont::ArrayHandle<vtkm::Vec3f>::ReadPortalType;
public:
  using ControlSignature = void(CellSetIn Cells,
                                WholeArrayIn Coords,
                                WholeArrayIn coeff_1D,
                                WholeArrayIn coeff_2D,
                                FieldOut Psi,
                                FieldOut B0,
                                FieldOut CurlB0,
                                FieldOut CurlNB0,
                                FieldOut GradPsi);

  using ExecutionSignature = void(InputIndex, PointCount, PointIndices, _2, _3, _4, _5, _6, _7, _8, _9);
  using InputDomain = _1;

  ~ComputeBCellWorklet()
  {
  }

  ComputeBCellWorklet()
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

  template <typename CoordType, typename VertexIndexType, typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC void operator()(const vtkm::Id& /*cellId*/,
                            const vtkm::IdComponent& numPoints,
                            const VertexIndexType& ptIndices,
                            const CoordType& coords,
                            const Coeff_1DType& Coeff_1D,
                            const Coeff_2DType& Coeff_2D,
                            vtkm::FloatDefault& Psi,
                            vtkm::Vec3f& B0,
                            vtkm::Vec3f& CurlB0,
                            vtkm::Vec3f& CurlNB0,
                            vtkm::Vec3f& GradPsi) const
  {
    //Get the center point.
    vtkm::FloatDefault R=0, Z=0;
    for (vtkm::IdComponent i = 0; i < numPoints; i++)
    {
      auto coordRZ = coords.Get(ptIndices[i]);
      R += coordRZ[0];
      Z += coordRZ[1];
    }
    auto div = 1/static_cast<vtkm::FloatDefault>(numPoints);
    R *= div;
    Z *= div;

    vtkm::Vec3f ptRPZ(R, 0, Z);
    vtkm::Vec<vtkm::Vec3f, 3> jacobian_rzp;
    this->HighOrderB(ptRPZ, Coeff_1D, Coeff_2D,
                     B0, jacobian_rzp, CurlB0, CurlNB0, Psi, GradPsi);

    /*
    if (cellId == 31174180) //5228208)
    {
      std::cout<<std::setprecision(20)
        <<"ComputeB: "<<ptRPZ<<" :: "<<Psi<<" "<<B0<<" "<<CurlB0<<" "<<CurlNB0<<" "<<GradPsi<<std::endl;
    }
    */

    /*
    auto divergence = vtkm::Abs(jacobian_rzp[0][0] + B0[0]/ptRPZ[0] + jacobian_rzp[1][1]);
    if (vtkm::Abs(divergence) > 1e-15)
      std::cout<<"ComputeBCell:: divergence is BAD: "<<divergence<<std::endl;
    */
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
  bool HighOrderB(const vtkm::Vec3f& ptRPZ,
                  const Coeff_1DType& Coeff_1D,
                  const Coeff_2DType& Coeff_2D,
                  vtkm::Vec3f& B0_rzp,
                  vtkm::Vec<vtkm::Vec3f, 3>& jacobian_rzp,
                  vtkm::Vec3f& curlB_rzp,
                  vtkm::Vec3f& curl_nb_rzp,
                  vtkm::FloatDefault& PSI,
                  vtkm::Vec3f& gradPsi_rzp) const
  {
    vtkm::FloatDefault R = ptRPZ[0], Z = ptRPZ[2], P = ptRPZ[1];
    vtkm::Vec3f ptRZ(R,Z,0);

    /*
    if (0)
    {
      vtkm::Vec3f particlePos_param;
      vtkm::Vec<vtkm::Id,3> particlePos_vids;
      if (!this->PtLoc(ptRZ, locator, cellSet, particlePos_param, particlePos_vids))
        return false;
      std::cout<<"B: ids= "<<particlePos_vids<<std::endl;
      std::cout<<"B: par= "<<particlePos_param<<std::endl;
      auto B0_rzp = this->EvalV(B_RZP, 0, particlePos_param, particlePos_vids);
      B0_rzp[2] /= R;
      auto res_rpz = vtkm::Vec3f(B0_rzp[0], B0_rzp[2], B0_rzp[1]);

      std::cout<<"B_rpz("<<R<<" "<<Z<<") = "<<res_rpz<<std::endl;
    }
    */

    int r_i = this->GetIndex(R, this->nr, this->rmin, this->dr_inv);
    int z_i = this->GetIndex(Z, this->nz, this->zmin, this->dz_inv);
    //std::cout<<"r_i, z_i = "<<r_i<<" "<<z_i<<std::endl;
    /*
    r_i = 76-1;
    z_i = 76-1;
    r_i = 88-1;
    z_i = 76-1;
    */

    // rc(i), zc(j)
    vtkm::FloatDefault Rc = rmin + (vtkm::FloatDefault)(r_i)*this->dr;
    vtkm::FloatDefault Zc = zmin + (vtkm::FloatDefault)(z_i)*this->dz;
    auto Rc_1 = Rc + this->dr;
    auto Zc_1 = Zc + this->dz;
    Rc = (Rc + Rc_1) * 0.5;
    Zc = (Zc + Zc_1) * 0.5;
    /*
    Rc = this->rmin + (this->rmax-this->rmin)/(this->nr-1) * vtkm::FloatDefault(r_i+1);
    Zc = this->zmin + (this->zmax-this->zmin)/(this->nr-1) * vtkm::FloatDefault(z_i+1);
    auto Rc_1 = this->rmin + (this->rmax-this->rmin)/(this->nr-1) * vtkm::FloatDefault(r_i+2);
    auto Zc_1 = this->zmin + (this->zmax-this->zmin)/(this->nr-1) * vtkm::FloatDefault(z_i+2);
    Rc = (Rc+Rc_1)/2.;
    Zc = (Zc+Zc_1)/2.;

    std::cout<<"Fix me: "<<__LINE__<<std::endl;
    std::cout<<"Rc, Zc= "<<Rc<<" "<<Zc<<std::endl;

    Rc = 2.999976000000000;
    Zc = 7.9990400000000666E-003;

    std::cout<<"SHOULD BE: Rc, Zc= "<<Rc<<" "<<Zc<<std::endl;
    std::cout<<"(r_i,z_i)= "<<r_i<<" "<<z_i<<"  Rc,Zc= "<<Rc<<" "<<Zc<<std::endl;
    */

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
    PSI = psi;
    gradPsi_rzp[0] = dpsi_dr;
    gradPsi_rzp[1] = dpsi_dz;
    gradPsi_rzp[2] = 0;
    /*
    std::cout<<" psi= "<<psi<<std::endl;
    std::cout<<" dpsi_dr = "<<dpsi_dr<<std::endl;
    std::cout<<" dpsi_dz = "<<dpsi_dz<<std::endl;
    */

    vtkm::FloatDefault fld_I = this->I_interpol(psi, 0, Coeff_1D);
    vtkm::FloatDefault fld_dIdpsi = this->I_interpol(psi, 1, Coeff_1D);

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
    curlB_rzp[0] = dBz_dp * over_r - dBp_dz;
    curlB_rzp[1] = Bp*over_r + dBp_dr - dBr_dp*over_r;
    curlB_rzp[2] = dBr_dz - dBz_dr;
    //std::cout<<"curl_B_rzp= "<<curlB_rzp<<std::endl;

    //calculate curl_nb
    /*
    !curl of norm b
    curl_nb(1) = curl_B(1)*over_B + (fld%bphi * dbdz                  ) * over_B2
    curl_nb(2) = curl_B(2)*over_B + (                - fld%bphi * dbdr) * over_B2
    curl_nb(3) = curl_B(3)*over_B + (fld%bz   * dbdr - fld%br   * dbdz) * over_B2
    */

    vtkm::FloatDefault Bmag = vtkm::Magnitude(B0_rzp);
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
    curl_nb_rzp[0] = curlB_rzp[0] * over_B + ( Bp * dBdz)*over_B2;
    curl_nb_rzp[1] = curlB_rzp[1] * over_B + (-Bp * dBdr)*over_B2;
    curl_nb_rzp[2] = curlB_rzp[2] * over_B + (Bz*dBdr - Br*dBdz)*over_B2;
    //std::cout<<"curl_nb_rzp= "<<curl_nb_rzp<<std::endl;


    /*
    std::cout<<"XXXXXX_B_rpz= "<<res<<std::endl;
    std::cout<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    std::cout<<"***************************************"<<std::endl;
    */
    //printf("  worklet %d\n", __LINE__);
    return true;
  }

  template <typename Coeff_1DType, typename Coeff_2DType>
  VTKM_EXEC
  void HighOrderB2(const vtkm::Vec3f& ptRPZ,
                  const Coeff_1DType& Coeff_1D,
                  const Coeff_2DType& Coeff_2D,
                  vtkm::Vec3f& B0_rzp,
                  vtkm::Vec<vtkm::Vec3f, 3>& jacobian_rzp,
                  vtkm::Vec3f& curlB_rzp,
                  vtkm::Vec3f& curl_nb_rzp,
                  vtkm::FloatDefault& PSI,
                  vtkm::Vec3f& gradPsi_rzp,
                  vtkm::FloatDefault& _dbr_dr,
                  vtkm::FloatDefault& _dbz_dz) const
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
    this->eval_bicub_2(R, Z, Rc, Zc, acoeff, psi, dpsi_dr, dpsi_dz, d2psi_drdz, d2psi_d2r, d2psi_d2z);
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
    _dbr_dr = dBr_dr;
    _dbz_dz = dBz_dz;
    auto divergence = dBr_dr + Br/R + dBz_dz;
    if (vtkm::Abs(divergence) > 1e-15)
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
  }


  int nr, nz;
  int ncoeff;
  vtkm::FloatDefault rmin, zmin, rmax, zmax;
  vtkm::FloatDefault dr, dz, dr_inv, dz_inv;
  vtkm::FloatDefault min_psi, max_psi;
  vtkm::FloatDefault one_d_cub_dpsi_inv;
  vtkm::FloatDefault sml_bp_sign = -1.0f;
};
