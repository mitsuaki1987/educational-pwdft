module rho_v
  !
  implicit none
  !
  real(8),allocatable :: &
  & Vks(:), & !< (g_rh%nr) Kohn-Sham potential [Htr]
  & Vps(:), & !< (g_rh%nr) Pseudopotential [Htr]
  & rho(:) !< (g_rh%nr) Charge density
  !
contains
  !>
  !! Initialize rho, Vps, Vks
  !!
  subroutine init_rho_V()
    !
    use gvec, only : g_rh
    use atm_spec, only : Vcell, nelec
    !
    allocate(Vps(g_rh%nr), Vks(g_rh%nr), rho(g_rh%nr))
    !
    rho(1:g_rh%nr) = nelec / Vcell
    !
    call generate_vps()
    Vks(1:g_rh%nr) = Vps(1:g_rh%nr)
    call hartree_pot(Vks)
    call xc_pot(Vks)
    !
  end subroutine init_rho_V
  !>
  !! Pseudopotential Vps is computed
  !!
  subroutine generate_vps()
    !
    use constant, only : pi
    use gvec, only : g_rh
    use atm_spec, only : bvec, spec, atm, nat, ntyp, Vcell
    use fftw_wrapper, only : fft_g2r
    !
    integer :: jtyp, iat, nat2, ipw
    real(8) :: prod(nat), pos2(3,nat), gv(3), ctail, glen, dr, &
    &          frmfac(g_rh%nr)
    complex(8) :: strfac(g_rh%nr), VpsG(g_rh%nr)
    real(8),allocatable :: sinc(:)
    !
    VpsG(1:g_rh%nr) = 0.0d0
    !
    do jtyp = 1, ntyp
       !
       nat2 = 0
       do iat = 1, nat
          if(atm(iat)%ityp == jtyp) then
             nat2 = nat2 + 1
             pos2(1:3,nat2) = atm(iat)%pos(1:3)
          end if
       end do
       !
       ! Structure factor
       !
       do ipw = 1, g_rh%nr
          prod(1:nat2) = 2.0d0 * pi * matmul(dble(g_rh%mill(1:3,ipw)),pos2(1:3,1:nat2))
          strfac(ipw) = sum(exp(cmplx(0.0d0, prod(1:nat), 8))) / Vcell
       end do
       !
       ! Form factor
       !
       allocate(sinc(spec(jtyp)%mmax))
       dr = spec(jtyp)%psr(spec(jtyp)%mmax) / dble(spec(jtyp)%mmax - 1)
       !
       ! G = 0 : Compensation to the Hartree potential 
       !
       ctail = 0.5d0 * spec(jtyp)%Zion * spec(jtyp)%psr(spec(jtyp)%mmax)**2
       !
       sinc(1:spec(jtyp)%mmax) = dr
       sinc(1              ) = dr * 0.5d0
       sinc(spec(jtyp)%mmax) = dr * 0.5d0
       !
       frmfac(1) = sum( sinc(1:spec(jtyp)%mmax) &
       &    * spec(jtyp)%psr(1:spec(jtyp)%mmax)**2 &
       &    * spec(jtyp)%psV(1:spec(jtyp)%mmax)) &
       &         + ctail
       frmfac(1) = 4.0d0 * pi * frmfac(1)
       !
       ! G /= 0
       !
       do ipw = 2, g_rh%nr
          !
          gv(1:3) = matmul(bvec(1:3,1:3), dble(g_rh%mill(1:3,ipw)))
          glen = sqrt(dot_product(gv(1:3),gv(1:3)))
          sinc(1:spec(jtyp)%mmax) = glen*spec(jtyp)%psr(1:spec(jtyp)%mmax)
          !
          ctail = - spec(jtyp)%Zion * cos(sinc(spec(jtyp)%mmax)) / glen**2
          !
          sinc(2:spec(jtyp)%mmax) = sin(sinc(2:spec(jtyp)%mmax)) &
          &                       / sinc(    2:spec(jtyp)%mmax) &
          &                       * dr
          sinc(1              ) = dr * 0.5d0
          sinc(spec(jtyp)%mmax) = sinc(spec(jtyp)%mmax) * 0.5d0
          !
          frmfac(ipw) = sum( sinc(1:spec(jtyp)%mmax) &
          &      * spec(jtyp)%psr(1:spec(jtyp)%mmax)**2 &
          &      * spec(jtyp)%psV(1:spec(jtyp)%mmax)) &
          &           + ctail
          frmfac(ipw) = 4.0d0 * pi * frmfac(ipw)
          !
       end do
       !
       deallocate(sinc)
       !
       VpsG(1:g_rh%nr) = VpsG(1:g_rh%nr) + strfac(1:g_rh%nr) * frmfac(1:g_rh%nr)
       !
    end do
    !
    call fft_g2r(VpsG, Vps)
    !
  end subroutine generate_vps
  !>
  !! Add Hartree potential
  !!
  subroutine hartree_pot(Vloc)
    !
    use constant, only : pi
    use gvec, only : g_rh
    use atm_spec, only : bvec
    use fftw_wrapper, only : fft_g2r, fft_r2g
    !
    real(8),intent(inout) :: Vloc(g_rh%nr)
    !
    integer :: ir
    real(8) :: g3(3), VH(g_rh%nr)
    complex(8) :: rhog(g_rh%nr)
    !
    call fft_r2G(rho, rhoG)
    !
    ! G = 0 : Compensation to the ionic potential 
    !
    rhog(1) = 0.0d0
    !
    do ir = 2, g_rh%nr
       g3(1:3) = matmul(bvec(1:3,1:3), dble(g_rh%mill(1:3,ir)))
       rhog(ir) = 4.0d0 * pi * rhog(ir) / dot_product(g3(1:3),g3(1:3))
    end do
    !
    call fft_G2r(rhoG, VH)
    !
    Vloc(1:g_rh%nr) = Vloc(1:g_rh%nr) + VH(1:g_rh%nr)
    !
  end subroutine hartree_pot
  !>
  !! Add XC potential (LDA)
  !!
  subroutine xc_pot(Vloc)
    !
    use constant, only : pi
    use gvec, only : g_rh
    !
    real(8),intent(inout) :: Vloc(g_rh%nr)
    !
    integer :: ir
    real(8) :: rs, exc, Vxc
    !
    do ir = 1, g_rh%nr
       !
       if(rho(ir) > 1.0d-10) then
          !
          rs = (3.0d0 / (4.0d0 * pi * rho(ir)))**(1.0d0/3.0d0)
          !
          if(rs < 1.0d0) then
             exc = -3.0d0 / (4.0d0*pi) * (9.0d0*pi/4)**(1.0d0/3.0d0) / rs &
             &   - 0.0480d0 + 0.031d0*log(rs) - 0.0116d0*rs + 0.0020d0*rs*log(rs)
             !
             vxc = exc -1.0d0 / (4.0d0*pi) * (9.0d0*pi/4)**(1.0d0/3.0d0) / rs**2 &
             & -0.0096d0*rs + 0.031d0 + 0.002d0*rs*log(rs)
          else
             exc = -3.0d0 / (4.0d0*pi) * (9.0d0*pi/4)**(1.0d0/3.0d0) / rs &
             &   - 0.1423d0 / (1.0d0 + 1.0529d0*sqrt(rs) + 0.3334d0*rs)
             !
             vxc = exc -1.0d0 / (4.0d0*pi) * (9.0d0*pi/4)**(1.0d0/3.0d0) / rs**2 &
             & - 0.0474333d0*(0.3334d0*rs + 0.52645*sqrt(rs)) &
             & /(1.0d0 + 1.0529d0*sqrt(rs) + 0.3334d0*rs)**2
          end if
          !
          Vloc(ir) = Vloc(ir) + vxc
          !
       end if
       !
    end do
    !
  end subroutine xc_pot
  !
end module rho_v
