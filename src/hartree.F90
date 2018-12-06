module hartree
  !
  implicit none
  !
contains
  !
  subroutine hartree_pot(Vloc)
    !
    use constant, only : pi
    use gvec, only : g_rh
    use atm_spec, only : bvec
    use io_vloc, only : rho
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
  !
  subroutine xc_pot(Vloc)
    !
    use constant, only : pi
    use io_vloc, only : rho
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
end module hartree
