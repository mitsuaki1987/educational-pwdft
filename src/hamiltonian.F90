module hamiltonian
  !
  implicit none
  !
contains
  !
  subroutine h_psi(kvec,psi, hpsi)
    !
    use gvec, only : g_wf, g_rh
    use solver, only : nbnd
    use fftw_wrapper, only : fft_g2r_w, fft_r2g_w
    use io_vloc, only : Vks
    use atm_spec, only : bvec
    !
    real(8),intent(in) :: kvec(3)
    complex(8),intent(in) :: psi(g_wf%npw,nbnd)
    complex(8),intent(out) :: hpsi(g_wf%npw,nbnd)
    !
    integer :: ipw, ibnd
    real(8) :: kgv(3)
    complex(8) :: psir(g_rh%nr)
    !
    ! Local potential
    !
    do ibnd = 1, nbnd
       call fft_g2r_w(psi(1:g_wf%npw,ibnd), psir)
       psir(1:g_rh%nr) = psir(1:g_rh%nr) * Vks(1:g_rh%nr)
       call fft_r2g_w(psir, hpsi(1:g_wf%npw,ibnd))
    end do
    !
    ! Kinetic energy term
    !
    do ipw = 1, g_wf%npw
       kgv(1:3) = kvec(1:3) + matmul(bvec(1:3,1:3), dble(g_wf%mill(1:3,g_wf%map(ipw))))
       hpsi(ipw,1:nbnd) = hpsi(ipw,1:nbnd) &
       &                + 0.5d0 * dot_product(kgv,kgv) * psi(ipw,1:nbnd)
    end do
    !
  end subroutine h_psi
  !
end module hamiltonian
