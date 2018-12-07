module k_point
  !
  implicit none
  !
  integer,save :: nk
  integer,save :: kgrd(3)
  integer,allocatable,save :: tetra(:,:)
  real(8),allocatable,save :: kvec(:,:)
  !
contains
  !
  subroutine ksum_rho()
    !
    use atm_spec, only : bvec, nelec
    use kohn_sham, only : nbnd, ef, eval, evec
    use fftw_wrapper, only : fft_g2r_w
    use rho_v, only : rho
    use gvec, only : g_rh, g_wf
    use libtetrabz, only : libtetrabz_fermieng
    !
    integer :: ik, ibnd
    real(8) :: occ(nbnd,nk)
    complex(8) :: psir(g_rh%nr)
    !
    CALL libtetrabz_fermieng(2,bvec,nbnd,kgrd,eval,kgrd,occ,ef,nelec*0.5d0)
    !
    rho(1:g_rh%nr) = 0.0d0
    do ik = 1, nk
       do ibnd = 1, nbnd
          call fft_g2r_w(evec(1:g_wf%npw,ibnd,ik), psir)
          rho(1:g_rh%nr) = rho(1:g_rh%nr) &
          &              + dble(conjg(psir(1:g_rh%nr))*psir(1:g_rh%nr))*occ(ibnd,ik)*2.0d0
       end do
    end do
    !
  end subroutine ksum_rho
  !
end module k_point
