module hamiltonian
  !
  implicit none
  !
contains
  !
  subroutine h_psi(kvec,psi, hpsi)
    !
    use gvec, only : wfft
    use solver, only : nbnd
    !
    real(8),intent(in) :: kvec(3)
    complex(8),intent(in) :: psi(wfft%npw,nbnd)
    complex(8),intent(out) :: hpsi(wfft%npw,nbnd)
    !
    integer :: ipw
    real(8) :: kgv(3)
    !
    hpsi(1:wfft%npw,1:nbnd) = 0.0d0
    !
    ! Kinetic energy term
    !
    do ipw = 1, wfft%npw
       kgv(1:3) = kvec(1:3) + wfft%gv(1:3,ipw)
       hpsi(ipw,1:nbnd) = 0.5d0 * dot_product(kgv,kgv) * psi(ipw,1:nbnd)
    end do
    !    
  end subroutine h_psi
  !
end module hamiltonian
