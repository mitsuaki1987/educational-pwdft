module diag_direct
  !
  implicit none
  !
contains
  !
  subroutine direct(kvec,evec,eval)
    !
    use solver, only : nbnd
    use gvec, only : rfft, wfft
    use io_vloc, only : Vloc
    !
    real(8),intent(in) :: kvec(3)
    complex(8),intent(out) :: evec(wfft%npw,nbnd)
    complex(8),intent(out) :: eval(nbnd)
    !
    integer :: info, ipw, jpw, lwork, dmill(3)
    real(8) :: rwork(3*wfft%npw-2), eval_full(wfft%npw), kgv(3)
    complex(8) :: ham(wfft%npw,wfft%npw), &
    &             VlocG(rfft%npw3(1), rfft%npw3(2), rfft%npw3(3))
    complex(8),allocatable :: work(:)
    !
    include 'fftw3.f'
    integer(8) :: plan
    !
    ! Local potential term
    !
    call dfftw_plan_dft_3d(plan, rfft%npw3(1), rfft%npw3(2), rfft%npw3(3), Vloc, VlocG, &
    &                      fftw_forward, fftw_estimate)
    call dfftw_execute_dft(plan, Vloc, VlocG)
    call dfftw_destroy_plan(plan)
    !
    do ipw = 1, wfft%npw
       do jpw = 1, wfft%npw
          dmill(1:3) = wfft%mill(1:3,jpw) - wfft%mill(1:3,ipw)
          dmill(1:3) = modulo(dmill(1:3), rfft%npw3(1:3)) + 1
          ham(jpw,ipw) = VlocG(dmill(1), dmill(2), dmill(3))
       end do
    end do
    !
    ! Kinetic energy term
    !
    do ipw = 1, wfft%npw
       kgv(1:3) = kvec(1:3) + wfft%gv(1:3,ipw)
       ham(ipw,ipw) = ham(ipw,ipw) + 0.5d0 * dot_product(kgv,kgv)
    end do
    !
    lwork = -1
    allocate(work(1))
    call zheev('V', 'U', wfft%npw, ham, wfft%npw, eval_full, work, lwork, rwork, info)
    lwork = nint(dble(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call zheev('V', 'U', wfft%npw, ham, wfft%npw, eval_full, work, lwork, rwork, info)
    deallocate(work)
    !
    eval(           1:nbnd) = eval_full(     1:nbnd)
    evec(1:wfft%npw,1:nbnd) = ham(1:wfft%npw,1:nbnd)
    !
  end subroutine direct
  !
end module diag_direct
