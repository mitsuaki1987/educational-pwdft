!
! Copyright (c) 2018 Mitsuaki Kawamura
!
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
module diag_direct
  !
  implicit none
  !
contains
  !
  subroutine direct(npw,kvec,evec,eval)
    !
    use atm_spec, only : bvec
    use kohn_sham, only : nbnd
    use gvec, only : g_rh, g_wf
    use rho_v, only : Vks
    use fftw_wrapper, only : fft_r2g
    !
    integer,intent(in) :: npw
    real(8),intent(in) :: kvec(3)
    complex(8),intent(out) :: evec(npw,nbnd)
    real(8),intent(out) :: eval(nbnd)
    !
    integer :: info, ipw, jpw, lwork, dmill(3)
    real(8) :: rwork(3*npw-2), eval_full(npw), kgv(3)
    complex(8) :: ham(npw,npw), &
    &             VksG(g_rh%nft(1), g_rh%nft(2), g_rh%nft(3))
    complex(8),allocatable :: work(:)
    !
    ! Local potential term
    !
    call fft_r2G(Vks, VksG)
    !
    do ipw = 1, npw
       do jpw = 1, npw
          dmill(1:3) = g_wf%mill(1:3,g_wf%map(jpw)) - g_wf%mill(1:3,g_wf%map(ipw))
          dmill(1:3) = modulo(dmill(1:3), g_rh%nft(1:3)) + 1
          ham(jpw,ipw) = VksG(dmill(1), dmill(2), dmill(3))
       end do
    end do
    !
    ! Kinetic energy term
    !
    do ipw = 1, npw
       kgv(1:3) = kvec(1:3) + matmul(bvec(1:3,1:3), dble(g_wf%mill(1:3,g_wf%map(ipw))))
       ham(ipw,ipw) = ham(ipw,ipw) + 0.5d0 * dot_product(kgv,kgv)
    end do
    !
    lwork = -1
    allocate(work(1))
    call zheev('V', 'U', npw, ham, npw, eval_full, work, lwork, rwork, info)
    lwork = nint(dble(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call zheev('V', 'U', npw, ham, npw, eval_full, work, lwork, rwork, info)
    deallocate(work)
    !
    eval(      1:nbnd) = eval_full(1:nbnd)
    evec(1:npw,1:nbnd) = ham(1:npw,1:nbnd)
    !
  end subroutine direct
  !
end module diag_direct
