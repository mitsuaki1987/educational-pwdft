module fftw_wrapper
  !
  implicit none
  !
  integer(8),save :: &
  & plan_g2r, & !< FFTW plan for F(G)e^{iGr} -> f(r) (backward)
  & plan_r2g !< FFTW plan for f(r)e^{-iGr} -> F(G) (forward)
  !
  integer,allocatable,save :: &
  & w2r(:) !< (g_wf%npw) g_wf%npw -> g_rh%nr
  !
  complex(8),allocatable,save :: &
  & fft_buffer(:) !< (g_rh%npw) FFT buffer
  !
  include "fftw3.f"
  !
  private plan_g2r, plan_r2g, fft_buffer
  !
contains
  !
  subroutine init_fft()
    !
    use gvec, only : g_rh, g_wf
    !
    integer :: ipw, igv(3)
    !
    allocate(fft_buffer(g_rh%npw))
    !
    call dfftw_plan_dft_3d(plan_G2r, g_rh%nft(1), g_rh%nft(2), g_rh%nft(3), &
    &                      fft_buffer, fft_buffer, &
    &                      fftw_backward, fftw_estimate)
    !
    call dfftw_plan_dft_3d(plan_r2G, g_rh%nft(1), g_rh%nft(2), g_rh%nft(3), &
    &                      fft_buffer, fft_buffer, &
    &                      fftw_forward, fftw_estimate)
    !
    allocate(w2r(g_wf%npw))
    do ipw = 1, g_wf%npw
       igv(1:3) = g_wf%mill(1:3, g_wf%map(ipw))
       igv(1:3) = modulo(igv(1:3), g_rh%nft(1:3))
       w2r(ipw) = 1 + igv(1) + g_rh%nft(1)*(igv(2) + g_rh%nft(2)*igv(3))
    end do
    !
  end subroutine init_fft
  !>
  !! v(r) e^{-iGr} -> V(G)
  !!
  subroutine fft_r2g(VlocR, VlocG)
    !
    use gvec, only : g_rh
    !
    real(8),intent(in) :: VlocR(g_rh%npw)
    complex(8),intent(out) :: VlocG(g_rh%npw)
    !
    fft_buffer(1:g_rh%npw) = VlocR(1:g_rh%npw)
    call dfftw_execute_dft(plan_r2g, fft_buffer, fft_buffer)
    VlocG(1:g_rh%npw) = fft_buffer(1:g_rh%npw) / dble(g_rh%nr)
    !
  end subroutine fft_r2g
  !>
  !! V(G) e^{iGr} -> v(r)
  !!
  subroutine fft_g2r(VlocG, VlocR)
    !
    use gvec, only : g_rh
    !
    complex(8),intent(in) :: VlocG(g_rh%npw)
    real(8),intent(out) :: VlocR(g_rh%npw)
    !
    fft_buffer(1:g_rh%npw) = VlocG(1:g_rh%npw)
    call dfftw_execute_dft(plan_g2r, fft_buffer, fft_buffer)
    VlocR(1:g_rh%npw) = dble(fft_buffer(1:g_rh%npw))
    !
  end subroutine fft_g2r
  !>
  !! w(r) e^{-iGr} -> W(G) -> compress
  !!
  subroutine fft_r2g_w(wfR, wfG)
    !
    use atm_spec, only : Vcell
    use gvec, only : g_rh, g_wf
    !
    complex(8),intent(in) :: wfR(g_rh%npw)
    complex(8),intent(out) :: wfG(g_wf%npw)
    !
    fft_buffer(1:g_rh%npw) = wfR(1:g_rh%npw)
    call dfftw_execute_dft(plan_r2g, fft_buffer, fft_buffer)
    wfG(1:g_wf%npw) = fft_buffer(w2r(1:g_wf%npw)) * sqrt(Vcell) / dble(g_rh%nr)
    !
  end subroutine fft_r2g_w
  !>
  !! Uncompress -> W(G) e^{iGr} -> w(r)
  !!
  subroutine fft_g2r_w(wfG, wfR)
    !
    use atm_spec, only : Vcell
    use gvec, only : g_rh, g_wf
    !
    complex(8),intent(in) :: wfG(g_wf%npw)
    complex(8),intent(out) :: wfR(g_rh%npw)
    !
    fft_buffer(1:g_rh%npw) = 0.0d0
    fft_buffer(w2r(1:g_wf%npw)) = wfG(1:g_wf%npw)
    call dfftw_execute_dft(plan_g2r, fft_buffer, fft_buffer)
    wfR(1:g_rh%npw) = fft_buffer(1:g_rh%npw) / sqrt(Vcell)
    !
  end subroutine fft_g2r_w
  !
end module fftw_wrapper
