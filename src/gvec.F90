module gvec
  !
  implicit none
  !
  type fft
     integer :: npw
     integer :: npw3(3)
     integer,allocatable :: mill(:,:)
     real(8),allocatable :: gv(:,:)
  end type fft
  !
  type(fft),save :: rfft
  type(fft),save :: wfft
  !
contains
  !
  subroutine setup_gvec(xfft,g2cut)
    !
    use atm_spec, only : avec, bvec
    !
    real(8),intent(in) :: g2cut
    type(fft),intent(out) :: xfft
    !
    integer :: ii, i1, i2, i3, nmax3(3), nmax
    real(8) :: pi = acos(-1.0d0), norm, gv0(3), g2
    !
    do ii = 1, 3
       norm = sqrt(dot_product(avec(1:3,ii), avec(1:3,ii)))
       nmax3(ii) = ceiling(sqrt(g2cut)*norm/(2.0d0*pi))
    end do
    !
    nmax = product(nmax3(1:3)*2+1)
    allocate(xfft%gv(3,nmax), xfft%mill(3,nmax))
    !
    xfft%npw = 0
    do i1 = -nmax3(1), nmax3(1)
       do i2 = -nmax3(2), nmax3(2)
          do i3 = -nmax3(3), nmax3(3)
             gv0(1:3) = matmul(bvec(1:3,1:3), dble((/i1,i2,i3/)))
             g2 = dot_product(gv0,gv0)
             if(g2 < g2cut) then
                xfft%npw = xfft%npw + 1
                xfft%mill(1:3,xfft%npw) = (/i1,i2,i3/)
                xfft%gv(1:3,xfft%npw) = gv0(1:3)
             end if
          end do
       end do
    end do
    !
    write(*,*) "    Numver of PW : ", xfft%npw
    !
    ! Find FFT grid
    !
    do ii = 1, 3
       xfft%npw3(ii) = maxval(xfft%mill(ii,1:xfft%npw))*2 + 1
       call base2357(xfft%npw3(ii))
    end do
    !
    write(*,*) "    FFT grid : ", xfft%npw3(1:3)
    !
  end subroutine setup_gvec
  !
  subroutine base2357(nfft)
    !
    integer,intent(inout) :: nfft
    !
    integer :: base(4) = (/2,3,5,7/), ibase, iexp, &
    &          start, last, nexp, nfft1, nfft2
    !
    start = nfft
    last = 2 * nfft
    !
    do nfft1 = start, last
       !
       nfft2 = nfft1
       !
       do ibase = 1, 4
          !
          nexp = ceiling(log(dble(nfft2))/log(dble(base(ibase))))
          !
          do iexp = 1, nexp
             if(mod(nfft2, base(ibase)) == 0) then
                nfft2 = nfft2 / base(ibase)
             else
                exit
             end if
          end do
          !
          if(nfft2 == 1) then
             nfft = nfft1
             return
          end if
          !
       end do
       !
    end do
    !
  end subroutine base2357
  !
end module gvec
