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
module gvec
  !
  implicit none
  !
  type fft
     integer :: &
     & nft(3), & !< FFT grid
     & nr, &  !< Number of r in whole grid, product(nft(1:3))
     & npw !< Number of PW inside sphere
     integer,allocatable :: &
     & map(:), & !< (npw) npw(sphere) -> nr(grid) mapper
     & mill(:,:) !< (3,npw) Miller index  
  end type fft
  !
  type(fft),save :: g_rh
  type(fft),save :: g_wf
  !
contains
  !
  subroutine setup_gvec(xfft,g2cut)
    !
    use constant, only : pi
    use atm_spec, only : avec, bvec
    !
    real(8),intent(in) :: g2cut
    type(fft),intent(out) :: xfft
    !
    integer :: ii, i1, i2, i3, nmax3(3), nmax
    real(8) :: norm, gv0(3), g2
    integer,allocatable :: mill0(:,:)
    !
    do ii = 1, 3
       norm = sqrt(dot_product(avec(1:3,ii), avec(1:3,ii)))
       nmax3(ii) = ceiling(sqrt(g2cut)*norm/(2.0d0*pi))
    end do
    !
    nmax = product(nmax3(1:3)*2+1)
    allocate(mill0(3,nmax))
    !
    xfft%npw = 0
    do i3 = -nmax3(3), nmax3(3)
       do i2 = -nmax3(2), nmax3(2)
          do i1 = -nmax3(1), nmax3(1)
             gv0(1:3) = matmul(bvec(1:3,1:3), dble((/i1,i2,i3/)))
             g2 = dot_product(gv0,gv0)
             if(g2 < g2cut) then
                xfft%npw = xfft%npw + 1
                mill0(1:3,xfft%npw) = (/i1,i2,i3/)
             end if
          end do
       end do
    end do
    !
    write(*,*) "    Numver of PW inside sphere : ", xfft%npw
    !
    ! Find FFT grid
    !
    do ii = 1, 3
       xfft%nft(ii) = maxval(mill0(ii,1:xfft%npw))*2 + 1
       call base2357(xfft%nft(ii))
    end do
    !
    write(*,*) "    FFT grid : ", xfft%nft(1:3)
    xfft%nr = product(xfft%nft(1:3))
    write(*,*) "    Numver of PW in grid : ", xfft%nr
    allocate(xfft%mill(3,xfft%nr), xfft%map(xfft%npw))
    !
    ! Whole grid : Miller index
    !
    ii = 0
    do i3 = 0, xfft%nft(3) - 1
       do i2 = 0, xfft%nft(2) - 1
          do i1 = 0, xfft%nft(1) - 1
             ii = ii + 1
             xfft%mill(1:3,ii) = modulo((/i1, i2, i3/) + xfft%nft(1:3)/2, xfft%nft(1:3)) &
             &                                         - xfft%nft(1:3)/2
          end do
       end do
    end do
    !
    ! Mapper
    !
    do ii = 1, xfft%npw
       mill0(1:3,ii) = modulo(mill0(1:3,ii), xfft%nft(1:3))
       xfft%map(ii) = 1 + mill0(1,ii) + xfft%nft(1)*(mill0(2,ii) + xfft%nft(2)*mill0(3,ii))
    end do
    !
    deallocate(mill0)
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
