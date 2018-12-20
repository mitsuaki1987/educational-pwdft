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
    call libtetrabz_fermieng(2,bvec,nbnd,kgrd,eval,kgrd,occ,ef,nelec*0.5d0)
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
  !>
  !! DOS
  !!
  subroutine ksum_dos()
    !
    use constant, only : htr2ev
    use atm_spec, only : bvec
    use kohn_sham, only : nbnd, ef, eval
    use libtetrabz, only : libtetrabz_dos, libtetrabz_intdos
    !
    integer,parameter :: ne = 500
    integer :: ie, fo = 20
    real(8) :: e0(ne), dos(ne), intdos(ne), wght(ne,nbnd,nk), &
    &          emax, emin, de
    !
    emax = maxval(eval(1:nbnd,1:nk))
    emin = minval(eval(1:nbnd,1:nk))
    de = (emax - emin) / dble(ne)
    do ie = 1, ne
       e0(ie) = emin + de * (ie -1)
    end do
    !
    call libtetrabz_dos(2,bvec,nbnd,kgrd,eval,kgrd,wght,ne,e0)
    !
    do ie = 1, ne
       dos(ie) = sum(wght(ie,1:nbnd,1:nk))
    end do
    !
    call libtetrabz_intdos(2,bvec,nbnd,kgrd,eval,kgrd,wght,ne,e0)
    !
    do ie = 1, ne
       intdos(ie) = sum(wght(ie,1:nbnd,1:nk))
    end do
    !
    write(*,*) "  Output dos.dat"
    open(fo, file = "dos.dat")
    !
    write(fo,*) "# Energy[eV]  DOS[eV^-1]  IntDOS"
    !
    do ie = 1, ne
       write(fo,*) (e0(ie)-ef)*htr2eV, dos(ie)*2.0d0/htr2ev, intdos(ie)*2.0d0
    end do
    !
    close(fo)
    !
  end subroutine ksum_dos
  !
end module k_point
