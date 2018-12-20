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
module plot
  !
  implicit none
  !
contains
  !>
  !! Output band dispersion
  !!
  subroutine band_plot()
    !
    use constant, only : htr2ev
    use k_point, only : kvec, nk
    use kohn_sham, only : nbnd, eval, ef
    !
    integer :: fo = 20, ibnd, ik
    real(8) :: xk(nk), dkvec(3), dxk
    !
    write(*,*) "  Output band.dat"
    open(fo, file="band.dat")
    !
    xk(1) = 0.0d0
    do ik = 2, nk
       dkvec(1:3) = kvec(1:3,ik) - kvec(1:3,ik-1)
       dxk = sqrt(dot_product(dkvec(1:3), dkvec(1:3)))
       xk(ik) = xk(ik-1) + dxk
    end do
    !
    write(fo,*) "# k energy[eV]" 
    do ibnd = 1, nbnd
       do ik = 1, nk
          write(fo,*) xk(ik), (eval(ibnd, ik)-ef)*htr2ev
       end do
       write(fo,*) ""
    end do
    !
    close(fo)
    !
  end subroutine band_plot
  !>
  !! Output FermiSurfer file
  !!
  subroutine fermi_plot()
    !
    use atm_spec, only : bvec
    use k_point, only : kgrd, nk, kvec
    use kohn_sham, only : eval, nbnd, evec, ef
    use gvec, only : g_wf
    use constant, only : htr2ev
    !
    integer :: ik, ibnd, fo = 20, ii, ipw
    real(8) :: kgv(g_wf%npw,3), avf(nbnd,nk), vf(3)
    !
    do ik = 1, nk
       do ipw = 1, g_wf%npw
          kgv(ipw,1:3) = kvec(1:3,ik) + matmul(bvec(1:3,1:3), dble(g_wf%mill(1:3,g_wf%map(ipw))))
       end do
       do ibnd = 1, nbnd
          do ii = 1, 3
             vf(ii) = dble(sum(conjg(evec(1:g_wf%npw,ibnd,ik)) &
             &                     * evec(1:g_wf%npw,ibnd,ik)  &
             &                      * kgv(1:g_wf%npw,ii     )  ))
          end do
          avf(ibnd,ik) = sqrt(dot_product(vf(1:3),vf(1:3)))
       end do
    end do
    !
    write(*,*) "  Output vf.frmsf"
    open(fo, file = "vf.frmsf")
    write(fo,*) kgrd(3), kgrd(2), kgrd(1)
    write(fo,*) 1
    write(fo,*) nbnd
    write(fo,*) real(bvec(1:3,3))
    write(fo,*) real(bvec(1:3,2))
    write(fo,*) real(bvec(1:3,1))
    do ibnd = 1, nbnd
       do ik = 1, nk
          write(fo,*) real((eval(ibnd,ik) - ef)*htr2ev)
       end do
    end do
    do ibnd = 1, nbnd
       do ik = 1, nk
          write(fo,*) real(avf(ibnd,ik))
       end do
    end do
    close(fo)    
    !
  end subroutine fermi_plot
  !
end module plot
