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
!>
!! @mainpage
!!
!! Plane-wave DFT code for education
!!
!! pwdft.F90 : Main program
!!
program pwdft
  !
  use kohn_sham, only : calculation, nbnd, eval, evec, ef
  use stdin, only : read_stdin
  use constant, only : htr2ev
  use gvec, only : g_rh, g_wf
  use rho_v, only : Vks, init_rho_v
  use k_point, only : nk, ksum_dos
  use pp, only : read_pp
  use scf, only : scf_loop, kohn_sham_eq
  use fftw_wrapper, only : init_fft
  use griddata, only : read_griddata, write_griddata
  use plot, only : fermi_plot, band_plot
  use energy, only : total_e
  !
  implicit none
  !
  real(8) :: t1, t2, t3, t4
  real(8),allocatable :: dummy(:)
  !
  call cpu_time(t1)
  !
  call read_stdin()
  !
  call read_pp()
  !
  call init_fft()
  call init_rho_v()
  !
  allocate(eval(nbnd,nk), evec(g_wf%npw,nbnd,nk))
  !
  if(calculation == "scf") then
     !
     call cpu_time(t3)
     call scf_loop()
     call cpu_time(t4)
     write(*,*) "  SCF time : ", t4 - t3, " sec"
     !
     call total_e()
     !
     Vks(1:g_rh%nr) = (Vks(1:g_rh%nr) - ef)*htr2ev
     call write_griddata("vks.xsf", Vks)
     Vks(1:g_rh%nr) = Vks(1:g_rh%nr)/htr2ev + ef
     !
  else
     !
     if(calculation == "bands" .or. calculation == "nscf") then
        call read_griddata("vks.xsf", Vks)
        write(*,*) "  Average potential [eV] : ", &
        &  sum(Vks(1:g_rh%nr)) / dble(g_rh%nr)
        Vks(1:g_rh%nr) = Vks(1:g_rh%nr) / htr2ev
     end if
     !
     allocate(dummy(g_rh%nr))
     call cpu_time(t3)
     call kohn_sham_eq(.true., dummy)
     call cpu_time(t4)
     write(*,*) "  Kohn-Sham time : ", t4 - t3, " sec"
     !
  end if
  !
  if(calculation == "scf" .or. calculation=="nscf") then
     call ksum_dos()
     call fermi_plot()
  else
     call band_plot()
  end if
  !
  call cpu_time(t2)
  !
  write(*,*) "  Total time : ", t2 - t1, " sec"
  write(*,*) "  End"
  !
end program pwdft
