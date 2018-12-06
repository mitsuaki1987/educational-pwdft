program pwdft
  !
  use solver, only : calculation, nbnd, eval, evec
  use stdin, only : read_stdin
  use gvec, only : g_rh, g_wf
  use io_vloc, only : Vks, Vps, rho, read_vks
  use diag_direct, only : direct
  use k_point, only : nk
  use pp, only : read_pp
  use scf, only : scf_loop, kohn_sham
  !
  implicit none
  !
  real(8) :: t1, t2
  real(8),allocatable :: rhs(:)
  !
  call cpu_time(t1)
  !
  call read_stdin()
  !
  call read_pp()
  !
  allocate(eval(nbnd,nk), evec(g_wf%npw,nbnd,nk))
  allocate(Vks(g_rh%nr), Vps(g_rh%nr), rho(g_rh%nr), rhs((g_rh%nr)))
  !
  if(calculation == "scf") then
     call scf_loop()
  else
     call read_vks()
     call kohn_sham(.true., rhs)
  end if
  !
  call cpu_time(t2)
  !
  write(*,*) "  End : ", t2 - t1, " sec"
  !
end program pwdft
