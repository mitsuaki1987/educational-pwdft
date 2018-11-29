program pwdft
  !
  use solver, only : calculation
  use stdin, only : read_stdin
  use gvec, only : rfft
  use io_vloc, only : vloc, read_vloc
  !
  implicit none
  !
  call read_stdin()
  !
  allocate(Vloc(rfft%npw3(1),rfft%npw3(2),rfft%npw3(3)))
  if(calculation /= "scf") call read_vloc()
  !
end program pwdft
