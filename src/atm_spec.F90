module atm_spec
  !
  implicit none
  !
  type atm_t
     integer :: ityp
     character(256) :: elem
     real(8) :: pos(3)
  end type atm_t
  !
  type spec_t
     character(256) :: elem
     character(256) :: ps_file
  end type spec_t
  !
  integer :: nat, ntyp
  real(8) :: avec(3,3)
  type(atm_t),allocatable :: atm(:)
  type(spec_t),allocatable :: spec(:)
  !
end module atm_spec
