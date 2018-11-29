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
  integer,save :: nat, ntyp
  real(8),save :: avec(3,3), bvec(3,3)
  type(atm_t),allocatable,save :: atm(:)
  type(spec_t),allocatable,save :: spec(:)
  !
end module atm_spec
