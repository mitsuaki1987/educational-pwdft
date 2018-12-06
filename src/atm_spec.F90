module atm_spec
  !
  implicit none
  !
  type atm_t
     integer :: ityp !< Element type index starts from 1
     character(256) :: elem !< Element name
     real(8) :: pos(3) !< Position in fractional coord.
  end type atm_t
  !
  type spec_t
     character(256) :: &
     & elem, & !< Element name
     & ps_file !< Pseudopotential file name
     real(8) :: &
     & zion, & !< Ion charge
     & zatom !< Nuclear charge
     integer :: &
     & mmax !< Number of radial grid for PseudoPot.
     real(8),allocatable :: &
     & psr(:), & !< Radial grid [Bohr]
     & psV(:) !< V_ps [Htr]
  end type spec_t
  !
  integer,save :: &
  & nat, &  !< Number of atoms
  & ntyp !< Number of species (elements)
  real(8),save :: &
  & Vcell, & !< Unit cell volume [Bohr^3]
  & nelec, & !< Number of electrons per u.c.
  & avec(3,3), & !< Unit lattice vector [Bohr]
  & bvec(3,3) !< Unit reciplocal lattice vector [Bohr^-1]
  type(atm_t),allocatable,save :: &
  & atm(:) !< (nat) Atom
  type(spec_t),allocatable,save :: &
  & spec(:) !< (ntyp) Species
  !
end module atm_spec
