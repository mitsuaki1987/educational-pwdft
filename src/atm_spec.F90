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
