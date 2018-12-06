module solver
  !
  implicit none
  !
  character(256),save :: &
  & calculation !< Calculation mode
  integer,save :: &
  & nbnd, & !< Number of bands
  & electron_maxstep !< Max number of iteration
  real(8),save :: &
  & mixing_beta, & !< Mixing for SCF
  & ef, & !< Fermi energy [Htr]
  & conv_thr !< Convergence threshold [Htr]
  real(8),allocatable,save :: &
  & eval(:,:) !< (nbnd,nk) Kohn-Sham eigenvalue (energy)  
  complex(8),allocatable,save :: &
  & evec(:,:,:) !< (g_wf%npw,nbnd,nk) Kohn-Sham eigenvector (orbital)
  !
end module solver
