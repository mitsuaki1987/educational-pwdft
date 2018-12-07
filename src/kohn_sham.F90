module kohn_sham
  !
  implicit none
  !
  character(256),save :: &
  & calculation !< Calculation mode
  integer,save :: &
  & nbnd !< Number of bands
  real(8),save :: &
  & ef = 0.0d0 !< Fermi energy [Htr]
  real(8),allocatable,save :: &
  & eval(:,:) !< (nbnd,nk) Kohn-Sham eigenvalue (energy)  
  complex(8),allocatable,save :: &
  & evec(:,:,:) !< (g_wf%npw,nbnd,nk) Kohn-Sham eigenvector (orbital)
  !
end module kohn_sham
