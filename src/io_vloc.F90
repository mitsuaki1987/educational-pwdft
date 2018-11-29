module io_vloc
  !
  implicit none
  !
  complex(8),allocatable :: Vloc(:,:,:)
  !
contains
  !
  subroutine read_vloc()
    !
    use atm_spec, only : nat, atm, bvec, avec
    use gvec, only : rfft
    !
    integer :: itemp(3), fi = 10, iat
    character(256) :: ctemp
    real(8) :: avec0(3,3), pi = acos(-1.0d0), &
    &          Vloc0(1:rfft%npw3(1)+1,1:rfft%npw3(2)+1,1:rfft%npw3(3)+1)
    !
    open(fi, file = "vloc.xsf")
    !
    read(fi,*) ctemp
    read(fi,*) ctemp
    read(fi,*) avec0(1:3,1:3)
    avec0(1:3,1:3) = avec0(1:3,1:3) / 0.529177249d0
    if(any(abs(avec0(1:3,1:3)-avec(1:3,1:3)) > 1.0d-3)) then
       write(*,*) "Error in read_vloc"
       write(*,*) "Direct lattice vector is different."
       stop 'error in read_vloc'
    end if
    read(fi,*) ctemp
    read(fi,*) itemp(1:2)
    if(nat /= itemp(1)) then
       write(*,*) "Error in read_vloc"
       write(*,*) "Number of atoms is different."
       stop 'error in read_vloc'
    end if
    do iat = 1, nat
       read(fi,*) ctemp, avec0(1:3,1)
       avec0(1:3,1) = avec0(1:3,1) / 0.529177249d0
       avec0(1:3,1) = matmul(avec0(1:3,1), bvec(1:3,1:3)) / (2.0d0*pi)
       if(any(abs(avec0(1:3,1)-atm(iat)%pos(1:3)) > 1.0d-3)) then
          write(*,*) "Error in read_vloc"
          write(*,*) "Position of atom ", iat, " is different."
          stop 'error in read_vloc'
       end if
    end do
    read(fi,*) ctemp
    read(fi,*) ctemp
    read(fi,*) ctemp 
    read(fi,*) itemp(1:3)
    if(any(itemp(1:3) /= rfft%npw3(1:3)+1)) then
       write(*,*) "Error in read_vloc"
       write(*,*) "FFT grid is different."
       stop 'error in read_vloc'
    end if
    read(fi,*) avec0(1:3,1)
    read(fi,*) avec0(1:3,1:3)
    read(fi,*) Vloc0(1:itemp(1),1:itemp(2),1:itemp(3))
    Vloc(     1:rfft%npw3(1),1:rfft%npw3(2),1:rfft%npw3(3)) &
    & = Vloc0(1:rfft%npw3(1),1:rfft%npw3(2),1:rfft%npw3(3))
    !
  end subroutine read_vloc
  !
end module io_vloc
