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
module griddata
  !
  implicit none
  !
contains
  !>
  !! Read unit-cell grid data such as Vks
  !!
  subroutine read_griddata(filename, gdata)
    !
    use constant, only : pi
    use atm_spec, only : nat, atm, bvec, avec
    use gvec, only : g_rh
    !
    character(*),intent(in) :: filename
    real(8),intent(out) :: gdata(g_rh%nft(1),g_rh%nft(2),g_rh%nft(3))
    !
    integer :: itemp(3), fi = 10, iat
    character(256) :: ctemp
    real(8) :: avec0(3,3), &
    &          gdata0(g_rh%nft(1)+1,g_rh%nft(2)+1,g_rh%nft(3)+1)
    !
    open(fi, file = trim(filename))
    !
    read(fi,*) ctemp
    read(fi,*) ctemp
    read(fi,*) avec0(1:3,1:3)
    avec0(1:3,1:3) = avec0(1:3,1:3) / 0.529177249d0
    if(any(abs(avec0(1:3,1:3)-avec(1:3,1:3)) > 1.0d-3)) then
       write(*,*) "Error in read_griddata"
       write(*,*) "Direct lattice vector is different."
       stop 'error in read_griddata'
    end if
    read(fi,*) ctemp
    read(fi,*) itemp(1:2)
    if(nat /= itemp(1)) then
       write(*,*) "Error in read_griddata"
       write(*,*) "Number of atoms is different."
       stop 'error in read_griddata'
    end if
    do iat = 1, nat
       read(fi,*) ctemp, avec0(1:3,1)
       avec0(1:3,1) = avec0(1:3,1) / 0.529177249d0
       avec0(1:3,1) = matmul(avec0(1:3,1), bvec(1:3,1:3)) / (2.0d0*pi)
       if(any(abs(avec0(1:3,1)-atm(iat)%pos(1:3)) > 1.0d-3)) then
          write(*,*) "Error in read_griddata"
          write(*,*) "Position of atom ", iat, " is different."
          stop 'error in read_griddata'
       end if
    end do
    read(fi,*) ctemp
    read(fi,*) ctemp
    read(fi,*) ctemp 
    read(fi,*) itemp(1:3)
    if(any(itemp(1:3) /= g_rh%nft(1:3)+1)) then
       write(*,*) "Error in read_griddata"
       write(*,*) "FFT grid is different."
       stop 'error in read_griddata'
    end if
    read(fi,*) avec0(1:3,1)
    read(fi,*) avec0(1:3,1:3)
    read(fi,*) gdata0(1:itemp(1),1:itemp(2),1:itemp(3))
    !
    close(fi)
    !
    gdata(     1:g_rh%nft(1),1:g_rh%nft(2),1:g_rh%nft(3)) &
    & = gdata0(1:g_rh%nft(1),1:g_rh%nft(2),1:g_rh%nft(3))
    !
  end subroutine read_griddata
  !>
  !! Write unit-cell grid data such as Vks
  !!
  subroutine write_griddata(filename, gdata)
    !
    use constant, only : bohr2ang
    use atm_spec, only : nat, atm, avec
    use gvec, only : g_rh
    !
    character(*),intent(in) :: filename
    real(8),intent(in) :: gdata(g_rh%nft(1),g_rh%nft(2),g_rh%nft(3))
    !    
    integer :: fo = 20, iat, i1, i2, i3, ii1, ii2, ii3
    real(8) :: gdata0(g_rh%nft(1)+1,g_rh%nft(2)+1,g_rh%nft(3)+1)
    !
    do i3 = 1, g_rh%nft(3)+1
       ii3 = modulo(i3 - 1, g_rh%nft(3)) + 1
       do i2 = 1, g_rh%nft(2)+1
          ii2 = modulo(i2 - 1, g_rh%nft(2)) + 1
          do i1 = 1, g_rh%nft(1)+1
             ii1 = modulo(i1 - 1, g_rh%nft(1)) + 1
             gdata0(i1, i2, i3) = gdata(ii1, ii2, ii3)
          end do
       end do
    end do
    !
    write(*,*) "  Output ", trim(filename)
    open(fo, file = trim(filename))
    !
    write(fo,'(a)') "CRYSTAL"
    write(fo,'(a)') "PRIMVEC"
    write(fo,'(3e20.10)') avec(1:3,1:3) * bohr2ang
    write(fo,'(a)') "PRIMCOORD"
    write(fo,'(2(2x,i0))') nat, 1
    do iat = 1, nat
       write(fo,'(a,2x,3e20.10)') trim(atm(iat)%elem), &
       &  matmul(avec(1:3,1:3), atm(iat)%pos(1:3)) * bohr2ang
    end do
    write(fo,'(a)') "BEGIN_BLOCK_DATAGRID_3D"
    write(fo,'(a)') "3D_PWSCF"
    write(fo,'(a)') "DATAGRID_3D_UNKNOWN"
    write(fo,'(3(2x,i0))') g_rh%nft(1:3)+1
    write(fo,'(3e20.10)') (/0.0d0, 0.0d0, 0.0d0/), avec(1:3,1:3) * bohr2ang
    write(fo,'(6e20.10)') gdata0(1:g_rh%nft(1)+1,1:g_rh%nft(2)+1,1:g_rh%nft(3)+1)
    write(fo,'(a)') "END_DATAGRID_3D"
    write(fo,'(a)') "END_BLOCK_DATAGRID_3D"
    !
    close(fo)
    !
  end subroutine write_griddata
  !
end module griddata
