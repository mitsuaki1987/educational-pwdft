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
module pp
  !
  implicit none
  !
contains
  !
  subroutine read_pp()
    !
    use atm_spec, only : spec, ntyp, nelec, atm, nat
    use kohn_sham, only : nbnd
    !
    integer :: ityp, fi=10, iat
    integer :: pspd, pspcode, pspxc, lmax, lloc, r2well, &
    &          rchrg, fchrg, qchrg, nproj(5), extension_switch, &
    &          ir, jr
    character(256) :: ctmp
    !
    do ityp = 1, ntyp
       !
       write(*,*) "  Reading ", trim(spec(ityp)%ps_file), &
       &          " for pseudopotential of ", trim(spec(ityp)%elem)
       !
       open(fi, file = trim(spec(ityp)%ps_file))
       !
       read(fi,*) ctmp
       !
       read(fi,*) spec(ityp)%zatom, spec(ityp)%zion, pspd
       write(*,*) "    Zatom : ", spec(ityp)%zatom
       write(*,*) "    Zion : ", spec(ityp)%Zion
       write(*,*) "    Gen. Date : ", pspd
       !
       read(fi,*) pspcode, pspxc, lmax, lloc, spec(ityp)%mmax, r2well
       write(*,*) "    PSP code : ", pspcode
       write(*,*) "    XC id : ", pspxc
       write(*,*) "    Max. L : ", lmax
       write(*,*) "    Localized L : ", lloc
       write(*,*) "    Radial grid : ", spec(ityp)%mmax
       !
       read(fi,*) rchrg, fchrg, qchrg
       write(*,*) "    R-charge : ", rchrg
       write(*,*) "    F-charge : ", fchrg
       write(*,*) "    Q-charge : ", qchrg
       !
       read(fi,*) nproj(1:5)
       read(fi,*) extension_switch
       read(fi,*) ir
       !
       allocate(spec(ityp)%psr(spec(ityp)%mmax), &
       &        spec(ityp)%psV(spec(ityp)%mmax)  )
       !
       do ir = 1, spec(ityp)%mmax
          read(fi,*) jr, spec(ityp)%psr(ir), &
          &              spec(ityp)%psV(ir)
       end do
       !
       close(fi)
       !
    end do
    !
    nelec = 0.0d0
    do iat = 1, nat
       nelec = nelec + spec(atm(iat)%ityp)%Zion
    end do
    write(*,*) "  Number of electrons : ", nelec
    if(nbnd == 0) nbnd = nint(nelec)
    write(*,*) "  Number of bands : ", nbnd
    !
  end subroutine read_pp
  !
end module pp
