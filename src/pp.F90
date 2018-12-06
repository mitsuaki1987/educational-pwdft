module pp
  !
  implicit none
  !
contains
  !
  subroutine read_pp()
    !
    use atm_spec, only : spec, ntyp, nelec, atm, nat
    use solver, only : nbnd
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
