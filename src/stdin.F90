module stdin
  !
  implicit none
  !
contains
  !
  subroutine read_stdin()
    !
    use solver, only : calculation, ecutwfc, nbnd
    use k_point, only : nk, kvec, kgrd
    use atm_spec, only : atm, spec, avec, nat, ntyp
    !
    integer :: iat, jtyp, ii, ik, jj
    character(256) :: key
    !
    namelist/control/ calculation
    namelist/system/ nat, ntyp, ecutwfc, nbnd
    !
    read(*,control)
    !
    write(*,*) "  calculation : ", trim(calculation)
    !
    read(*,system)
    !
    write(*,*) "  Number of atoms : ", nat
    write(*,*) "  Number of species : ", ntyp
    write(*,*) "  Plane-wave cutoff [Ry] : ", ecutwfc
    write(*,*) "  Number of bands : ", nbnd
    !
    allocate(atm(nat), spec(ntyp))
    !
    do jj = 1, 4
       !
       read(*,*) key
       !
       if(key == "ATOMIC_SPECIES") then
          !
          write(*,*) "  Atomic species :"
          !
          do jtyp = 1, ntyp
             read(*,*) spec(jtyp)%elem, spec(jtyp)%ps_file
             write(*,*) jtyp, trim(spec(jtyp)%elem), " ", trim(spec(jtyp)%ps_file)
          end do
          !
       else if(key == "ATOMIC_POSITIONS") then
          !
          do iat = 1, nat
             read(*,*) atm(iat)%elem, atm(iat)%pos(1:3)
          end do
          !
       else if(key == "K_POINTS") then
          !
          if(calculation == "scf" .or. calculation == "nscf") then
             read(*,*) kgrd(1:3)
          else
             read(*,*) nk
             allocate(kvec(3,nk))
             do ik = 1, nk
                read(*,*) kvec(1:3,ik)
             end do
          end if
          !
       else if(key == "CELL_PARAMETERS") then
          !
          read(*,*) avec(1:3,1:3)
          write(*,*) "  Cell parameters:"
          !
          do ii = 1, 3
             write(*,*) avec(1:3,ii)
          end do
          !
       end if
       !
    end do
    !
    write(*,*) "  Atomic position"
    !
    do iat = 1, nat
       !
       do jtyp = 1, ntyp
          if(trim(spec(jtyp)%elem) == trim(atm(iat)%elem)) then
             atm(iat)%ityp = jtyp
             exit
          end if
       end do
       !
       write(*,*) iat, trim(atm(iat)%elem), atm(iat)%ityp, atm(iat)%pos(1:3)
       !
    end do
    !
  end subroutine read_stdin
  !
end module stdin
