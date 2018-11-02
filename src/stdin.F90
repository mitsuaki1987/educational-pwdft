module stdin
  !
  implicit none
  !
contains
  !
  subroutine read_stdin()
    !
    use control, only : calculation
    use system, only : nat, ntype, ecutwfc
    use electrons, only : eectron_maxstep
    use k_point, only : nk, kvec, kgrd
    use atm_spec, only : atm, spec
    !
    integer :: iat, jtpe, ii, ik
    namelist/control/ calculation
    namelist/system/ nat, ntype, ecutwfc
    !
    read(*,control)
    !
    write(*,*) "  calculation : ", trim(calculation)
    !
    read(*,system)
    !
    write(*,*) "  Number o atoms : ", nat
    write(*.*) "  Number of species : ", ntype
    write(*,*) "  Plane-wave cutoff : ", ecutwc
    !
    allocate(atm(nat), spec(ntype))
    !
    read(*,*) key
    if(key == "ATOMIC_SPECIES") THEN
       write(*,*) " Atomic species :"
       do jtype = 1, ntype
          read(*,*) spec(jtype)%elem, spec(jtype)%ps_file
          write(*,*) jtype, spec(jtype)%elem, spec(jtype)%ps_file
       end do
    else if(key == "ATOMIC_POSITIONS") THEN
       do iat = 1, nat
          read(*,*) atm(iat)%elem, atm(iat)%pos(1:3)
       end do
    else if(key == "K_POINTS")
       if(calculation == "scf" .or. calculation == "nscf") then
          read(*,*) kgrd(1:3)
       else
          read(*,*) nk
          allocate(kvec(3,nk))
          do ik = 1, nk
             read(*,*) kvec(1:3,ik)
          end do
       end if
    else if(key == "CELL_PARAMETERS")
       read(*,*) avec(1:3,1:3)
       write(*,*) "  Cell parameters:"
       do ii = 1, 3
          write(*,*) avec(1:3,ii)
       end do
    end if
    !
    write(*,*) "  Atomic position & specy"
    do iat = 1, nat
       do jtype = 1, ntype
          if(trim(spec(itype)%elem) == trim(atm(iat)%elem)) &
          & atm(iat)%itype = jtype
          exit
       end do
       write(*,*) iat, trim(atm(iat)%elem)), atm(iat)%itype, atm(iat)%pos(1:3)
    end do
    !
  end subroutine read_stdin
  !
end module stdin
