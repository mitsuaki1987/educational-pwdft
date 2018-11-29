module stdin
  !
  implicit none
  !
contains
  !
  subroutine read_stdin()
    !
    use solver, only : calculation, nbnd, electron_maxstep, conv_thr
    use k_point, only : nk, kvec, kgrd
    use atm_spec, only : atm, spec, avec, nat, ntyp, bvec
    use gvec, only : rfft, wfft, setup_gvec
    !
    integer :: iat, jtyp, ii, ik, jj, lwork = 3, ipiv(3), info
    real(8) :: work(3), pi = acos(-1.0d0), ecutrho, ecutwfc
    character(256) :: key
    !
    namelist/control/ calculation
    namelist/system/ nat, ntyp, ecutwfc, ecutrho, nbnd
    namelist/electrons/ electron_maxstep, conv_thr
    !
    read(*,control)
    !
    write(*,*) "                     calculation : ", trim(calculation)
    !
    read(*,system)
    !
    write(*,*) "                 Number of atoms : ", nat
    write(*,*) "               Number of species : ", ntyp
    write(*,*) "    Plane-wave cutoff (wfc) [Ry] : ", ecutwfc
    write(*,*) "  Plane-wave cutoff (rho,V) [Ry] : ", ecutrho
    write(*,*) "                 Number of bands : ", nbnd
    !
    allocate(atm(nat), spec(ntyp))
    !
    electron_maxstep = 100
    conv_thr = 1.0d-10
    !
    read(*,electrons)
    !
    write(*,*) "                   Max iteration : ", electron_maxstep
    write(*,*) "           Convergence threshold : ", conv_thr
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
          avec(1:3,1:3) = avec(1:3,1:3) / 0.529177249d0
          write(*,*) "  Cell parameters [Bohr] :"
          !
          do ii = 1, 3
             write(*,*) avec(1:3,ii)
          end do
          !
          ! Reciprocal lattice vectors
          !
          bvec(1:3,1:3) = transpose(avec(1:3,1:3))
          call dgetrf(3, 3, bvec, 3, ipiv, info) 
          call dgetri(3, bvec, 3, ipiv, work, lwork, info)
          if(info /= 0) stop 'read_stdin : inverce of avec.'
          bvec(1:3,1:3) = 2.0d0 * pi * bvec(1:3,1:3)
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
    write(*,*) "  FFT and G-vector for wfc"
    call setup_gvec(wfft, ecutwfc)
    write(*,*) "  FFT and G-vector for rho and V"
    call setup_gvec(rfft, ecutrho)
    !
  end subroutine read_stdin
  !
end module stdin
