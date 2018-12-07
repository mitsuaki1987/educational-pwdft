module stdin
  !
  implicit none
  !
contains
  !
  subroutine read_stdin()
    !
    use constant, only : bohr2ang, pi, htr2ev
    use kohn_sham, only : calculation, nbnd
    use scf, only : electron_maxstep, conv_thr, mixing_beta
    use k_point, only : nk, kvec, kgrd
    use atm_spec, only : atm, spec, avec, nat, ntyp, bvec, Vcell
    use gvec, only : g_rh, g_wf, setup_gvec
    !
    integer :: iat, jtyp, ii, jj, lwork = 3, ipiv(3), &
    &          info, i1, i2, i3, nn, ik
    real(8) :: work(3), ecutrho, ecutwfc
    character(256) :: key
    !
    namelist/control/ calculation
    namelist/system/ nat, ntyp, ecutwfc, ecutrho, nbnd
    namelist/electrons/ electron_maxstep, conv_thr, mixing_beta
    !
    read(*,control)
    !
    write(*,*) "                     calculation : ", trim(calculation)
    !
    nbnd = 0
    !
    read(*,system)
    !
    write(*,*) "                 Number of atoms : ", nat
    write(*,*) "               Number of species : ", ntyp
    write(*,*) "    Plane-wave cutoff (wfc) [Ry] : ", ecutwfc
    write(*,*) "  Plane-wave cutoff (rho,V) [Ry] : ", ecutrho
    !
    allocate(atm(nat), spec(ntyp))
    !
    electron_maxstep = 100
    conv_thr = 1.0d-10
    mixing_beta = 0.3d0
    !
    read(*,electrons)
    !
    write(*,*) "                   Max iteration : ", electron_maxstep
    write(*,*) "           Convergence threshold : ", conv_thr
    write(*,*) "                  Initial mixing : ", mixing_beta
    !
    conv_thr = conv_thr / htr2ev
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
          if(calculation == "bands") then
             read(*,*) nn
             nk = 10*(nn-1)+1
             allocate(kvec(3,nk))
             read(*,*) kvec(1:3,1)
             do ii = 1, nn - 1
                read(*,*) kvec(1:3,10*ii+1)
                do ik = 2, 10
                   kvec(1:3,10*(ii-1)+ik) &
                   & = kvec(1:3,10*(ii-1)+1) * dble(11-ik)/10.0d0 &
                   & + kvec(1:3,10* ii   +1) * dble(ik-1)/10.0d0
                end do
             end do
          else if (calculation == "direct" .or. calculation == "iterative") then
             read(*,*) nk
             allocate(kvec(3,nk))
             read(*,*) kvec(1:3,1:nk)
          else
             read(*,*) kgrd(1:3)
             nk = product(kgrd(1:3))
             allocate(kvec(3,nk))
             nk = 0
             do i3 = 0, kgrd(3) - 1
                do i2 = 0, kgrd(2) - 1
                   do i1 = 0, kgrd(1) - 1
                      nk = nk + 1
                      kvec(1:3,nk) = dble(modulo((/i1,i2,i3/)+kgrd(1:3)/2,kgrd(1:3)) &
                      &                                      -kgrd(1:3)/2) &
                      &            / dble(kgrd(1:3))
                   end do
                end do
             end do
          end if
          !
          write(*,*) "  Number of k : ", nk
          !
       else if(key == "CELL_PARAMETERS") then
          !
          read(*,*) avec(1:3,1:3)
          avec(1:3,1:3) = avec(1:3,1:3) / bohr2ang
          write(*,*) "  Cell parameters [Bohr] :"
          !
          do ii = 1, 3
             write(*,*) "    ", avec(1:3,ii)
          end do
          Vcell = (avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)) * avec(1,3) &
          &     + (avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)) * avec(2,3) &
          &     + (avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)) * avec(3,3)
          Vcell = abs(Vcell)
          write(*,*) "  Cell volume [Bohr^3] :", Vcell
          !
          ! Reciprocal lattice vectors
          !
          bvec(1:3,1:3) = transpose(avec(1:3,1:3))
          call dgetrf(3, 3, bvec, 3, ipiv, info) 
          call dgetri(3, bvec, 3, ipiv, work, lwork, info)
          if(info /= 0) stop 'read_stdin : inverce of avec.'
          bvec(1:3,1:3) = 2.0d0 * pi * bvec(1:3,1:3)
          write(*,*) "  Reciplocal lattice vector [Bohr^-1] :"
          do ii = 1, 3
             write(*,*) "    ", bvec(1:3,ii)
          end do
          !
       end if
       !
    end do
    kvec(1:3,1:nk) = matmul(bvec(1:3,1:3), kvec(1:3,1:nk))
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
    call setup_gvec(g_wf, ecutwfc)
    write(*,*) "  FFT and G-vector for rho and V"
    call setup_gvec(g_rh, ecutrho)
    !
  end subroutine read_stdin
  !
end module stdin
