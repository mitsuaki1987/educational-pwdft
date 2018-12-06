module scf
  !
  implicit none
  !
contains
  !
  subroutine scf_loop()
    !
    use gvec, only : g_rh
    use solver, only : electron_maxstep, mixing_beta, conv_thr
    use io_vloc, only : Vks
    use constant, only : htr2ev
    !
    integer :: istep, jstep
    real(8) :: res, rhs(g_rh%nr), drhs(g_rh%nr), rhs0(g_rh%nr), dVks(g_rh%nr), &
    &          jacob1(g_rh%nr,electron_maxstep), &
    &          jacob2(g_rh%nr,electron_maxstep), &
    &          alpha
    !
    call init_rho_V()
    !
    istep = 0
    write(*,*) "Iteration ", istep
    !
    call kohn_sham(.true., rhs)
    res = sqrt(dot_product(rhs(1:g_rh%nr), rhs(1:g_rh%nr))) / dble(g_rh%nr)
    write(*,*) "  delta Vks [eV] : ", res * htr2eV
    !
    if(res < conv_thr) GOTO 5
    !
    dVks(1:g_rh%nr) = - mixing_beta * rhs(1:g_rh%nr)
    !
    do istep = 1, electron_maxstep
       !
       write(*,*) "Iteration ", istep
       !
       Vks(1:g_rh%nr) = Vks(1:g_rh%nr) + dVks(1:g_rh%nr)
       !
       rhs0(1:g_rh%nr) = rhs(1:g_rh%nr)
       call kohn_sham(.false., rhs)
       res = sqrt(dot_product(rhs(1:g_rh%nr), rhs(1:g_rh%nr))) / dble(g_rh%nr)
       write(*,*) "  delta Vks [eV] : ", res * htr2eV
       !
       if(res < conv_thr) then
          !       
          GOTO 5
          !
       end if
       !
       ! Update Jacobian with drhs
       !
       drhs(1:g_rh%nr) = rhs(1:g_rh%nr) - rhs0(1:g_rh%nr)
       !
       jacob1(1:g_rh%nr,istep) = - mixing_beta * drhs(1:g_rh%nr)
       do jstep = 1, istep - 1
          alpha = dot_product(jacob2(1:g_rh%nr,jstep), drhs(1:g_rh%nr))
          jacob1(1:g_rh%nr,istep) = jacob1(1:g_rh%nr,istep) - jacob1(1:g_rh%nr,jstep) * alpha
       end do
       jacob1(1:g_rh%nr,istep) = dVks(1:g_rh%nr) + jacob1(1:g_rh%nr,istep)
       alpha = dot_product(drhs(1:g_rh%nr), drhs(1:g_rh%nr))
       jacob2(1:g_rh%nr,istep) = drhs(1:g_rh%nr) / alpha
       !
       ! Compute dVks with new Jacobian & rhs
       !
       dVks(1:g_rh%nr) = - mixing_beta * rhs(1:g_rh%nr)
       do jstep = 1, istep
          alpha = dot_product(jacob2(1:g_rh%nr,jstep), rhs(1:g_rh%nr))
          dVks(1:g_rh%nr) = dVks(1:g_rh%nr) - jacob1(1:g_rh%nr,jstep) * alpha 
       end do
       !
    end do ! istep
    !
    write(*,'(/,7x,"Not converged ! res = ",e12.5,/)') res
    RETURN
    !
5   CONTINUE
    !
    write(*,'(/,7x,"Converged ! iter = ",i0,/)') istep
    !
  end subroutine scf_loop
  !
  !
  !
  subroutine kohn_sham(linit,rhs)
    !
    use gvec, only : g_wf, g_rh
    use solver, only : nbnd, eval, evec, calculation
    use k_point, only : kvec, nk, ksum_rho
    use hartree, only : hartree_pot, xc_pot
    use io_vloc, only : Vks, Vps
    use diag_direct, only : direct
    use lobpcg, only : lobpcg_main
    !
    logical,intent(in) :: linit
    real(8),intent(out) :: rhs(g_rh%nr)
    !
    integer :: ik
    !
    if(calculation == "direct") then
       call direct(g_wf%npw,kvec(1:3,1),evec(1:g_wf%npw,1:nbnd,1),eval(1:nbnd,1))
    else
       do ik = 1, nk
          call lobpcg_main(linit, g_wf%npw, kvec(1:3,ik),&
          &                evec(1:g_wf%npw,1:nbnd,ik),eval(1:nbnd,ik))
       end do
    end if
    !
    if(calculation == "scf") then
       !
       call ksum_rho()
       !
       rhs(1:g_rh%nr) = Vps(1:g_rh%nr)
       call hartree_pot(rhs)
       call xc_pot(rhs)
       !
       rhs(1:g_rh%nr) = Vks(1:g_rh%nr) - rhs(1:g_rh%nr)
       !
    else
       rhs(1:g_rh%nr) = 0.0d0
    end if
    !
  end subroutine kohn_sham
  !
  !
  !
  subroutine init_rho_V()
    !
    use gvec, only : g_rh
    use io_vloc, only : rho, Vks, Vps
    use hartree, only : hartree_pot, xc_pot
    use atm_spec, only : Vcell, nelec
    !
    rho(1:g_rh%nr) = nelec / Vcell
    !
    Vks(1:g_rh%nr) = Vps(1:g_rh%nr)
    call hartree_pot(Vks)
    call xc_pot(Vks)
    !
  end subroutine init_rho_V
  !
end module scf
