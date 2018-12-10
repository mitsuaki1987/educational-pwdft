module scf
  !
  implicit none
  !
  integer,save :: &
  & electron_maxstep !< Max number of iteration
  real(8),save :: &
  & mixing_beta, & !< Mixing for SCF
  & conv_thr !< Convergence threshold [Htr]
  !
contains
  !
  subroutine scf_loop()
    !
    use gvec, only : g_rh
    use rho_v, only : Vks
    use constant, only : htr2ev
    !
    integer :: istep, jstep
    real(8) :: alpha, res
    real(8),allocatable :: jacob1(:,:), jacob2(:,:), rhs0(:), dVks(:), rhs(:), drhs(:)
    !
    allocate(rhs0(g_rh%nr),dVks(g_rh%nr),rhs(g_rh%nr), drhs(g_rh%nr), &
    &        jacob1(g_rh%nr,electron_maxstep), jacob2(g_rh%nr,electron_maxstep))
    !
    istep = 0
    write(*,*) "  Iteration ", istep
    !
    call kohn_sham_eq(.true., rhs)
    res = sqrt(dot_product(rhs(1:g_rh%nr), rhs(1:g_rh%nr))) / dble(g_rh%nr)
    write(*,*) "    delta Vks [eV] : ", res * htr2eV
    !
    if(res < conv_thr) electron_maxstep = 0
    !
    dVks(1:g_rh%nr) = - mixing_beta * rhs(1:g_rh%nr)
    !
    do istep = 1, electron_maxstep
       !
       write(*,*) "  Iteration ", istep
       !
       Vks(1:g_rh%nr) = Vks(1:g_rh%nr) + dVks(1:g_rh%nr)
       !
       rhs0(1:g_rh%nr) = rhs(1:g_rh%nr)
       call kohn_sham_eq(.false., rhs)
       res = sqrt(dot_product(rhs(1:g_rh%nr), rhs(1:g_rh%nr))) / dble(g_rh%nr)
       write(*,*) "    delta Vks [eV] : ", res * htr2eV
       !
       if(res < conv_thr) exit
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
    if(istep >= electron_maxstep) then
       write(*,'(/,7x,"Not converged ! res = ",e12.5,/)') res
    else
       write(*,'(/,7x,"Converged ! iter = ",i0,/)') istep
    end if
    !
  end subroutine scf_loop
  !
  !
  !
  subroutine kohn_sham_eq(linit,rhs)
    !
    use gvec, only : g_wf, g_rh
    use kohn_sham, only : nbnd, eval, evec, calculation
    use k_point, only : kvec, nk, ksum_rho
    use rho_v, only : Vks, Vps, hartree_pot, xc_pot
    use diag_direct, only : direct
    use lobpcg, only : lobpcg_main
    !
    logical,intent(in) :: linit
    real(8),intent(out) :: rhs(g_rh%nr)
    !
    integer :: ik, istep, avestep
    !
    if(calculation == "direct") then
       do ik = 1, nk
          call direct(g_wf%npw,kvec(1:3,ik),evec(1:g_wf%npw,1:nbnd,ik),eval(1:nbnd,ik))
       end do
    else
       avestep = 0
       if(linit) call initialize_wf(evec(1:g_wf%npw,1:nbnd,1))
       do ik = 1, nk
          call lobpcg_main(g_wf%npw, kvec(1:3,ik),&
          &                evec(1:g_wf%npw,1:nbnd,ik),eval(1:nbnd,ik), istep)
          avestep = avestep + istep
          if(linit .and. ik /= nk) evec(1:g_wf%npw,1:nbnd,ik+1) &
          &                      = evec(1:g_wf%npw,1:nbnd,ik)
       end do
       write(*,*) "    Average LOBPCG steps : ", avestep / nk
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
  end subroutine kohn_sham_eq
  !>
  !! Initialize wave function with random number
  !!
  subroutine initialize_wf(psi)
    !
    use kohn_sham, only : nbnd
    use gvec, only : g_wf
    !
    complex(8),intent(out) :: psi(g_wf%npw,nbnd)
    !
    integer :: ibnd, nseed
    integer :: seed(256)
    real(8) :: rpsi(g_wf%npw,nbnd), ipsi(g_wf%npw,nbnd), norm
    !
    call random_seed(size=nseed)
    seed(1:nseed)=2
    call random_seed(put=seed)
    call random_number(rpsi)
    call random_number(ipsi)
    !
    psi(1:g_wf%npw,1:nbnd) = cmplx(rpsi(1:g_wf%npw,1:nbnd), &
    &                              ipsi(1:g_wf%npw,1:nbnd), 8)
    !
    do ibnd = 1, nbnd
       norm = sqrt(dble(dot_product(psi(1:g_wf%npw,ibnd), psi(1:g_wf%npw,ibnd))))
       psi(1:g_wf%npw,ibnd) = psi(1:g_wf%npw,ibnd) / norm
    end do
    !
  end subroutine initialize_wf
  !
end module scf
