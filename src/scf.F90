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
    logical :: linit
    real(8) :: alpha, res, &
    &          Fvec0(g_rh%nr), dVks(g_rh%nr), Fvec(g_rh%nr), dFvec(g_rh%nr), &
    &          jacob1(g_rh%nr,2:electron_maxstep), jacob2(g_rh%nr,2:electron_maxstep)
    !
    linit = .true.
    !
    do istep = 1, electron_maxstep
       !
       write(*,*) "  Iteration ", istep
       !
       call kohn_sham_eq(linit, Fvec)
       linit = .false.
       res = sqrt(dot_product(Fvec(1:g_rh%nr), Fvec(1:g_rh%nr))) / dble(g_rh%nr)
       write(*,*) "    delta Vks [eV] : ", res * htr2eV
       if(res < conv_thr) exit
       !
       if(istep > 1) then
          !
          ! Update Jacobian with dFvec
          !
          dFvec(1:g_rh%nr) = Fvec(1:g_rh%nr) - Fvec0(1:g_rh%nr)
          !
          jacob1(1:g_rh%nr,istep) = mixing_beta * dFvec(1:g_rh%nr) + dVks(1:g_rh%nr)
          do jstep = 2, istep - 1
             alpha = dot_product(jacob2(1:g_rh%nr,jstep), dFvec(1:g_rh%nr))
             jacob1(1:g_rh%nr,istep) = jacob1(1:g_rh%nr,istep) - jacob1(1:g_rh%nr,jstep) * alpha
          end do
          !
          alpha = dot_product(dFvec(1:g_rh%nr), dFvec(1:g_rh%nr))
          jacob2(1:g_rh%nr,istep) = dFvec(1:g_rh%nr) / alpha
          !
       end if
       !
       ! Compute dVks with new Jacobian & Fvec
       !
       dVks(1:g_rh%nr) = mixing_beta * Fvec(1:g_rh%nr)
       do jstep = 2, istep
          alpha = dot_product(jacob2(1:g_rh%nr,jstep), Fvec(1:g_rh%nr))
          dVks(1:g_rh%nr) = dVks(1:g_rh%nr) - jacob1(1:g_rh%nr,jstep) * alpha 
       end do
       Fvec0(1:g_rh%nr) = Fvec(1:g_rh%nr)
       !
       Vks(1:g_rh%nr) = Vks(1:g_rh%nr) + dVks(1:g_rh%nr)
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
  subroutine kohn_sham_eq(linit,Fvec)
    !
    use gvec, only : g_wf, g_rh
    use kohn_sham, only : nbnd, eval, evec, calculation
    use k_point, only : kvec, nk, ksum_rho
    use rho_v, only : Vks, Vps, hartree_pot, xc_pot
    use diag_direct, only : direct
    use lobpcg, only : lobpcg_main
    !
    logical,intent(in) :: linit
    real(8),intent(out) :: Fvec(g_rh%nr)
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
       Fvec(1:g_rh%nr) = Vps(1:g_rh%nr)
       call hartree_pot(Fvec)
       call xc_pot(Fvec)
       !
       Fvec(1:g_rh%nr) = Fvec(1:g_rh%nr) - Vks(1:g_rh%nr)
       !
    else
       Fvec(1:g_rh%nr) = 0.0d0
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
