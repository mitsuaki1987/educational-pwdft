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
module energy
  !
  implicit none
  !
contains
  !
  subroutine total_e()
    !
    use constant, only : htr2ev
    !
    real(8) :: Etot
    !
    Etot = 0.0d0
    !
    write(*,*) "  Energy per u.c. [eV]"
    !
    call kinetic(Etot)
    call hartree(Etot)
    call atomic(Etot)
    call ewald(Etot)
    call xc(Etot)
    !
    write(*,*) "    Total energy : ", Etot*htr2ev
    !
  end subroutine total_e
  !
  subroutine kinetic(Etot)
    !
    use constant, only : htr2ev
    use atm_spec, only : bvec
    use k_point, only : nk, kvec, kgrd
    use kohn_sham, only : nbnd, ef, eval, evec
    use gvec, only : g_wf
    use libtetrabz, only : libtetrabz_occ
    !
    real(8),intent(inout) :: Etot
    !
    integer :: ik, ibnd, ipw
    real(8) :: occ(nbnd,nk), kgv(3), gg05(g_wf%npw), Ekin
    !
    eval(1:nbnd,1:nk) = eval(1:nbnd,1:nk) - ef
    call libtetrabz_occ(2,bvec,nbnd,kgrd,eval,kgrd,occ)
    eval(1:nbnd,1:nk) = eval(1:nbnd,1:nk) + ef
    !
    Ekin = 0.0d0
    do ik = 1, nk
       !
       do ipw = 1, g_wf%npw
          kgv(1:3) = kvec(1:3,ik) + matmul(bvec(1:3,1:3), dble(g_wf%mill(1:3,g_wf%map(ipw))))
          gg05(ipw) = 0.5d0 * dot_product(kgv,kgv)
       end do
       !
       do ibnd = 1, nbnd
          Ekin = Ekin + occ(ibnd,ik) &
          &           * dble(sum(conjg(evec(1:g_wf%npw,ibnd,ik)) &
          &                          * evec(1:g_wf%npw,ibnd,ik)  &
          &                          * gg05(1:g_wf%npw        )  ))
       end do
    end do
    !
    ! Spin
    !
    Ekin = Ekin * 2.0d0
    !
    write(*,*) "    Kinetic energy : ", Ekin*htr2ev
    Etot = Etot + Ekin
    !
  end subroutine kinetic
  !
  subroutine hartree(Etot)
    !
    use constant, only : pi, htr2ev
    use gvec, only : g_rh
    use atm_spec, only : bvec, Vcell
    use fftw_wrapper, only : fft_r2g
    use rho_v, only : rho
    !
    real(8),intent(inout) :: Etot
    !
    integer :: ir
    real(8) :: g3(3), EH
    complex(8) :: rhog(g_rh%nr)
    !
    call fft_r2G(rho, rhoG)
    !
    ! G = 0 : Compensation to the ionic potential 
    !
    EH = 0.0d0
    do ir = 2, g_rh%nr
       g3(1:3) = matmul(bvec(1:3,1:3), dble(g_rh%mill(1:3,ir)))
       EH = EH + 4.0d0 * pi * dble(conjg(rhog(ir))*rhog(ir)) &
       &       / dot_product(g3(1:3),g3(1:3))
    end do
    EH = EH * 0.5d0 * Vcell
    !
    write(*,*) "    Hartree energy : ", EH*htr2ev
    Etot = Etot + EH
    !
  end subroutine hartree
  !
  subroutine atomic(Etot)
    !
    use constant, only : htr2ev
    use rho_v, only : rho, Vps
    use gvec, only : g_rh
    use atm_spec, only : Vcell
    !
    real(8),intent(inout) :: Etot
    !
    real(8) :: Eps
    !
    Eps = dot_product(rho(1:g_rh%nr), Vps(1:g_rh%nr)) * Vcell / dble(g_rh%nr)
    !
    write(*,*) "    rho*V : ", Eps*htr2ev
    Etot = Etot + Eps
    !
  end subroutine atomic
  !
  subroutine ewald(Etot)
    !
    use constant, only : htr2ev, pi
    use atm_spec, only : avec, bvec, Vcell, atm, spec, nat, nelec
    !
    real(8),intent(inout) :: Etot
    !
    integer :: nmax3(3), iat, jat, i1, i2, i3, ii
    real(8) :: Eew, eta, cutoff, dtau(3), ZZ, gv(3), g2, phase, norm, rv(3)
    !
    Eew = 0.0d0
    !
    eta = 0.0d0
    do ii = 1, 3
       norm = sqrt(dot_product(bvec(1:3,ii), bvec(1:3,ii)))
       eta = eta + norm
    end do
    eta = eta / 3.0d0 / sqrt(2.0d0) / sqrt(2.0d0*pi)
    write(*,*) "      Ewald eta [Bohr^-1] : ", eta
    !
    ! Reciprocal space
    !
    cutoff = 10.0d0 * eta * 2.0d0
    do ii = 1, 3
       norm = sqrt(dot_product(avec(1:3,ii), avec(1:3,ii)))
       nmax3(ii) = ceiling(cutoff*norm/(2.0d0*pi))
    end do
    !
    write(*,*) "      Ewald grid (G) : ", nmax3(1:3)
    !
    do iat = 1, nat
       do jat = 1, nat
          !
          dtau(1:3) = atm(iat)%pos(1:3) - atm(jat)%pos(1:3)
          ZZ = spec(atm(iat)%ityp)%Zion * spec(atm(jat)%ityp)%Zion
          !
          do i3 = -nmax3(3), nmax3(3)
             do i2 = -nmax3(2), nmax3(2)
                do i1 = -nmax3(1), nmax3(1)
                   !
                   if(all((/i1,i2,i3/)==0)) cycle
                   !
                   gv(1:3) = matmul(bvec(1:3,1:3), dble((/i1,i2,i3/)))
                   g2 = dot_product(gv,gv)
                   phase = 2.0d0*pi*dot_product(dtau(1:3), dble((/i1,i2,i3/)))
                   !
                   Eew = Eew &
                   &   + 0.5d0 * ZZ * cos(phase) * 4.0d0 * pi * exp(-g2/(4.0d0*eta**2)) / (Vcell*g2)
                   !
                end do
             end do
          end do
       end do
    end do
    !
    Eew = Eew - 0.5d0 * pi * nelec**2 / (Vcell*eta**2)
    !
    ! Real space
    !
    cutoff = 10.0d0 / eta
    do ii = 1, 3
       norm = sqrt(dot_product(bvec(1:3,ii), bvec(1:3,ii)))
       nmax3(ii) = ceiling(cutoff*norm/(2.0d0*pi))
    end do
    !
    write(*,*) "      Ewald grid (R) : ", nmax3(1:3)
    !
    do iat = 1, nat
       do jat = 1, nat
          !
          dtau(1:3) = atm(iat)%pos(1:3) - atm(jat)%pos(1:3)
          ZZ = spec(atm(iat)%ityp)%Zion * spec(atm(jat)%ityp)%Zion
          !
          do i3 = -nmax3(3), nmax3(3)
             do i2 = -nmax3(2), nmax3(2)
                do i1 = -nmax3(1), nmax3(1)
                   !
                   if(all((/i1,i2,i3/)==0) .and. iat == jat) cycle
                   !
                   rv(1:3) = matmul(avec(1:3,1:3), dtau(1:3) + dble((/i1,i2,i3/)))
                   norm = sqrt(dot_product(rv(1:3), rv(1:3)))
                   !
                   Eew = Eew &
                   &   + 0.5d0 * ZZ * erfc(norm*eta) / norm
                   !
                end do
             end do
          end do
       end do
       !
       Eew = Eew - eta * spec(atm(iat)%ityp)%Zion**2 / sqrt(pi)
       !
    end do
    !
    write(*,*) "    Ewald energy : ", Eew*htr2ev
    Etot = Etot + Eew
    !
  end subroutine ewald
  !
  subroutine xc(Etot)
    !
    use constant, only : htr2ev, pi
    use gvec, only : g_rh
    use rho_v, only : rho
    use atm_spec, only : Vcell
    !
    real(8),intent(inout) :: Etot
    !
    integer :: ir
    real(8) :: rs, exc0, Exc
    !
    Exc = 0.0d0
    do ir = 1, g_rh%nr
       !
       if(rho(ir) > 1.0d-10) then
          !
          rs = (3.0d0 / (4.0d0 * pi * rho(ir)))**(1.0d0/3.0d0)
          !
          if(rs < 1.0d0) then
             exc0 = -3.0d0 / (4.0d0*pi) * (9.0d0*pi/4)**(1.0d0/3.0d0) / rs &
             &    - 0.0480d0 + 0.031d0*log(rs) - 0.0116d0*rs + 0.0020d0*rs*log(rs)
          else
             exc0 = -3.0d0 / (4.0d0*pi) * (9.0d0*pi/4)**(1.0d0/3.0d0) / rs &
             &   -  0.1423d0 / (1.0d0 + 1.0529d0*sqrt(rs) + 0.3334d0*rs)
          end if
          !
          Exc = Exc + exc0 * rho(ir)
          !
       end if
       !
    end do
    Exc = Exc * Vcell / dble(g_rh%nr)
    !
    write(*,*) "    XC energy : ", Exc*htr2ev
    Etot = Etot + Exc
    !
  end subroutine xc
  !
end module energy
