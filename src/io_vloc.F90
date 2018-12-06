module io_vloc
  !
  implicit none
  !
  real(8),allocatable :: &
  & Vks(:), & !< (g_rh%nr) Kohn-Sham potential [Htr]
  & Vps(:), & !< (g_rh%nr) Pseudopotential [Htr]
  & rho(:) !< (g_rh%nr) Charge density
  !
contains
  !
  subroutine read_vks()
    !
    use constant, only : pi, htr2ev
    use atm_spec, only : nat, atm, bvec, avec
    use gvec, only : g_rh
    !
    integer :: itemp(3), fi = 10, iat
    character(256) :: ctemp
    real(8) :: avec0(3,3), &
    &          Vks0(1:g_rh%nft(1)+1,1:g_rh%nft(2)+1,1:g_rh%nft(3)+1)
    !
    open(fi, file = "vks.xsf")
    !
    read(fi,*) ctemp
    read(fi,*) ctemp
    read(fi,*) avec0(1:3,1:3)
    avec0(1:3,1:3) = avec0(1:3,1:3) / 0.529177249d0
    if(any(abs(avec0(1:3,1:3)-avec(1:3,1:3)) > 1.0d-3)) then
       write(*,*) "Error in read_vks"
       write(*,*) "Direct lattice vector is different."
       stop 'error in read_vks'
    end if
    read(fi,*) ctemp
    read(fi,*) itemp(1:2)
    if(nat /= itemp(1)) then
       write(*,*) "Error in read_vks"
       write(*,*) "Number of atoms is different."
       stop 'error in read_vks'
    end if
    do iat = 1, nat
       read(fi,*) ctemp, avec0(1:3,1)
       avec0(1:3,1) = avec0(1:3,1) / 0.529177249d0
       avec0(1:3,1) = matmul(avec0(1:3,1), bvec(1:3,1:3)) / (2.0d0*pi)
       if(any(abs(avec0(1:3,1)-atm(iat)%pos(1:3)) > 1.0d-3)) then
          write(*,*) "Error in read_vks"
          write(*,*) "Position of atom ", iat, " is different."
          stop 'error in read_vks'
       end if
    end do
    read(fi,*) ctemp
    read(fi,*) ctemp
    read(fi,*) ctemp 
    read(fi,*) itemp(1:3)
    if(any(itemp(1:3) /= g_rh%nft(1:3)+1)) then
       write(*,*) "Error in read_vks"
       write(*,*) "FFT grid is different."
       stop 'error in read_vks'
    end if
    read(fi,*) avec0(1:3,1)
    read(fi,*) avec0(1:3,1:3)
    read(fi,*) Vks0(1:itemp(1),1:itemp(2),1:itemp(3))
    Vks(1:g_rh%nr) = reshape(Vks0(1:g_rh%nft(1),1:g_rh%nft(2),1:g_rh%nft(3)), (/g_rh%nr/))
    !
    write(*,*) "  Average potential [eV] : ", &
    &  sum(Vks(1:g_rh%nr)) / dble(g_rh%nr)
    Vks(1:g_rh%nr) = Vks(1:g_rh%nr) / htr2ev
    !
  end subroutine read_vks
  !
  subroutine generate_vps()
    !
    use constant, only : pi
    use gvec, only : g_rh
    use atm_spec, only : bvec, spec, atm, nat, ntyp, Vcell
    use fftw_wrapper, only : fft_g2r
    !
    integer :: jtyp, iat, nat2, ipw
    real(8) :: prod(nat), pos2(3,nat), gv(3), ctail, glen, dr, &
    &          frmfac(g_rh%nr)
    complex(8) :: strfac(g_rh%nr), VpsG(g_rh%nr)
    real(8),allocatable :: sinc(:)
    !
    VpsG(1:g_rh%nr) = 0.0d0
    !
    do jtyp = 1, ntyp
       !
       nat2 = 0
       do iat = 1, nat
          if(atm(iat)%ityp == jtyp) then
             nat2 = nat2 + 1
             pos2(1:3,nat2) = atm(iat)%pos(1:3)
          end if
       end do
       !
       ! Structure factor
       !
       do ipw = 1, g_rh%nr
          prod(1:nat2) = 2.0d0 * pi * matmul(dble(g_rh%mill(1:3,ipw)),pos2(1:3,1:nat2))
          strfac(ipw) = sum(exp(cmplx(0.0d0, prod(1:nat), 8))) / Vcell
       end do
       !
       ! Form factor
       !
       allocate(sinc(spec(jtyp)%mmax))
       dr = spec(jtyp)%psr(spec(jtyp)%mmax) / dble(spec(jtyp)%mmax - 1)
       !
       ! G = 0 : Compensation to the Hartree potential 
       !
       ctail = 0.5d0 * spec(jtyp)%Zion * spec(jtyp)%psr(spec(jtyp)%mmax)**2
       !
       sinc(1:spec(jtyp)%mmax) = dr
       sinc(1              ) = dr * 0.5d0
       sinc(spec(jtyp)%mmax) = dr * 0.5d0
       !
       frmfac(1) = sum( sinc(1:spec(jtyp)%mmax) &
       &    * spec(jtyp)%psr(1:spec(jtyp)%mmax)**2 &
       &    * spec(jtyp)%psV(1:spec(jtyp)%mmax)) &
       &         + ctail
       frmfac(1) = 4.0d0 * pi * frmfac(1)
       !
       ! G /= 0
       !
       do ipw = 2, g_rh%nr
          !
          gv(1:3) = matmul(bvec(1:3,1:3), dble(g_rh%mill(1:3,ipw)))
          glen = sqrt(dot_product(gv(1:3),gv(1:3)))
          sinc(1:spec(jtyp)%mmax) = glen*spec(jtyp)%psr(1:spec(jtyp)%mmax)
          !
          ctail = - spec(jtyp)%Zion * cos(sinc(spec(jtyp)%mmax)) / glen**2
          !
          sinc(2:spec(jtyp)%mmax) = sin(sinc(1:spec(jtyp)%mmax)) &
          &                       / sinc(    1:spec(jtyp)%mmax) &
          &                       * dr
          sinc(1              ) = dr * 0.5d0
          sinc(spec(jtyp)%mmax) = sinc(spec(jtyp)%mmax) * 0.5d0
          !
          frmfac(ipw) = sum( sinc(1:spec(jtyp)%mmax) &
          &      * spec(jtyp)%psr(1:spec(jtyp)%mmax)**2 &
          &      * spec(jtyp)%psV(1:spec(jtyp)%mmax)) &
          &           + ctail
          frmfac(ipw) = 4.0d0 * pi * frmfac(ipw)
          !
       end do
       !
       deallocate(sinc)
       !
       VpsG(1:g_rh%nr) = VpsG(1:g_rh%nr) + strfac(1:g_rh%nr) * frmfac(1:g_rh%nr)
       !
    end do
    !
    call fft_g2r(VpsG, Vps)
    !
  end subroutine generate_vps
  !
end module io_vloc
