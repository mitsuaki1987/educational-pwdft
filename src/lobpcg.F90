module lobpcg
  !
  implicit none
  !
contains
  !
  subroutine diag_ovrp(nsub,hsub,ovlp,esub)
    !
    integer,intent(in) :: nsub
    complex(8),intent(inout) :: hsub(nsub,nsub), ovlp(nsub,nsub)
    real(8),intent(out) :: esub(nsub)
    !
    integer :: lwork, isub, nsub2, info
    real(8) :: rwork(3*nsub-2)
    complex(8),allocatable :: work(:)
    !
    lwork = -1
    allocate(work(1))
    call zheev('V', 'U', nsub, ovlp, nsub, esub, work, lwork, rwork, info)
    lwork = nint(dble(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call zheev('V', 'U', nsub, ovlp, nsub, esub, work, lwork, rwork, info) 
    deallocate(work)
    !
    nsub2 = 0
    do isub = 1, nsub
       if(esub(isub) > 1.0d-14) then
          nsub2 = nsub2 + 1
          ovlp(1:nsub,nsub2) = ovlp(1:nsub,isub) / sqrt(esub(isub))
       end if
    end do
    ovlp(1:nsub,nsub2+1:nsub) = 0.0d0
    !
    hsub(1:nsub2,1:nsub2) = matmul(conjg(transpose(ovlp(1:nsub,1:nsub2))), &
    &                              matmul(hsub(1:nsub,1:nsub), ovlp(1:nsub,1:nsub2)))
    !
    lwork = -1
    allocate(work(1))
    call zheev('V', 'U', nsub2, hsub, nsub, esub, work, lwork, rwork, info)
    lwork = nint(dble(work(1)))
    deallocate(work)
    allocate(work(lwork))    
    call zheev('V', 'U', nsub2, hsub, nsub, esub, work, lwork, rwork, info)
    deallocate(work)
    !
    hsub(1:nsub, 1:nsub2) = matmul(ovlp(1:nsub,1:nsub2), hsub(1:nsub2,1:nsub2))
    hsub(1:nsub, nsub2+1:nsub) = 0.0d0
    !
  end subroutine diag_ovrp
  !
  subroutine lobpcg_main(npw,kvec,evec,eval,istep)
    !
    use kohn_sham, only : nbnd, calculation
    use hamiltonian, only : h_psi
    use atm_spec, only : bvec
    use gvec, only : g_wf
    !
    integer,intent(in) :: npw
    real(8),intent(in) :: kvec(3)
    real(8),intent(out) :: eval(nbnd)
    complex(8),intent(inout) :: evec(npw,nbnd)
    integer,intent(out) :: istep
    !
    integer :: ii, ibnd, nsub, ipw, cg_maxstep = 100
    real(8) :: norm, maxnorm, ekin(npw), pre(npw), gv(3), ekin0, &
    &          cg_thr = 1.0d-4, esub(3*nbnd)
    complex(8) :: wxp(npw,nbnd,3), hwxp(npw,nbnd,3), xp(npw,nbnd,2:3), &
    &             hsub(nbnd,3,3*nbnd), ovlp(3*nbnd,3*nbnd), rotmat(nbnd,3,nbnd,2:3)
    !
    do ipw = 1, g_wf%npw
       gv(1:3) = matmul(bvec(1:3,1:3), dble(g_wf%mill(1:3,g_wf%map(ipw))))
       ekin(ipw) = dot_product(gv,gv)
    end do
    !             
    nsub = 3 * nbnd
    xp(  1:npw,1:nbnd,2:3) = 0.0d0
    wxp( 1:npw,1:nbnd,1:3) = 0.0d0
    hwxp(1:npw,1:nbnd,1:3) = 0.0d0
    !
    wxp(1:npw,1:nbnd,2) = evec(1:npw,1:nbnd)
    !
    call h_psi(kvec,wxp(1:npw,1:nbnd,2), hwxp(1:npw,1:nbnd,2))
    !
    do ibnd = 1, nbnd
       esub(ibnd) = dble(dot_product(wxp(1:npw,ibnd,2), hwxp(1:npw,ibnd,2)))
    end do
    !
    if(calculation=="iterative") write(*,*) "Step  Residual"
    !
    do istep = 1, cg_maxstep
       !
       maxnorm = 0.0d0
       do ibnd = 1, nbnd
          !
          wxp(1:npw,ibnd,1) = hwxp(1:npw,ibnd,2) - esub(ibnd) * wxp(1:npw,ibnd,2)
          norm = sqrt(dble(dot_product(wxp(1:npw,ibnd,1), wxp(1:npw,ibnd,1))))
          maxnorm = max(norm, maxnorm)
          !
          ! Preconditioning
          !
          ekin0 = sum(dble(conjg(wxp(1:npw,ibnd,2))*(wxp(1:npw,ibnd,2))*ekin(1:npw)))
          pre(1:npw) = ekin(1:npw) / ekin0
          pre(1:npw) = (27.0d0 + pre(1:npw)*(18.0d0 + pre(1:npw)*(12.0d0 * pre(1:npw)*8.0d0))) &
          &          / (27.0d0 + pre(1:npw)*(18.0d0 + pre(1:npw)*(12.0d0 + pre(1:npw)*(8.0d0 + pre(1:npw)*16.0d0))))
          !
          wxp(1:npw,ibnd,1) = wxp(1:npw,ibnd,1) * pre(1:npw)
          !
          ! Normalize
          !
          norm = sqrt(dble(dot_product(wxp(1:npw,ibnd,1), wxp(1:npw,ibnd,1))))
          wxp(1:npw,ibnd,1) = wxp(1:npw,ibnd,1) / norm
          !
       end do
       if(calculation == "iterative") write(*,*) "    ", istep, maxnorm
       if(maxnorm < cg_thr) exit
       !
       call h_psi(kvec,wxp(1:npw,1:nbnd,1), hwxp(1:npw,1:nbnd,1))
       !
       call zgemm("C", "N", 3*nbnd, 3*nbnd, npw, &
       &          (1.0d0,0.0d0), wxp, npw, hwxp, npw, (0.0d0,0.0d0), hsub, 3*nbnd) 
       call zgemm("C", "N", 3*nbnd, 3*nbnd, npw, &
       &          (1.0d0,0.0d0), wxp, npw,  wxp, npw, (0.0d0,0.0d0), ovlp, 3*nbnd) 
       !
       call diag_ovrp(nsub,hsub,ovlp,esub)
       !
       rotmat(1:nbnd,1:3,1:nbnd,2) = hsub(1:nbnd,1:3,1:nbnd)
       rotmat(1:nbnd,  1,1:nbnd,3) = hsub(1:nbnd,  1,1:nbnd)
       rotmat(1:nbnd,  2,1:nbnd,3) = 0.0d0
       rotmat(1:nbnd,  3,1:nbnd,3) = hsub(1:nbnd,  3,1:nbnd)
       !
       call zgemm("N", "N", npw, 2*nbnd, 3*nbnd, &
       &          (1.0d0,0.0d0), wxp, npw, rotmat, 3*nbnd, (0.0d0,0.0d0), xp, npw)
       wxp(1:npw,1:nbnd,2:3) = xp(1:npw,1:nbnd,2:3)
       call zgemm("N", "N", npw, 2*nbnd, 3*nbnd, &
       &          (1.0d0,0.0d0), hwxp, npw, rotmat, 3*nbnd, (0.0d0,0.0d0), xp, npw)
       hwxp(1:npw,1:nbnd,2:3) = xp(1:npw,1:nbnd,2:3)
       !
       do ii = 2, 3
          do ibnd = 1, nbnd
             norm = sqrt(dble(dot_product(wxp(1:npw,ibnd,ii), wxp(1:npw,ibnd,ii))))
             wxp( 1:npw,ibnd,ii) =  wxp(1:npw,ibnd,ii) / norm
             hwxp(1:npw,ibnd,ii) = hwxp(1:npw,ibnd,ii) / norm
          end do
       end do
       !
    end do
    !
    if(istep >= cg_maxstep) then
       write(*,*) "      Not converged at kvec : ", kvec(1:3), ", norm :  ", maxnorm
    end if
    !
    evec(1:npw,1:nbnd) = wxp(1:npw,1:nbnd,2)
    eval(1:nbnd) = esub(1:nbnd)
    !
  end subroutine lobpcg_main
  !
end module lobpcg
