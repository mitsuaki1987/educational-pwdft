module lobpcg
  !
  implicit none
  !
contains
  !
  subroutine initialize(psi)
    !
    use solver, only : nbnd
    use gvec, only : wfft
    !
    complex(8),intent(out) :: psi(wfft%npw,nbnd)
    !
    real(8) :: rpsi(wfft%npw,nbnd), ipsi(wfft%npw,nbnd)
    !
    call random_number(rpsi)
    call random_number(ipsi)
    !
    psi(1:wfft%npw,1:nbnd) = cmplx(rpsi(1:wfft%npw,1:nbnd), &
    &                              ipsi(1:wfft%npw,1:nbnd), 8)
    !
  end subroutine initialize
  !
  subroutine diag_ovrp(nsub,hsub,ovlp,eval)
    !
    integer,intent(in) :: nsub
    complex(8),intent(inout) :: hsub(nsub,nsub), ovlp(nsub,nsub)
    real(8),intent(out) :: eval(nsub)
    !
    integer :: lwork, isub, nsub2, info
    real(8) :: rwork(3*nsub-2)
    complex(8),allocatable :: work(:)
    !
    lwork = -1
    allocate(work(1))
    call zheev('V', 'U', nsub, ovlp, nsub, eval, work, lwork, rwork, info)
    lwork = nint(dble(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call zheev('V', 'U', nsub, ovlp, nsub, eval, work, lwork, rwork, info) 
    deallocate(work)
    !
    nsub2 = 0
    do isub = 1, nsub
       if(eval(isub) > 1.0d-14) then
          nsub2 = nsub2 + 1
          ovlp(1:nsub,nsub2) = ovlp(1:nsub,isub) / sqrt(eval(isub))
       end if
    end do
    ovlp(1:nsub,nsub2+1:nsub) = 0.0d0
    !
    hsub(1:nsub2,1:nsub2) = matmul(conjg(transpose(ovlp(1:nsub,1:nsub2))), &
    &                              matmul(hsub(1:nsub,1:nsub), ovlp(1:nsub,1:nsub2)))
    !
    lwork = -1
    allocate(work(1))
    call zheev('V', 'U', nsub2, hsub, nsub, eval, work, lwork, rwork, info)
    lwork = nint(dble(work(1)))
    deallocate(work)
    allocate(work(lwork))    
    call zheev('V', 'U', nsub2, hsub, nsub, eval, work, lwork, rwork, info)
    deallocate(work)
    !
    hsub(1:nsub, 1:nsub2) = matmul(ovlp(1:nsub,1:nsub2), hsub(1:nsub2,1:nsub2))
    hsub(1:nsub, nsub2+1:nsub) = 0.0d0
    !
  end subroutine diag_ovrp
  !
  subroutine lobpcg_main(linit,npw,nbnd,kvec,evec,eval)
    !
    use solver, only : electron_maxstep
    use hamiltonian, only : h_psi
    !
    logical,intent(in) :: linit
    integer,intent(in) :: npw, nbnd
    real(8),intent(in) :: kvec(3)
    real(8),intent(out) :: eval(nbnd)
    complex(8),intent(out) :: evec(npw,nbnd)
    !
    integer :: ii, ibnd, iter, nsub
    real(8) :: norm
    complex(8) :: wxp(npw,nbnd,3), hwxp(npw,nbnd,3), xp(npw,nbnd,2:3), &
    &             hsub(nbnd,3,3*nbnd), ovlp(3*nbnd,3*nbnd), rotmat(nbnd,3,nbnd,2:3)
    !
    nsub = 3 * nbnd
    xp(  1:npw,1:nbnd,2:3) = 0.0d0
    wxp( 1:npw,1:nbnd,1:3) = 0.0d0
    hwxp(1:npw,1:nbnd,1:3) = 0.0d0
    !
    if(linit) call initialize(wxp(1:npw,1:nbnd,2))
    !
    call h_psi(kvec,wxp(1:npw,1:nbnd,2), hwxp(1:npw,1:nbnd,2))
    !
    do ibnd = 1, nbnd
       eval(ibnd) = dble(dot_product(wxp(1:npw,ibnd,2), hwxp(1:npw,ibnd,2)))
    end do
    !
    do iter = 1, electron_maxstep
       !
       do ibnd = 1, nbnd
          !
          wxp(1:npw,ibnd,1) = hwxp(1:npw,ibnd,2) - eval(ibnd) * wxp(1:npw,ibnd,2)
          ! todo precondition
          norm = sqrt(dble(dot_product(wxp(1:npw,ibnd,1), wxp(1:npw,ibnd,1))))
          wxp(1:npw,ibnd,1) = wxp(1:npw,ibnd,1) / norm
          !
       end do
       !
       call h_psi(kvec,wxp(1:npw,1:nbnd,1), hwxp(1:npw,1:nbnd,1))
       !
       call zgemm("C", "N", 3*nbnd, 3*nbnd, npw, &
       &          (1.0d0,0.0d0), wxp, npw, hwxp, npw, (0.0d0,0.0d0), hsub, 3*nbnd) 
       call zgemm("C", "N", 3*nbnd, 3*nbnd, npw, &
       &          (1.0d0,0.0d0), wxp, npw,  wxp, npw, (0.0d0,0.0d0), ovlp, 3*nbnd) 
       !
       call diag_ovrp(nsub,hsub,ovlp,eval)
       !
       rotmat(1:nbnd,1:3,1:nbnd,2) = hsub(1:nbnd,1:3,1:nbnd)
       rotmat(1:nbnd,  1,1:nbnd,3) = hsub(1:nbnd,  1,1:nbnd)
       rotmat(1:nbnd,  2,1:nbnd,3) = 0.0d0
       rotmat(1:nbnd,  3,1:nbnd,3) = hsub(1:nbnd,  3,1:nbnd)
       !
       call zgemm("C", "N", npw, 2*nbnd, 3*nbnd, &
       &          (1.0d0,0.0d0), wxp, npw, rotmat, 3*nbnd, (0.0d0,0.0d0), xp, npw)
       wxp(1:npw,1:nbnd,2:3) = xp(1:npw,1:nbnd,2:3)
       call zgemm("C", "N", npw, 2*nbnd, 3*nbnd, &
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
    evec(1:npw,1:nbnd) = wxp(1:npw,1:nbnd,2)
    !
  end subroutine lobpcg_main
  !
end module lobpcg
