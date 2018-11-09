module lobpcg
  !
  implicit none
  !
contains
  !
  subroutine diag_ovrp()
    !
    integer,intent(in) :: nsub
    complex(8),intent(inout) :: hsub(nsub,nsub), ovlp(nsub,nsub)
    complex(8),intent(out) :: eval(nsub)
    !
    liwork = 5 * nsub + 3
    lwork = nsub*nsub + 2 * nsub
    lrwork = 3 * nsub*nsub + (4 + (int)log2(nsub) + 1) * nsub + 1
    !
    call zheev('V', 'U', nsub, ovlp, nsub, eval, work, lwork, rwork, info) 
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
    call zheev('V', 'U', nsub2, hsub, nsub, eval, work, lwork, rwork, info)
    !
    hsub(1:nsub, 1:nsub2) = matmul(ovlp(1:nsub,1:nsub2), hsub(1:nsub2,1:nsub2))
    hsub(1:nsub, nsub2+1:nsub) = 0.0d0
    !
  end subroutine diag_ovrp
  !
  subroutine lobpcg_main()
    !
    nsub = 3 * nbnd
    !
    call initialize(wxp(1:npw,1:nbnd,1))
    !
    hwxp(1:npw,1:nbnd,1) = 0.0d0
    !
    call h_psi(wxp(1:npw,1:nbnd,1), hwxp(1:npw,1:nbnd,1))
    !


    
  end subroutine lobpcg_main
  !
end module lobpcg
