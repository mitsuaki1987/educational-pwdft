module plot
  !
  implicit none
  !
contains
  !
  subroutine band_plot()
    !
    use constant, only : htr2ev
    use k_point, only : kvec, nk
    use kohn_sham, only : nbnd, eval
    !
    integer :: fo = 20, ibnd, ik
    real(8) :: xk(nk), dkvec(3), dxk
    !
    open(fo, file="band.dat")
    !
    xk(1) = 0.0d0
    do ik = 2, nk
       dkvec(1:3) = kvec(1:3,ik) - kvec(1:3,ik-1)
       dxk = sqrt(dot_product(dkvec(1:3), dkvec(1:3)))
       xk(ik) = xk(ik-1) + dxk
    end do
    !
    write(fo,*) "# k energy[eV]" 
    do ibnd = 1, nbnd
       do ik = 1, nk
          write(fo,*) xk(ik), eval(ibnd, ik)*htr2ev
       end do
       write(fo,*) ""
    end do
    !
    close(fo)
    !
  end subroutine band_plot
  !
  subroutine fermi_plot()
  end subroutine fermi_plot
  !
end module plot
