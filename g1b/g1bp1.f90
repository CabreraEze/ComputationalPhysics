program g1bp1
    use precision, only: pr=>dp
    use, intrinsic  :: iso_c_binding

    implicit none
    include "/usr/include/fftw3.f03"
    type(C_ptr)                   :: plan_rc, plan_cr
    character(len=20), parameter  :: fmt = "(10(e28.20,5x))"
    real(pr)                      :: x, factor, xmin, xmax
    real(pr), allocatable         :: in(:), out_cr(:)
    complex(pr), allocatable      :: out(:), in_cr(:)      
    integer                       :: i, n

    !######################
    !     COND INI
    !######################

    n = 1024
    xmin = 0._pr
    xmax = 4._pr

    !######################

    allocate(in(n), out(n/2+1), in_cr(n/2+1), out_cr(n))

    plan_rc = fftw_plan_dft_r2c_1d(n , in , out , FFTW_MEASURE)
    plan_cr = fftw_plan_dft_c2r_1d(n , in_cr , out_cr , FFTW_MEASURE)

    do i=1,n
        x = xmin + (real(i,pr)-1._pr)*(xmax-xmin)/(real(n,pr)-1._pr)
        in(i) = f(x)
    enddo

    !######################
    !     PLAN
    !######################

    call fftw_execute_dft_r2c(plan_rc, in , out)
    in_cr(:) = out(:)
    call fftw_execute_dft_c2r(plan_cr, in_cr , out_cr)

    open(33,file="g1bp1.d",status='replace')
    factor = 1._pr/(xmax-xmin)
    
    do i = -n/2,n/2
      if ( i < 0) then
        write(33,fmt) factor*real(i,pr),real( conjg( out(-i+1) )/real(n,pr) ),aimag(conjg( out(-i+1) )/real(n,pr) )
      else if (i == 0) then
        write(33,fmt) factor*real(i,pr),real(out(1)/real(n,pr)),aimag(out(1)/real(n,pr))
      else
        write(33,fmt) factor*real(i,pr),real( out(i+1)/real(n,pr) ),aimag( out(i+1)/real(n,pr) )
      endif
    enddo
    close(33)

    call fftw_destroy_plan(plan_cr)
    call fftw_destroy_plan(plan_rc)
    
    !######################

    contains
    function f(t)
        real(pr), intent(in)    :: t
        real(pr)                :: f, pi=4._pr*atan(1._pr)
            f=sin(0.5_pr*pi*t)+cos(20*pi*t)
    end function
end program g1bp1
