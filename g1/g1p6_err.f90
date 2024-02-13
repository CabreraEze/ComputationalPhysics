program g1p6_err
    use precision, only: pr=>dp
    use ode_solver

    implicit none  
    real(pr)                             :: x0(0:1), t0, t1, h, t, err_eu, err_rk2, err_rk4, x_val
    integer(pr)                          :: j, i, tmp
    integer(pr), parameter               :: n_max=7
    real(pr), dimension(0:10**n_max,0:1) :: x_eu, x_rk2, x_rk4
    character(len=20), parameter         :: fmt="(15(e28.20,5x))"

!#####################################
!           Parametros iniciles
!                    e
!         intervalos de evaluacion
!#####################################

    x0(0) = 1._pr
    x0(1) = 1._pr

    t0 = 0._pr
    t1 = 10._pr
 
    x_val = cos(10._pr) + sin(10._pr)

!#####################################
!
!#####################################

    open(11,file="ode_err.dat")

    do i=1,n_max

        x_rk2(:,0) = x0(0)
        x_rk2(:,1) = x0(1)
        x_rk4 = x_rk2
        x_eu = x_rk2

        tmp = 10**i
        h = (t1-t0)/tmp

        do j=0,tmp

            t = t0 + j*h

            if (j<=tmp-1) then
                call euler(h,dx,t,x_eu(j,:),x_eu(j+1,:))
                call rk2(h,dx,t,x_rk2(j,:),x_rk2(j+1,:))
                call rk4(h,dx,t,x_rk4(j,:),x_rk4(j+1,:))  
            else if (j==tmp) then
                err_eu  = abs((x_eu(j,0)-x_val)/x_val)
                err_rk2 = abs((x_rk2(j,0)-x_val)/x_val)
                err_rk4 = abs((x_rk4(j,0)-x_val)/x_val)
            endif
        enddo

        write(11,fmt) 1/h, err_eu, err_rk2, err_rk4

    enddo

    close(11)

contains

    function dx(t,x)
        use precision, only: pr=>dp
        implicit none
        real(pr), intent(in)                :: t, x(0:)
        real(pr), dimension(0:size(x)-1)    :: dx
        real(pr)                            :: k=1._pr

        dx(0) = x(1)
        dx(1) = -k*x(0)

    end function dx

end program g1p6_err