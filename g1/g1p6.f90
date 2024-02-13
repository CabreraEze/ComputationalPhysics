program g1p6
    use precision, only: pr=>dp
    use ode_solver

    implicit none  
    real(pr)                             :: x0(0:1), t0, t1, h, t, ener_eu, ener_rk2, ener_rk4
    integer(pr)                          :: j
    integer(pr), parameter               :: n_max=10**2
    real(pr), dimension(0:n_max,0:1)     :: x_eu, x_rk2, x_rk4
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

!#####################################
!
!#####################################

    open(11,file="ode.dat")
    open(12,file="ode_ener.dat")

    x_rk2(:,0) = x0(0)
    x_rk2(:,1) = x0(1)
    x_rk4 = x_rk2
    x_eu = x_rk2

    do j=0,n_max

        h = (t1-t0)/n_max
        t = t0 + j*h

        if (j<=n_max-1) then
            call euler(h,dx,t,x_eu(j,:),x_eu(j+1,:))
            call rk2(h,dx,t,x_rk2(j,:),x_rk2(j+1,:))
            call rk4(h,dx,t,x_rk4(j,:),x_rk4(j+1,:))  
        endif

        ener_eu = ener(t,x_eu(j,:))
        ener_rk2 = ener(t,x_rk2(j,:))
        ener_rk4 = ener(t,x_rk4(j,:))
            
        write(11,fmt) t, x_eu(j,0), x_eu(j,1), x_rk2(j,0), x_rk2(j,1), x_rk4(j,0), x_rk4(j,1)
        write(12,fmt) t, ener_eu, ener_rk2, ener_rk4

    enddo

    close(11)
    close(12)

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

    function ener(t,x)
        use precision, only: pr=>dp
        implicit none
        real(pr), intent(in)    :: t, x(0:)
        real(pr)                :: ener
        real(pr)                :: k=1._pr, m=1._pr

        ener = 0.5_pr*m*x(1)*x(1)+ 0.5_pr*k*x(0)*x(0)

    end function ener
end program g1p6