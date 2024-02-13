program g1bp3_frac
    use precision, only: pr=>dp
    use ode_solver

    implicit none
    character(len=50), parameter    :: fmt="(7(e28.20,5x))", loc="p3/"

    real(pr), parameter             :: g=9.8_pr, alp=1._pr, bet=1._pr, gam=1._pr, m1=1._pr, pi=4._pr*atan(1._pr)
    real(pr), dimension(0:3)        :: th0, th
    real(pr)                        :: m2, l1, l2, t0, t1, t, h, h1, h2, th1(0:1), th2(0:1)
    integer                         :: i, j, k, n, n1, n2

    !###############################
    !           FRACTAL
    !###############################
    
    m2=m1*alp
    l1=g/gam
    l2=l1*bet

    t0 = 0._pr
    t1 = 360._pr
    n = 24000 !int(t1*2._pr)
    
    th1(0) = -3._pr
    th1(1) = 3._pr
    n1 = 600
    
    th2(0) = 0._pr
    th2(1) = 3._pr
    n2 = 300
    
    open(33,file=trim(adjustl(loc))//"frac.d",status="replace")
    write(33,"(7(a15,15x))") "#theta1", "theta2", "t_flip"
        h = (t1-t0)/n
        h1 = (th1(1)-th1(0))/n1
        h2 = (th2(1)-th2(0))/n2
        do k=0,n1
            do j=0,n2
                th0(0) = th1(0)+k*h1
                th0(1) = th2(0)+j*h2
                th0(2) = 0._pr
                th0(3) = 0._pr
                
                if ((2._pr*cos(th0(0))+cos(th0(1)))>1._pr) then
                    ! write(33,fmt) th0(0), th0(1), t1
                    ! write(33,fmt) -th0(0), -th0(1), t1
                    cycle
                endif
                th(:) = th0(:)
                do i=0,n
                    t = t0 + i*h
                    call rk4(h,df,t,th(:),th(:))
                    if ((abs(th(0))>pi).or.(abs(th(1))>pi)) then
                        write(33,fmt) th0(0), th0(1), t
                        write(33,fmt) -th0(0), -th0(1), t
                        exit
                    endif
                enddo
            enddo
        enddo
    close(33)

contains

    !###############################
    !          FUNCIONES
    !###############################

function df(t,ang)
    use precision, only: pr=>dp
    implicit none
    real(pr), intent(in)                :: t, ang(0:)
    real(pr), dimension(0:size(ang)-1)  :: df
    real(pr)                            :: sn, cs
    !real(pr)                            :: denom, omega1, omega2

    sn = sin(ang(0)-ang(1)) 
    cs = cos(ang(0)-ang(1))

    df(0) = ang(2)
    df(1) = ang(3)

    df(2) = -(1._pr+alp)*gam*sin(ang(0)) - &
            alp*bet*ang(3)*ang(3)*sn - &
            alp*cs*(ang(2)*ang(2)*sn-gam*sin(ang(1)))
    df(2) = df(2)/(1._pr+alp*sn*sn)
    
    df(3) = (1._pr+alp)*(ang(2)*ang(2)*sn-gam*sin(ang(1))) + &
            cs*((1._pr+alp)*gam*sin(ang(0))+alp*bet*ang(3)*ang(3)*sn)
    df(3) = df(3)/(bet*(1._pr+alp*sn*sn))

    !#####################################

    ! sn = sin(ang(1)-ang(0)) 
    ! cs = cos(ang(1)-ang(0))

    ! df(0) = ang(2)
    ! df(1) = ang(3)

    ! denom = l1*(m1 + m2*sn*sn)
    ! omega1 = m2*l2*sn*ang(3)*ang(3) - (m1 + m2)*g*sin(ang(0))
    ! omega1 = omega1 + (m2*l1*sn*ang(2)*ang(2) + m2*g*sin(ang(1)))*cs
    ! omega1 = omega1/denom
    
    ! df(2) = omega1

    ! denom = l2*(-m1*m2-m2*m2*sn*sn)
    ! omega2 = (m2*l2*sn*ang(3)*ang(3)-(m1+m2)*g*sin(ang(0)))*m2*cs
    ! omega2 = omega2 + (m2*l1*sn*ang(2)*ang(2)+m2*g*sin(ang(1)))*(m1+m2)
    ! omega2 = omega2/denom

    ! df(3) = omega2

end function df

end program g1bp3_frac
