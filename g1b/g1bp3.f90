program g1bp3
    use precision, only: pr=>dp
    use ode_solver
    use, intrinsic :: iso_c_binding

    implicit none
    include "/usr/include/fftw3.f03"
    type(C_ptr)                     :: plan1, plan2
    character(len=50), parameter    :: fmt="(7(e28.20,5x))", loc="p3/"
    character(len=50)               :: tmp

    real(pr), parameter             :: g=9.8_pr, alp=1._pr/3._pr, bet=0.5_pr, gam=0.5_pr, m1=1._pr, pi2=8._pr*atan(1._pr)
    real(pr), dimension(0:3)        :: th0, th, fase
    real(pr)                        :: m2, l1, l2, t0, t1, t, h, x1, y1, x2, y2, sn1, sn2, cs1, cs2, e(0:4), e1, factor
    integer                         :: i, j, k, nmax, test, n_ft

    real(pr), allocatable           :: in1(:), in2(:)
    complex(pr), allocatable        :: out1(:), out2(:)

    !###############################
    !           PEND_1
    !###############################

    m2=m1*alp
    l1=g/gam
    l2=l1*bet

    th0(0) = 0._pr
    th0(1) = 0._pr
    th0(2) = 0.332_pr
    th0(3) = 0.845_pr

    t0 = 0._pr
    t1 = 150._pr

    nmax=1000

    open(33,file=trim(adjustl(loc))//"pend0.d",status="replace")
    write(33,"(7(a15,15x))") "#theta1", "theta2", "x1","y1","x2","y2"
        h = (t1-t0)/nmax
        th(:)=th0(:)
        
        do i=0,nmax
            t = t0 + i*h
            call rk4(h,df,t,th(:),th(:))
            
            do while (abs(th(0))>=pi2)
                if (th(0)<0) then
                    th(0) = th(0)+pi2
                else if (th(0)>0) then
                    th(0) = th(0)-pi2
                end if
            enddo
            do while (abs(th(1))>=pi2)
                if (th(1)<0) then
                    th(1) = th(1)+pi2
                else if (th(1)>0) then
                    th(1) = th(1)-pi2
                end if
            enddo
            
            sn1=sin(th(0))
            sn2=sin(th(1))
            cs1=cos(th(0))
            cs2=cos(th(1))

            y1 = -l1*cs1
            y2 = -l1*cs1-l2*cs2
            x1 = l1*sn1
            x2 = l1*sn1+l2*sn2

            if (t<=150) then
                write(33,fmt) th(1), th(2), x1, y1, x2, y2
            endif
        enddo
    close(33)

    !###############################
    !          PEND_2
    !###############################

    th0(0) = 0._pr
    th0(1) = 0._pr
    th0(2) = sqrt(1.125_pr)
    th0(3) = 0._pr

    t0 = 0._pr
    t1 = 150._pr

    open(34,file=trim(adjustl(loc))//"pend1.d",status="replace")
    write(34,"(7(a15,15x))") "#theta1", "theta2", "x1","y1","x2","y2"
        h = (t1-t0)/nmax
        th(:)=th0(:)
    
        do i=0,nmax
            t = t0 + i*h
            call rk4(h,df,t,th(:),th(:))

            do while (abs(th(0))>=pi2)
                if (th(0)<0) then
                    th(0) = th(0)+pi2
                else if (th(0)>0) then
                    th(0) = th(0)-pi2
                end if
            enddo
            do while (abs(th(1))>=pi2)
                if (th(1)<0) then
                    th(1) = th(1)+pi2
                else if (th(1)>0) then
                    th(1) = th(1)-pi2
                end if
            enddo

            sn1=sin(th(0))
            sn2=sin(th(1))
            cs1=cos(th(0))
            cs2=cos(th(1))

            y1 = -l1*cs1
            y2 = -l1*cs1-l2*cs2
            x1 = l1*sn1
            x2 = l1*sn1+l2*sn2

            if (t<=150) then
                write(34,fmt) th(1), th(2), x1, y1, x2, y2
            endif
        enddo
    close(34)
    
    !###############################
    !          POINCARE
    !###############################

    t0 = 0._pr
    t1 = 1000._pr

    th0(0) = 0._pr
    th0(2) = 0._pr
    th0(3) = 0._pr

    e = (/-0.745_pr,-0.6_pr,-0.58_pr,-0.575_pr,0._pr/)

    nmax=1000
    h = (t1-t0)/nmax

    do k=0,size(e)-1
        write(tmp,*) k
        open(35,file=trim(adjustl(loc))//"poincare_e"//trim(adjustl(tmp))//".d",status="replace")
        write(35,"(2(a15,15x))") "#theta2","p_theta2"
            do j=0,50
                th0(1) = acos(0.5_pr*sqrt(2._pr)) + j*0.01_pr
                th(:)=th0(:)
                do i=0,nmax
                    t = t0 + i*h
                    test = int(sign(1.1_pr,th(0)))
                    call rk4(h,df,t,th(:),th(:))
                    if ((test.ne.int(sign(1.1_pr,th(0)))).and.(fase(2)>0._pr)) then
                        write(35,fmt) fase(1), fase(3)
                    endif
                enddo
            enddo
        close(35)
    enddo

    !###############################
    !           FFT
    !###############################

    t0 = 0._pr
    t1 = 150._pr

    nmax=1000
    n_ft=nmax+1

    e1=0.18
    e1=gam*bet*alp*cos(e1)
    e1=0.745_pr-e1
    e1=e1/((alp+1)*gam)

    h = (t1-t0)/nmax
    factor = pi2
    factor = 1._pr/factor

    allocate(in1(n_ft),in2(n_ft))
    allocate(out1(n_ft/2+1),out2(n_ft/2+1))

    plan1 = fftw_plan_dft_r2c_1d(nmax,in1(:),out1(:),FFTW_MEASURE)
    plan2 = fftw_plan_dft_r2c_1d(nmax,in2(:),out2(:),FFTW_MEASURE)

        open(36,file=trim(adjustl(loc))//"fft_e0.d",status="replace")
        write(36,"(7(a15,15x))") "#frec1", "real1", "img1","real2","img2"

        open(37,file=trim(adjustl(loc))//"fft_e1.d",status="replace")
        write(37,"(7(a15,15x))") "#frec1", "real1", "img1","real2","img2"
            do k=0,1
                if (k==0) then 
                    th0(0) = 0._pr
                    th0(1) = 0._pr
                    th0(2) = sqrt(1.125_pr)
                    th0(3) = 0._pr
                else
                    th0(0) = acos(e1)
                    th0(1) = 0.18
                    th0(2) = 0._pr
                    th0(3) = 0._pr
                endif

                th(:)=th0(:)
                
                do i=0,nmax
                    t = t0 + i*h
                    call rk4(h,df,t,th(:),th(:))
                    
                    do while (abs(th(0))>=pi2)
                        if (th(0)<0) then
                            th(0) = th(0)+pi2
                        else if (th(0)>0) then
                            th(0) = th(0)-pi2
                        end if
                    enddo
                    do while (abs(th(1))>=pi2)
                        if (th(1)<0) then
                            th(1) = th(1)+pi2
                        else if (th(1)>0) then
                            th(1) = th(1)-pi2
                        end if
                    enddo
                    in1(i+1) = th(0)
                    in2(i+1) = th(1)
                enddo
                call fftw_execute_dft_r2c(plan1,in1(:),out1(:))
                call fftw_execute_dft_r2c(plan2,in2(:),out2(:))

                if (k==0) then 
                    do i = -n_ft/2,n_ft/2
                        if ( i < 0) then
                            write(36,fmt) factor*real(i,pr), &
                                        real( conjg( out1(-i+1) )/real(n_ft,pr) ), &
                                        aimag(conjg( out1(-i+1) )/real(n_ft,pr) ), &
                                        real( conjg( out2(-i+1) )/real(n_ft,pr) ), &
                                        aimag(conjg( out2(-i+1) )/real(n_ft,pr) )
                        else if (i == 0) then
                            write(36,fmt) factor*real(i,pr), &
                                        real(out1(1)/real(n_ft,pr)), &
                                        aimag(out1(1)/real(n_ft,pr)), &
                                        real(out2(1)/real(n_ft,pr)), &
                                        aimag(out2(1)/real(n_ft,pr))
                        else
                            write(36,fmt) factor*real(i,pr), &
                                        real( out1(i+1)/real(n_ft,pr) ), &
                                        aimag( out1(i+1)/real(n_ft,pr) ), &
                                        real( out2(i+1)/real(n_ft,pr) ), &
                                        aimag( out2(i+1)/real(n_ft,pr) )
                        endif
                    enddo
                else
                    do i = -n_ft/2,n_ft/2
                        if ( i < 0) then
                            write(37,fmt) factor*real(i,pr), &
                                        real( conjg( out1(-i+1) )/real(n_ft,pr) ), &
                                        aimag(conjg( out1(-i+1) )/real(n_ft,pr) ), &
                                        real( conjg( out2(-i+1) )/real(n_ft,pr) ), &
                                        aimag(conjg( out2(-i+1) )/real(n_ft,pr) )
                        else if (i == 0) then
                            write(37,fmt) factor*real(i,pr), &
                                        real(out1(1)/real(n_ft,pr)), &
                                        aimag(out1(1)/real(n_ft,pr)), &
                                        real(out2(1)/real(n_ft,pr)), &
                                        aimag(out2(1)/real(n_ft,pr))
                        else
                            write(37,fmt) factor*real(i,pr), &
                                        real( out1(i+1)/real(n_ft,pr) ), &
                                        aimag( out1(i+1)/real(n_ft,pr) ), &
                                        real( out2(i+1)/real(n_ft,pr) ), &
                                        aimag( out2(i+1)/real(n_ft,pr) )
                        endif
                    enddo
                endif
            enddo
        close(36)
        close(37)
    call fftw_destroy_plan(plan1)
    call fftw_destroy_plan(plan2)

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

        sn = sin(ang(0)-ang(1)) 
        cs = cos(ang(0)-ang(1))

        df(0) = ang(2)
        df(1) = ang(3)

        df(2) = (-1._pr-alp)*gam*sin(ang(0)) - &
                alp*bet*ang(3)*ang(3)*sn - &
                alp*cs*(ang(2)*ang(2)*sn-gam*sin(ang(1)))
        df(2) = df(2)/(1._pr+alp*sn*sn)
        
        df(3) = (1._pr+alp)*(ang(2)*ang(2)*sn-gam*sin(ang(1))) + &
                cs*((1._pr+alp)*gam*sin(ang(0))+alp*bet*ang(3)*ang(3)*sn)
        df(3) = df(3)/(bet*(1._pr+alp*sn*sn))

    end function df

    function p(ener,ang)
        use precision, only: pr=>dp
        implicit none
        real(pr), intent(in)                :: ener, ang(0:)
        real(pr), dimension(0:size(ang)-1)  :: p
        real(pr)                            :: sn, cs

        cs = cos(ang(0)-ang(1))
        sn = sin(ang(0)-ang(1))

        p(0) = ang(0)
        p(1) = ang(1)

        p(3) = m2*l2*l2*ang(3) + m2*l1*l2*ang(2)*cs
        !p(2) = (m1+m2)*l1*l1*ang(2) + m2*l1*l2*ang(3)*cs
        p(2) = p(3)/(bet*m1*l1*l1) + &
               sqrt((2._pr+alp*(1._pr-cs*cs-sn*sn))* &
               (2._pr*bet*bet*alp*(ener+gam*((1._pr+alp)*cos(ang(0))+bet*alp*cos(ang(1))))- &
               (p(3)/(l1*l1*m1))**2))/(bet*sqrt(2._pr*alp))
        p(2) = p(2)*m1*l1*l1

    end function p

end program g1bp3
