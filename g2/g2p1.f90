program g2p1
    use precision, only: pr=>dp
    use pde_solver

    implicit none
    character(len=50), parameter            :: loc='p1/', fmt='(10(e16.8,5x))'

    real(pr), parameter                     :: k=237._pr, c=900._pr, p=2700._pr, l=1._pr, dt_step=0.3_pr
    real(pr)                                :: const, dt, dx, nu, t_real , x_real, tmp, tol
    real(pr)                                :: solucion, start, finish, err_explicita, err_implicita, err_cranknic
    real(pr), allocatable, dimension(:,:)   :: x_explicita, x_implicita, x_cranknic

    integer(pr)                             :: i, j, m, q, r, nx, nt
    
    !####################################
    !           CONSTANTES
    !####################################

    const = c*p*l*l
    const = k/const

    nx = 101
    nt = 3000

    dt = dt_step*const
    dx = 1._pr/real(nx-1,pr)
    dx = dx/l

    nu = dt/(dx*dx) 

    !####################################
    !           INC VAL           
    !####################################

    allocate(x_explicita(nx,2), x_implicita(nx,2), x_cranknic(nx,2))

    x_explicita(1,1) = tempta(0._pr)
    x_implicita(1,1) = tempta(0._pr)
    x_cranknic(1,1)  = tempta(0._pr)

    x_explicita(nx,1) = temptb(0._pr)
    x_implicita(nx,1) = temptb(0._pr)
    x_cranknic(nx,1)  = temptb(0._pr)

    do i = 2,nx-1
        x_explicita(i,1) = tempx(real(i-1,pr)*dx*l)
        x_implicita(i,1) = tempx(real(i-1,pr)*dx*l)
        x_cranknic(i,1)  = tempx(real(i-1,pr)*dx*l)
    enddo   

    open(33,file=trim(adjustl(loc))//"temp.d",status="replace")
    write(33,'(6(a10,18x))') '#t', 'x', 'explicita', 'implicita', 'cranknic'
        do i= 1, nx
            x_real = real(i-1,pr)*dx*l
            write(33,fmt) 0._pr, x_real, x_explicita(i,1), x_implicita(i,1), x_cranknic(i,1)
        enddo
        do i = 1, nt
            t_real = real(i-1,pr)*dt/const

            x_explicita(:,2) = 0._pr
            x_implicita(:,2) = 0._pr
            x_cranknic(:,2)  = 0._pr

            x_explicita(1,2) = tempta(t_real)
            x_implicita(1,2) = tempta(t_real)
            x_cranknic(1,2)  = tempta(t_real)

            x_explicita(nx,2) = temptb(t_real)
            x_implicita(nx,2) = temptb(t_real)
            x_cranknic(nx,2)  = temptb(t_real)

            call explicita(nu,x_explicita)
            call implicita(nu,x_implicita)
            call cranknic(nu,x_cranknic)

            if (mod(i,300)==0) then
                write(33,'(a2)') ''
                do j = 1,nx
                    x_real = real(j-1,pr)*dx*l
                    write(33,fmt) t_real, x_real, x_explicita(j,2), x_implicita(j,2), x_cranknic(j,2)
                enddo
            endif

            x_explicita(:,1) = x_explicita(:,2)
            x_implicita(:,1) = x_implicita(:,2)
            x_cranknic(:,1)  = x_cranknic(:,2)

        enddo
    close(33)

    deallocate(x_explicita, x_implicita, x_cranknic)

    !####################################
    !           CONSTANTES
    !####################################

    dt = 45e-4_pr
    nt = int(const*1800._pr/dt) + 1

    nx = 11 
    dx = 1._pr/real(nx-1,pr)

    !revision-----
    dt = const*0.3_pr
    dx = 0.01_pr/l
    nt = int(1800._pr/0.3_pr) + 1 
    nx = 101
    !-------------
    
    nu = dt/(dx*dx) 

    !####################################
    !           INC VAL           
    !####################################

    allocate(x_explicita(nx,2), x_implicita(nx,2), x_cranknic(nx,2))

    open(43,file=trim(adjustl(loc))//"evo_temp.d",status="replace")
    write(43,'(6(a10,18x))') '#t', 'x', 'explicita', 'implicita', 'cranknic'

    open(45,file=trim(adjustl(loc))//"temp_comp.d",status="replace")
    write(45,'(6(a10,18x))') '#t', 'x', 'explicita', 'implicita', 'cranknic', 'analitica'

        x_explicita(1,1) = tempta(0._pr)
        x_implicita(1,1) = tempta(0._pr)
        x_cranknic(1,1)  = tempta(0._pr)

        x_explicita(nx,1) = temptb(0._pr)
        x_implicita(nx,1) = temptb(0._pr)
        x_cranknic(nx,1)  = temptb(0._pr)

        do i = 2,nx-1
            x_explicita(i,1) = tempx(real(i-1,pr)*dx*l)
            x_implicita(i,1) = tempx(real(i-1,pr)*dx*l)
            x_cranknic(i,1)  = tempx(real(i-1,pr)*dx*l)
        enddo   

        write(43,'(a2)') ''
        do j = 1,nx
            x_real = real(j-1,pr)*dx*l
            write(43,fmt) 0._pr, x_real, x_explicita(j,1), x_implicita(j,1), x_cranknic(j,1)
        enddo

        open(34,file=trim(adjustl(loc))//"error_evo.d",status="replace")
        write(34,'(7(a15,10x))') '#t', 'x', 'err_explicita', 'err_implicita', 'err_cranknic'
            open(35,file=trim(adjustl(loc))//"cpu_time.d",status="replace")
            write(35,'(4(a10,18x))') '#time', 'cpu_time'

                call cpu_time(start)
                do i = 1, nt
                    t_real = real(i,pr)*dt/const

                    x_explicita(:,2) = 0._pr
                    x_implicita(:,2) = 0._pr
                    x_cranknic(:,2)  = 0._pr

                    x_explicita(1,2) = tempta(t_real)
                    x_implicita(1,2) = tempta(t_real)
                    x_cranknic(1,2)  = tempta(t_real)
        
                    x_explicita(nx,2) = temptb(t_real)
                    x_implicita(nx,2) = temptb(t_real)
                    x_cranknic(nx,2)  = temptb(t_real)

                    call explicita(nu,x_explicita)
                    call implicita(nu,x_implicita)
                    call cranknic(nu,x_cranknic)

                    if ((i==nt).or.(i==nt/10)) then
                        call cpu_time(finish)
                        write(35,fmt) t_real, (finish-start)

                        write(45,fmt) t_real, 0._pr, x_explicita(1,2), x_implicita(1,2), x_cranknic(1,2), analitica(0._pr,t_real)
                        do j = 2,nx-1
                            x_real = real(j-1,pr)*dx*l
                            solucion = analitica(x_real,t_real)

                            err_explicita = abs((solucion-x_explicita(j,2))/solucion)
                            err_implicita = abs((solucion-x_implicita(j,2))/solucion)
                            err_cranknic  = abs((solucion-x_cranknic(j,2))/solucion)

                            write(45,fmt) t_real, x_real, x_explicita(j,2), x_implicita(j,2), x_cranknic(j,2), solucion
                            write(34,fmt) t_real, x_real, err_explicita, err_implicita, err_cranknic
                        enddo
                        write(45,fmt) t_real, 1._pr, x_explicita(nx,2), x_implicita(nx,2), x_cranknic(nx,2), analitica(1._pr,t_real)
                        write(45,'(a4)') ''
                        write(45,'(a4)') ''

                        write(34,'(a4)') ''
                        write(34,'(a4)') ''
                    endif

                    x_explicita(:,1) = x_explicita(:,2)
                    x_implicita(:,1) = x_implicita(:,2)
                    x_cranknic(:,1)  = x_cranknic(:,2)

                    if ((mod(i,300)==0).or.(i==nt)) then
                        write(43,'(a2)') ''
                        do j = 1,nx
                            x_real = real(j-1,pr)*dx*l
                            write(43,fmt) t_real, x_real, x_explicita(j,2), x_implicita(j,2), x_cranknic(j,2)
                        enddo
                    endif  

                enddo
            close(35)
        !write(34,fmt) nu
        close(34)

    close(45)
    close(43)

    deallocate(x_explicita, x_implicita, x_cranknic)

    !####################################
    !           CONSTANTES
    !####################################
    nx = 101 
    dx = 1._pr/real(nx-1,pr)

    m =  1

    tol = 1e-3_pr

    !####################################
    !           INC VAL           
    !####################################

    allocate(x_explicita(nx,2), x_implicita(nx,2), x_cranknic(nx,2))

    open(36,file=trim(adjustl(loc))//"cpu_time_comparado.d",status="replace")
    write(36,'(4(a10,18x))') '#method', 'cpu_time', 'dt', 'error'
    do q = 1,3
        do r = 1, m
            call cpu_time(start)
            dt = 5*10**(-real(3+r,pr))
            nt = int(const*180._pr/dt,pr) + 1
            nu = dt/(dx*dx) 
            !write(*,*) start, dt, r, nt, q

            select case (q)
                case (1)
                    x_explicita(1,1) = tempta(0._pr)
                    x_explicita(nx,1) = temptb(0._pr)
                    do i = 2,nx-1
                    x_explicita(i,1) = tempx(real(i-1,pr)*dx*l)
                    enddo

                    do i = 1, nt
                        t_real = real(i,pr)*dt/const
                        x_explicita(:,2) = 0._pr

                        x_explicita(1,2) = tempta(t_real)
                        x_explicita(nx,2) = temptb(t_real)
                        call explicita(nu,x_explicita)
                        x_explicita(:,1) = x_explicita(:,2)

                        if (i==nt) then
                            tmp = 0._pr
                            do j = 2,nx-1
                                x_real = real(j-1,pr)*dx*l
                                solucion = analitica(x_real,t_real)
                                err_explicita = abs((solucion-x_explicita(j,2))/solucion)
                                if (err_explicita>=tmp) then
                                    tmp = err_explicita
                                endif
                            enddo
                        endif
                    enddo
                    if (tmp<=tol) then
                        exit
                    endif
                case (2)
                    x_implicita(1,1) = tempta(0._pr)
                    x_implicita(nx,1) = temptb(0._pr)
                    do i = 2,nx-1
                    x_implicita(i,1) = tempx(real(i-1,pr)*dx*l)
                    enddo

                    do i = 1, nt
                        t_real = real(i,pr)*dt/const
                        x_implicita(:,2) = 0._pr

                        x_implicita(1,2) = tempta(t_real)
                        x_implicita(nx,2) = temptb(t_real)
                        call implicita(nu,x_implicita)
                        x_implicita(:,1) = x_implicita(:,2)

                        if (i==nt) then
                            tmp = 0._pr
                            do j = 2,nx-1
                                x_real = real(j-1,pr)*dx*l
                                solucion = analitica(x_real,t_real)
                                err_implicita = abs((solucion-x_implicita(j,2))/solucion)
                                if (err_implicita>=tmp) then
                                    tmp = err_implicita
                                endif
                            enddo
                        endif
                    enddo
                    if (tmp<=tol) then
                        exit
                    endif

                    write(*,*) "estoy trabajando aqui"
                case (3)
                    x_cranknic(1,1)  = tempta(0._pr)
                    x_cranknic(nx,1)  = temptb(0._pr)
                    do i = 2,nx-1
                    x_cranknic(i,1)  = tempx(real(i-1,pr)*dx*l)
                    enddo

                    do i = 1, nt
                        t_real = real(i,pr)*dt/const
                        x_cranknic(:,2) = 0._pr

                        x_cranknic(1,2) = tempta(t_real)
                        x_cranknic(nx,2) = temptb(t_real)
                        call cranknic(nu,x_cranknic)
                        x_cranknic(:,1) = x_cranknic(:,2)

                        if (i==nt) then
                            tmp = 0._pr
                            do j = 2,nx-1
                                x_real = real(j-1,pr)*dx*l
                                solucion = analitica(x_real,t_real)
                                err_cranknic = abs((solucion-x_cranknic(j,2))/solucion)
                                if (err_cranknic>=tmp) then
                                    tmp = err_cranknic
                                endif
                            enddo
                        endif
                    enddo
                    if (tmp<=tol) then
                        exit
                    endif
            end select    
        enddo
        call cpu_time(finish)
        write(36,fmt) real(q,pr), (finish-start), dt, tmp
    enddo
    close(36)     
    
    deallocate(x_explicita, x_implicita, x_cranknic)

contains
    function tempx(x)
        use precision, only: pr=>dp
        implicit none
        real(pr), intent(in)    :: x
        real(pr)                :: tempx, temp0=100._pr

        tempx = temp0
    end function tempx

    function tempta(t)
        use precision, only: pr=>dp
        implicit none
        real(pr), intent(in)    :: t
        real(pr)                :: tempta, tempa=0._pr

        tempta = tempa
    end function tempta

    function temptb(t)
        use precision, only: pr=>dp
        implicit none
        real(pr), intent(in)    :: t
        real(pr)                :: temptb, tempb=0._pr

        temptb = tempb
    end function temptb

    function analitica(x,t)
        use precision, only: pr=>dp
        implicit none
        real(pr), intent(in)    :: x, t
        real(pr)                :: analitica, temp0=100._pr, kn, pi=4._pr*atan(1._pr), adim, constante
        integer                 :: s
        
        adim = c*p
        adim = k/adim

        constante = 4._pr*temp0/pi

        analitica = 0._pr
        do s = 17,1,-2
            kn = real(s,pr)*pi/l
            analitica = analitica + constante*sin(kn*x)/(exp(kn*kn*adim*t)*real(s,pr))
        enddo 
    end function analitica

end program g2p1