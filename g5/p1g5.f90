program p1g5
    use precision, only: pr=>qp
    use molecular_dyn
    implicit none
    
    real(pr), allocatable        :: r(:,:), v(:,:), f(:,:), req(:,:), veq(:,:)
    real(pr)                     :: ro, dt, l, a, rcut, t0, t, lambda, temp_eq
    real(pr), allocatable        :: e_cin(:), e_pot(:), e_tot(:), press(:), temp(:)
    real(pr)                     :: de_tot, de_cin, de_pot, dtemp, dpress, timeini, timefinal
    real(pr), parameter          :: pi = 4._pr*atan(1._pr)

    integer                      :: base, n, i, j, t_iter, t_run, t_eq, t_scal, varmx
    !integer, allocatable         :: head(:), list(:)

    character(len=50), parameter :: fmt='(10(e18.10,2x))', loc='p1/', pfmt='(2(a3,e18.10))'

    logical                      :: gauss, dt_variable, rcut_variable, base_variable

    !#########################################
    !       INITIAL VALUE
    !#########################################

    !ro: densidad del sistema
    !base: cantidad de celdas unitarias que abarca un lado de la caja de simulacion
    !rcut: distancia de corte para primeros vecinos
    !dt: intervalo de tiempo
    !t0: temperatura inicial

    !t_eq: numero de pasos de simulacion para equilibrar el sistema
    !t_scal: cada t_scal durante la equilibracion se reescalaran las velocidades para mantener la temperatura constante
    !trun: pasos de corrida de la simulacion

    ro = 0.8_pr
    base = 4
    rcut = 2.5_pr 
    dt = 0.005_pr
    t0 = 1.1_pr

    t_eq = 1000
    t_scal = 50
    t_run = 1000

    !·····DISTRIBUCION GAUSSIANA DE VELOCIDADES COMO CONDICION INICIAL·····
    gauss = .true.

    !·····VARIABILIDAD DE RCUT, DT Y BASE·····
    dt_variable = .true.
    rcut_variable = .false.
    base_variable = .false.
    varmx = 7

    !#########################################

    if ( (dt_variable.and.(rcut_variable.or.base_variable)).or.(rcut_variable.and.base_variable) ) then
        write(*,*) '··only one variable parameter··'
        stop
    else if (base_variable.and.(varmx>8)) then
        write(*,*) '··maximum value for the base exceeded··'
        stop
    endif 

    if (.not.base_variable) then
        write(*,*) 'Setting initial condition···'
        l = 0._pr
        a = 0._pr
        call initialize(r, v, ro ,base, l ,a)
        n=size(r,dim=2)
        if (gauss) then
            do i=1,3
                do j=1,n
                    v(i,j) = gasdev()/(2._pr*pi*t0)
                enddo
            enddo
        endif

        allocate(f(3,n), req(3,n), veq(3,n))
        f(:,:)    = 0._pr
        temp_eq = 0._pr

        !termalizacion del sistema
            do t_iter = 1,t_eq
                call forces(0, r, f, rcut, l)
                call dynamic_eq(r, v, f, rcut, l ,dt)
                call pbc(r, l)

                !cada t_escal se reescalea las velocidades para 
                !mantener la temperatura constante
                if (mod(t_iter, t_scal)==0) then
                    temp_eq = temperature(v)
                    lambda = sqrt(t0/temp_eq)
                    v = v*lambda
                endif

            enddo
            req = r
            veq = v
        !·························
        
    endif

    if (dt_variable.or.rcut_variable.or.base_variable) then

        !···FILES···
            if (dt_variable) then
                open(33,file=trim(adjustl(loc))//'energy_vs_dt.d',status='replace')
                write(33,'(8(a20,2x))') '#dt', 'Etot/Ekin'
            elseif (base_variable) then
                open(34,file=trim(adjustl(loc))//'time_vs_dim.d',status='replace')
                write(34,'(8(a20,2x))') '#time', 'dimension'
            elseif (rcut_variable) then   
                open(35,file=trim(adjustl(loc))//'energy_vs_rcut.d',status='replace')
                write(35,'(8(a20,2x))') '#rcut', 'Etot/Ekin', 'Etot/Epot'
            endif
        !············
        
        do i=1,varmx
            
            !···DINAMIC PARAMETERS···
                if (dt_variable) then
                    dt = 0.04_pr/real(i*2,pr)
                    t_run = int(5._pr/dt)
                    write(*,*) 'implementing dt:'
                    write(*,'(e18.10)') dt
                elseif (base_variable) then
                    base = 2+i

                    write(*,*) 'implementing base:'
                    write(*,'(i3)') base

                    l = 0._pr
                    a = 0._pr
                    call initialize(r, v, ro ,base, l ,a)
                    n=size(r,dim=2)
                    
                    allocate(f(3,n), req(3,n), veq(3,n))
                    f(:,:)    = 0._pr
                    temp_eq = 0._pr

                    !termalizacion del sistema
                        do t_iter = 1,t_eq
                            call forces(0, r, f, rcut, l)
                            call dynamic_eq(r, v, f, rcut, l ,dt)
                            call pbc(r, l)

                            !cada t_escal se reescalea las velocidades para 
                            !mantener la temperatura constante
                            if (mod(t_iter, t_scal)==0) then
                                temp_eq = temperature(v)
                                lambda = sqrt(t0/temp_eq)
                                v = v*lambda
                            endif

                        enddo
                        req = r
                        veq = v
                    !·························

                elseif (rcut_variable) then   
                    rcut = 1._pr + 4._pr*real(i-1,pr)/real(varmx,pr)
                    dt = 0.0005_pr

                    write(*,*) 'implementing rcut:'
                    write(*,'(e18.10)') rcut
                endif
            !························
            !
            !
            ! 
            !···SIMULATION···
                r = req
                v = veq
                allocate(e_tot(t_run), e_cin(t_run), e_pot(t_run), press(t_run))
                call cpu_time(timeini)
                e_cin = 0._pr   
                e_pot = 0._pr
                press = 0._pr
                do t_iter = 1,t_run

                    t = dt*real(t_iter,pr)

                    if (dt_variable.or.rcut_variable) then
                        call forces(1, r, f, rcut, l, e_pot(t_iter), press(t_iter))
                        e_cin(t_iter) = kinetic(v) 
                        e_tot(t_iter) = e_cin(t_iter) + e_pot(t_iter)
                    else
                        call forces(0, r, f, rcut, l)
                    endif

                    call dynamic_eq(r, v, f, rcut, l ,dt)
                    call pbc(r, l)

                enddo
                call cpu_time(timefinal)
            !················
            !
            !
            ! 
            !···RESULTS···
                if (dt_variable) then
                    de_tot = 0._pr
                    de_cin = 0._pr
                    do j=1,t_run
                        de_cin  = de_cin + (e_cin(j) - sum(e_cin)/real(t_run,pr))**2
                        de_tot  = de_tot + (e_tot(j) - sum(e_tot)/real(t_run,pr))**2
                    enddo
                    write(33,fmt) dt, sqrt(de_tot/de_cin)
                elseif (base_variable) then
                    write(34,fmt) real(base,pr), timeini-timefinal
                elseif (rcut_variable) then
                    de_tot = 0._pr
                    de_cin = 0._pr
                    de_pot = 0._pr
                    do j=1,t_run
                        de_cin  = de_cin + (e_cin(j) - sum(e_cin)/real(t_run,pr))**2
                        de_pot  = de_pot + (e_pot(j) - sum(e_pot)/real(t_run,pr))**2
                        de_tot  = de_tot + (e_tot(j) - sum(e_tot)/real(t_run,pr))**2
                    enddo   
                    write(35,fmt) rcut, sqrt(de_tot/de_cin), sqrt(de_tot/de_pot)
                endif
            !·············

            deallocate(e_tot, e_cin, e_pot, press)
            if (base_variable) then
                deallocate(r,v,f,req,veq)
            endif
        enddo

    else

        !···FILES···
            open(33,file=trim(adjustl(loc))//'obser.d',status='replace')
            write(33,'(8(a20,2x))') '#t', 'etot', 'ecin', 'epot', 'temp', 'pres'
            open(34,file=trim(adjustl(loc))//'dyn.d',status='replace')
            write(34,'(8(a20,2x))') '#Rx', 'Ry', 'Rz'
            open(35,file=trim(adjustl(loc))//'distro.d',status='replace')
            write(35,'(8(a20,2x))') '#V', 'Vx', 'Vy', 'Vz', 'Veq', 'Veqx', 'Veqy', 'Veqz'
        !············
        !
        !
        ! 
        !···SIMULATION···
            write(*,*) 'Simulation in progress···'
            allocate(e_tot(t_run), e_cin(t_run), e_pot(t_run), temp(t_run), press(t_run))
            
            r = req
            v = veq
            e_cin = 0._pr   
            e_pot = 0._pr
            press = 0._pr
            temp  = 0._pr

            ! do i=1,n
            !     write(35,fmt) sqrt(sum(v(:,i)*v(:,i))), v(1,i), v(2,i), v(3,i)
            ! enddo

            do t_iter = 1,t_run

                t = dt*real(t_iter,pr)
                call forces(1, r, f, rcut, l, e_pot(t_iter), press(t_iter))

                e_cin(t_iter) = kinetic(v) 
                temp(t_iter) = temperature(v)
                press(t_iter) = press(t_iter)/(3._pr*real(n,pr)*l**3) + ro*temp(t_iter)
                e_pot(t_iter) = e_pot(t_iter)
                e_tot(t_iter) = e_cin(t_iter) + e_pot(t_iter)
    
                write(33,fmt) t, e_tot(t_iter), e_cin(t_iter), e_pot(t_iter), temp(t_iter), press(t_iter)

                call dynamic_eq(r, v, f, rcut, l ,dt)
                call pbc(r, l)

                do i = 1,n
                    write(34,fmt) r(1,i), r(2,i), r(3,i)
                enddo
                write(34,*) ''
                write(34,*) ''

            enddo
        !················
        !
        !
        ! 
        !···RESULTS···
            write(*,*) 'Getting results···'
            do i=1,n
                write(35,fmt) sqrt(sum(v(:,i)*v(:,i))), v(1,i), v(2,i), v(3,i), &
                            sqrt(sum(veq(:,i)*veq(:,i))), veq(1,i), veq(2,i), veq(3,i)
            enddo

            de_tot = 0._pr
            de_cin = 0._pr
            de_pot = 0._pr
            dpress = 0._pr
            dtemp = 0._pr
            do i=1,t_run
                de_cin  = de_cin + (e_cin(i) - sum(e_cin)/real(t_run,pr))*(e_cin(i) - sum(e_cin)/real(t_run,pr))
                de_pot  = de_pot + (e_pot(i) - sum(e_pot)/real(t_run,pr))*(e_pot(i) - sum(e_pot)/real(t_run,pr))
                de_tot  = de_tot + (e_tot(i) - sum(e_tot)/real(t_run,pr))*(e_tot(i) - sum(e_tot)/real(t_run,pr))
                dpress  = dpress + (press(i) - sum(press)/real(t_run,pr))*(press(i) - sum(press)/real(t_run,pr))
                dtemp   = dtemp + (temp(i) - sum(temp)/real(t_run,pr))*(temp(i) - sum(temp)/real(t_run,pr))
            enddo

            call cpu_time(timefinal)

            write(33,*) '#RESULTADOS'
            write(33,*) '#Energia total'
            write(33,pfmt) '#',sum(e_tot)/real(t_run,pr), ' +-', sqrt(de_tot/real(t_run,pr))
            write(33,pfmt) '#'
            write(33,*) '#Energia cinetica'
            write(33,pfmt) '#',sum(e_cin)/real(t_run,pr), ' +-', sqrt(de_cin/real(t_run,pr))
            write(33,pfmt) '#'
            write(33,*) '#Energia potencial'
            write(33,pfmt) '#',sum(e_pot)/real(t_run,pr), ' +-', sqrt(de_pot/real(t_run,pr))
            write(33,pfmt) '#'
            write(33,*) '#Presion'
            write(33,pfmt) '#',sum(press)/real(t_run,pr), ' +-', sqrt(dpress/real(t_run,pr))
            write(33,pfmt) '#'
            write(33,*) '#Temperatura'
            write(33,pfmt) '#',sum(temp)/real(t_run,pr), ' +-', sqrt(dtemp/real(t_run,pr))
            write(33,*) '#Tail correc. presion'
            write(33,pfmt) '#',16._pr*pi*ro*ro*(2._pr*((1._pr/rcut)**9)/3._pr-(1._pr/rcut)**3)/3._pr
            write(33,*) '#Tail correc. pot'
            write(33,pfmt) '#',8._pr*pi*ro*(((1._pr/rcut)**9)/3._pr-(1._pr/rcut)**3)/3._pr
        !·············

    endif

    if (.not.base_variable) then
        deallocate(r,v,f,req,veq)
    endif

    close(33)
    close(34)
    close(35)
    

end program p1g5