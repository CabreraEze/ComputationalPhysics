program final
    use precision, only: pr=>qp
    use molecular_dyn
    implicit none

    real(pr), allocatable       :: r(:,:), f(:,:), dr(:,:), g(:), gtot(:), pd(:,:)
    real(pr)                    :: rcut, ro, nu, D0, T0, l, a, dg, rmx, t, dt, mds, dt_scale
    real(pr), parameter         :: pi = 4._pr*atan(1._pr)

    integer                     :: base, t_iter, t_eq, t_run, m, n, i, j

    character(len=50), parameter:: fmt='(10(e18.10,2x))', loc='data/'

    logical                     :: ro_variable, evo

    !#########################################
    !       INITIAL VALUE
    !######################################### 
    
    base = 4
    rcut = 2.5_pr
    ro   = 0.8_pr
    T0   = 1.1_pr
    nu   = 2.87_pr

    t_eq   = 1
    t_run = 1000
    dt_scale = 10e-5_pr

    rmx = 3.5_pr
    dg = rmx/100._pr


    D0 = T0/(3._pr*pi*nu)

    !....OPTIONAL....
    ro_variable = .true.
    evo = .false.

    !######################################### 

    open(33,file=trim(adjustl(loc))//'mediciones.d',status='replace')
    if (evo) then
        open(22,file=trim(adjustl(loc))//'evolution.d',status='replace')
    endif

    if (ro_variable) then 

        do j=0,15
            ro = 0.2_pr + 0.85_pr*real(j,pr)/15._pr
            write(*,*) 'ro simulado:', ro

            write(33,'(a13,e11.3,a2)') '"ro=', ro, '"'

            l = 0._pr
            a = 0._pr
            call FCC(r, ro, base, l ,a)
            n=size(r,dim=2)
        
            dt = (rcut/2.5_pr)**2/D0
            dt = dt*dt_scale
        
            write(*,*) 'dt determinado:', dt
        
            m = int(rmx/dg)
            allocate(f(3,n), dr(3,n), pd(3,n))
            f(:,:) = 0._pr
            dr(:,:)= 0._pr
        
        
            !termalizacion del sistema
            do t_iter=1,t_eq
                t = dt*real(t_iter,pr)
        
                f(:,:) = 0._pr
                call forces(0, r, f, rcut, l)
                call BM(r, f, nu, D0, dt)
                call pbc(r, l)
        
            enddo
            !·························
   
            !corrida de simulacion
            do t_iter=1,t_run
                mds = 0._pr
                t = dt*real(t_iter,pr)
        
                f(:,:) = 0._pr
                dr(:,:) = r(:,:)
                call forces(0, r, f, rcut, l)
                call BM(r, f, nu, D0, dt)
                dr(:,:) = r(:,:) - dr(:,:)
        
                call pbc(r, l)
        
                pd(:,:) = pd(:,:) + dr(:,:)
                do i=1,n
                    mds =  mds + sum(pd(:,i)*pd(:,i))
                enddo

                write(33,fmt) t, mds/real(n,pr)
                
            enddo
            !·························

            write(33,*) ''
            write(33,*) ''

            deallocate(f, dr, pd, r)
        enddo


    else

        l = 0._pr
        a = 0._pr
        call FCC(r, ro, base, l ,a)
        n=size(r,dim=2)

        dt = (rcut/2.5_pr)**2/D0
        dt = dt*dt_scale

        write(*,*) 'dt determinado:', dt

        m = int(rmx/dg)
        allocate(f(3,n), dr(3,n), pd(3,n), g(m), gtot(m))
        f(:,:) = 0._pr
        dr(:,:)= 0._pr


        !termalizacion del sistema
        do t_iter=1,t_eq
            t = dt*real(t_iter,pr)

            f(:,:) = 0._pr
            call forces(0, r, f, rcut, l)
            call BM(r, f, nu, D0, dt)
            call pbc(r, l)

        enddo
        !·························

        write(33,'(8(a20,2x))') '#t', '<dr(t)**2>(t)'
        !corrida de simulacion
        do t_iter=1,t_run
            mds = 0._pr
            t = dt*real(t_iter,pr)

            f(:,:) = 0._pr
            dr(:,:) = r(:,:)
            call forces(0, r, f, rcut, l)
            call BM(r, f, nu, D0, dt)
            dr(:,:) = r(:,:) - dr(:,:)

            call pbc(r, l)
            call correlation(g, r, dg, l, ro)
            gtot = gtot + g

            pd(:,:) = pd(:,:) + dr(:,:)
            do i=1,n
                mds =  mds + sum(pd(:,i)*pd(:,i))
            enddo

            write(33,fmt) t, mds/real(n,pr)
            if (evo) then
                do i = 1,n
                    write(22,fmt) r(1,i), r(2,i), r(3,i)
                enddo
                write(22,*) ''
                write(22,*) ''
            endif
        enddo
        !·························
        
        write(33,*) ''
        write(33,*) ''
        write(33,'(8(a20,2x))') '#r', 'g(r)'
        
        gtot = gtot/real(t_run,pr)
        do i=1,m
            write(33,fmt) (real(i-1,pr) + 0.5)*dg, gtot(i)
        enddo

        deallocate(f, dr, pd, g, gtot)
    endif

    close(33)
    if (evo) then
        close(22)
    endif
end program final