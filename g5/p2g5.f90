program p2g5
    use precision, only: pr=>qp
    use molecular_dyn
    implicit none
    
    real(pr), allocatable        :: r1(:,:), v1(:,:), f1(:,:), g1(:), g1tot(:), r1d(:,:), pd1(:,:)
    real(pr), allocatable        :: r2(:,:), v2(:,:), f2(:,:), g2(:), g2tot(:), r2d(:,:), pd2(:,:)
    real(pr)                     :: dt, rcut, t, t0, lambda, temp_eq, rmx, dr
    real(pr)                     :: l1 ,a1 ,ro1 ,s1(2), k1(3), mds1
    real(pr)                     :: l2 ,a2 ,ro2 ,s2(2), k2(3), mds2
    real(pr), parameter          :: pi = 4._pr*atan(1._pr)

    integer                      :: base, n, m, i, t_iter, t_run, t_eq

    character(len=50), parameter :: fmt='(10(e18.10,2x))', loc='p2/', pfmt='(2(a3,e18.10))'

    !#########################################
    !       INITIAL VALUE
    !#########################################

    !ro: densidad del sistema
    !base: cantidad de celdas unitarias que abarca un lado de la caja de simulacion
    !rcut: distancia de corte para primeros vecinos
    !dt: intervalo de tiempo
    !t0: temperatura inicial
    !dr: distancia para la correlacion de pares
    !rmx: distancia maxima a evaluar la correlacion de pares

    !t_eq: numero de pasos de simulacion para equilibrar el sistema
    !t_scal: cada t_scal durante la equilibracion se reescalaran las velocidades para mantener la temperatura constante
    !trun: pasos de corrida de la simulacion

    ro1 = 0.8_pr
    ro2 = 1.2_pr

    base = 5
    rcut = 2.5_pr 
    dt = 0.005_pr
    t0 = 0.01_pr
    rmx = 3.5_pr
    dr = rmx/100._pr

    t_eq = 500
    t_run = 1000

    !#########################################
    
    
    m = int(rmx/dr)
    allocate(g1(m), g2(m))
    allocate(g1tot(m), g2tot(m))
    g1tot = 0._pr
    g2tot = 0._pr

    l1 = 0._pr
    a1 = 0._pr
    call initialize(r1, v1, ro1 ,base ,l1 ,a1)
    k1(1) = -2._pr*pi/a1
    k1(2) = 2._pr*pi/a1
    k1(3) = -2._pr*pi/a1

    l2 = 0._pr
    a2 = 0._pr
    call initialize(r2, v2, ro2 ,base ,l2 ,a2)
    k2(1) = -2._pr*pi/a2
    k2(2) = 2._pr*pi/a2
    k2(3) = -2._pr*pi/a2

    n=size(r1,dim=2)
    allocate(f1(3,n),f2(3,n),r1d(3,n),r2d(3,n))
    allocate(pd1(3,n),pd2(3,n))
    pd1(:,:) = 0._pr
    pd2(:,:) = 0._pr
    f1(:,:) = 0._pr
    f2(:,:) = 0._pr
    temp_eq = 0._pr

    !termalizacion del sistema
    do t_iter = 1,t_eq

        !···DINAMICA 1···
            call forces(0, r1, f1, rcut, l1)
            call dynamic_eq(r1, v1, f1, rcut, l1, dt)
            call pbc(r1, l1)
            lambda = sqrt(t0/temperature(v1))
            v1 = v1*lambda

        !···DINAMICA 2···
            call forces(0, r2, f2, rcut, l2)
            call dynamic_eq(r2, v2, f2, rcut, l2, dt)
            call pbc(r2, l2)
            lambda = sqrt(t0/temperature(v2))
            v2 = v2*lambda

    enddo
    !·························

    write(*,*) "starting simulation..."

    open(33,file=trim(adjustl(loc))//'temp.d',status='replace')
    write(33,'(8(a20,2x))') '#t', 'S(k,t)_1', '<dr(t)**2>(t)_1', 'S(k,t)_2', '<dr(t)**2>(t)_2'
    open(34,file=trim(adjustl(loc))//'corre.d',status='replace')
    write(34,'(8(a20,2x))') '#r', 'g(r)_1', 'g(r)_2'

        !corrida de simulacion
        do t_iter = 1,t_run
            s1(:) = 0._pr
            s2(:) = 0._pr
            mds1 = 0._pr
            mds2 = 0._pr
            t = dt*real(t_iter,pr)

            !···CONDITIONS 1···
                r1d(:,:) = r1(:,:)
                call forces(0, r1, f1, rcut, l1)
                call dynamic_eq(r1, v1, f1, rcut, l1, dt)
                r1d(:,:) = r1(:,:) - r1d(:,:)
                temp_eq = temperature(v1)
                lambda = sqrt(t0/temp_eq)
                v1 = v1*lambda
                call pbc(r1, l1)
                call correlation(g1, r1, dr, l1, ro1)
                g1tot = g1tot + g1

            !···CONDITIONS 2···
                r2d(:,:) = r2(:,:)
                call forces(0, r2, f2, rcut, l2)
                call dynamic_eq(r2, v2, f2, rcut, l2, dt)
                r2d(:,:) = r2(:,:) - r2d(:,:)
                temp_eq = temperature(v2)
                lambda = sqrt(t0/temp_eq)
                v2 = v2*lambda
                call pbc(r2, l2)
                call correlation(g2, r2, dr, l2, ro2)
                g2tot = g2tot + g2

            pd1(:,:) = pd1(:,:) + r1d(:,:)
            pd2(:,:) = pd2(:,:) + r2d(:,:)
            do i=1,n
                s1(1) = s1(1) + cos(sum(k1(:)*r1(:,i)))
                s1(2) = s1(2) + sin(sum(k1(:)*r1(:,i)))

                s2(1) = s2(1) + cos(sum(k2(:)*r2(:,i)))
                s2(2) = s2(2) + sin(sum(k2(:)*r2(:,i)))

                mds1 =  mds1 + sum(pd1(:,i)*pd1(:,i))
                mds2 =  mds2 + sum(pd2(:,i)*pd2(:,i))
            enddo

            write(33,fmt) t ,sum(s1*s1)/real(n*n,pr) ,mds1/real(n,pr) ,sum(s2*s2)/real(n*n,pr) ,mds2/real(n,pr)

        enddo
        g1tot = g1tot/real(t_run,pr)
        g2tot = g2tot/real(t_run,pr)

        do i=1,m
            write(34,fmt) (real(i-1,pr) + 0.5)*dr, g1tot(i), g2tot(i)
        enddo   

    close(33)
    close(34)


    deallocate(g1,g2,g1tot,g2tot)
    deallocate(r1,r2,v1,v2,f1,f2)

end program p2g5