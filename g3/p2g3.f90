!###################################
!#  El programa realiza una        #
!#  caminata al azar en 2 y 3      #
!#  dimensiones utilizando 3       #
!#  rng distintos, es              #
!#  necesario tener un             #
!#  directorio con el nombre       #
!#  'p2' donde el programa         #
!#  almacene los datos generados   #
!#  para su posterior analisis.    #
!###################################

program p2g3
    use precision, only: pr=>dp
    use mtmod, only: mtran=>grnd !, mtseed=>sgrnd
    use mzranmod, only: mzran=>rmzran
    use ranmod
    implicit none

    character(len=50), parameter    :: loc='p2/', fmt='(i5,10(e20.12,2x))'

    real(pr)                :: t_r2, t_mt, t_mz, x(3), y(3), z(3), despl(3)
    real(pr), parameter     :: r=1._pr, x0=0._pr, y0=0._pr, z0=0._pr, pi=4._pr*atan(1._pr) 

    integer                :: r0seed, r2seed, i, cuad(3), k, nstep
    integer, parameter     :: nmax=10**6, initial_seed=2516 !Seed inicial para todo el programa

    !##########################
    !           2D
    !##########################

    open(33,file=trim(adjustl(loc))//'caminata_sample.d',status='replace')
    open(34,file=trim(adjustl(loc))//'promediaciones.d',status='replace')
        write(34,'(a100)') '#promedios de ingreso al 1° cuadrante y desplazamiento cuadratico en funcion de n'
        write(34,'(10(a20,5x))') '#nstep', 'ran2 cuad1', 'mtran cuad1', 'mzran cuad1', &
                                'ran2 despl medio', 'mtran despl medio', 'mzran despl medio'

        write(33,'(a50)') '#ejemplo de caminata con n=500'
        write(33,'(5(a20,5x))') '#nstep', 'ran2 x-y', 'mtran x-y', 'mzran x-y'
        write(33,fmt) 500, x0 ,y0, x0 ,y0, x0 ,y0

        do nstep=50,1000,50
            r0seed=initial_seed
            r2seed=initial_seed
            !call mtseed(initial_seed)

            despl(:) = 0._pr
            cuad(:) = 0
            do k=1,nmax
                x(:)=x0
                y(:)=y0

                do i=1,nstep



                    !###    RAN2    ###
                    t_r2=ran2(r2seed)
                    if (t_r2<=0.25_pr) then
                        x(1) = x(1) + r
                    else if ((t_r2>0.25_pr).and.(t_r2<=0.5_pr)) then
                        y(1) = y(1) + r
                    else if ((t_r2>0.5_pr).and.(t_r2<=0.75_pr)) then
                        x(1) = x(1) - r
                    else if (t_r2>=0.75_pr) then
                        y(1) = y(1) - r
                    endif

                    if ((x(1)>=0._pr).and.(y(1)>=0._pr).and.(i==nstep)) then
                        cuad(1) = cuad(1) + 1
                    endif
                    !###################



                    !###    MTRAN    ###
                    t_mt=mtran()
                    if (t_mt<=0.25_pr) then
                        x(2) = x(2) + r
                    else if ((t_mt>0.25_pr).and.(t_mt<=0.5_pr)) then
                        y(2) = y(2) + r
                    else if ((t_mt>0.5_pr).and.(t_mt<=0.75_pr)) then
                        x(2) = x(2) - r
                    else if (t_mt>=0.75_pr) then
                        y(2) = y(2) - r
                    endif

                    if ((x(2)>=0._pr).and.(y(2)>=0._pr).and.(i==nstep)) then
                        cuad(2) = cuad(2) + 1
                    endif
                    !###################



                    !###    MZRAN    ###
                    t_mz=mzran()
                    if (t_mz<=0.25_pr) then
                        x(3) = x(3) + r
                    else if ((t_mz>0.25_pr).and.(t_mz<=0.5_pr)) then
                        y(3) = y(3) + r
                    else if ((t_mz>0.5_pr).and.(t_mz<=0.75_pr)) then
                        x(3) = x(3) - r
                    else if (t_mz>=0.75_pr) then
                        y(3) = y(3) - r
                    endif

                    if ((x(3)>=0._pr).and.(y(3)>=0._pr).and.(i==nstep)) then
                        cuad(3) = cuad(3) + 1
                    endif
                    !###################



                    if ((k==1).and.(nstep==500)) then
                        write(33,fmt) nstep, x(1), y(1), x(2), y(2), x(3), y(3)
                    endif
                enddo
                despl(:) = despl(:) + x(:)*x(:) + y(:)*y(:)
            enddo
            write(34,fmt) nstep, real(cuad(1),pr)/real(nmax,pr), real(cuad(2),pr)/real(nmax,pr), &
                        real(cuad(3),pr)/real(nmax,pr), despl(1)/real(nmax,pr), despl(2)/real(nmax,pr), despl(3)/real(nmax,pr) 
        enddo
    close(34)
    close(33)

    !##########################
    !           3D
    !##########################

    open(35,file=trim(adjustl(loc))//'caminata_sample3D.d',status='replace')
    open(36,file=trim(adjustl(loc))//'promediaciones3D.d',status='replace')
        write(36,'(a100)') '#promedios de ingreso al 1° cuadrante y desplazamiento cuadratico en funcion de n'
        write(36,'(10(a20,5x))') '#nstep', 'ran2 cuad1', 'mtran cuad1', 'mzran cuad1', &
                                'ran2 despl medio', 'mtran despl medio', 'mzran despl medio'

        write(35,'(a50)') '#ejemplo de caminata con n=500'
        write(35,'(5(a20,5x))') '#nstep', 'ran2 x-y-z', 'mtran x-y-z', 'mzran x-y-z'
        write(35,fmt) 500, x0 ,y0 ,z0, x0 ,y0 ,z0 ,x0 ,y0 ,z0

        do nstep=50,1000,50
            r0seed=initial_seed
            r2seed=initial_seed
            !call mtseed(initial_seed)

            despl(:) = 0._pr
            cuad(:) = 0
            do k=1,nmax
                x(:)=x0
                y(:)=y0
                z(:)=z0
                
                do i=1,nstep



                    !###    RAN2    ###
                    t_r2=ran2(r2seed)
                    if (t_r2<=1._pr/6._pr) then
                        x(1) = x(1) + r
                    else if ((t_r2>1._pr/6._pr).and.(t_r2<=1._pr/3._pr)) then
                        y(1) = y(1) + r
                    else if ((t_r2>1._pr/3._pr).and.(t_r2<=0.5_pr)) then
                        x(1) = x(1) - r
                    else if ((t_r2>0.5_pr).and.(t_r2<=2._pr/3._pr)) then
                        y(1) = y(1) - r
                    else if ((t_r2>2._pr/3._pr).and.(t_r2<=5._pr/6._pr)) then
                        z(1) = z(1) + r
                    else if ((t_r2>5._pr/6._pr).and.(t_r2<=1._pr)) then
                        z(1) = z(1) - r
                    endif

                    if ((x(1)>=0._pr).and.(y(1)>=0._pr).and.(z(1)>=0._pr).and.(i==nstep)) then
                        cuad(1) = cuad(1) + 1
                    endif
                    !###################



                    !###    MTRAN    ###
                    t_mt=mtran()
                    if (t_mt<=1._pr/6._pr) then
                        x(2) = x(2) + r
                    else if ((t_mt>1._pr/6._pr).and.(t_mt<=1._pr/3._pr)) then
                        y(2) = y(2) + r
                    else if ((t_mt>1._pr/3._pr).and.(t_mt<=0.5_pr)) then
                        x(2) = x(2) - r
                    else if ((t_mt>0.5_pr).and.(t_mt<=2._pr/3._pr)) then
                        y(2) = y(2) - r
                    else if ((t_mt>2._pr/3._pr).and.(t_mt<=5._pr/6._pr)) then
                        z(2) = z(2) + r
                    else if ((t_mt>5._pr/6._pr).and.(t_mt<=1._pr)) then
                        z(2) = z(2) - r
                    endif

                    if ((x(2)>=0._pr).and.(y(2)>=0._pr).and.(z(2)>=0._pr).and.(i==nstep)) then
                        cuad(2) = cuad(2) + 1
                    endif
                    !###################



                    !###    MZRAN    ###
                    t_mz=mzran()
                    if (t_mz<=1._pr/6._pr) then
                        x(3) = x(3) + r
                    else if ((t_mz>1._pr/6._pr).and.(t_mz<=1._pr/3._pr)) then
                        y(3) = y(3) + r
                    else if ((t_mz>1._pr/3._pr).and.(t_mz<=0.5_pr)) then
                        x(3) = x(3) - r
                    else if ((t_mz>0.5_pr).and.(t_mz<=2._pr/3._pr)) then
                        y(3) = y(3) - r
                    else if ((t_mz>2._pr/3._pr).and.(t_mz<=5._pr/6._pr)) then
                        z(3) = z(3) + r
                    else if ((t_mz>5._pr/6._pr).and.(t_mz<=1._pr)) then
                        z(3) = z(3) - r
                    endif

                    if ((x(3)>=0._pr).and.(y(3)>=0._pr).and.(z(3)>=0._pr).and.(i==nstep)) then
                        cuad(3) = cuad(3) + 1
                    endif
                    !###################



                    if ((k==1).and.(nstep==500)) then
                        write(35,fmt) nstep, x(1), y(1), z(1), x(2), y(2), z(2), x(3), y(3), z(3)
                    endif

                enddo
                despl(:) = despl(:) + x(:)*x(:) + y(:)*y(:) + z(:)*z(:)
            enddo
            write(36,fmt) nstep, real(cuad(1),pr)/real(nmax,pr), real(cuad(2),pr)/real(nmax,pr), &
                        real(cuad(3),pr)/real(nmax,pr), despl(1)/real(nmax,pr), despl(2)/real(nmax,pr), despl(3)/real(nmax,pr) 
        enddo
    close(36)
    close(35)

end program p2g3