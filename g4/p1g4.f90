!######################################
!#  Es necesario tener un directorio  #
!#  con el nombre 'p1' donde el       #
!#  programa almacene los datos       #
!#  generados para su posterior       #
!#  analisis.                         #
!######################################

program p1g4
    use precision, only: pr=>dp
    use mzranmod, only: mzran=>rmzran
    use ising2D
    implicit none

    character(len=50), parameter    :: loc='p1/', fmt1='(i5,2x,10(e20.12,2x))', fmt2='(10(e20.12,2x))'

    integer, allocatable            :: s1(:,:), s2(:,:), s1mc(:,:), s2mc(:,:)
    integer                         :: n, i, j, k, l, nmc, mcs, mct, n1, n2, nt, ndt

    real(pr)                        :: mag(2), mag_abs(2), ene(2)
    real(pr)                        :: t, dt, tinf, tsup, ene2, mag2, mag4, temp1, temp2
    real(pr)                        :: mce, mcmag, mce2, mcmag2, mcmag4
 
    !#######################################################
    !           seteando red a T=0 y T->inf
    !#######################################################

    n=10    
    allocate(s1(n,n), s2(n,n))

    s1(:,:) = 1
    do i=1,n
        do j=1,n
            if (mzran().ge.0.5_pr) then
                s2(i,j) = 1
            else
                s2(i,j) = -1
            endif
        enddo
    enddo
    
    !#######################################################
    !           inciso a)
    !#######################################################

    n1=size(s1,dim=1)
    n2=size(s1,dim=2)
    allocate(s1mc(n1,n2), s2mc(n1,n2))

    open(33,file=trim(adjustl(loc))//'incisoa.d',status='replace')
    write(33,'(8(a20,2x))') '#MCstep', 'energia (T=0)', 'magnetizacion (T=0)', '|magnetizacion| (T=0)', &
                            'energia (T->inf)', 'magnetizacion (T->inf)', '|magnetizacion| (T->inf)'

        nmc = 10000

        t = 2._pr
        s1mc = s1
        s2mc = s2
        do i=1,nmc
            call flip(s1mc,t)
            call flip(s2mc,t)

            write(33,fmt1) i, real(energia(s1mc),pr), magnetizacion(s1mc), magnetizacion_abs(s1mc), &
                        real(energia(s2mc),pr), magnetizacion(s2mc), magnetizacion_abs(s2mc)
        enddo
        write(33,*) ''
        write(33,*) ''

        t = 2.2676_pr
        s1mc = s1
        s2mc = s2
        do i=1,nmc
            call flip(s1mc,t)
            call flip(s2mc,t)

            write(33,fmt1) i, real(energia(s1mc),pr), magnetizacion(s1mc), magnetizacion_abs(s1mc), &
                        real(energia(s2mc),pr), magnetizacion(s2mc), magnetizacion_abs(s2mc)
        enddo
        write(33,*) ''
        write(33,*) ''

        t = 3.3_pr
        s1mc = s1
        s2mc = s2
        do i=1,nmc
            call flip(s1mc,t)
            call flip(s2mc,t)

            write(33,fmt1) i, real(energia(s1mc),pr), magnetizacion(s1mc), magnetizacion_abs(s1mc), &
                        real(energia(s2mc),pr), magnetizacion(s2mc), magnetizacion_abs(s2mc)
        enddo
        write(33,*) ''
        write(33,*) ''
 
    close(33)
    
    deallocate(s1 ,s2)
    deallocate(s1mc ,s2mc)

    !#######################################################
    !           inciso b) c)
    !#######################################################

    ! nmc: numero total de pasos de MC
    ! mcs: largo de un bloque
    ! mct: tiempo de termalizacion
    ! dt: delta de temperatura
    ! nt: numeros de temperaturas a evaluar
    !  
    ! El programa adicionalmente agrega ndt-1 puntos entre las temperaturas
    ! en el intervalo [tinf:tsup]

    nmc = 100000
    mcs = 2000
    mct = 4000
    dt = 0.033_pr
    nt = 100
    tinf = 1.8_pr
    tsup = 2.8_pr
    ndt = 10

    open(34,file=trim(adjustl(loc))//'incisob.d',status='replace')
    write(34,'(8(a25,5x))') '#t', 'energia', 'magnetizacion', 'calor esp.', 'suceptibilidad', 'cumulante'

        do j=0,2
            n=10*2**j

            write(*,*) 'implementando L=', n, '...'

            allocate(s2(n,n))
            do i=1,n
                do k=1,n
                    if (mzran().ge.0.8_pr) then
                        s2(i,k) = 1
                    else
                        s2(i,k) = -1
                    endif
                enddo
            enddo

            allocate(s2mc(n,n))
            s2mc=s2

            do k=1,nt
                t = real(k,pr)*dt

                mag(:) = 0._pr
                ene(:) = 0._pr
                ene2 = 0._pr
                mag2 = 0._pr
                mag4 = 0._pr
                mce = 0._pr
                mce2 = 0._pr
                mcmag = 0._pr
                mcmag2 = 0._pr
                mcmag4 = 0._pr

                if ((t<tsup).and.(t>tinf)) then
                    do l=0,ndt-1
                        t = (real(k,pr)+ real(l,pr)/10._pr)*dt

                        mag(:) = 0._pr
                        ene(:) = 0._pr
                        ene2 = 0._pr
                        mag2 = 0._pr
                        mag4 = 0._pr
                        mce = 0._pr
                        mce2 = 0._pr
                        mcmag = 0._pr
                        mcmag2 = 0._pr
                        mcmag4 = 0._pr

                        do i=1,nmc
                            call flip(s2mc,t)

                            if (mod(i,mcs)==0) then
                                mag(:) = mag(:)/real(mcs,pr)
                                ene(:) = ene(:)/real(mcs,pr)
                                ene2 = ene2/real(mcs,pr)
                                mag2 = mag2/real(mcs,pr)
                                mag4 = mag4/real(mcs,pr)

                                mce = mce + ene(2)
                                mce2 = mce2 + ene2
                                mcmag = mcmag + mag(2)
                                mcmag2 = mcmag2 + mag2
                                mcmag4 = mcmag4 + mag4
                            
                                mag(:) = 0._pr
                                ene(:) = 0._pr
                                ene2 = 0._pr
                                mag2 = 0._pr
                                mag4 = 0._pr
                            endif
                            if (i>mct) then
                                temp1 = magnetizacion_abs(s2mc)
                                temp2 = real(energia(s2mc),pr)

                                mag(2) = mag(2) + temp1
                                ene(2) = ene(2) + temp2

                                mag2 = mag2 + temp1*temp1
                                ene2 = ene2 + temp2*temp2

                                mag4 = mag4 + temp1*temp1*temp1*temp1
                            endif
                        enddo
                        temp1 = real(nmc-mct,pr)/real(mcs,pr) 

                        mce = mce/temp1
                        mce2 = mce2/temp1
                        mcmag = mcmag/temp1
                        mcmag2 = mcmag2/temp1
                        mcmag4 = mcmag4/temp1
                        write(34,fmt2) t, mce, mcmag, (mce2 - mce*mce)/(t*t), (mcmag2 - mcmag*mcmag)/t, &
                                    (1._pr - mcmag4/(3._pr*mcmag2*mcmag2))
                    enddo
                else
                    do i=1,nmc
                        call flip(s2mc,t)

                        if (mod(i,mcs)==0) then
                            mag(:) = mag(:)/real(mcs,pr)
                            ene(:) = ene(:)/real(mcs,pr)
                            ene2 = ene2/real(mcs,pr)
                            mag2 = mag2/real(mcs,pr)
                            mag4 = mag4/real(mcs,pr)

                            mce = mce + ene(2)
                            mce2 = mce2 + ene2
                            mcmag = mcmag + mag(2)
                            mcmag2 = mcmag2 + mag2
                            mcmag4 = mcmag4 + mag4
                        
                            mag(:) = 0._pr
                            ene(:) = 0._pr
                            ene2 = 0._pr
                            mag2 = 0._pr
                            mag4 = 0._pr
                        endif
                        if (i>mct) then
                            temp1 = magnetizacion_abs(s2mc)
                            temp2 = real(energia(s2mc),pr)

                            mag(2) = mag(2) + temp1
                            ene(2) = ene(2) + temp2

                            mag2 = mag2 + temp1*temp1
                            ene2 = ene2 + temp2*temp2

                            mag4 = mag4 + temp1*temp1*temp1*temp1
                        endif
                    enddo
                    temp1 = real(nmc-mct,pr)/real(mcs,pr) 

                    mce = mce/temp1
                    mce2 = mce2/temp1
                    mcmag = mcmag/temp1
                    mcmag2 = mcmag2/temp1
                    mcmag4 = mcmag4/temp1
                    write(34,fmt2) t, mce, mcmag, (mce2 - mce*mce)/(t*t), (mcmag2 - mcmag*mcmag)/t, &
                                (1._pr - mcmag4/(3._pr*mcmag2*mcmag2))
                endif
            enddo

            write(34,*) ''
            write(34,*) ''
            deallocate(s2mc,s2)
        enddo

    close(34)
contains

    

end program p1g4