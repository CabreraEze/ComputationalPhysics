module molecular_dyn
    use precision, only: pr=>qp
    use mzranmod, only: mzran=>rmzran
    
contains

    !#########################################
    !       FUNCIONES
    !#########################################

        !potencial de leonard jones
        real(pr) function lju(r2)
            implicit none
            real(pr), intent(in)    :: r2
            real(pr)                :: r
            
            r = 1._pr/r2
            lju = 4._pr*(r**6 - r**3)

        end function lju



        !vector de fuerza para la particula i
        function force(dr, r2)
            implicit none
            real(pr), intent(in)    :: dr(3) ,r2
            real(pr)                :: force(3), r

            r = 1._pr/r2
            force(:) = 48._pr*dr(:)*r*((r**6) - 0.5_pr*(r**3))

        end function force



        !distancia entre dos particulas, considerando pbc
        reaL(pr) function distance2_pbc(dr, l)
            implicit none
            real(pr), intent(in)       :: l
            real(pr), intent(inout)    :: dr(3)

            dr(1) = dr(1) - anint((dr(1)/l), pr)*l
            dr(2) = dr(2) - anint((dr(2)/l), pr)*l
            dr(3) = dr(3) - anint((dr(3)/l), pr)*l

            distance2_pbc = sum(dr*dr)

        end function distance2_pbc


        
        real(pr) function gasdev()
            implicit none
            real(pr)        :: rsq,v1,v2
            real(pr), save  :: g
            logical, save   :: gaus_stored=.false.
            if (gaus_stored) then
                gasdev=g
                gaus_stored=.false.
            else
                do
                    v1=2.0_pr*mzran()-1.0_pr
                    v2=2.0_pr*mzran()-1.0_pr
                    rsq=v1**2+v2**2
                    if (rsq > 0.0_pr .and. rsq < 1.0_pr) exit
                end do
                rsq=sqrt(-2.0_pr*log(rsq)/rsq)
                gasdev=v1*rsq
                g=v2*rsq
                gaus_stored=.true.
            end if
        end function gasdev

    !#########################################
    !       OBSERVABLES
    !#########################################

        !energia cinetica
        real(pr) function kinetic(v)
            implicit none
            real(pr), intent(in)    :: v(:,:)
            integer                 :: n, i

            n = size(v,dim=2)
            kinetic = 0._pr
            do i=1,n
                kinetic = kinetic + sum(v(:,i)*v(:,i))
            enddo 
            kinetic = 0.5_pr*kinetic

        end function kinetic



        !temperatura
        real(pr) function temperature(v)
            implicit none
            real(pr), intent(in)    :: v(:,:)
            real(pr)                :: v2
            integer                 :: n, i

            n = size(v,dim=2)
            temperature = 0._pr
            do i=1,n
                v2 = v(1,i)*v(1,i)
                temperature = temperature + v2
            enddo
            temperature = temperature/real(3*n,pr)

        end function temperature



        subroutine correlation(g ,r ,dr ,l ,ro)
            implicit none
            real(pr), intent(in)    :: r(:,:), dr, l ,ro
            real(pr), intent(inout) :: g(:)
            real(pr)                :: drij(3) ,rij ,rlw ,rup ,gfac
            real(pr), parameter     :: pi = 4._pr*atan(1._pr)

            integer                 :: n, m, i, j, ir, index, ig(size(g))

            n = size(r,dim=2)
            m = size(g)

            ig = 0
            do i=1,n-1
                do j=i+1,n
                    drij = r(:,i) - r(:,j)
                    rij = sqrt(distance2_pbc(drij, l))

                    ir = int(rij/dr) + 1
                    !index = iabs(ir)
                    if (ir<=m) then
                        ig(ir) = ig(ir) + 1
                    endif
                enddo
            enddo

            do i=1,m
                rlw = real(i-1,pr)*dr
                rup = rlw + dr

                gfac = real(n,pr)*pi*4._pr*ro*(rup**3 - rlw**3)/3._pr
                gfac = 2._pr/gfac
                gfac = gfac/real(m,pr)

                g(i) = real(ig(i),pr)*gfac
            enddo
        end subroutine correlation



    !#########################################
    !       SUBRUTINAS
    !#########################################



        !subrutina para la evolucion de la dinamica del sistema a traves de Velocity-Verlet
        !manteniendo en 0 la posicion y velocidad del centro de masa
        subroutine dynamic_eq(r, v, f, rcut, l, dt)
            implicit none
            real(pr), intent(inout) :: r(:,:), v(:,:)
            real(pr), intent(in)    :: dt, rcut, l,  f(:,:)
            real(pr), allocatable   :: rnew(:,:), df(:,:), dv(:), dr(:)

            integer                 :: i, n

            n = size(r,dim=2)
            allocate(rnew(3,n), df(3,n), dv(3), dr(3))
            rnew = 0._pr
            df = 0._pr
            dv = 0._pr
            dr = 0._pr

            do i=1,n
                rnew(:,i) = r(:,i) + v(:,i)*dt + f(:,i)*0.5_pr*dt**2
                dr(:) = dr(:) + rnew(:,i)
            enddo

            ! move the center of mass
                dr(:) = dr(:)/real(n, pr) 
                do i=1,n
                    rnew(:,i) = rnew(:,i) - dr(:)
                enddo
            !------------------------
            r = rnew

            call pbc (rnew, l)
            call forces(0, rnew ,df ,rcut, l)

            do i=1,n
                v(:,i) = v(:,i) + (df(:,i) + f(:,i))*dt*0.5_pr
                dv(:) = dv(:) + v(:,i)
            enddo

            ! move the velocity of the center of mass
                dv(:) = dv(:)/real(n, pr) 
                do i=1,n
                    v(:,i) = v(:,i) - dv(:)
                enddo
            !------------------------

            deallocate(rnew, df, dv, dr)

        end subroutine dynamic_eq



        !subrutina para la evolucion de la dinamica del sistema a traves de Velocity-Verlet
        subroutine dynamic(r, v, f, rcut, l, dt)
            implicit none
            real(pr), intent(inout) :: r(:,:), v(:,:)
            real(pr), intent(in)    :: dt, rcut, l,  f(:,:)
            real(pr), allocatable   :: rnew(:,:), df(:,:)

            integer                 :: i, n

            n = size(r,dim=2)
            allocate(rnew(3,n), df(3,n))
            rnew = 0._pr
            df = 0._pr

            do i=1,n
                rnew(:,i) = r(:,i) + v(:,i)*dt + f(:,i)*dt*dt*0.5_pr
            enddo
            r = rnew

            call pbc (rnew, l)
            call forces(0, rnew ,df ,rcut, l)

            do i=1,n
                v(:,i) = v(:,i) + (df(:,i) + f(:,i))*dt*0.5_pr
            enddo

            deallocate(rnew, df)

        end subroutine dynamic



        !organiza el sistema en FCC para las condiciones iniciales
        !y proporciona velocidades al azar para cada particula
        subroutine initialize(r, v, ro, base, l, a)
            implicit none
            real(pr), intent(in)                     :: ro
            real(pr), intent(inout)                  :: l, a
            real(pr), allocatable, intent(inout)     :: r(:,:), v(:,:)
            real(pr)                                 :: r0x, r0y, r0z

            !integer, optional, allocatable, intent(inout)      :: head(:), list(:)
            integer, intent(in)                                :: base
            integer                                            :: n, cell, i ,j ,k, cx, cy ,cz, pos_c

            n = 4*base**3
            l = real(n,pr)/ro
            l = l**(1._pr/3._pr)
            a = l/real(base,pr)

            allocate(r(3,n), v(3,n))
            
            r(:,:) = 0._pr
            v(:,:) = 0._pr

            ! if (present(head).and.present(list)) then
            !     allocate(head(base**3), list(n))
            !     head(:) = 0
            !     list(:) = 0
            ! endif

            i = 0
            j = 0
            k = 0

            r0x = -l/2._pr
            r0y = -l/2._pr
            r0z = -l/2._pr

            do cell = 1,base**3
                cx = mod(cell-1,base) + 1
                if (cx==1) then
                    j = j + 1
                    cy = mod(j-1,base) + 1
                    if (cy==1) then
                        k = k + 1
                        cz = mod(k-1,base) + 1
                    endif
                endif
 
                pos_c = ((cz-1)*base*base + (cy-1)*base + (cx-1))*4

                r(1,pos_c+1) = (cx-1)*a
                r(2,pos_c+1) = (cy-1)*a
                r(3,pos_c+1) = (cz-1)*a

                r(1,pos_c+2) = (cx-1)*a + a/2._pr
                r(2,pos_c+2) = (cy-1)*a + a/2._pr
                r(3,pos_c+2) = (cz-1)*a

                r(1,pos_c+3) = (cx-1)*a + a/2._pr
                r(2,pos_c+3) = (cy-1)*a
                r(3,pos_c+3) = (cz-1)*a + a/2._pr

                r(1,pos_c+4) = (cx-1)*a
                r(2,pos_c+4) = (cy-1)*a + a/2._pr
                r(3,pos_c+4) = (cz-1)*a + a/2._pr

                ! if (present(head).and.present(list)) then
                !     head(cell) = pos_c + 1
                !     list(pos_c+1) = pos_c + 2
                !     list(pos_c+2) = pos_c + 3
                !     list(pos_c+3) = pos_c + 4
                ! endif
            enddo

            r(1,:) = r(1,:) + r0x
            r(2,:) = r(2,:) + r0y
            r(3,:) = r(3,:) + r0z

            do i=1,n
                do j=1,3
                    v(j,i) = (mzran()-0.5_pr)
                enddo
            enddo

        end subroutine initialize



        !subrutina encargada de mantener las condiciones periodicas de contorno
        subroutine pbc(r, l)
            implicit none
            real(pr), intent(inout)    :: r(:,:)
            real(pr), intent(in)       :: l
            
            integer                    :: i, n

            n = size(r,dim=2)
            do i=1,n
                r(1,i) = r(1,i) - anint((r(1,i)/l), pr)*l
                r(2,i) = r(2,i) - anint((r(2,i)/l), pr)*l
                r(3,i) = r(3,i) - anint((r(3,i)/l), pr)*l
            enddo

        end subroutine pbc



        !subrutina para determinar las fuerzas que experimenta cada particula, teniendo
        !en cuenta solo primeros vecinos

        !obs: 0 si no es necesario calcular los observables, 1 si es necesario
        !r: matriz con las posiciones de todas las particulas
        !rcut: distancia de corte para primeros vecinos
        !l: largo de la caja
        !e_pot: (OPCIONAL) observable de la energia potencial
        !press: (OPCIONAL) observable de la presion
        subroutine forces(obs, r, f, rcut, l, e_pot, press)
            implicit none
            real(pr), intent(inout)                  :: f(:,:)
            real(pr), optional, intent(out)          :: e_pot, press
            real(pr), intent(in)                     :: r(:,:), rcut, l
            real(pr)                                 :: df(3), dr(3), r2, obs_press, obs_pot, tail

            integer, intent(in)     :: obs
            integer                 :: n, i, j, k

            n = size(r,dim=2)
            f(:,:) = 0._pr
            k = 0

            obs_press = 0._pr
            obs_pot = 0._pr

            do i = 1,n-1
                do j = i+1,n
                    dr(:) = r(:,i) - r(:,j)
                    r2 = distance2_pbc(dr(:),l)

                    !check para primeros vecinos
                    if (r2<=(rcut*rcut)) then
                        df = force(dr,r2)
                        f(:,i) = f(:,i) + df
                        f(:,j) = f(:,j) - df
                        
                        if (obs==1) then
                            k = k + 1
                            obs_pot = obs_pot + lju(r2)
                            obs_press = obs_press + sum(f(:,i)*dr(:))
                        endif
                    endif
                enddo
            enddo

            if (present(e_pot).and.present(press)) then
                tail = lju(rcut*rcut)

                press = obs_press/real(k,pr)
                e_pot = obs_pot - real(k,pr)*tail
            endif

        end subroutine forces

    !#########################################
    !       BROWNIAN MOTION
    !#########################################
        
        subroutine BM(r, f, nu ,D0 ,dt)
            implicit none
            real(pr), intent(in)        :: f(:,:), nu, D0, dt
            real(pr), intent(inout)     :: r(:,:)
            real(pr)                    :: dr(3)
            integer                     :: i, n
            real(pr), parameter         :: pi = 4._pr*atan(1._pr)
            
            n = size(r,dim=2)
            do i=1,n
                dr(1) = gasdev()
                dr(2) = gasdev()
                dr(3) = gasdev()
                r(:,i) = r(:,i) + dt*f(:,i)/(3._pr*pi*nu) + sqrt(2._pr*D0*dt)*dr(:)
            enddo
        
        end subroutine BM
    
        subroutine FCC(r, ro, base, l, a)
            implicit none
            real(pr), intent(in)                     :: ro
            real(pr), intent(inout)                  :: l, a
            real(pr), allocatable, intent(inout)     :: r(:,:)
            real(pr)                                 :: r0x, r0y, r0z
        
            integer, intent(in)                                :: base
            integer                                            :: n, cell, i ,j ,k, cx, cy ,cz, pos_c
        
            n = 4*base**3
            l = real(n,pr)/ro
            l = l**(1._pr/3._pr)
            a = l/real(base,pr)
        
            allocate(r(3,n))
            
            r(:,:) = 0._pr
        
            i = 0
            j = 0
            k = 0
        
            r0x = -l/2._pr
            r0y = -l/2._pr
            r0z = -l/2._pr
        
            do cell = 1,base**3
                cx = mod(cell-1,base) + 1
                if (cx==1) then
                    j = j + 1
                    cy = mod(j-1,base) + 1
                    if (cy==1) then
                        k = k + 1
                        cz = mod(k-1,base) + 1
                    endif
                endif
        
                pos_c = ((cz-1)*base*base + (cy-1)*base + (cx-1))*4
        
                r(1,pos_c+1) = (cx-1)*a
                r(2,pos_c+1) = (cy-1)*a
                r(3,pos_c+1) = (cz-1)*a
        
                r(1,pos_c+2) = (cx-1)*a + a/2._pr
                r(2,pos_c+2) = (cy-1)*a + a/2._pr
                r(3,pos_c+2) = (cz-1)*a
        
                r(1,pos_c+3) = (cx-1)*a + a/2._pr
                r(2,pos_c+3) = (cy-1)*a
                r(3,pos_c+3) = (cz-1)*a + a/2._pr
        
                r(1,pos_c+4) = (cx-1)*a
                r(2,pos_c+4) = (cy-1)*a + a/2._pr
                r(3,pos_c+4) = (cz-1)*a + a/2._pr
        
            enddo
        
            r(1,:) = r(1,:) + r0x
            r(2,:) = r(2,:) + r0y
            r(3,:) = r(3,:) + r0z
        
        end subroutine FCC
        

end module molecular_dyn