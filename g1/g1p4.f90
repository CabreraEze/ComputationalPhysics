program g1p4
    use precision !modulo de presicion, dp doble presicion
    use gaussmod  !modulo de metodo de gauss
    implicit none

    integer, parameter             :: k_max=5, n_max=5*2**k_max
    integer                        :: k, n ,i
    real(dp), dimension(7,n_max)   :: x_n, inte, inte_simp
    real(dp), dimension(7)         :: h, suma, suma_simp, suma_gauss
    real(dp), dimension(7)         :: err, err_simp, err_gauss
    real(dp), dimension(n_max)     :: x, w

    x_n = 0._dp
    inte = 0._dp
    suma = 0._dp
    inte_simp = 0._dp
    suma_simp = 0._dp
    err = 0._dp
    err_simp = 0._dp

    suma_gauss = 0._dp
    !seteando matrices en 0
    
    !open(10,file="g1p4_ini.dat")
    !read

    open(11,file="g1p4.dat")

    do k=1,k_max
        n=5*2**k+1
        h(k)=1._dp/n    !alojo todos los h(n)

        do i=1,n
            x_n(k,i) = h(k)*(i-1)   !genero los x_n puntos, pero solo uso hasta n+1

            if (i==1 .or. i==n) then                    !identifico los casos extremos
                inte(k,i) = f(x_n(k,i))*h(k)*0.5_dp
                inte_simp(k,i) = f(x_n(k,i))*h(k)/3._dp
            else
                inte(k,i) = f(x_n(k,i))*h(k)

                if (mod(i,2)==0) then                       
                    inte_simp(k,i) = f(x_n(k,i))*h(k)*4._dp/3._dp   !casos pares para simpson
                else
                    inte_simp(k,i) = f(x_n(k,i))*h(k)*2._dp/3._dp   !casos impares para simpson
                endif

            endif

        enddo

        call gauss(n,0,0._dp,1._dp,x,w) 

        do i=1,size(x)
            suma_gauss(k) = suma_gauss(k)+f(x(i))*w(i) 
        enddo

        suma(k) = sum(inte(k,:))                        !alojo los resultados de la integral en un vector
        suma_simp(k) = sum(inte_simp(k,:))              !lo mismo pero para simpson 
        

        err(k) = abs(suma(k)+exp(-1._dp)-1._dp)              !error trapecio
        err_simp(k) = abs(suma_simp(k)+exp(-1._dp)-1._dp)    !error simpson
        err_gauss(k) = abs(suma_gauss(k)+exp(-1._dp)-1._dp) 

        write(11,'(i3,5x,3(e8.2,5x))') n, err(k), err_simp(k), err_gauss(k)

    enddo

    contains

    real(dp) function f(xx)
        implicit none
        real(dp), intent(in) :: xx

        f = exp(-xx)
    end function f

    end program g1p4
