program g1p3
    use precision !modulo de precision, dp precision doble
    implicit none

    integer(dp), parameter         :: n_max=2**8+1
    integer(dp)                    :: k, n ,i
    real(dp), dimension(7,n_max)   :: x_n, inte, inte_simp
    real(dp), dimension(7)         :: h, suma, suma_simp, err, err_simp

    x_n = 0._dp
    inte = 0._dp
    suma = 0._dp
    inte_simp = 0._dp
    suma_simp = 0._dp
    err = 0._dp
    err_simp = 0._dp
    !seteando matrices en 0

    open(11,file="g1p3.dat")    !datos de los errores

    do k=1,7
        n=2**(k+1)
        h(k)=1._dp/n    !alojo todos los h(n)

        do i=1,(n+1)
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

        suma(k) = sum(inte(k,:))                        !alojo los resultados de la integral en un vector
        suma_simp(k) = sum(inte_simp(k,:))              !lo mismo pero para simpson 

        err(k) = abs(suma(k)-exp(1._dp)+1)              !error trapecio
        err_simp(k) = abs(suma_simp(k)-exp(1._dp)+1)    !error simpson

        write(11,'(i3,5x,2(e8.2,5x))') n, err(k), err_simp(k)

    enddo

    close(11)

    !comprobando output
    ! write(*,*) suma
    ! write(*,*) ""
    ! write(*,*) suma_simp

contains

function f(xx)      !funcion a integrar
    implicit none
    real(dp)                :: f
    real(dp), intent(in)    :: xx
        f = exp(xx)
end function f

end program
