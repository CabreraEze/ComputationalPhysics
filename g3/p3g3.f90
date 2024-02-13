!###################################
!#  El programa realiza una        #
!#  integracion de un polinomio    #
!#  de la forma x**n, utilizando   #
!#  el rng "Mersenne Twister", es  #
!#  necesario tener un directorio  #
!#  con el nombre 'p3' donde el    #
!#  programa almacene los datos    #
!#  generados para su posterior    #
!#  analisis.                      #
!###################################

program p3g3
    use precision, only: pr=>dp
    use mtmod, only: mtran=>grnd !, mtseed=>sgrnd
    use mzranmod, only: mzran=>rmzran
    use ranmod 
    implicit none

    character(len=50), parameter    :: loc='p3/', fmt='(i5,10(e20.12,2x))'

    real(pr)                :: inte(3), err(3), x
    real(pr), parameter     :: a=0._pr, b=1._pr

    integer                :: i, j, ninte
    integer, parameter     :: n=3

    open(33,file=trim(adjustl(loc))//'mcinte.d',status='replace')
    write(33,'(10(a20,5x))') '#n', 'mc_inte', 'is_mc_g2', 'is_mc_g3', 'err mc_inte', 'err is_mc_g2', 'err is_mc_g3'
        do i=1,1000
            ninte = 10*i
            inte(:) = 0._pr 

            do j=1,ninte
                
                !#####################
                x = mtran()*(b-a) + a
                inte(1) = inte(1) + f(x,n)

                !#####################
                x = u(2)
                inte(2) = inte(2) + f(x,n)/g(x,2)

                !#####################
                x = u(3)
                inte(3) = inte(3) + f(x,n)/g(x,3)

            enddo
            inte(1) = inte(1)*(b-a)
            inte(:) = inte(:)/real(ninte,pr)
            err(:) = abs(inte(:)-(1._pr/real(n+1,pr)))

            write(33,fmt) ninte, inte(1), inte(2), inte(3), err(1), err(2), err(3)
        enddo
    close(33)

contains

    ! function fint(x,n)
    !     real(pr), intent(in)    :: x
    !     integer, intent(in)     :: n
    !     real(pr)                :: fint

    !     fint = 1._pr/real(n+1,pr)
    !     fint = fint*x**(n+1)

    ! end function fint

    function f(x,n)
        real(pr), intent(in)    :: x
        integer, intent(in)     :: n
        real(pr)                :: f

        f = x**n

    end function f

    function g(x,k)
        real(pr), intent(in)    :: x
        integer, intent(in)     :: k
        real(pr)                :: g

        g = real(1+k,pr)*x**k

    end function

    function u(k)
        integer, intent(in)     :: k
        real(pr)                :: u

        u = mtran()
        u = u**(1._pr/real(k+1,pr))

    end function u

end program p3g3