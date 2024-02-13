program g1p4_1
    use precision, only: pc=>dp
    use gaussmod

    implicit none
    character(len=20), parameter   :: fmt="(i5,5x,6(e28.20,5x))"
    integer, parameter     :: k_max = 10
    integer                :: n, i, j
    real(pc), allocatable  :: x(:), w(:)
    real(pc), parameter    :: inte_val = 1._pc - 1._pc/exp(1._pc), a = 0., b = 1.
    real(pc)               :: h, er(0:k_max,3), inte(0:k_max,3)

    open(11,file="g1p4_1.dat")
    !write(11,"(7(a20,5x))") "n", "erTRAPECIO", "erSIMPSON", "er GAUSS", "intTRAPECIO", "intSIMPSON", "intGAUSS"

    do i=0,k_max

        if (i == 0) then
            n = 2
            h = (b-a)/real(n,pc)
            allocate(x(0:n))
            allocate(w(0:n))

            do j=0,n
                x(j) = a + h*j
            enddo

            !##trapecio##
            w(:) = 0.5_pc ; w(1) = 1._pc

            inte(i,1) = sum(w*f(x)*h)
            er(i,1) = abs((inte(i,1)-inte_val)/inte_val)

            !##simpson##
            w(:) = 1._pc/3._pc ; w(1) = 4._pc/3._pc

            inte(i,2) = sum(w*f(x)*h)
            er(i,2) = abs((inte(i,2)-inte_val)/inte_val)
            
            !##gauss##
            w(:) = 0._pc
            x(:) = 0._pc       
            call gauss(n,0,a,b,x,w)

            inte(i,3) = sum(w*f(x))
            er(i,3) = abs((inte(i,3)-inte_val)/inte_val)
            
            deallocate(x)
            deallocate(w)

            write(11,fmt) n, er(i,1), er(i,2), er(i,3), inte(i,1), inte(i,2), inte(i,3)

        else
            n = 5*2**i
            h = (b-a)/real(n,pc)
            allocate(x(0:n))
            allocate(w(0:n))
            
            do j=0,n
                x(j) = a + h*j
            enddo

            !##trapecio##
            w(:) = 0.5_pc ; w(1:(n-1)) = 1._pc

            inte(i,1) = sum(w*f(x)*h)
            er(i,1) = abs((inte(i,1)-inte_val)/inte_val)

            !##simpson##
            w(:) = 1._pc/3._pc
            w(1:(n-1):2) = 4._pc/3._pc ; w(2:(n-1):2) = 2._pc/3._pc

            inte(i,2) = sum(w*f(x)*h)
            er(i,2) = abs((inte(i,2)-inte_val)/inte_val)
            
            !##gauss##
            w(:) = 0._pc
            x(:) = 0._pc 
            call gauss(n,0,a,b,x,w)

            inte(i,3) = sum(w*f(x))
            er(i,3) = abs((inte(i,3)-inte_val)/inte_val)
            
            deallocate(x)
            deallocate(w)

            write(11,fmt) n, er(i,1), er(i,2), er(i,3), inte(i,1), inte(i,2), inte(i,3)

        endif

    enddo

contains

function f(y)
    implicit none
    real(pc), dimension(:), intent(in) :: y
    real(pc), dimension(size(y))       :: f 
    f(:) = 1._pc/exp(y(:))
end function f

end program g1p4_1