module ode_solver
    use precision, only: pr=>dp

contains

!#####################################
!           EULER
!#####################################

subroutine euler(h,g,t,y0,y)
    implicit none
    real(pr), intent(in)               :: h, y0(0:), t
    real(pr), intent(out)              :: y(0:size(y0)-1)
  
    interface
        function g(t,x)
            use precision, only: pr=>dp
            implicit none
            real(pr), intent(in)                :: t, x(0:)
            real(pr), dimension(0:size(x)-1)    :: g
        end function g
    end interface 

    y(:) = y0(:) + h*g(t,y0(:))

end subroutine euler

!#####################################
!           RUNGE KUTTA 2
!#####################################

subroutine rk2(h,g,t,y0,y)
    implicit none
    real(pr), intent(in)                :: h, y0(0:), t
    real(pr), intent(out)               :: y(0:size(y0)-1)
    real(pr), dimension(0:size(y0)-1)   :: k1, k2
  
    interface
        function g(t,x)
            use precision, only: pr=>dp
            implicit none
            real(pr), intent(in)                :: t, x(0:)
            real(pr), dimension(0:size(x)-1)    :: g
        end function g
    end interface 

    k1(:) = h*g(t,y0(:))
    k2(:) = h*g(t+0.5_pr*h,y0(:)+0.5_pr*k1(:))

    y(:) = y0(:) + k2(:)

end subroutine rk2

!#####################################
!           RUNGE KUTTA 4
!#####################################

subroutine rk4(h,g,t,y0,y)
    implicit none
    real(pr), intent(in)                :: h, y0(0:), t
    real(pr), intent(out)               :: y(0:size(y0)-1)
    real(pr), dimension(0:size(y0)-1)   :: k1, k2, k3, k4
  
    interface
        function g(t,x)
            use precision, only: pr=>dp
            implicit none
            real(pr), intent(in)                :: t, x(0:)
            real(pr), dimension(0:size(x)-1)    :: g
        end function g
    end interface 
    
    k1(:) = h*g(t,y0(:))
    k2(:) = h*g(t+0.5_pr*h,y0(:)+0.5_pr*k1(:))
    k3(:) = h*g(t+0.5_pr*h,y0(:)+0.5_pr*k2(:))
    k4(:) = h*g(t+h,y0(:)+k3(:))

    y(:) = y0(:) + (k1(:)+2._pr*k2(:)+2._pr*k3(:)+k4(:))/6._pr 

end subroutine rk4

!#####################################
!
!#####################################

end module ode_solver

