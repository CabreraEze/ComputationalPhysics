module pde_solver
    use precision, only: pr=>dp

contains

subroutine triangularizar(a,b,c,d,x)
    implicit none
    real(pr), intent(in), dimension(:)     :: a,b,c,d
    real(pr), intent(out)                  :: x(:)
    real(pr), dimension(size(x))           :: h, p
    integer                                :: i, n

    n=size(x)
    h(1) = c(1)/d(1)
    p(1) = b(1)/d(1)
    do i = 2, n
        h(i) = c(i)/(d(i)-a(i)*h(i-1))
        p(i) = (b(i)-a(i)*p(i-1))/(d(i)-a(i)*h(i-1))
    enddo

    x(:) = p(:)
    do i = n-1,1,-1
        x(i) = p(i) - h(i)*x(i+1)
    enddo

end subroutine triangularizar

subroutine explicita(nu,x)
    implicit none
    real(pr), intent(inout) :: x(:,:)
    real(pr), intent(in)    :: nu
    integer                 :: i, n

    n=size(x,1)
    do i = 2, n-1
        x(i,2) = x(i,1) + nu*(x((i+1),1) + x((i-1),1) - 2._pr*x(i,1))
    enddo

end subroutine explicita

subroutine implicita(nu,x)
    implicit none
    real(pr), intent(inout) :: x(:,:)
    real(pr), intent(in)    :: nu
    integer                 :: n
    real(pr), allocatable   :: a(:), b(:), c(:), d(:), tmp(:)

    n=size(x,1)-2
    allocate(a(n),b(n),c(n),d(n),tmp(n))
    a(:) = -nu
    b(:) = x(2:(n+1),1)
        b(1) = b(1) + nu*x(1,2)
        b(n) = b(n) + nu*x((n+2),2)
    c(:) = -nu
    d(:) = 1._pr + 2._pr*nu

    call triangularizar(a,b,c,d,tmp)
    x(2:(n+1),2) = tmp(:)

end subroutine implicita

subroutine cranknic(nu,x)
    implicit none
    real(pr), intent(inout) :: x(:,:)
    real(pr), intent(in)    :: nu
    real(pr)                :: const
    integer                 :: i, n
    real(pr), allocatable   :: a(:), b(:), c(:), d(:), tmp(:)

    n=size(x,1)-2
    allocate(a(n),b(n),c(n),d(n),tmp(n))

    const = 2._pr/nu
    b(1) = x(1,2) + x(1,1) + (const-2._pr)*x(2,1) + x(3,1)
    b(n) = x(n,1) + (const-2._pr)*x((n+1),1) + x((n+2),1) + x((n+2),2)
    do i = 2, n-1
        b(i) = x(i,1) + (const-2._pr)*x((i+1),1) + x((i+2),1)
    enddo

    a(:) = -1._pr
    d(:) = 2._pr + const
    c(:) = -1._pr

    call triangularizar(a,b,c,d,tmp)
    x(2:(n+1),2) = tmp(:)
    
end subroutine cranknic

end module pde_solver