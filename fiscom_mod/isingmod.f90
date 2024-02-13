module ising2D
    use precision, only: pr=>dp
    use mzranmod, only: mzran=>rmzran
    implicit none

contains

    !#######################################################
    !           M :Magnetizacion inst.
    !#######################################################
    function magnetizacion(ss)
        integer, intent(in)     :: ss(:,:)
        real(pr)                :: magnetizacion
        integer                 :: n1, n2, n

        n1 = size(ss,dim=1)
        n2 = size(ss,dim=2)
        n = n1*n2

        magnetizacion = real(sum(ss),pr)/real(n,pr)

    end function magnetizacion

    !#######################################################
    !           |M| :Valor abs. magnetizacion inst.
    !#######################################################
    function magnetizacion_abs(ss)
        integer, intent(in)     :: ss(:,:)
        real(pr)                :: magnetizacion_abs
        integer                 :: n1, n2, n

        n1 = size(ss,dim=1)
        n2 = size(ss,dim=2)
        n = n1*n2

        magnetizacion_abs = real(abs(sum(ss)),pr)/real(n,pr)

    end function magnetizacion_abs

    !#######################################################
    !           H :Energia inst.
    !#######################################################
    function energia(ss)
        implicit none
        integer, intent(in)     :: ss(:,:)
        integer                 :: energia
        integer                 :: n1, n2, i, j 

        n1 = size(ss,dim=1)
        n2 = size(ss,dim=2)
        
        energia = 0
        do i=1,n1
            do j=1,n2
                energia = energia - ss(i,j)*ss((mod(i,n1)+1),j) - ss(i,j)*ss(i,(mod(j,n2)+1)) 
            enddo
        enddo

    end function energia

    !#######################################################
    !           p(dE) = exp(-dE/kT) 
    !#######################################################
    function boltzmann(e,t)
        implicit none 
        integer, intent(in)     :: e
        real(pr), intent(in)    :: t
        real(pr)                :: boltzmann

        boltzmann = exp(-real(e,pr)/t)

    end function boltzmann

    !#######################################################
    !           Subroutina de flipeo
    !#######################################################
    subroutine flip(ss,t)
        implicit none
        integer, intent(inout)  :: ss(:,:)
        real(pr), intent(in)    :: t
        integer                 :: n1 , n2, n, de, i, j, k

        n1 = size(ss,dim=1)
        n2 = size(ss,dim=2)
        n = n1*n2

        do k=1,n
            i = int(mzran()*n1) + 1
            j = int(mzran()*n2) + 1

            de = 2*ss(i,j)*(ss((mod(i,n1))+1,j) + ss(i,(mod(j,n2)+1)) + ss((mod(n1+i-2,n1)+1),j) + ss(i,(mod(n2+j-2,n2)+1)))

            if (de > 0) then
                if (boltzmann(de,t) > mzran()) then
                    ss(i,j) = -ss(i,j)
                endif
            else
                ss(i,j) = -ss(i,j)
            endif
        enddo

    end subroutine flip

end module ising2D