    program p3
    use integrales, only: pr,simpson
    implicit none
!    integer, parameter    :: pr=selected_real_kind(p=12)
    integer               :: i,nn
    real(pr)              :: a,b,s
    
    
    do
      write(*,*) 'Ingrese el nro de puntos (ojo que no es el nro de intervalos)'
      read(*,*) nn
      if ( mod(nn,2) == 1 ) then
        exit
      else
        write(*,*) 'n debe ser impar!!!  '
      endif
    enddo
    write(*,*) 'a?'
    read(*,*) a
    write(*,*) 'b?'
    read(*,*) b
    
    write(*,*) a,b
    
    i = 1
    
    call simpson(a,b,nn,ff,s)
  
    
    write(*,*) ' Integral = ',s
    write(*,*) ' error absoluto = ',abs(s - exp(1._pr) + 1._pr)
     
    CONTAINS
  

      function ff(xx)
      implicit none
      real(pr), intent(in)   :: xx
      real(pr)               :: ff
      ff = exp(xx)
      end function ff

    end program
    
    
    
    
