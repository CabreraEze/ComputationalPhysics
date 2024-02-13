       module mzranmod
         INTEGER,PARAMETER :: K4B=selected_int_kind(9)
         INTEGER,PARAMETER :: DP =KIND(1.0D0)
       contains
!***************************************************************************C
!***************************************************************************C
!!
!*********************************************************************
!**                Random Number Generator  (period > 2^94)         **
!**  See  G. Marsaglia and A. Zaman, Comp. Phys. 8, 117-121 (1994)  **
!*********************************************************************
        function rmzran()
        implicit none
        integer (k4b) ::  mzran,n,i,j,k,is,js,ks,ns,mzranset
        real (DP)     ::  rmzran
        save i,j,k,n
        data i,j,k,n/521288629,362436069,16163801,1131199299/
        mzran = i-k
        if(mzran.lt.0) mzran = mzran + 2147483579
        i = j
        j = k
        k = mzran
        n = 69069*n + 1013904243
!        n = ishft(3533*ishft(n,-16)+iand(n,65535),16)+3533*iand(n,65535))
!        n = n + 1013904243
        mzran = mzran + n
!       For random reals on  (0,1): UNI() = .5+.2328306E-9*mzran()
!       For random reals on (-1,1): VNI() = .4656613E-9*mzran()
!
        rmzran = real(mzran,kind=DP)*0.2328306E-9_DP + 0.5_DP 
!
        return
        entry mzranset(is,js,ks,ns)
        i = 1+iabs(is)
        j = 1+iabs(js)
        k = 1+iabs(ks)
        n = ns
        mzranset = n
        return
        end function rmzran
!
!*********************************************************************
!*********************************************************************
!
!
!
!*********************************************************************
!**     subrutina de para asegurar que los enteros sean de          **
!**     bits y cumplan los requisitos necesarios para el generador. **
!**     
!**       Extraida del Numerical Recipes                            ** 
!*********************************************************************
!
         SUBROUTINE mzran_init(is,js,ks,n)
           IMPLICIT NONE
           INTEGER(K4B), intent(in), optional :: n,is,js,ks
           INTEGER(K4B),PARAMETER ::hg=huge(1_K4B),hgm=-hg,hgng=hgm-1
           INTEGER(K4B)::new,j,hgt,nl
             hgt=hg
!           The following lines check that kind value K4B is in fact a 
!           32-bit integer with the usual properties that we expect it 
!           to have (under negation and wrap-around addition).If all of 
!           these tests are satisfied,then the routines that use this 
!           module are portable,even though they go eyond Fortran 90 's 
!           integer model.
           if (hg /=2147483647)call nrerror('ran_init:arith assump 1 fails ')
           if (hgng >=0)call nrerror('ran_init:arith assump 2 fails ')
           if (hgt+1 /=hgng)call nrerror('ran_init:arith assump 3 fails ')
           if (not(hg)>=0)call nrerror('ran_init:arith assump 4 fails ')
           if (not(hgng)<0)call nrerror('ran_init:arith assump 5 fails ')
           if (hg+hgng >=0)call nrerror('ran_init:arith assump 6 fails ')
           if (not(-1_k4b)<0)call nrerror('ran_init:arith assump 7 fails ')
           if (not(0_k4b)>=0)call nrerror('ran_init:arith assump 8 fails ')
           if (not(1_k4b)>=0)call nrerror('ran_init:arith assump 9 fails ')
!
!          si se se quiere iniciar la secuencia de mzran con otras semaillas
!          (no usando las definidas el la funcion mzran)
!
           if (present(is) .and. present(js) .and. present(ks) .and. present(n)) then
             nl = mzranset(is,js,ks,n)
           endif
!
           return
         end subroutine mzran_init
!
!
         SUBROUTINE nrerror(char)
           character, intent(in) :: char
           write(*,*) char
           stop
         end subroutine nrerror
         end module mzranmod
