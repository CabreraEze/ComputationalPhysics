module dftmod

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 5.1   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
SUBROUTINE DFT (FR,FI,GR,GI,N)
!
! Subroutine to perform the discrete Fourier transform with
! FR and FI as the real and imaginary parts of the input and
! GR and GI as the corresponding  output.
!
!f2py intent(in) n
!f2py intent(in) fr
!f2py intent(in) fi
!f2py intent(out) gr
!f2py intent(out) gi
!f2py depend(n) gr
!f2py depend(n) gi
  IMPLICIT NONE
  integer, parameter                    :: pr = selected_real_kind(p=13)
  INTEGER, INTENT (IN)                  :: N
  INTEGER                               :: I,J
  REAL(pr)                              :: PI,X,Q
  REAL(pr), INTENT (IN), DIMENSION (N)  :: FR,FI
  REAL(pr), INTENT (OUT), DIMENSION (N) :: GR,GI
!
  PI = 4.0_pr*ATAN(1.0_pr)
  X  = 2._pr*PI/real(N,pr)
!
  DO I = 1, N
    GR(I) = 0.0_pr
    GI(I) = 0.0_pr
    DO J = 1, N
      Q = X*(J-1)*(I-1)
      GR(I) = GR(I) + FR(J)*COS(Q)+FI(J)*SIN(Q)
      GI(I) = GI(I) + FI(J)*COS(Q)-FR(J)*SIN(Q)
    END DO
  END DO
END SUBROUTINE DFT


end module
