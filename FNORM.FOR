c================================================================
      SUBROUTINE fnorm(x,y,f,xpr,kp,nn)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION F[reference](1),Y[reference](1),xpr[reference](1)
      n=xpr(1)
      kv=xpr(3)
      call fnormvn(x,y,f,xpr,n,kp,kv)
      return
      END
c=================================================================
