      subroutine pechgr(ir,xp,kp,eps)
      IMPLICIT REAL*8(A-B,D-H,O-Z)
      DIMENSION XP[reference](1)
      if (ir.eq.0) then
      WRITE(3,1) eps
   1  FORMAT(3X,'Точность расчета eps=',e12.5)
c      WRITE(3,4) eps
c   4  FORMAT(3X,'печать параметров ,которые фиксированы',e12.5)
      
      WRITE(3,2) 
   2  FORMAT(3X,' alfa   betta   om   Gr   Pr   Prd   epsilon')
      else
      WRITE(3,3) (XP(i),i=2,kp)
   3  FORMAT(3X,e12.5,',',e12.5,',',e12.5,',',e12.5,',',e12.5,','
     * ,e12.5,',',e12.5)
      endif
      return
      end
