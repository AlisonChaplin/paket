      subroutine pechpar(xpm,kp,ip,ip1)
      IMPLICIT REAL*8(A-B,D-H,O-Z)
      DIMENSION XPM[reference](1)
c ip1=0- печать в файл результатов, ip1=1-печать на экране
c ip=1- печать в файл промежуточных результатов
      if(ip1.eq.0) then
      WRITE(6,1) (XPM(i),i=3,kp)
  1   FORMAT(3X, 'betta=',f7.5,'   om=',f7.4, '   Gr='
     * ,f7.2, '   Pr=',f7.3,3x,'Prd=',f9.4/3x,'eps=',f8.5)
      endif
      if(ip1.eq.1) then
      WRITE(*,1) (XPM(i),i=3,kp)
      endif
      if(ip.eq.1) then
      WRITE(4,1) (XPM(i),i=3,kp)
      endif
      return
      end
