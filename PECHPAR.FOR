      subroutine pechpar(xpm,kp,ip,ip1)
      IMPLICIT REAL*8(A-B,D-H,O-Z)
      DIMENSION XPM[reference](1)
c ip1=0- печать в файл результатов, ip1=1-печать на экране
c ip=1- печать в файл промежуточных результатов
      if(ip1.eq.0) then
      WRITE(6,1) (XPM(i),i=2,kp)
  1   FORMAT(3X,' alfa=',e12.5,' betta=',e12.5,' om=',e12.5/ ' Gr='
     * ,e12.5,3x,' Pr=',e12.5,3x,' Prd=',e12.5,3x,' epsilon=',e15.8)
      endif
      if(ip1.eq.1) then
      WRITE(*,1) (XPM(i),i=2,kp)
      endif
      if(ip.eq.1) then
      WRITE(4,1) (XPM(i),i=2,kp)
      endif
      return
      end
