      subroutine fcy0(x,y,f,xpm,kp,nn)
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      dimension cy(8),cf(8)
      dimension f[reference](1),y[reference](1),xpm[reference](1)  
c alf=0,betta=0 , c=omega,G,Pr,Prd,eps 
      do 1 i=1,8
      k=i*2
   1  cy(i)=dcmplx(y(k-1),y(k))
      ci=dcmplx(0.d0,1.d0)
      om=xpm(4)
      Gr=xpm(5)
      Pr=xpm(6)
      Prd =xpm(7)
      eps =xpm(8)
      cf(1)=cy(5)
      cf(2)=cy(6)
      cf(3)=cy(7)
      cf(4)=Prd*ci*om*eps*cy(3)+Prd*ci*om*cy(8)
      cf(5)=ci*om*cy(1)
      cf(6)=ci*om*cy(2)-(1.0+eps)*cy(3)-cy(8)
      cf(7)=Pr*ci*om*cy(3)
      cf(8)=cy(4)
c      if (abs(y-0.5).le.1.e-06) then
c      write (4,1) y,(cy(i),i=1,5),(cf(i),i=1,5)
c  1   format(3x,e12.5/5x,'cy'/5(7x,2(e12.5,2x)/),5x,
c     *'cf',/5(7x,2(e12.5,2x)/))
c      endif
      do 2 i=1,8
      f(i*2-1)=dreal(cf(i))
   2  f(i*2)=dimag(cf(i))
      return
      end
