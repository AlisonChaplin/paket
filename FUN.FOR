   	 subroutine fun(y,cy,cf,xpm)
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      dimension cf[reference](1),cy[reference](1),xpm[reference](1)
      ci=dcmplx(0.d0,1.d0)
       al =xpm(2)
       bet =xpm(3)
       om =xpm(4)
       Gr =xpm(5)
       Pr =xpm(6)
       Prd =xpm(7)
       eps =xpm(8)
       v0 =(1.0+eps)*(y*y*y-y)/6
       dv0 =(1.0+eps)*(3*y*y-1.0)/6
       ca =Gr*v0*ci*bet
       al2 =al*al+bet*bet
       ca2 =ci*om+ca+al2
       ca3 =al2+Prd*ca+ci*om*Prd
       cf(1)=cy(7)
       cf(2)=-ci*al*cy(1)-ci*bet*cy(3)
       cf(3)=cy(8)
       cf(4)=cy(9)
       cf(5)=-Gr*eps*Prd*cy(2)+eps*Prd*(ca+ci*om)*cy(4)+ca3*cy(10)
       cf(6)=-ca2*cy(2)-ci*al*cy(7)-ci*bet*cy(8)
       cf(7)=ca2*cy(1)+ci*al*cy(6)
       cf(8)=Gr*dv0*cy(2)+ca2*cy(3)-(1.0+eps)*cy(4)+ci*bet*cy(6)-cy(10)
       cf(9)=-Pr*Gr*cy(2)+(Pr*(ci*om+ca)+al2)*cy(4)
       cf(10)=cy(5)		 
      return
      end
