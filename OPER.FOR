      SUBROUTINE operb(y,cu,cf,xpr)
      IMPLICIT REAL*8(A-b,d-H,O-Z),complex*16(c)
      DIMENSION cu[reference](1),xpr[reference](1),cf[reference](1)
c
      do i=1,8
      cf(i)=cu(i)
      enddo      
      RETURN
      END
c=====================================================================
      SUBROUTINE operj(y,cu,cf,xpr)
      IMPLICIT REAL*8(A-b,d-H,O-Z),complex*16(c)
      DIMENSION cu[reference](1),xpr[reference](1),cf[reference](1)
c cf=j(cu), а xpr-параметры, соответств. cf
      do i=1,10
      cf(i)=cu(i)
      enddo
      xpr(2)=-xpr(2)
      cf(1)=-cu(1)
      cf(7)=-cu(7)      
      RETURN
      END
c=====================================================================
      SUBROUTINE operk(y,cy,cf,xpr)
      IMPLICIT REAL*8(A-b,d-H,O-Z),complex*16(c)
      DIMENSION cy[reference](1),xpr[reference](1),cf[reference](1)
c 
      ci=dcmplx(0.0,1.0)
      eps=xpr(8)
      bet =xpr(3)
      w=(1.0+eps)*(y*y*y-y)/6
      dw=(1.0+eps)*(3*y*y-1.0)/6
      do i=1,4
      cf(i)=cy(i)*ci*bet*w
      enddo 
      cf(3)=cf(3)+cy(2)*dw
      cf(4)=cf(4)-cy(2)
      cf(5)=-eps*cy(2)+eps*ci*bet*w*cy(4)+ci*bet*w*cy(10)
      RETURN
      END
c=====================================================================
      SUBROUTINE operko(y,u,v,xpr,ipr,n)
      IMPLICIT REAL*8(A-b,d-H,O-Z),complex*16(c)
      DIMENSION u[reference](1),xpr[reference](1),v[reference](1),
     *cf[allocatable](:),cy[allocatable](:)
      allocate(cf(n/2),stat=ier)
      allocate(cy(n/2),stat=ier)
      if(ier.ne.0) go to 1 
      m=n/2  
      do i=1,m
      cf(i)=dcmplx(0.0,0.0)
      enddo
      do i=1,m
      k=i*2
      cy(i)=dcmplx(u(k-1),u(k))
      enddo
      if (ipr.eq.0) then
      call operk(y,cy,cf,xpr)
      else
      if (ipr.eq.1) then
      call operk(y,cy,cf,xpr)
      else
      if (ipr.eq.2) then
      call operj(y,cy,cf,xpr)
      else
      if (ipr.gt.2) then
      call operb(y,cy,cf,xpr)
      endif
      endif
      endif
      endif
      do i=1,m
      k=i*2
      v(k-1)=dreal(cf(i))
      v(k)=dimag(cf(i))
      enddo
      deallocate(cy,cf,stat=ier)
      if(ier.ne.0) go to 1 
      return
   1  write(6,2) ier
   2  format(3x,'авост в operko ier=',i3)
      stop
      END
c=====================================================================
      subroutine operl(y,cu,cv,cf,xpr,ipr)
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      dimension cu[reference](1),cv[reference](1),cf[reference](1),
     *xpr[reference](1)
c 
      al=xpr(2)
      bet=xpr(3)
      Gr=xpr(5)
      Pr=xpr(6)
      Prd =xpr(7)

      if(ipr.eq.1) then
      j=6
      else
      j=1
      endif
      ci=dcmplx(0.d0,1.d0)
      do i=1,5  
      cf(j-1+i)=cv(i)*ci*(al*cu(1)+bet*cu(3))
      enddo
      cf(j)=cf(j)+cu(2)*cv(6) 
      cf(j+1)=cf(j+1)-cu(2)*ci*(al*cv(1)+bet*cv(3)) 
      cf(j+2)=cf(j+2)+cu(2)*cv(7) 
      cf(j+3)=cf(j+3)+cu(2)*cv(8)
      cf(j+4)=cf(j+4)+cu(2)*cv(9)	  
      if(ipr.eq.1) then
	  cr=cf(6)
      cr1=cf(7)
      cr2=cf(10)
      cf(5)=-Prd*cr2
      cf(6)=cr1
      cf(7)=-cr
      cf(8)=-cf(8)
      cf(9)=-Pr*cf(9)
      cf(10)=0	   
      endif  
      do 6 i=1,10
   6  cf(i)=-cf(i)*Gr
      return
      end
c=====================================================================
      subroutine scal(y,cu,cv,xpr,csc)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 cu(5),cv(5),csc,cr
      DIMENSION xpr[reference](1)
c 
      csc=dcmplx(0.0,0.0)
      do 2 i=1,5
      cr=dconjg(cv(i))
c   2  csc=csc+cv(i)*cu(i) 
   2  csc=csc+cr*cu(i) 
c      csc=csc*y 
      RETURN
      END
c=====================================================================