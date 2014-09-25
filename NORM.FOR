c нормировка интеграл от u*u-сопряженное*dy равен 1
c 3.141592654
      subroutine norm(xpm,prmt,cc,n,kv,kp,cbb)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION prmt(30),xpm[reference](1),
     *cc[reference](1),xpr[allocatable](:),
     *y[allocatable](:),y0[allocatable](:)
      external fnorm,outc  
      allocate(y(n+2),stat=ier)
      allocate(y0(n+2),stat=ier)
      allocate(xpr(kp+4),stat=ierr)
      if(ierr.ne.0) go to 4
      write(4,*)' norm14 '   
      ky=prmt(20)
      Pi=datan(1.0)*4
      if(dabs(xpm(2)).lt.0.1e-5) then
      anorm=Pi*2/xpm(3)
      else
      anorm=Pi*Pi*4/xpm(3)/xpm(2)
      endif
      do  i=1,kp+1
      xpr(3+i)=xpm(i)
      enddo
      xpr(1)=n
      xpr(2)=kp
      xpr(3)=kv
      j=n+2
      do 1 i=1,j
  1   y0(i)=0.0
      jj=n-ky*2
      do  i=1,ky
      y0(jj+2*i-1)=dreal(cc(i))
      y0(jj+2*i)=dimag(cc(i))
      enddo
      call drkf(prmt,y0,y,j,ier,xpr,fnorm,outc)
      if(ier-10) 2,2,3
   2  y0(n+1)=y0(n+1)*anorm 
      y0(n+2)=y0(n+2)*anorm 
      cbb=dcmplx(y0(n+1),y0(n+2))
      deallocate(xpr,y,y0,STAT=IERR)
      IF(IERR.NE.0) GO TO 4
      write(6,*) 'нормировка интеграл от u*u-сопряженное*dy равен 1'
      return
  4   write(6,11) ierr
  11  format(3x,'ошибка размещения или удаления массивов в koef,ier=',
     *i6)
      stop
  3   write(6,15) ier
  15  format(3x,'авост в norm (drkf),ier=',i2)
      stop
      END
