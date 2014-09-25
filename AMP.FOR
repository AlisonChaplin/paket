c koefa,fkfa,fnzd,fnzvn,nzd,oll,scal1
      subroutine koefa(xpm,prmt,cc,n,kp,cal)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION prmt(30),xpm[reference](1),
     *cc[reference](1),xpr[allocatable](:),
     *y[allocatable](:),y0[allocatable](:)
      external fkfa,outc  
      allocate(y(n*2+2),stat=ier)
      allocate(y0(n*2+2),stat=ier)
      allocate(xpr(kp+1),stat=ierr)
      if(ierr.ne.0) go to 4
      write(4,*)' koefa13 ' 	  
      ky=prmt(20)
      do  i=1,kp+1
      xpr(i)=xpm(i)
      enddo
      do  i=1,2*n+2
      y0(i)=0.0
      enddo
      j=n-2*ky 
      do  i=1,ky
      y0(j+2*i-1)=dreal(cc(i))
      y0(j+2*i)=dimag(cc(i))
      y0(n+j+2*i-1)=dreal(cc(ky*2+i))
      y0(n+j+2*i)=dimag(cc(ky*2+i))
      enddo
      j=2*n+2 
      write(4,*) 'koefa xpr',xpr
c      prmt(4)=1
      if(prmt(23).eq.10) print *,'вычисляются коэф. при лин. слагаемых'
      call drkf(prmt,y0,y,j,ier,xpr,fkfa,outc)
      if(ier-10) 2,2,3
  2   cal=dcmplx(y0(2*n+1),y0(2*n+2))
      write(4,*)' koefa35 '  
      deallocate(xpr,y,y0,STAT=IERR)
      IF(IERR.NE.0) GO TO 4
      return
  4   write(6,11) ierr
  11  format(3x,'ошибка размещения или удаления массивов в koef,ier=',
     *i6)
      stop
  3   write(6,15) ier
  15  format(3x,'авост в koef2 (drkf),ier=',i2)
      stop
      END
c=================================================================
      SUBROUTINE fkfa(X,Y,F,xpr,kp,nn)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION F[reference](1),Y[reference](1),xpr[reference](1),
     *yr[allocatable](:),v[allocatable](:)
      allocate(yr(nn-1),stat=ier)
      allocate(v(nn-1),stat=ier)
      if(ier.ne.0) go to 30
      n1=nn-1
      n2=n1/2
      xpr(kp+1)=0
      call fcy(x,y,f,xpr,kp,n2)
c правые части для линейной задачи
      xpr(kp+1)=2
      do  i=1,n1
      yr(i)=y(n1+i)
      enddo
      call fcy(x,yr,v,xpr,kp,n2)
c правые части для сопряженной задачи
      do  i=1,n1
      f(n1+i)=v(i)
      enddo
      call operko(x,y,v,xpr,0,n1)
c v=B1(fi0)
      cs=(0.0,0.0)
      do i=1,4
      cr=dcmplx(v(2*i-1),v(2*i))
      cr1=dcmplx(yr(2*i-1),yr(2*i))
      cs=cs+cr*cr1
      enddo
      f(2*n1+1)=dreal(cs) 
      f(2*n1+2)=dimag(cs) 
      deallocate(yr,v,stat=ier)
      if(ier.ne.0) go to 30
      RETURN
   30 write(6,31) ier
   31 format(3x,'Ошибка в размещении и удалении динамических',
     *'массивов в fkfavn',i3)
      stop
      END
c================================================================
      SUBROUTINE fnzd(X,Y,F,xpr,kp,nn)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION F[reference](1),Y[reference](48),xpr[reference](21)
      kv=xpr(1)
      n=xpr(kp+kv*3+4)
      nn=n/2
      k=xpr(kp+2)
      l=xpr(kp+kv+2)
      DO 1 I=1,n
      Y(n*3+I)=Y(I)
   1  y(i)=y(n*2+i)
      DO 2 I=1,kv
   2  xpr(i)=xpr(kp+kv*2+1+I)
      call fcy(x,y,f,xpr,kp,nn)
      DO 3 I=1,n
      f(n*2+i)=f(i)
   3  Y(I)=Y(n+I)
      DO 4 I=1,kv
   4  xpr(i)=xpr(kp+kv+1+I)
      call fcy(x,y,f,xpr,kp,nn)
      DO 5 I=1,n
      f(n+i)=f(i)
   5  Y(I)=Y(3*n+I)
      if(k.ne.l) then
      DO 6 I=1,kv
   6  xpr(i)=xpr(kp+1+I)
      call fcy(x,y,f,xpr,kp,nn)
      endif
      call fnzvn(x,y,f,xpr,n,kv,kp)
      xpr(1)=kv
      return
      END
c==================================================================
      SUBROUTINE fnzvn(x,y,f,xpr,n,kv,kp)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION u[allocatable](:),v[allocatable](:),w[allocatable](:),
     *w1[allocatable](:),y[reference](1),xpr[reference](21),
     *f[reference](1)
      allocate (u(n),stat=ierr)
      allocate (v(n),stat=ierr)
      allocate (w(n),stat=ierr)
      allocate (w1(n),stat=ierr)
      if(ierr.ne.0) go to 7
      ih=xpr(kp+kv*3+3)
      is=xpr(kp+kv*3+2)
      do 1 i=1,n
      u(i)=y(i)
   1  v(i)=y(n+i)
      if(ih.eq.2) then
      do i=2,n,2
      v(i)=-v(i)
      enddo
      do i=1,kv
      xpr(i)=-xpr(kp+kv+1+i)
      enddo
      else
      do  i=1,kv
      xpr(i)=xpr(kp+kv+1+i)
      enddo
      endif
      if(is.eq.1) call operko(x,v,v,xpr,2,n)
      call oll(x,u,v,w,xpr,1,n)
      if(ih.ne.1) go to 6
      do 4 i=1,kv
   4  xpr(i)=xpr(kp+1+i)
      call oll(x,v,u,w1,xpr,1,n)
      do 5 i=1,n
   5  w(i)=w(i)+w1(i)
   6  do 9 i=1,n
   9  f(2*n+i)=f(2*n+i)+w(i)
      deallocate(u,v,w,w1,stat=ierr)
      if(ierr.ne.0) go to 7
      return
   7  write(6,8) ierr
      write(*,8) ierr
      stop
   8  format(3x,'Ошибка в fnzvn при размещении или удалении',
     *'массивов u,v,w,w1',i6)
      end
c===========================================================
      SUBROUTINE nzd(prmt,xpr,cc,kv,n,ky)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      dimension prmt(30),xpr[reference](1),
     *cc[reference](1),y0[allocatable](:),
     *y[allocatable](:),c1[allocatable](:,:),ck[allocatable](:),
     *yr[allocatable](:)
      external fcy,fnzd,outc
      allocate(y0(4*n),stat=ierr)
      allocate(y(4*n),stat=ierr)
      allocate(yr(n),stat=ierr)
      allocate(ck(ky),stat=ierr)
      allocate(c1(ky,ky),stat=ierr)
      if(ierr.ne.0) go to 15
      lz=prmt(14)
      ipp=prmt(30)
      prmt(14)=5.0
      call detc(prmt,xpr,cd,n,ky,ir,ck,fcy,c1)
      if(ir.le.10) go to 3
      write(6,12) ir
   12 format(2x,'Авост в nzd(detc) ir=',i3)
      stop
    3 do 5 i=1,n*4
    5 y0(i)=0.0
      prmt(14)=lz
      ik=prmt(29)
      il=prmt(25)
      do 6 i=1,ky
      j=n-ky*2
      y0(j+(2*i-1))=dreal(cc(i+ik*ky))
      y0(j+2*i)=dimag(cc(i+ik*ky))
      y0(j+(2*i-1)+n)=dreal(cc(il*ky+i))
   6  y0(j+i*2+n)=dimag(cc(il*ky+i))
c  В массиве y0 от 1 до n начальные значения для Фok,
c В массиве y0 от n до 2n начальные значения для Фol
      xpr(1)=kv
c Элемент массива xpr временно испортила
      call drkf(prmt,y0,y,n*3,ier,xpr,fnzd,outc)
      if (ier.le.10) go to 13
      write(6,14) ier
  14  format(3x,'авост в nzd (drkf) ier=',i2)
      stop
  13  do 11 i=1,ky
      K=I*2
  11  ck(i)=-dcmplx(y0(2*n+k-1),y0(2*n+k))
c      write(4,*) ' nzd ck ' 
      jpr=0
      do 8 k=1,ky
      do  i=1,ky
      if(cdabs(c1(i,k)).gt.0.00001) go to 8
      enddo
      if(jpr.eq.0) then
      jpr=k
      if(cdabs(ck(k)).gt.0.00001) then
      write(6,*) 'в nzd неоднородная система несовместна'
      stop
      endif
      else
      write(6,*) 'в nzd две строки матрицы c1  нулевые'
      stop
      endif
   8  continue
      if (jpr.eq.0) then
      call cmlin(c1,ck,ky,ky,1,cd,1,ipp)
      do 2 i=1,ky
   2  cc(ky*4+I)=ck(i)
      else
      jk=1
      do 17 k=1,ky
      if(k.eq.jpr) go to 17
      ji=1
      do 21 i=1,ky
      if(i.eq.jpr) go to 21
      c1(jk,ji)=c1(i,k)
      ji=ji+1
  21  continue
      ck(jk)=ck(k)
      jk=jk+1
  17  continue
      jk=ky-1
      call cmlin(c1,ck,ky,jk,1,cd,1,ipp)
      ji=1
      do 20 i=1,ky
      if(i.eq.jpr) then
      cc(ky*4+i)=1.0
      else
      cc(ky*4+I)=ck(ji)
      ji=ji+1
      endif
  20  continue
      endif
      deallocate(y,y0,c1,ck,yr,stat=ierr)
      if(ierr.ne.0) go to 15
      RETURN
  15  write(6,16) ierr
  16  format(3x,'Ошибка при размещении или удалении массивов в nzd ier=',
     *i3)
      stop
      END
c============================================================
      subroutine oll(x,u,v,w,xpr,ipr,n)
      IMPLICIT REAL*8(a-b,d-h,o-z),complex*16(c)
      dimension u[reference](1),v[reference](1),w[reference](1),
     *xpr[reference](1),cu[allocatable](:),cv[allocatable](:),
     *cf[allocatable](:)
      allocate(cu(n/2),stat=ier) 
      allocate(cv(n/2),stat=ier) 
      allocate(cf(n/2),stat=ier)
      if(ier.ne.0) go to 2
      j=n/2
      DO 1 I=1,j
      K=I*2
      cf(i)=DCMPLX(0.0,0.0)
      CU(I)=DCMPLX(U(K-1),U(K))
   1  CV(I)=DCMPLX(V(K-1),V(K))
      call operl(x,cu,cv,cf,xpr,ipr)
      DO  I=1,j
      K=I*2
      w(k-1)=dreal(cf(i))
      w(k)=dimag(cf(i))
      enddo
      deallocate(cu,cv,cf,stat=ier) 
      if(ier.ne.0) go to 2   
      RETURN
   2  write(6,3) ier
   3  format(3x,'ошибка в размещении массивов в oll, ier=',i3)   
      stop
      END
c===============================================================
      subroutine scal1(x,u,v,xpr,csc,n)
      IMPLICIT REAL*8(a-b,d-h,O-Z),COMPLEX*16(c)
      DIMENSION U[reference](1),v[reference](1),xpr[reference](1),
     *cu[allocatable](:),cv[allocatable](:)
      allocate(cu(n/2),stat=ier) 
      allocate(cv(n/2),stat=ier)
      if(ier.ne.0) goto 2
      m=n/2  
      do 1 i=1,m
      k=i*2
      cu(i)=dcmplx(u(k-1),u(k))
   1  cv(i)=dcmplx(v(k-1),v(k))
      csc=dcmplx(0.0,0.0)
      call scal(x,cu,cv,xpr,csc)
      deallocate(cv,cu,stat=ier)
      if(ier.ne.0) goto 2
      RETURN
   2  write(6,3) ier
   3  format(3x,'авост в scal1',i3)
      stop  
      END
