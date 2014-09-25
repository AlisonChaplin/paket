c==================================================================
      SUBROUTINE fnzd0(X,Y,F,xpr,kp,nn)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION F[reference](1),Y[reference](48),xpr[reference](21)
      kv=xpr(1)
      n=xpr(kp+kv*3+4)
      nn=n/2
      k=xpr(kp+2)
      l=xpr(kp+kv+2)
      n0=n-4
      DO  I=1,n
      Y(n*3+I)=Y(I)
      if(i.le.n0) y(i)=y(n*2+i)
      enddo
      DO 2 I=1,kv
   2  xpr(i)=xpr(kp+kv*2+1+I)
      call fcy0(x,y,f,xpr,kp,nn)
      DO 3 I=1,n
      if(i.le.n0) f(n*2+i)=f(i)
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
      call fnzvn0(x,y,f,xpr,n,kv,kp)
      xpr(1)=kv
      return
      END
c===================================================================
      SUBROUTINE fnzvn0(x,y,f,xpr,n,kv,kp)
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
      do  i=1,n
      w(i)=w(i)+w1(i)
      enddo
c 6    if(abs(x-1).le.0.00000001) then
c      write(4,*) 'oll w ',w
c      write(4,*) 'oll  f',(f(2*n+i),i=1,n)
c      endif
   6  do 9 i=1,n-4
   9  f(2*n+i)=f(2*n+i)+w(i+2)
c      if(abs(x-1).le.0.00000001) then
c      write(4,*) 'oll  f',(f(2*n+i),i=1,n)
c      endif
      deallocate(u,v,w,w1,stat=ierr)
      if(ierr.ne.0) go to 7
      return
   7  write(6,8) ierr
      write(*,8) ierr
      stop
   8  format(3x,'Ошибка в fnzvn при размещении или удалении',
     *'массивов u,v,w,w1',i6)
      end
c===============================================================
      SUBROUTINE nzd0(prmt,xpr,cc,kv,n,ky)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      dimension prmt(30),xpr[reference](1),
     *cc[reference](1),y0[allocatable](:),
     *y[allocatable](:),c1[allocatable](:,:),ck[allocatable](:),
     *yr[allocatable](:)
      external fcy0,fnzd0,outc
      allocate(y0(4*n),stat=ierr)
      allocate(y(4*n),stat=ierr)
      allocate(yr(n),stat=ierr)
      allocate(ck(ky),stat=ierr)
      allocate(c1(ky-1,ky-1),stat=ierr)
      if(ierr.ne.0) go to 15 
      lz=prmt(14)
      ipp=prmt(30)
      prmt(14)=5.0
      n0=n-4
      ky0=ky-1
      call detc(prmt,xpr,cd,n0,ky0,ir,ck,fcy0,c1)
      write(4,*)' nzd0114 '	  
      if(ir.le.10) go to 3
      write(6,12) ir
   12 format(2x,'Авост в nzd(detc) ir=',i3)
      stop
c посчитана матрица С1 размера (6 Х 6)
    3 do 5 i=1,n*4-4
    5 y0(i)=0.0
      prmt(14)=lz
      ik=prmt(29)
      il=prmt(25)
c      write(4,*) 'ik il ',ik,il
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
      call drkf(prmt,y0,y,n*3-4,ier,xpr,fnzd0,outc)
      if (ier.le.10) go to 13
      write(6,14) ier
  14  format(3x,'авост в nzd (drkf) ier=',i2)
      stop
  13  do 11 i=1,ky-1
      K=I*2
  11  ck(i)=-dcmplx(y0(2*n+k-1),y0(2*n+k))
c      write(4,*) 'ck  ',ck
      jpr=0
      do 8 k=1,ky-1
      do  i=1,ky-1
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
      write(4,*)' nzd0165 ' 	  
      call cmlin(c1,ck,ky-1,ky-1,1,cd,1,ipp)
      write(4,*)' nzd0167 '   	  
      do 2 i=1,ky-1
   2  cc(ky*4+I)=ck(i)
      cc(ky*5)=dcmplx(1.0,0.0)     
      else
      jk=1
      do 17 k=1,ky-1
      if(k.eq.jpr) go to 17
      ji=1
      do 21 i=1,ky-1
      if(i.eq.jpr) go to 21
      c1(jk,ji)=c1(i,k)
      ji=ji+1
  21  continue
      ck(jk)=ck(k)
      jk=jk+1
  17  continue
      jk=ky-2
      call cmlin(c1,ck,ky-2,jk,1,cd,1,ipp)  
      do 20 i=1,ky-1
      if(i.eq.jpr) then
      cc((ky-1)*4+i)=1.0
      else
      cc((ky-1)*4+I)=ck(ji)
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
c=================================================================
