c условие нормировки (u,u*)=1 ( не забыть, что для уравнения 
c разветвления не надо делить на (fi0,fi0s)
      SUBROUTINE fnormvn(x,y,f,xpr,n,kp,kv)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION F[reference](1),Y[reference](1),xpr[reference](1),
     *yr[allocatable](:),fr[allocatable](:),xp1[allocatable](:)
      allocate(yr(n),stat=ier)
      allocate(fr(n),stat=ier)
      allocate(xp1(kp+1),stat=ier)
      if(ier.ne.0) go to 30
c 
      nn=n/2 
      do 1 i=1,kp+1
    1 xp1(i)=xpr(3+i)
      xp1(kp+1)=1.0
      call fcy(x,y,f,xp1,kp,nn)
c в массиве f(от 1 до n) находятся правые части для Фоm
      cs=dcmplx(0.0,0.0)
      do i=1,5
      cu=dcmplx(y(2*i-1),y(2*i))
      cs=cs+cu*dconjg(cu)
      enddo 
      f(n+1)=dreal(cs)
      f(n+2)=dimag(cs)
      deallocate(yr,fr,xp1,stat=ier)
      if(ier.ne.0) go to 30
      RETURN
   30 write(6,31) ier
   31 format(3x,'Ошибка в размещении и удалении динамических',
     *'массивов в fnorm',i3)
      stop
      END
