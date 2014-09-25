      subroutine koef(xpm,xpn,prmt,cc,ipar,n,kv,kp,ii,ca,infkf)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION prmt(30),ipar(5),xpm[reference](1),xpn[reference](1),
     *ca[reference](1),cc[reference](1),infkf[reference](7,1),
     *xpr[allocatable](:),
     *y[allocatable](:),y0[allocatable](:),xp[allocatable](:)
      external fkf,outc
      allocate(y(n*5+ii*2),stat=ier)
      allocate(y0(n*5+ii*2),stat=ier)
      allocate(xp(kp+12+kv*2+ii*6),stat=ier)
      allocate(xpr(kp+kv*3+5),stat=ierr)
      if(ierr.ne.0) go to 13
      write(4,*)' koef14 ' 
      do 12 i=1,kp+1
      xp(7+i)=xpm(i)
   12 xpr(i)=xpm(i)
      write(4,*)' koef18 '
      ipri=prmt(23)
      ky=prmt(20)
      xpr(kp+kv*3+4)=n
      ih=ipar(1)
      is=ipar(2)
      xpr(kp+kv*3+3)=ih
      xpr(kp+kv*3+2)=is
      iak=ipar(3)
      ial=ipar(4)
      if(iak.eq.1) then
      do  i=1,kv
      xpr(kp+1+i)=xpm(i)
      xp(7+i)=xpm(i)
      xp(8+i+kp)=xpn(i)
      enddo
      write(4,*)' koef34 '
      ik=0
      else
      do  i=1,kv
      xpr(kp+1+i)=xpn(i)
      xp(7+i)=xpn(i)
      xp(8+i+kp)=xpm(i)
      enddo
      ik=1
      endif
      if(ial.eq.1) then
      do  i=1,kv
      xpr(kp+1+kv+i)=xpm(i)
      enddo
      write(4,*)' koef48 '
      il=0
      else
      do  i=1,kv
      xpr(kp+1+kv+i)=xpn(i)
      enddo
      il=1
      endif
      do  i=1,kv
      xpr(i)=xpr(kp+1+kv+i)
      enddo
c      write(4,*) 'koef xpr',xpr 
      do i=1,n
      y0(i)=0.0
      enddo
      write(4,*)' koef63 '
      if(is.eq.1) call operko(xaxa,y0,y0,xpr,2,n)
c xaxa и у0 здесь можно портить, надо только изменить знак волн. числа. 
      do i=1,kv
      if(ih.eq.2) xpr(i)=-xpr(i)
      xpr(i)=xpr(i)+xpr(kp+1+i)
      end do
      prmt(29)=ik
      prmt(25)=il
      xp(5)=ik
      xp(7)=il
      DO 1 i=1,kv
    1 xpr(kp+1+kv*2+i)=xpr(i)
c В массиве xpr от 1 до kv находятся волновые числа для (Фskl)h
c В массиве xpr от kp+2 до kp+1+kv находятся волновые числа для Фok
c В массиве xpr от kp+2+kv до kp+1+2kv находятся волновые числа для Фol
c В массиве xpr от kp+2+2kv до kp+1+3kv находятся волновые
c  числа для (Фskl)h
      if(ipri.eq.10) then
      print *,'  разыскивается решение неоднородной задачи для'
      print 19,ih,is,iak,ial
  19  format(6x,'h=',i3,' j=',i3,' p=',i3,' q=',i3)
      endif
      write(4,*)' koef86 '
c если xpr(2)=alfa=0 и xpr(3)=betta=0,то есть вырождение
      if(abs(xpr(2)).gt.0.00001.or.abs(xpr(3)).gt.0.00001) then 
      call nzd(prmt,xpr,cc,kv,n,ky)
      else
      call nzd0(prmt,xpr,cc,kv,n,ky)	  
      endif
      write(4,*)' koef92 '	  
      if(iak.ne.ial) then
      ipr=1  
      else
      iak1=infkf(4,1)
      write(4,*)' koef97 '	  
      if(iak.ne.iak1) then
      ipr=1
      else
      ipr=3
      write(4,*)' koef102 '	  
      if(ii.ne.1) then
      do i=2,ii
      if(infkf(4,i).ne.iak1) then
      ipr=1
      exit
      endif
      end do
      endif
      endif
      endif
      if(ii.eq.1) then
      ipr=ipr+1
      else
      ial1=infkf(6,1)
      do i=2,ii
      if(ial1.ne.infkf(6,i)) go to 3
      end do
      write(4,*)' koef118 '
      ipr=ipr+1
      endif
   3  prmt(17)=ipr
      j1=n*4
      if(ipr.eq.1) j1=j1+n
      if(ipr.eq.4) j1=j1-n
      j=j1+ii*2
      do 16 i=1,j
  16  y0(i)=0.0
      jj=n-ky*2
      do  i=1,ky
      y0(jj+2*i-1)=dreal(cc(ky*ik+i))
      y0(jj+2*i)=dimag(cc(ky*ik+i))
      if(ipr.le.2) then
      y0(jj+n+2*i-1)=dreal(cc(ky*(1-ik)+i))
      y0(jj+n+2*i)=dimag(cc(ky*(1-ik)+i))
      endif
      end do
      write(4,*)' koef137 '
      jj=jj+n
      if(infkf(6,1).eq.1) il=2
      if(infkf(6,1).ne.1) il=3
      if(ipr.le.2) jj=jj+n
      do 17 i=1,ky
      y0(jj+2*i-1)=dreal(cc(ky*il+i))
 17   y0(jj+2*i)=dimag(cc(ky*il+i))
      if((ipr.eq.2).or.(ipr.eq.4)) go to 5
      jj=jj+n
      do 4 i=1,ky
      y0(jj+2*i-1)=dreal(cc(ky*(5-il)+i))
   4  y0(jj+2*i)=dimag(cc(ky*(5-il)+i))
   5  jj=jj+n
      do 6 i=1,ky
      y0(jj+2*i-1)=dreal(cc(ky*4+i))
   6  y0(jj+2*i)=dimag(cc(ky*4+i))
      xp(1)=n
      xp(2)=kp
      xp(3)=kv
      xp(4)=ii
      xp(6)=ipr
      i1=kp+8+kv
      do 18 i=1,kv
   18 xp(i1+i)=xpr(kp+kv*2+1+i)
      i1=kp+8+2*kv
      do 9 i=1,4
  9   xp(i1+i)=ipar(i)
      i1=i1+3
      do 8 i=1,ii
      do 8 i2=2,7
   8  xp(i1+6*(i-1)+i2)=infkf(i2,i)
c      write(4,*) 'koef y0 ',y0 
c      write(4,*) 'koef xp ',xp 
      if(ipri.eq.10) then
      print *,'вычисляются слагаемые коэффициентов амр. ур-ний'
      endif
      write(4,*)' koef174 '
      call drkf(prmt,y0,y,j,ier,xp,fkf,outc)
      if(ier-10) 10,10,14
  10  do 7 i=1,ii
      i1=infkf(1,i)
c вставили для печати слагаемых 3 коэфф
c      if(i1.eq.3) then
c      write(6,*)' часть коэф. для решения неоднородной задачи '
c      write(6,19) ih,is,iak,ial
c      write(6,111) y0(j1+i*2-1),y0(j1+i*2)
c  111 format(3x,e14.7,2x,e14.7)
c      endif 
c вставили для печати слагаемых 3 коэфф
   7  ca(i1)=ca(i1)+dcmplx(y0(j1+i*2-1),y0(j1+i*2))
      deallocate(xpr,y,y0,xp,STAT=IERR)
      IF(IERR.NE.0) GO TO 13
      return
  13  write(6,11) ierr
  11  format(3x,'ошибка размещения или удаления массивов в koef,ier=',
     *i6)
      stop
  14  write(6,15) ier
  15  format(3x,'авост в koef (drkf),ier=',i2)
      stop
      END
