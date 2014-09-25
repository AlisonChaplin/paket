c здесь файлы:drkf,bisek,cmlin cdot,outc,fcy
c=====================================================================
      subroutine drkf(prmt,y0,y,n,ier,par,fcy,outc)
      implicit real*8(a-h,o-z)
      dimension prmt(1),y0[reference](1),y[reference](1),
     *par[reference](1),aux[allocatable](:,:)
      allocate(aux(6,n),stat=ierror)
      if(ierror.ne.0) go to 17
      e=dsqrt(prmt(3))
      h=dsqrt(e)
      x=prmt(1)
      xend=prmt(2)
      nn=n/2
      hp=prmt(5)
      kp=prmt(21)
      n1=(xend-x)/h
      if(n1.eq.0) go to 16
      if(n1.lt.0) n1=-n1
      h=(xend-x)/n1
      hmax=h*2.0
      prmt(6)=0.0
      prmt(7)=x
      if(dabs(prmt(5)).lt.dabs(h)) prmt(5)=h
      if(h*prmt(5).lt.0) prmt(5)=-prmt(5)
      call outc(x,y0,0,n,prmt)
    1 ier=0
   14 do 2 i=1,n
    2 aux(6,i)=y0(i)
      call fcy(x,y0,y,par,kp,nn)
      do 3 i=1,n
      aux(1,i)=y(i)
    3 y0(i)=aux(6,i)+y(i)/4.0*h
      x1=x+h/4.0
      call fcy(x1,y0,y,par,kp,nn)
      do 4 i=1,n
      aux(2,i)=y(i)
    4 y0(i)=aux(6,i)+(aux(1,i)+y(i)*3.0)*0.09375*h
      x1=x+0.375*h
      call fcy(x1,y0,y,par,kp,nn)
      do 5 i=1,n
      aux(3,i)=y(i)
    5 y0(i)=aux(6,i)+(aux(1,i)*1932.0-aux(2,i)*7200.0
     *      +y(i)*7296.0)/2197.0*h
      x1=x+12.0/13.0*h
      call fcy(x1,y0,y,par,kp,nn)
      do 6 i=1,n
      aux(4,i)=y(i)
    6 y0(i)=aux(6,i)+(aux(1,i)*439.0/216.0-aux(2,i)*8.0
     *      +(aux(3,i)*3680.0-y(i)*105.625)/513.0)*h
      x1=x+h
      call fcy(x1,y0,y,par,kp,nn)
      do 7 i=1,n
      aux(5,i)=y(i)
    7 y0(i)=aux(6,i)+(aux(2,i)*2.0-aux(1,i)*8.0/27.0
     *      -aux(3,i)*3544.0/2565.0+aux(4,i)*1859.0/4104.0
     *      -y(i)*0.275)*h
      x1=x+h/2.0
      call fcy(x1,y0,y,par,kp,nn)
      do 8 i=1,n
    8 y0(i)=aux(6,i)
      ir=1
      do 9 i=1,n
      e=dabs(y0(i))
      if (e.le.1.0) e=1.0
      r=dabs(aux(1,i)/360.0-aux(3,i)*384.0/12825.0
     *-aux(4,i)*2197.0/75240.0+aux(5,i)/50.0+y(i)*2.0/55.0)*h
      if(r.gt.e*prmt(3)) go to 12
      if(r.gt.e*prmt(3)/32.0) ir=0
    9 continue
      do 10 i=1,n
      y(i)=y0(i)+(aux(1,i)*16.0/135.0+aux(3,i)*6656.0/12825.0
     *+aux(4,i)*28561.0/56430.0-aux(5,i)*0.18+y(i)*2.0/55.0)*h
   10 y0(i)=y(i)
      x=x+h
      if(dabs(hp).le.dabs(h)) prmt(5)=h
      call outc(x,y0,ier,n,prmt)
      if(prmt(6).ne.0) go to 13
      if(dabs(xend-x).le.prmt(3)*1.0e-04) go to 13
      if(ir.eq.1) h=h+h
      if(dabs(h).gt.dabs(hmax)) h=hmax
      if(x+h-xend) 1,1,11
   11 h=xend-x
      go to 1
   12 h=h/2.0
      ier=ier+1
      if(dabs(h).le.prmt(3)*prmt(3)) go to 19
      if(ier.le.10) go to 14
      call outc(x,y0,ier,n,prmt)
      prmt(5)=hp
   13 deallocate(aux,stat=ierror)
      if(ierror.ne.0) ier=-101
   18 return
   16 ier=13
      go to 13
   19 ier=12
      go to 13
   17 ier=-100
      go to 18
      end
c==========================================================
c февраль 1996.
      SUBROUTINE BISEK(xpm,prmt,CD,IBX,AX,AY,IZN,AC,U,V)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX *16 CD
      DIMENSION IZN(4),U(4),V(4),AC(4),
     *prmt(1),xpm[reference](1)
      ipri=prmt(23)
      IF(DREAL(CD).LT.0.0) IZ=-1
      IF(DREAL(CD).GT.0.0) IZ=1
      IF(DIMAG(CD).LT.0.0) IZI=-1
      IF(DIMAG(CD).GT.0.0) IZI=1
      IBX=0
      IF(IZ) 4,12,1
    1 IF(IZI) 3,15,2
    2 K=1
      GO TO 7
    3 K=2
      GO TO 7
    4 IF(IZI) 5,16,6
    5 K=3
      GO TO 7
    6 K=4
    7 IF(IZN(K).NE.0) GO TO 11
      IZN(K)=IZ
    8 AC(K)=CDABS(CD)
      ix1=prmt(15)
      U(K)=xpm(ix1)
      ix2=prmt(16)
      V(K)=xpm(ix2)
    9 IF(K.GT.2) J=K-2
      IF(K.le.2) J=K+2
      IF(IZN(J).EQ.0) GO TO 10
      ibx=1
      AX=(U(K)+U(J))/2.0
      AY=(V(K)+V(J))/2.0
      R=DABS(U(K)-U(J))/2.0
      R1=DABS(V(K)-V(J))/2.0
      IF(DABS(U(K)).GT.1.0) R=R/DABS(U(K))
      IF(DABS(V(K)).GT.1.0) R1=R1/DABS(V(K))
      IF((R.LE.prmt(3)).AND.(R1.LE.prmt(3))) IBX=-1
   10 RETURN
   11 IF(CDABS(CD)-AC(K)) 8,14,14
   14 IBX=2
      GO TO 10
   12 if(ipri.eq.10) PRINT 13,CD
   13 FORMAT(2X,'REAL(CD)=0   CD=',2E12.5)
      GO TO 1
   15 if(ipri.eq.10) PRINT 17,CD
   17 FORMAT(2X,'IMAG(CD)=0   CD=',2E12.5)
      GO TO 2
   16 if(ipri.eq.10) PRINT 17,CD
      GO TO 5
      END
c=========================================================
      SUBROUTINE CMLIN( A, B, IDIM1, N, M, DET, MODE,ipp)
      implicit real*8(a-h,o-z)
      LOGICAL DSCALE, NOTINV
      COMPLEX*16 A, B, DET, CTEMP, CDOT
      DIMENSION A[reference](IDIM1,N), B[reference](IDIM1,M),
     *V[allocatable](:)
      DATA CONST/0.90949470242E-15/
      allocate(v(n),stat=ierr)
      write(4,*)' cmlin163 '  
 1100 NMINUS=N-1
      DSCALE = (MODE.GE.2)
      NOTINV = (MODE.EQ.0).OR.(MODE.EQ.2)
      PTEST = (4.0*CONST*dFLOAT(N))**2
      write(4,*)' cmlin168 '  	  
      IF((N.LT.1).OR.(M.LT.0).OR.(MODE.LT.0).OR.(MODE.GT.3)) GO TO 5100
 1200 IF(N.NE.1) GO TO 1300
      write(4,*)' cmlin171 '  
      DET=A(1,1)
      IF(dABS(dREAL(DET))+dABS(dIMAG(DET)) .EQ. 0. ) GO TO 5201
      A(1,1) = (1.0,0.0)/A(1,1)
      write(4,*)' cmlin175 '  
      IF(M.EQ.0) go to 6002
 1210 DO 1219 K=1,M
      B(1,K)=A(1,1)*B(1,K)
 1219 CONTINUE
      write(4,*)' cmlin180 '  
      go to 6002
 1300 DO 1309 I=1,N
      SUM=0.
 1310 DO 1319 J=1,N
      SUM = SUM +dREAL(A(I,J))**2 + dIMAG(A(I,J))**2
      write(4,30) SUM
   30 format(2x,'Sum=', e12.5)	  
      write(4,*)' cmlin186 '
 1319 CONTINUE
      if (sum.eq.0) then
      write (4,*) '---------SUM=0!!!!--------'
      if (i.eq.4) then
      write (4,31) A(4,1), A(4,2), A(4,3), A(4,4)
   31 format(2x, 'A(4,1)=', e12.5, 2x, 'A(4,2)=', e12.5, 2x, 'A(4,3)='
     * , e12.5, 2x, 'A(4,4)=', e12.5, 2x)
      endif	  
      endif
      if (sum.eq.0) SUM=0.001
      V(I)=1./SUM
      write(4,*)' cmlin190 '	  
 1309 CONTINUE
      write(4,6007) v(1),v(2),v(3),v(4)
 6007 format(2x,'v(i): ',e12.5,1x,e12.5,1x,e12.5,1x,e12.5)
      write(4,*)' cmlin192 ' 
 2000 DET=1.0
 2010 DO 2019 K=1,N
      KPLUS=K+1
      KMINUS=K-1
      L=K
      PSQMAX=0.
 2020 DO 2029 I=K,N
      CTEMP = -CDOT(-A(I,K),A(I,1),A(1,K),KMINUS,IDIM1,1)
      A(I,K)=CTEMP
      PSQ = V(I)*(dREAL(CTEMP)**2 + dIMAG(CTEMP)**2 )
      IF(PSQ.LE.PSQMAX) GO TO 2029
      PSQMAX=PSQ
      L=I
      write(4,*)' cmlin201 ' 	  
 2029 CONTINUE
      VTEMP=V(K)
      IF(L.EQ.K) GO TO 2011
 2040 DO 2049 J=1,N
      CTEMP=A(K,J)
      A(K,J)=A(L,J)
      A(L,J)=CTEMP
 2049 CONTINUE
      VTEMP=V(L)
      V(L)=V(K)
      DET=-DET
 2011 DET = A(K,K)*DET
       write(4,*)' cmlin214 ' 
c      write(4,6006) k,a(k,k)
c 6006 format(2x,'Ведущий злемент ',i5,' столбца ',e12.5,1x,e12.5)
      IF(DSCALE) DET=dSQRT(VTEMP)*DET
      V(K)=dFLOAT(L)
      IF(PSQMAX.LE.PTEST) GO TO 5200
      CTEMP = 1.0/A(K,K)
      A(K,K)=CTEMP
      IF(KPLUS.GT.N) GO TO 2019
 2050 DO 2059 J=KPLUS,N
      A(K,J) = -CTEMP*CDOT(-A(K,J),A(K,1),A(1,J),KMINUS,IDIM1,1)
 2059 CONTINUE
 2019 CONTINUE
 3000 IF(M.EQ.0) GO TO 4000
 3100 DO 3109 I=1,N
      IROW=INT(V(I))
      IF(IROW.EQ.I) GO TO 3109
 3110 DO 3119 K=1,M
      CTEMP=B(I,K)
      B(I,K)=B(IROW,K)
      B(IROW,K)=CTEMP
      write(4,*)' cmlin235 ' 	  
 3119 CONTINUE
 3109 CONTINUE
 3200 DO 3209 K=1,M
 3210 DO 3219 I=1,N
      B(I,K) = -A(I,I)*CDOT(-B(I,K),A(I,1),B(1,K),I-1,IDIM1,1)
 3219 CONTINUE
 3300 DO 3309 L=1,NMINUS
      I=N-L
      B(I,K) = -CDOT(-B(I,K),A(I,I+1),B(I+1,K),L,IDIM1,1)
 3309 CONTINUE
 3209 CONTINUE
 4000 IF(NOTINV) go to 6002
 4100 DO 4109 K=1,NMINUS
      KPLUS=K+1
 4110 DO 4119 I=KPLUS,N
      A(I,K) = -A(I,I)*CDOT((0.0d00,0.0d00),A(I,K),A(K,K),I-K,IDIM1,1)
      A(K,I) = -CDOT(A(K,I),A(K,KPLUS),A(KPLUS,I),I-K-1,IDIM1,1)
 4119 CONTINUE
 4109 CONTINUE
 4200 DO 4209 K=1,N
 4210 DO 4219 I=1,N
      IF(I.GE.K) GO TO 4212
 4211 A(I,K) = CDOT((0.0d00,0.0d00),A(I,K),A(K,K),N-K+1,IDIM1,1)
      GO TO 4219
 4212 A(I,K) = CDOT(A(I,K),A(I,I+1),A(I+1,K),N-I,IDIM1,1)
 4219 CONTINUE
 4209 CONTINUE
 4300 DO 4309 L=1,N
      K=N-L+1
      KCOL=INT(V(K))
      IF(KCOL.EQ.K) GO TO 4309
 4310 DO 4319 I=1,N
      CTEMP=A(I,K)
      A(I,K)=A(I,KCOL)
      A(I,KCOL)=CTEMP
 4319 CONTINUE
 4309 CONTINUE
      go to 6002
 5100 write(6,5190)N,M,MODE
 5190 FORMAT( / 1X 'CMLIN ..... NON-PERMISSIBLE VALUE OF N, M, OR MODE.'
     1        /  6X 'N    =' I15
     2        /  6X 'M    =' I15
     3        /  6X 'MODE =' I15 / )
      go to 6002
 5200 DET=(0.0,0.0)
 5201 IF( (M.EQ.0).AND.(NOTINV) ) go to 6004
       write(6,5290)
 5290 FORMAT( / 1X 'CMLIN ..... MATRIX IS SINGULAR.' / )
      go to 6002
 6004 if(ipp.eq.1) write(4,6005) l,psqmax
 6005 format(2x,'При преобразовании ',i15,' столбца его максимальный',
     */2x,'по модулю нормированный элемент равен',e12.5)
      go to 6002
 6000 write(6,6001) ierr
 6001 format(3x,'ошибка при размещении раб. массива в CMLIN ier=',i5)
      return
 6002 deallocate(v,stat=ier)
      if(ier.eq.0) return
      write(6,6003) ier
 6003 format(3x,'ошибка закрытии раб. массива в CMLIN ier=',i5)
      return
      END
      COMPLEX*16 FUNCTION CDOT(Z,A,B,N,IASTEP,IBSTEP)
      implicit real*8(a-h,o-z)
      COMPLEX*16 Z,W,A[reference](2),b[reference](2)
      W=(0.0d00,0.0d00)
      IF(N.LE.0) GOTO 2
      DO 1 J=1,N
      W=W+A(1+(J-1)*IASTEP)*B(1+(J-1)*IBSTEP)
 1    CONTINUE
 2    CDOT=Z+W
      RETURN
      END
c================================================================
      SUBROUTINE OUTC(X,Y,IER,N,PRMT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION prmt(1),Y[reference](12)
      IF (PRMT(4).EQ.0.0) GO TO 3
      IF (DABS( X-prmt(2)).LT.1.D-04) GO TO 1
      IF (X.LT.PRMT(7)) GO TO 3
      PRMT(7)=PRMT(7)+PRMT(5)
   1  ipp=prmt(30)
      if(ipp.eq.1) then
      write(4,2) X
      write(4,4) (Y(i),i=1,n)
   2  FORMAT(4X,'X=',E12.5)
   4  FORMAT(7X,2E12.5,2X,2E12.5)
      endif
   3  RETURN
      END
c============================================================
      subroutine fcy(x,y,f,xpm,kp,nn)
      implicit real*8(a-b,d-h,o-z),complex*16(c)
      dimension cy[allocatable](:),cf[allocatable](:)
      dimension f[reference](1),y[reference](1),xpm[reference](1)
      allocate(cy(nn),stat=ier)
      allocate(cf(nn),stat=ier)
      if(ier.ne.0) goto 3
      do 1 i=1,nn
      k=i*2
   1  cy(i)=dcmplx(y(k-1),y(k))
      i=kp+1
      if (xpm(i).lt.2) then
      call fun(x,cy,cf,xpm)
      else
      call funs(x,cy,cf,xpm)
      endif
      do 2 i=1,nn
      f(i*2-1)=dreal(cf(i))
   2  f(i*2)=dimag(cf(i))
      deallocate(cy,cf,stat=ier)
      if(ier.ne.0) goto 3
      return
  3   write(*,*) 'ошибка в размещении массивов в fcy' 
      write(6,*) 'ошибка в размещении массивов в fcy' 
      stop 
      end
