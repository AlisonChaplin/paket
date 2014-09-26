c здесь файлы: frn1,lnzdn,detc
c===================================================================
      subroutine frn1(xpm,prmt,ck,n,ky,kp,ip)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      COMPLEX*16 CK[reference](1),c1[allocatable](:,:)
      DIMENSION PRMT(33),XPM[extern](1)
      EXTERNAL FCY
      allocate(c1(ky,ky),stat=ier)
      if(ier.ne.0) go to 5
      lz=prmt(14)
      xpm(kp+1)=lz
      IF((LZ.EQ.1).OR.(LZ.EQ.3)) GO TO 6
      CALL LNZDN(IR,PRMT,xpm,N,KY,CK,FCY)
      IF (IR) 1,2,1
   1  ipp=prmt(30)
      if(ipp.eq.1) write(4,3) IR
      ipp=prmt(23)
      if(ipp.eq.1) write(*,3) IR
      write(6,3) ir
   3  FORMAT('ЛИHEЙHAЯ ЗAДAЧA (frn1) ДAET ABOCT IR ',I2)
      STOP
   2  if(prmt(23).eq.10) call pechpar(xpm,kp,ip,1)
      if(prmt(31).ne.0) call pechpar(xpm,kp,ip,2)
      IF(PRMT(4).ne.0) GO TO 6
   4  deallocate(c1,stat=ier)
      if(ier.ne.0) go to 5
      return
   6  CALL DETC(prmt,xpm,cd,n,KY,IR,CK,FCY,c1)
      IF (IR) 1,4,1
   5  write(6,7) ier
   7  format(3x,'Ошибка в размещении или удалении с1 в frn1',i3)
      stop
      END
c===================================================================
c февраль 1996 г.
      subroutine LNZDN(IR,prmt,xpm,N,M,CK,FCY)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 CK[reference](1),c[allocatable](:,:),
     *CD,CD1,CD2,CD3
      DIMENSION prmt(33),IZN(4),V(4),
     *U(4),VC(4),xpm[reference](1)
      EXTERNAL FCY
      allocate(c(m,m),stat=ier)
      KOR=0
  3   if(ier.ne.0) go to 8
      DO 17 I=1,4
  17  IZN(I)=0
      ITER=0
      IR=0
      IT=prmt(13)
      ipri=prmt(23)
      ipp=prmt(30)
      EPS=prmt(3)
      ix1=prmt(15)
      ix2=prmt(16)
      RE0=xpm(ix1)
      C0=xpm(ix2)
      R1=dsqrt(EPS)
      IF(C0.GE.1.0) THEN
      dc=R1*C0
      else
      dc=r1
      endif
      IF(RE0.GE.1.0)THEN
      DRE=R1*RE0
      ELSE
      DRE=R1
      endif
      call detc(prmt,xpm,CD,N,M,IR,CK,FCY,c)
      if(ipp.eq.1) write(4,31) re0,c0
      if(ipri.eq.10) WRITE(*,31) re0,c0
  31  FORMAT(4X,'начальные значения ',F12.5,4X,F12.5)
      anev=cdabs(cd)
      if(anev.le.eps) go to 19
      IF(IR.NE.0) GO TO 2
      call bisek(xpm,prmt,cd,ibx,ax,ay,izn,vc,u,v)
      if(ibx) 15,1,15
   15 if(ipp.eq.1) write(4,18)
      if(ipri.eq.10) write(*,18)
   18 format(3x,'Сбой в bisek')
      stop  
    1 if(ipp.eq.1) WRITE(4,23) ANEV
      if(ipri.eq.10) WRITE(*,23) ANEV
   23 FORMAT(4X,'соответствующая им невязка NEV=    ',E12.5)
      xpm(ix1)=RE0+DRE
      call detc(prmt,xpm,CD1,N,M,IR,CK,FCY,c)
      IF(IR.NE.0) GO TO 2
      xpm(ix1)=RE0
      xpm(ix2)=C0+dc
      call detc(prmt,xpm,CD2,N,M,IR,CK,FCY,c)
      IF(IR.NE.0) GO TO 2
      F1=DREAL(CD)
      F2=DIMAG(CD)
      R1=(DREAL(CD1)-F1)/DRE
      R2=(DIMAG(CD1)-F2)/DRE
      F3=(DREAL(CD2)-F1)/DC
      F4=(DIMAG(CD2)-F2)/DC
      D=R1*F4-R2*F3
      F5=F2*F3-F4*F1
      F6=R2*F1-R1*F2
      ITER1=0
      R1=F5/D
      R2=F6/D
      F7=DABS(R1)
      F8=DABS(R2)
      IF(xpm(ix1).GT.1.0) F7=F7/DABS(xpm(ix1))
      IF(xpm(ix2).GT.1.0) F8=F8/DABS(xpm(ix2))
      IF(DABS(D).GE.EPS) GO TO 9
      if(ipp.eq.1) WRITE(4,24) D
      if(ipri.eq.10) WRITE(*,24) d
  24  format(4x,'Система и методе Ньютона (lnzd) плохообусловлена',
     *'определитель d=',e14.7)
      stop
   9  xpm(ix1)=RE0+R1
      xpm(ix2)=C0+R2
      IF(F7.LE.EPS.AND.F8.LE.EPS) GOTO 12
      call detc(prmt,xpm,CD3,N,M,IR,CK,FCY,c)
      IF(IR.NE.0) GO TO 2
      ANEV1=CDABS(CD3)
      if(dabs(cd3).lt.eps) go to 19
      AR1=xpm(ix1)
      AC1=xpm(ix2)
      call bisek(xpm,prmt,cd3,ibx,ax,ay,izn,vc,u,v)
      if(ibx) 27,28,30
   27 if(ipp.eq.1) write(4,29)
      if(ipri.eq.10) write(*,29)
   29 format(3x,'Выход из LNZD с помощью bisek')
      go to 2
   30 if(ibx.eq.2) go to 28
      xpm(ix1)=AX
      xpm(ix2)=AY
      call detc(prmt,xpm,CD1,N,M,IR,CK,FCY,c)
      IF(IR.NE.0) GO TO 2
      if(cdabs(cd1).le.eps) go to 19
      if(cdabs(cd1).le.anev1) then
      ANEV1=CDABS(CD1)
      if(ipp.eq.1) WRITE(4,32) xpm(ix1),xpm(ix2),anev1
      if(ipri.eq.10) WRITE(*,32) xpm(ix1),xpm(ix2),anev1
   32 FORMAT(10X,' bisek собст.знач. ',F12.5,4X,F12.5,' nev=',e12.5)
      cd3=cd1
      go to 4
      else
      xpm(ix1)=ar1
      xpm(ix2)=ac1
      endif
  28  IF(ITER1.EQ.0) ANEV2=ANEV1
      IF (ANEV1-ANEV) 4,4,5
   5  ITER1=ITER1+1
      IF(ANEV1.GT.ANEV2) GO TO 11
      AR=xpm(ix1)
      AC=xpm(ix2)
      ANEV2=ANEV1
C     BЫЧИCЛEHHЫЙ ПO HЬЮTOHУ ШAГ ДEЛИTCЯ ПOПОЛAM 5 PAЗ
  11  IF(ITER1-4) 6,6,7
   7  ANEV1=ANEV2
      xpm(ix1)=AR
      xpm(ix2)=AC
      GO TO 4
   6  D=D*2.0
      R1=F5/D
      R2=F6/D
      if(ipp.eq.1) WRITE(4,21)ANEV,ANEV1,CD3,xpm(ix1),xpm(ix2)
      if(ipri.eq.10) WRITE(*,21)ANEV,ANEV1,CD3,xpm(ix1),xpm(ix2)
   21 FORMAT(4X,'NEV И NEV1: ',2E12.5,4X/' DET=',E12.5,'+I*',E12.5/
     *14X,'1-OE COБCT. ЗH.=',E12.7/14X,'2-OE COБCT. ЗH.=',E12.7)
      GOTO 9
   4  ANEV=ANEV1
      CD=CD3
      RE0=xpm(ix1)
      C0=xpm(ix2)
      if(ipp.eq.1) WRITE(4,25) xpm(ix1),xpm(ix2)
      if(ipri.eq.10) WRITE(*,25) xpm(ix1),xpm(ix2)
   25 FORMAT(10X,'COБCT.ЗHAЧ. ',F12.5,4X,F12.5)
      ITER=ITER+1
      IF(ITER-IT) 1,1,10
  10  ipc=prmt(32)
      hpc=prmt(33)
c     корректировка
      xpm(ipc)=xpm(ipc)+hpc/5.0
      KOR=KOR+1
      write(*,'(a)') '++++++++++++'
      write(3,'(a)') '++++++++++++'
      write(7,16) KOR, xpm(ipc) 
  16  FORMAT(3x, 'KOR=',i3, 3x, 'xpm(ip)=',F14.7, )
      if (KOR.EQ.30) then
      write(7,'(a)') 'KOR=30'
      STOP
      endif	  
      go to 3
   2  if(ipp.eq.1) write(4,14)xpm(ix1),xpm(ix2),ANEV
      if(ipri.eq.10) WRITE(*,14) xpm(ix1),xpm(ix2),ANEV
   14 FORMAT(6X,'1-OE COБCT. ЗH. =',F14.7/
     *6X,'2-OE COБCT. ЗH.=',F14.7,'     NEB=',E10.3)
      deallocate(c,stat=ier)
      if(ier.ne.0) go to 8
      RETURN
  12  if(ipp.eq.1) WRITE(4,13)
      if(ipri.eq.10) WRITE(*,13)
  13  FORMAT(2X,'BЫXOД ИЗ LNZD ПO ШAГУ ИTEPAЦИИ')
      xpm(ix1)=RE0
      xpm(ix2)=C0
      GO TO 2
  19  if(ipp.eq.1) WRITE(4,20)
      if(ipri.eq.10) WRITE(*,20)
  20  FORMAT(2X,'BЫXOД ИЗ LNZD ПO НЕВЯЗКЕ')
      go to 2
   8  write(6,26) ier
  26  format(3x,'Ошибка в размещении или удалении массива с в lnzdn',
     *'  ier=',i3)
      stop
      END
c======================================================================
      subroutine DETC(prmt,xpm,CD,N,M,IR,CK,FCY,c)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 CK[reference](1),C[reference](m,1),
     *CD,C1[allocatable](:,:),CR
      DIMENSION prmt(30),xpm[reference](8),Y[allocatable](:),
     *Y0[allocatable](:)
      EXTERNAL FCY,OUTC
      ALLOCATE(C1(M,M),STAT=IERR)
      ALLOCATE(Y(n),STAT=IERR)
      ALLOCATE(Y0(n),STAT=IERR)
      IF(IERR.NE.0) GO TO 9
      IR=0
      EPS=prmt(3)
      LZ=prmt(14)
      ipri=prmt(23) 
      ipech=prmt(4)
      ipp=prmt(30)
      prmt(4)=0.0
      DO 3 K=1,M
      DO 2 I=1,N
    2 Y0(I)=0.0
      Y0(2*M-1+K*2)=1.0
      CALL drkf(PRMT,Y0,Y,N,IER,xpm,FCY,OUTC)
      IF (IER-10) 5,5,4
    5 DO 8 I=1,M
      C(i,k)=DCMPLX(Y0(I*2-1),Y0(I*2))
    8 C1(i,k)=C(i,k)
    3 CONTINUE
      prmt(4)=ipech
      IF(lz.eq.5) GO TO 7
      CALL cmlin(C1,CK,m,M,0,CD,0,ipp)
      IF((LZ.EQ.1).OR.(LZ.EQ.3)) GO TO 1
      GO TO 7
  1   kk=m-1
      DO 6 K=1,kk
  6   CK(K)=-C(K,M)
      CALL cmlin(C,CK,m,kk,1,cr,1,ipp)
      CK(m)=DCMPLX(1.0,0.0)
      IF(CDABS(CR).gt.prmt(3)) GO TO 10
      if(ipp.eq.10) write(4,101)
  101 FORMAT(2X,'Минор СЛАУ для CK мал')
      GO TO 7
   10 IF(prmt(4).EQ.0) GO TO 7
      DO 12 I=1,N
   12 Y0(I)=0.0
      DO 13 I=1,M
      K=I+M
      Y0(2*K-1)=DREAL(CK(I))
   13 Y0(2*K)=DIMAG(CK(I))
      CALL drkf(PRMT,Y0,Y,N,IER,xpm,FCY,OUTC)
      IF (IER.gt.10) then
      if(ipri.eq.10) PRINT 16
   16 FORMAT(2X,'OШИБKA ПPИ BЫЧИCЛEHИИ COБCTB. Ф-ЦИИ')
      IR=IER
      endif
    7 deallocate(C1,Y,Y0,stat=ierr)
      IF(IERR.NE.0) IR=-101
      return
   9  IR=-100
      GO TO 7
   4  IR=IER
      go to 7
      END
c=====================================================================
