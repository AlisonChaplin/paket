c ���� kritzn ��� ���᫥��� ����᪨� �ᥫ(�ணࠬ�� frn1).
c ����� �࣠������� �뤠� १���⮢ � 䠩� fname,�஬�������
c १���⮢ � 䠩� fname1 � १���⮢, ����室���� ��� ����p�����
c ��䨪�� � ������� ����� Maple, � fname2. 
c ����� �࣠������� 横� �� ������ �� ��ࠬ��஢
c ��ࠢ�﫠 9.12.98
      SUBROUTINE kritzn(n,ky,kp,prmt)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION PRMT(33),XPM[ALLOCATABLE](:),ck[ALLOCATABLE](:)
      dimension ig[ALLOCATABLE](:)
      character*12 fname,fname1,fname2,fname4,fname5
      ALLOCATE(XPM(KP+1),STAT=IER)
      ALLOCATE(ck(ky),STAT=IER)
      ALLOCATE(ig(kp),STAT=IER)
      IF(IER.NE.0) GO TO 5
      read(5,*) (xpm(i),i=1,kp)
      READ(5,*) IX1,IX2,IT,IP,HP,NP
      close(5)
c new
      write(*,1) (xpm(i),i=1,kp)
   1  format(2x,7(f10.5,2x))	  
      write(*,80) IX1,IX2,IT
      write(*,81) IP,HP,NP	  
  80  format(3X,'ix1=',i3,', ix2=',i3,', it=',i3)	
  81  format(3X,',ip=',i3,',hp=',e12.5,',np=',i3)

      xpm(kp+1)=prmt(14)
      ipri=10
      fname='rez.txt'
      write(*,74) fname
  74  format(3x,'��� 䠩�� १���⮢:',a)      
      open(6,access="append",file=fname)
      ipp=1
      fname1='prez.txt'
      write(*,75) fname1
  75  format(3x,'��� 䠩�� �஬������� १���⮢:',a)
      open(4,access="append",file=fname1)

      igr=1
      fname2='gr.txt'
      write(*,76) fname2
  76  format(3x,'��� 䠩�� ��� ����஥��� ��䨪��:',a)
      open(3,access="append",file=fname2)
	  
c      fname4='test.txt'
c      open(7,access="append",file=fname4)	 

      fname5='grcl.txt'
      open(8,access="append",file=fname5) 

      iip=0
      PRMT(30)=ipp
      PRMT(13)=IT
      lz=PRMT(14)
      PRMT(15)=IX1
      PRMT(16)=IX2
      prmt(20)=ky
      prmt(22)=n
      prmt(21)=kp
      prmt(32)=ip
      prmt(33)=hp
      WRITE(6,9) prmt(3)
      write(*,9) prmt(3)
      call pechpar(xpm,kp,1,0)
      call pechpar(xpm,kp,1,1)
      write(3,29) prmt(3)
  29  format(3x,'eps=',e12.5)
c  31  format(3x,'����쪮 ��ࠬ��஢ ���ᨢ� XPM �ॡ����� ��� ����஥�',
c     *'�� ��䨪�� (������ �᫮)')
      kk=3
c  30  format(3x,'������ �� �� ���浪� (��',i3,' ���)')
      ig(1)=3
      ig(2)=4
      ig(3)=5
      call pechpar(xpm,kp,1,3)      

      WRITE(6,8) IX1,IX2
      write(*,8) ix1,ix2
   9  FORMAT(3X,'��筮��� ���� eps=',e12.5,'. ���A�����:')
   8  format(3X,'CO�CT. �HA�. O���B�EH� XP(',i3,') � XP(',i3,')')
      write(6,3)
      write(*,3)
      do 4 k=1,np
      open(6,access="append",file=fname)
      open(4,access="append",file=fname1)
      open(3,access="append",file=fname2)
c      open(7,access="append",file=fname4)
      open(8,access="append",file=fname5)	  
      CALL frn1(xpm,prmt,ck,n,ky,kp,1)
	  
c	  alison
      if((xpm(3).lt.0).or.(xpm(4).lt.0).or.(xpm(5).lt.0)) then
      write(3,2)
   2  format('�訡��! ����᪨� �᫠ �⠫� ����⥫�묨')
      STOP   
      endif
      if((xpm(3).gt.3).or.(xpm(4).gt.50).or.(xpm(5).gt.10000)) then
      write(3,11)
  11  format('�訡��! ����᪨� �᫠ �⠫� ᫨誮� ����訬�')
      STOP
      endif	  
	  
      call pechpar(xpm,kp,1,0)
      write(3,32) (xpm(i), i=3,8)
   32 format('[',f7.5,', ',f7.4,', ',f7.2,', ',f7.3,', ',f9.4,', ',f8.5,
     *'],')

      write(8,12) (xpm(i), i=3,5)
   12 format('[',f7.5,', ',f7.4,', ',f7.2,'],')
	 
      if(prmt(4).ne.0) then
      if(lz.eq.0) prmt(14)=1.0
      if(lz.eq.2) prmt(14)=3.0
      CALL frn1(xpm,prmt,ck,n,ky,kp,1)
      prmt(14)=lz
      endif
      write(6,3)
      if(ipri.eq.10) write(*,3)
   3  format(63('-'))
  20  format('� ���ᨢ� ',i3,'  ������⮢')
	  
	  xpm(ip)=xpm(ip)+hp
      close(6)
      close(4)
      close(3)
c      close(7)
      close(8)	  
   4  continue
      open(3,access="append",file=fname2)
      write(3,20) np
	  
  7   DEALLOCATE(XPM,ck,ig,STAT=IERR)
      IF(IERR.NE.0) GO TO 55
      close(6)
      close(3)
      close(4)
c      close(7)
      close(8)	  
 10   STOP
   5  write(6,6) ier
   6  format(3x,'kritzn �訡�� ࠧ��饭�� ���ᨢ�� xpm ��� ck ier=',i6)
      go to 7
  55  write(6,56) ierr
  56  format(3x,'kritzn �訡�� 㤠����� ���ᨢ�� xpm ��� ck ier=',i6)
      go to 7
      END

c===================================================================
