c � �⮩ �ணࠬ�� �࣠������� ���� ������, � ���饭�� � ����� ��  
c �ணࠬ�: kritzn ��� yrv ��� ���. 
c 
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION PRMT(33),iaa(5,20),inf(20,7,9),xp(20),kof(10)
      character*12 fname
      write(*,'(a)') '  ������ ����� �� �㤥� ���� ������.  '
      write(*,'(a)') '           ��� �⮣� ������ '
      write(*,'(a)') '  0-�᫨ ��������� ⮫쪮 ����᪨� �᫠ '
      write(*,'(a)') '  1-�᫨ ��������� �����樥��� �ࠢ�����  '
      write(*,'(a)') '     ࠧ��⢫���� '
      read *,lz
      write(*,'(a)') '                               '
      write(*,'(a)') '  ����� ������⥫�� �⢥� �ॡ�� ����� 1,'
      write(*,'(a)') '          � ����⥫�� - 0.'
      write(*,'(a)') '  -------------------------------------------'
      write(*,'(a)') '                               '
c      write(*,'(a)') '  �� ��� �뤠���� १����� �� ��ᯫ��?'
c      read *,i
      i=1
      if(i.eq.1) then
      ipri=10
      else
      ipri=0
      endif
      prmt(23)=ipri
      if(lz.eq.0) then
c ���� ��� ���. �ᥫ - ����� � ��࠭� 
c      write(*,'(a)') '  �� �㤥� ������� ����� � ��࠭�? '
c      read *,ipr
      ipr=0
      else
c      write(*,'(a)') '                               '
c      write(*,'(a)') '  ��������! �� ���� �����樥�⮢ �������-'
c      write(*,'(a)') '  ��� �ࠢ����� ���� �஬������, ���⮬� �������'
c      write(*,'(a)') '  � ��࠭� �����!'
c      write(*,'(a)') '                               '
      ipr=0
      endif
      if(ipr.eq.0) then
c ��� ������ � ��࠭�
c      write(*,'(a)') '  �� �㤥� ᮧ������ 䠩� ������? '
c      read *,ipr1
      ipr1=0
      if(ipr1.eq.0) go to 4
      else
      ipr1=0
      endif
c ᮧ������ 䠩� ������
      write(*,'(a)') '  ������ ���冷� ��⥬�. '
      read *,n
      n=n*2
      write(*,'(a)') '  ������ ��砫� � ����� ��१�� [a,b] � '
      read *,prmt(1),prmt(2)
      write(*,'(a)') '            �筮��� ����. '
      read *,prmt(3)
      write(*,'(a)') '  ������ ������⢮ ��ࠬ��஢ �����. '
      read *,kp
      if(ipr1.eq.1) then
      write(*,'(a)') '  ������ �� ���浪� ��ࠬ���� �����   '
      write(*,'(a)') '  (� ���ᨢ ��ࠬ��஢  XP(1:��) )  '
      read *,(xp(i),i=1,kp)
      write(*,'(a)') '  ������ ������� ������⮢ ���ᨢ� ��ࠬ��஢ '
      write(*,'(a)') '  XP, ����� �������� �᪮�묨 ᮡ�⢥��묨'
      write(*,'(a)') '  ���祭�ﬨ.'
      read *,ix1,ix2
      write(*,'(a)') '  ������ ���ᨬ��쭮 �����⨬�� �᫮ ���権  '
      write(*,'(a)') '  �� ����  '
      read *,it
      write(*,'(a)') '  �� �㤥� �࣠�����뢠�� 横� �� ������ ��  '
      write(*,'(a)') '  ��ࠬ��஢?   '
      read *,i
      if(i.eq.0) then
      ip=3
      hp=0.0
      np=1
      else
      write(*,'(a)') '  ������ ������ �⮣� �������, ����稭� 蠣�  '
      write(*,'(a)') '  �� ��ࠬ����, ������⢮ ����७�� 横��.    '
      read *,ip,hp,np
      endif
      endif
      write(*,'(a)') '  ������ ������⢮ �ࠥ��� �᫮���, �������� '
      write(*,'(a)') '  � ��砫� ��१��.    '
      read *,ky
      write(*,'(a)') '  �� �㤥� �ᯥ��뢠�� ���祭�� ᮡ�⢥���� '
      write(*,'(a)') '  �襭��?'
      read *,prmt(4)
      if(prmt(4).eq.0) then
      prmt(5)=0.0
      else
      write(*,'(a)') '  ������ ���祭�� 蠣� ����.'
      read *,prmt(5)
      endif
   4  if((ipr1.eq.1).or.(ipr.eq.0)) then
c      write(*,'(a)') '  ������ ��� 䠩�� ������:'
c      read(*,'(a)') fname
      if(lz.eq.0) fname='krit.dat'
      if(lz.eq.1) fname='alf0.dat'
      write(*,76) fname
  76  format(3x,'��� 䠩�� ������:',a)
      endif
      if(lz.eq.0) then
      if(ipr.ne.0) then
      if(ipr1.eq.1) then
      open(6,file=fname)
      write(6,15) n,(prmt(i),i=1,5),kp,ky,lz
  15  format(i4,2(f10.5),e12.5,2(f10.5),3(i4))
      write(6,8) (xp(i),i=1,kp)
      write(6,14) ix1,ix2,it,ip,hp,np
  14  format(4(i4),f10.5,i4)
      close(6)
      endif 
      else
      open(5,file=fname)
      READ(5,*) N,(PRMT(I),I=1,5),KP,KY,lz
      endif
      prmt(14)=0.0
      call kritzn(n,ky,kp,prmt)
      STOP
      else
c lz=1
      if(ipr1.eq.1) then
      open(6,file=fname)
      write(6,16) n,(prmt(i),i=1,5),kp,ky
  16  format(i4,2(f10.5),e12.5,2(f10.5),2(i4))
      endif
      prmt(14)=0.0
      print *,'lz=1'
      if(ipr1.eq.1) then
      m=xp(1)
      write(*,20)
   20 format('   ����쪮 �ࠢ����� ࠧ��⢫���� �ॡ���� ���᫨��? (��
     *yt ����� 10')
      read *,kkyrv
      if(kkyrv.gt.10) write(*,21)
   21 format(3x,'��ࠢ��� � ���ᠭ�� osn.for ࠧ��� ���ᨢ� kof')
      do i=1,kkyrv
      write(*,19) i
   19 format('  ����쪮 �����樥�⮢ � ',i3,' �ࠢ����� ࠧ��⢫����')
      read *,kof(i)
      enddo
      kobr=1
      imax=1
  2   write(*,'(a)') '  ������ ���ଠ�� � ��㯯� ��࠭��� ᫠������'
      write(*,'(a)') '  ���砫� ������� �襭�� ������த��� �����:'
      write(*,'(a)') '  ���孨�, � ��⥬ �� ������ �� ���浪� '
      read *,iaa(1,kobr),iaa(2,kobr),i,ij
      ijk=0
      do j=3,4
  3   if (i.eq.m) then
      iaa(j,kobr)=1
      else
      if (i.ne.m) then
      iaa(j,kobr)=2
      else
      print *,'  �� ���ࠢ��쭮 ����� ������.������ ᭮��'
      read *,i
      ijk=ijk+1
      if (ijk.gt.3) then
      print *,'  �� �������⥫��, ����ﭭ� �訡�����.'
      stop
      endif
      go to 3
      endif
      endif
      i=ij
      enddo
      write(*,'(a)') '  �� ���ࠫ� ���ଠ�� � �襭�� ������த��� '
      write(*,'(a)') '  ����� '
      write(*,17) (iaa(j,kobr),j=1,4)
  17  format(3x,'h=',i4,' j=',i4,' p=',i4,' q=',i4) 
      write(*,'(a)') '  �� ��� ��ࠢ��� ��� ���ଠ��,'
      write(*,'(a)') '  � ����� ᭮��'
      read *,i
      if(i.eq.1) go to 2
      ii=1
  6   write(*,'(a)') '  ������ ����� �����樥��, ���஬� �ਭ����-'
      write(*,'(a)') '  ��� ����塞�� ᫠������ '
      read *,inf(kobr,1,ii)
      do ij=1,2
      write(*,'(a)') '  ������ ���� � ��ன ������� ᮡ�⢥����� '
      if(ij.eq.1) then
      jj=3
      write(*,'(a)') '  �襭�� �������� �����'
      else
      write(*,'(a)') '  �襭�� ᮯ�殮���� �����'
      jj=5
      endif
      read *,inf(kobr,jj,ii),i
      ijk=0 
  5   if (i.eq.m) then
      inf(kobr,jj+1,ii)=1
      else
      if (i.ne.m) then
      inf(kobr,jj+1,ii)=2
      else
      print *,'  �� ���ࠢ��쭮 ����� ',ij,' ������.������ ᭮��'
      read *,i
      ijk=ijk+1
      if (ijk.gt.3) then
      print *,'  �� �������⥫��, ����ﭭ� �訡�����.'
      stop
      endif
      go to 5
      endif
      endif
      if(ij.eq.1) then
      write(*,'(a)') '  �ᯮ������ �������᭮-ᮯ�殮���� '
      write(*,'(a)') '  ᮡ�⢥���� �襭�� �������� �����?'
      read *,inf(kobr,2,ii)
      endif
      enddo
      write(*,'(a)') '  ��ࠦ����, ������饥 �� �襭�� ������த���'
      write(*,'(a)') '  �����  U  �  �室�饥 � ��࠭��� ᫠������,'
      write(*,'(a)') '  ��।���� ⨯ ᫠�������. ������ ��� ⨯ ,'
      write(*,'(a)') '  �.�. ������   '
      write(*,'(a)') '    1 - �᫨ �ᯮ������ U   '
      write(*,'(a)') '    2 - �᫨ �ᯮ������ U* (�����. ᮯ�殮����)'
      write(*,'(a)') '    3 - �᫨ �ᯮ������ J(U)   '
      write(*,'(a)') '    4 - �᫨ �ᯮ������ J(U*)   '
      read *,inf(kobr,7,ii)
      write(*,18)  ii
  18  format(3x,'�� ���ࠫ� ���ଠ�� � ',i4,' ᫠������ ')
      write(*,11) (inf(kobr,i,ii),i=1,7)
      write(*,'(a)') '  �� ��� ��ࠢ��� ��� ���ଠ��,'
      write(*,'(a)') '  � ����� ᭮��'
      read *,i
      if(i.eq.1) go to 6
      write(*,'(a)') '  �� �㤥� ������� �� ���� �� ��㯯� '
      write(*,'(a)') '   ᫠������?'
      read *,i
      if(i.eq.1) then
      ii=ii+1
      if(ii.gt.9) then
      write(*,'(a)') '  �� ��ࠫ� �祭� ������� ��㯯� ᫠������, ࠧ-'
      write(*,'(a)') '  ����� ��, ��稭�� � ᫥���饣� ᫠�������.'
      i=1
      go to 7
      endif
      go to 6
      else
      if(imax.lt.ii) imax=ii
      iaa(5,kobr)=ii
      endif
      write(*,'(a)') '   �� �㤥� �த������ ���᫥���? '
      read *,i
  7   if(i.eq.1) then
      kobr=kobr+1
      if(kobr.gt.20) then
      write(*,'(a)') '  � ��� �祭� ����� ᫠������. �ॡ���� � 䠩��'
      write(*,'(a)') '  OSN.FOR 㢥����� �᫮ 20 � ���ᠭ�� '
      write(*,'(a)') '   ���ᨢ�� iaa � inf.'
      stop
      endif
      go to 2
      endif
      write(6,1) kkoef,imax,kobr,ipt,kkyrv
   1  format(5(i4))
      write(6,1) (kof(i),i=1,kkyrv)
      write(6,8) (xp(i),i=1,kp)
c �᫨ ��ࠬ��஢ ����� 7,� ��� ���� �������� ����� 祬 � ���� �����
   8  format(7(f12.5))
      do j=1,kobr
      write(6,10) (iaa(i,j),i=1,5)
  10  format(5(i4))
      jj=iaa(5,j)
      do ij=1,jj
      write(6,11) (inf(j,i,ij),i=1,7)
  11  format(7(i4))
      enddo
      enddo
c      write(6,12) iper,n1,h1
  12  format(2(i4),f10.5)
      write(6,13) ix1,ix2,it,ip,hp,np,kkyrv
  13  format(4(i4),f10.5,2(i4))
      write(6,10) (kof(i),i=1,kkyrv)
      close(6)
      endif
      open(5,file=fname)
      READ(5,*) N,(PRMT(I),I=1,5),KP,KY
      READ(5,*) kkoef,imax,kobr,ipt,kkyrv
      READ(5,*) (kof(i),i=1,kkyrv)
      prmt(26)=kkoef
      prmt(27)=kobr
      prmt(28)=imax
      call yrv(n,ky,kp,prmt,kkoef,kobr,imax,kkyrv,kof,ipt)
      STOP
      endif
      write(*,'(a)') '  �� ���ࠢ��쭮 㪠���� ������, ������ ᮡ�- '
      write(*,'(a)') '  ࠥ��� ����. '
      stop
      END
