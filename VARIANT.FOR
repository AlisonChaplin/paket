c  ��ࠦ����, ������饥 �� �襭�� ������த��� ����� U � �室�饥
c  � ᫠������ ������㤭��� �����樥��, ����� ࠧ���� ⨯. ��� ⨯ 
c  ��।������ ��ࠬ��஬  nvar: 
c      nvar=1 - �᫨ �ᯮ������ U,
c      nvar=2 - �᫨ �ᯮ������ U* (�����. ᮯ�殮����),
c      nvar=3 - �᫨ �ᯮ������ J(U),
c      nvar=4 - �᫨ �ᯮ������ J(U*). 4.12.98   
      SUBROUTINE variant(u,kv,xp1,nvar,n)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION u[reference](1),xp1[reference](1)
c  nvar=1   u,xp1
      if(nvar.ne.1) then
c nvar=2,3,4
      if(nvar.eq.3) then
      call operko(x,u,u,xp1,2,n)
c  nvar=3 � u ⥯��� J(u), ᮮ⢥�����饥 xp1
      else  
c nvar=2,4
      do 2 i=2,n,2
   2  u(i)=-u(i) 
      do 16 i=1,kv
   16 xp1(i)=-xp1(i) 
c  � u ⥯��� u*, ᮮ⢥�����饥 xp1
      if(nvar.eq.4) then
      call operko(x,u,u,xp1,2,n)
c  nvar=4 � u ⥯��� J(u*), ᮮ⢥�����饥 xp1
      endif
c  nvar=2 � u ⥯��� u*, ᮮ⢥�����饥 xp1
      endif
      endif
      return
      END
