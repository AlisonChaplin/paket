c
c
c
c
c==================================================================
      SUBROUTINE fkf(X,Y,F,xp,kp,nn)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION F[reference](1),Y[reference](1),xp[reference](1)
      n=xp(1)
      kv=xp(3)
      call fkfvn(x,y,f,xp,n,kp,kv)
      return
      END
c==================================================================
      SUBROUTINE fkfvn(X,Y,F,xp,n,kp,kv)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION F[reference](1),Y[reference](1),xp[reference](1),
     *yr[allocatable](:),fr[allocatable](:),xp2[allocatable](:),
     *u[allocatable](:),v[allocatable](:),xp1[allocatable](:),
     *xp3[allocatable](:),w[allocatable](:)
      allocate(yr(n),stat=ier)
      allocate(u(n),stat=ier)
      allocate(v(n),stat=ier)
      allocate(w(n),stat=ier)
      allocate(fr(n),stat=ier)
      allocate(xp1(kp+1),stat=ier)
      allocate(xp3(kp+1),stat=ier)
      allocate(xp2(kp+1),stat=ier)
      if(ier.ne.0) go to 30
      ii=xp(4)
      ipr=xp(6)
      ikk=xp(5)
      nn=n/2
      do 1 i=1,kp+1
      xp2(i)=xp(7+i)
      xp3(i)=xp(7+i)
    1 xp1(i)=xp(7+i)
      j1=kp+8
      call fcy(x,y,f,xp1,kp,nn)
c � ���ᨢ� f(�� 1 �� n) ��室���� �ࠢ� ��� ��� ��k
      jj=n
      if(ipr.le.2) then
      do 2 i=1,kv
    2 xp1(i)=xp(i+j1)
      do 3 i=1,n
    3 yr(i)=y(jj+i)
      call fcy(x,yr,fr,xp1,kp,nn)
      do 4 i=1,n
    4 f(jj+i)=fr(i)
c  � ���ᨢ� f(�� n �� 2n) ��室���� �ࠢ� ��� ��� ��n(n.ne.k),
c            �᫨ ipr<=2.
      jj=jj+n
      endif
      jr=j1+kv*2
      if(ikk.eq.0) then
      iamn=1
      ianm=2
      else
      iamn=2
      ianm=1
      endif
      iaks=xp(jr+9)
      if(iaks.eq.iamn) iks=0
      if(iaks.ne.iamn) iks=1
      do 5 i=1,kv
    5 xp1(i)=xp(7+i+iks*(kp+1))
      xp1(kp+1)=3
      do 6 i=1,n
    6 yr(i)=y(jj+i)
      call fcy(x,yr,fr,xp1,kp,nn)
      do 7 i=1,n
    7 f(jj+i)=fr(i)
c � ���ᨢ� f ��室���� �ࠢ� ��� ��� ��ks
      if((ipr.eq.2).or.(ipr.eq.4)) go to 9
      jj=jj+n
      do 8 i=1,kv
    8 xp1(i)=xp(7+i+(1-iks)*(kp+1))
      do 16 i=1,n
   16 yr(i)=y(jj+i)
      call fcy(x,yr,fr,xp1,kp,nn)
      do 13 i=1,n
   13 f(jj+i)=fr(i)
c � ���ᨢ� f ��室���� �ࠢ� ��� ��� ��ns(n.ne.k)
    9 jj=jj+n
      do 10 i=1,kv
   10 xp3(i)=xp(i+8+kp+kv)
      do 12 i=1,n
   12 yr(i)=y(jj+i)
      call fcy(x,yr,fr,xp3,kp,nn)
      do 23 i=1,n
   23 f(jj+i)=fr(i)
c � ���ᨢ� f ��室���� ���. ᫠����� ��. ��⥩ ��� (�skl)h
      do  i=1,kv
      xp1(i)=xp(7+i)
      enddo
      xp1(kp+1)=1.0
      do  i=1,n
      u(i)=y(i)
      enddo
c � ���ᨢ� u - �ok, ��� ᮮ⢥����� xp1
      ial=xp(jr+4)
      if(ial.eq.iamn) then
      do i=1,kv
      xp2(i)=xp(7+i)
      enddo
      do i=1,n
      v(i)=y(i)
      enddo
      else
      do i=1,kv
      xp2(i)=xp(kp+8+i)
      enddo
      do i=1,n
      v(i)=y(n+i)
      end do
      endif
c ���ᨢ 'u'- �� �ok(k=m ��� k=n), ��� ᮮ⢥����� xp1
c ���ᨢ 'v'- �� �ol(k=n ��� k=m), ��� ᮮ⢥����� xp2
      ih=xp(jr+1)
      is=xp(jr+2)
      if(ih.eq.2) then
      do 17 i=2,n,2
  17  v(i)=-v(i)
      do 15 i=1,kv
  15  xp2(i)=-xp2(i)
c  ���ᨢ 'v'-�������᭮-ᮯ�殮� �ol, �᫨ h=2
      endif
      if(is.eq.1) call operko(x,v,v,xp2,2,n)
c  ���ᨢ 'v'-J(�ol), �᫨ s=1.���ᨢ xp2 ᮮ⢥����� 'v'
      call oll(x,u,v,fr,xp2,1,n)
c  ����⠭� L(�ok,�sl(h))
      do 11 i=1,n
  11  f(jj+i)=fr(i)+f(jj+i)
      if(ih.eq.1) then
c   ���ᨢ xp1  ᮮ⢥����� �sl � ��࠭���� h
      call oll(x,v,u,fr,xp1,1,n)
c  ����⠭� L(�sl(h),�ok)
      do 19 i=1,n
   19 f(jj+i)=fr(i)+f(jj+i)
      endif
c � ���ᨢ� f ��室���� �ࠢ� ��� ��� (�skl)h
      jr=jr+4
      ir=jj-n
      if((ipr.eq.1).or.(ipr.eq.3)) ir=ir-n
      iajr=xp(jr+5)
      do 20 i=1,ii
      nvar=xp(jr+6*i)
      iaj=xp(jr+(i-1)*6+3)
      if(iamn.eq.iaj) then
      imk=0
      else
      if(ianm.eq.iaj) then
      imk=1
      else
      print *,'�訡�� � ���ᨢ� infor'
      write(6,*) '�訡�� � ���ᨢ� infor'
      stop
      endif
      endif
      do 21 i1=1,n
      if(ipr.gt.2) then
      yr(i1)=y(i1)
      else
      yr(i1)=y(i1+imk*n)
      endif
   21 fr(i1)=y(i1+jj)
      do 22 i1=1,kv
      xp1(i1)=xp3(i1)
   22 xp2(i1)=xp(7+i1+imk*(kp+1))
      iaj=xp(jr+(i-1)*6+1)
      if(iaj.eq.1) then
      do 14 i1=2,n,2
   14 yr(i1)=-yr(i1)
      do 27 i1=1,kv
   27 xp2(i1)=-xp2(i1)
      endif
      iaj=xp(jr+(i-1)*6+2)
      if(iaj.eq.1) call operko(x,yr,yr,xp2,2,n)
c ���� ���࠭� � L(�os,�  )  ��⮢
      call variant(fr,kv,xp1,nvar,n)
c ��ன ���࠭� � L(�os,(�skl)h)  ��⮢
      call oll(x,yr,fr,v,xp1,2,n)
      call oll(x,fr,yr,w,xp2,2,n)
      do 26 i1=1,n
   26 v(i1)=v(i1)+w(i1)
      iaj=xp(jr+(i-1)*6+5)
      if(iaks.eq.iaj) then
      imk=0
      else
      imk=1
      endif
      do  i1=1,kv
      xp2(i1)=-xp(7+i1+imk*(kp+1))
      enddo
      if(iajr.eq.iaj) then
      do i1=1,n
      fr(i1)=y(ir+i1)
      enddo
      else
      do i1=1,n
      fr(i1)=y(ir+i1+n)
      enddo
      endif
      iaj=xp(jr+(i-1)*6+4)
      if(iaj.eq.1) call operko(x,fr,fr,xp2,2,n)
      do i1=1,n,2
      fr(i1+1)=-fr(i1+1)
      enddo
      call scal1(x,v,fr,xp2,csc,n)
      f(jj+n+(2*i-1))=dreal(csc)
   20 f(jj+n+2*i)=dimag(csc)
      deallocate(yr,u,v,fr,xp3,xp1,xp2,w,stat=ier)
      if(ier.ne.0) go to 30
      RETURN
   30 write(6,31) ier
   31 format(3x,'�訡�� � ࠧ��饭�� � 㤠����� �������᪨�',
     *'���ᨢ�� � fkfvn',i3)
      stop
      END
