c Файл yrv для вычисления коэффициентов уравнения разветвления
c Здесь организуется выдача результатов в файл fname,промежуточных
c результатов в файл fname1 и результатов, необходимых для постpоения
c графиков с помощью пакета Maple, в fname2. 
c Также организуется цикл по одному из параметров
c Исправляла 1.06.2000
      SUBROUTINE yrv(n,ky,kp,prmt,kkoef,kobr,imax,kkyrv,kof,ipt)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION prmt(32),ipar(5),infkf[ALLOCATABLE](:,:),
     *cr[ALLOCATABLE](:),is[allocatable](:,:),kof(1),
     *xpn[ALLOCATABLE](:),XPM[ALLOCATABLE](:),cc[ALLOCATABLE](:),
     *ca[ALLOCATABLE](:),infor[ALLOCATABLE](:,:,:),ig[ALLOCATABLE](:)
      character*12 fname,fname1,fname2,fname3
      ALLOCATE(xpm(KP+1),STAT=IER)
      ALLOCATE(xpn(KP+1),STAT=IER)
      ALLOCATE(ca(kkoef),STAT=IER)
      ALLOCATE(cc(ky*5),STAT=IER)
      ALLOCATE(cr(ky),STAT=IER)
      ALLOCATE(infor(kobr,7,imax),STAT=IER)
      ALLOCATE(infkf(7,imax),STAT=IER)
      ALLOCATE(is(5,kobr),STAT=IER)
      ALLOCATE(ig(kp),STAT=IER)
      IF(IER.NE.0) GO TO 5
      READ(5,*) (XPM(I),I=1,KP)
c      write(*,'(a)') '  Вы будете изменять параметры ?'
c      read *,ippm
      ippm=0
      if(ippm.eq.1) then
c      write(*,'(a)') '  Вы будете изменять параметры массива XPM(kp) ?'
c      read *,ipp
      ipp=0
      if(ipp.eq.1) then
      write(*,25) (xpm(i),i=1,kp)
  25  format(2x,7(f10.5,2x))
  26  write(*,'(a)') '  Введите индекс и новое значение элемента'
      write(*,'(a)') '  массива XPM , который изменяется.    '
      read *,j,r
      xpm(j)=r
      write(*,'(a)') '  Вы будете еще изменять параметры XPM ?   '
      read *,jj
      if(jj.eq.1) go to 26
      write(*,25) (xpm(i),i=1,kp)
      endif
      endif
      do 16 i=1,kobr
      read(5,*) (is(k,i),k=1,5)
      ii= is(5,i)
      read(5,*) ((infor(i,k,j),k=1,7),j=1,ii)
   16 continue
      READ(5,*) IX1,IX2,IT,IP,HP,NP
      if(ippm.ne.0) then  
c      write(*,'(a)') '  Вы хотите изменить паpаметpы?'
c      read *,ippm
      ippm=0
      endif  
      if(ippm.eq.1) then
      write(*,'(a)') '  Вы измените индексы паpаметpов, обьявленных'
      write(*,'(a)') '  искомыми собственными значениями?'
      read *,i
      write(*,'(a)') '  Вы измените индекс паpаметpа, по котоpому оpга-'
      write(*,'(a)') '  зован цикл, шаг цикла и число повтоpений ?'
      read *,k
      if(i.eq.1) write(*,'(a)') '  Введите индексы IX1 иIX2'
      read *,ix1,ix2
      if(k.eq.1) write(*,'(a)') '  Введите индекс IP,шаг HP '
      if(k.eq.1) write(*,'(a)') '  и число повтоpений NP'
      read *,ip,hp,np
      endif
      close(5)
      ipri=prmt(23)
c      write(*,'(a)') '  Ввeдите имя файла результатов:'
c      read(*,'(a)') fname
      fname='rez1.txt'
      write(*,'(a)') '  Имя файла результатов: rez1.txt'	  
      open(6,access="append",file=fname)
      if(ipri.ne.10) then
c если нет печати на экран
      write(6,*) 'Промежуточныe результаты будут храниться в файле 
     *prez1.txt, если он не нужен, удалите его после расчета.' 
      ipp=1
      fname1='prez1.txt'
      open(4,access="append",file=fname1)
      else
c      write(*,'(a)') '  Вы хотите создать файл промежуточных '
c      write(*,'(a)') '          результатов?'
c      read *,ipp
      ipp=1	  
      if(ipp.eq.1) then
c      write(*,'(a)') '  Введите имя этого файла :'
c      read(*,'(a)') fname1
      fname1='prez1.txt'
      write(*,'(a)') '  Имя файла промежуточных результатов: prez1.txt'	  
      open(4,access="append",file=fname1)
      endif
      endif
c      write(*,'(a)') '  Вы хотите создать файл для построения графиков 
c     * зависимости критических чисел от параметров ?'
c      read *,igr
      igr=1
      if(igr.eq.1) then
c      write(*,'(a)') '  Введите имя этого файла :'
c      read(*,'(a)') fname2
      fname2='gr1.txt'
      write(*,'(a)') '  Имя файла для построения графиков 
     *зависимости критических чисел от параметров: gr1.txt'
      open(3,access="append",file=fname2)
      endif
c      write(*,'(a)') '  Вы хотите создать файл для построения графиков
c     * зависимости коэфф. ур-ний разветвления от параметров ?'
c      read *,igr1
      igr1=1
      if(igr1.eq.1) then
c      write(*,'(a)') '  Введите имя этого файла :'
c      read(*,'(a)') fname3
      fname3='grkoef1.txt'
      write(*,'(a)') '  Имя файла для построения графиков зависимости
     * коэфф. ур-ний разветвления от параметров: grkoef1.txt'
      open(7,access="append",file=fname3)
      endif
      kv=4
      PRMT(30)=ipp
      PRMT(13)=IT
      PRMT(15)=IX1
      PRMT(16)=IX2
c  PRMT(19)=n1
      PRMT(20)=ky
      PRMT(21)=kp
      PRMT(22)=n
      PRMT(24)=kv
c      PRMT(18)=h1
      PRMT(8)=igr
      ipech=prmt(4)
      WRITE(6,9) prmt(3)
      call pechpar(xpm,kp,1,0)
      WRITE(6,8) IX1,IX2
      WRITE(6,3) 
      if(ipri.eq.10) then
      WRITE(*,9) prmt(3)
      call pechpar(xpm,kp,1,1)
      WRITE(*,8) IX1,IX2
      WRITE(*,3)
      endif
   9  FORMAT(3X,'Точность расчета eps=',e12.5,'. ПАРAМЕТРЫ: ')
   3  FORMAT(63('='))
   8  format(3X,'COБCT. ЗHAЧ. OБЬЯBЛEHЫ XP(',i3,') И XP(',i3,')')
      if(igr.eq.1) then
      write(3,29) prmt(3)
  29  format(3x,'eps=',e12.5)
c      WRITE(*,31)
c  31  format(3x,'Сколько параметров массива XPM требуются для построен',
c     *'ия графиков'/' (введите число)')
c      read *,kk
      kk=3
c      write(*,30) kk 
c  30  format(3x,'Введите их по порядку (их',i3,' штук)')
c      read *,(ig(k),k=1,kk)
      ig(1)=3
      ig(2)=4
      ig(3)=5
      call pechpar(xpm,kp,1,3)      
      endif
      do 4 k=1,np
      open(6,access="append",file=fname)
      open(4,access="append",file=fname1)
      if(igr.eq.1)  open(3,access="append",file=fname2)
      if(igr1.eq.1)  open(7,access="append",file=fname3)
      if(ipt.eq.0) then
      prmt(4)=0.0
      prmt(14)=0.0
      call frn1(xpm,prmt,cc,n,ky,kp,1)
      endif 
      if(igr.eq.1) then
      write(3,12) (xpm(ig(ii)),ii=1,kk)
   12 format('[',5(e12.5,', '),'],')
      endif 
      call pechpar(xpm,kp,1,0)
      prmt(4)=ipech
   17 do i=1,kkoef
      ca(i)=dcmplx(0.0,0.0)
      enddo
      prmt(14)=1.0
      do i=1,kp+1
      xpn(i)=xpm(i)
      enddo
      xpn(1)=2
      xpn(3)=-xpn(3)      
      call frn1(xpm,prmt,cc,n,ky,kp,1)
      call frn1(xpn,prmt,cr,n,ky,kp,1)
      do i=1,ky
      cc(ky+i)=cr(i)
      enddo
      prmt(14)=3.0
      call frn1(xpm,prmt,cr,n,ky,kp,1)
      do i=1,ky
      cc(2*ky+i)=cr(i)
      enddo
      call frn1(xpn,prmt,cr,n,ky,kp,1)
      do i=1,ky
      cc(3*ky+i)=cr(i)
      enddo
      prmt(14)=0.0
c заполнять эти ячейки обязательно
      xpm(kp+1)=1
      xpn(kp+1)=1
      do 24 i=1,kobr
      ii=is(5,i)
      if(ipri.eq.10) then
      print *,' вычисляются ',ii,' слагаемых ',i,'-ой  группы.'
      write(4,*)' вычисляются ',ii,' слагаемых ',i,'-ой  группы.'
      endif
      write(4,*)' YRV211 ' 
      do 19 j=1,5
  19  ipar(j)=is(j,i)
      write(4,*)' YRV214 ' 
      do 18 j=1,7
      do 18 k1=1,ii
  18  infkf(j,k1)=infor(i,j,k1)
      write(4,*)' YRV218 '   
      call koef(xpm,xpn,prmt,cc,ipar,n,kv,kp,ii,ca,infkf)
      write(4,*)' YRV220 ' 
  24  continue
      call koefa(xpm,prmt,cc,n,kp,cal)
      write(4,*)' YRV223 ' 
      call norm(xpm,prmt,cc,n,kv,kp,cbb)
      write(4,*)' YRV225 ' 
      do ii=1,kkyrv
      ik=kof(ii)
      write(4,*)' YRV228 ' 
      write(6,2) ii
      write(*,3)
   2  format(3x,'   Значение коэффициентов для ',i2,' уравнения',
     *' разветвления')
      if(ipri.eq.10) write(*,2) ii
      if(ipri.eq.10) write(*,10) cal
      if(igr1.eq.1) write(7,10) cal
      write(6,10) cal
      do i=1,ik
      ca(i)=ca(i)/cbb
      if(ipri.eq.10) write(*,1) i,ca(i)
      write(6,1) i,ca(i)
      if(igr1.eq.1) write(7,1) i,ca(i)
      enddo
      if(igr1.eq.1) write(7,3)
   1  format(5x,i2,'  коэффициент =',e12.5,'+i(',e12.5,')')
  10  format(5x,'сигма=',e12.5,'+i(',e12.5,')')
      enddo
      call analise(kkoef,cal,ca,nkoef)
      WRITE(6,3)
      write(*,15)
      close(6)
      close(4)
  15  format(5x,'Если надо прекратить счет введите 1')
      read *,ii
      if(ii.eq.1) then
      if(igr.eq.1) then
      open(3,access="append",file=fname2)
      write(3,20) k
  20  format('у массива ',i2,'  элементов')
      close(3)
      close(7)
      endif 
      stop
      endif
      xpm(ip)=xpm(ip)+hp
      if(igr.eq.1) close(3)
      if(igr1.eq.1) close(7)
   4  continue
      if(igr.eq.1) then
      open(3,access="append",file=fname2)
      write(3,20) np
      endif
   7  DEALLOCATE(xpm,xpn,ca,cc,cr,infor,is,infkf,ig,STAT=IER)
      IF(IER.NE.0) GO TO 5
  11  stop
   5  write(6,6) ier
   6  format(3x,'ошибка размещения массивов в yrv,ier=',i6)
      stop
      END
