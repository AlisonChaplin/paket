c В этой программе организован ввод данных, и обращение к одной из  
c программ: kritzn или yrv для Иры. 
c 
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION PRMT(33),iaa(5,20),inf(20,7,9),xp(20),kof(10)
      character*12 fname
      write(*,'(a)') '  Укажите какую Вы будете решать задачу.  '
      write(*,'(a)') '           Для этого введите '
      write(*,'(a)') '  0-если вычисляются только критические числа '
      write(*,'(a)') '  1-если вычисляются коэффициенты уравнений  '
      write(*,'(a)') '     разветвлений '
      read *,lz
      write(*,'(a)') '                               '
      write(*,'(a)') '  Здесь положительный ответ требует ввода 1,'
      write(*,'(a)') '          а отрицательный - 0.'
      write(*,'(a)') '  -------------------------------------------'
      write(*,'(a)') '                               '
c      write(*,'(a)') '  Вы хотите выдавать результаты на дисплей?'
c      read *,i
      i=1
      if(i.eq.1) then
      ipri=10
      else
      ipri=0
      endif
      prmt(23)=ipri
      if(lz.eq.0) then
c ввод для крит. чисел - можно с экрана 
c      write(*,'(a)') '  Вы будете вводить данные с экрана? '
c      read *,ipr
      ipr=0
      else
c      write(*,'(a)') '                               '
c      write(*,'(a)') '  Внимание! При расчете коэффициентов амплитуд-'
c      write(*,'(a)') '  ных уравнений ввод громоздкий, поэтому вводить'
c      write(*,'(a)') '  с экрана нельзя!'
c      write(*,'(a)') '                               '
      ipr=0
      endif
      if(ipr.eq.0) then
c нет вввода с экрана
c      write(*,'(a)') '  Вы будете создавать файл данных? '
c      read *,ipr1
      ipr1=0
      if(ipr1.eq.0) go to 4
      else
      ipr1=0
      endif
c создается файл данных
      write(*,'(a)') '  Введите порядок системы. '
      read *,n
      n=n*2
      write(*,'(a)') '  Введите начало и конец отрезка [a,b] и '
      read *,prmt(1),prmt(2)
      write(*,'(a)') '            точность расчета. '
      read *,prmt(3)
      write(*,'(a)') '  Введите количество параметров задачи. '
      read *,kp
      if(ipr1.eq.1) then
      write(*,'(a)') '  Введите по порядку параметры задачи   '
      write(*,'(a)') '  (в массив параметров  XP(1:кр) )  '
      read *,(xp(i),i=1,kp)
      write(*,'(a)') '  Введите индексы элементов массива параметров '
      write(*,'(a)') '  XP, которые обьявляются искомыми собственными'
      write(*,'(a)') '  значениями.'
      read *,ix1,ix2
      write(*,'(a)') '  Введите максимально допустимое число итераций  '
      write(*,'(a)') '  при расчете  '
      read *,it
      write(*,'(a)') '  Вы будете организовывать цикл по одному из  '
      write(*,'(a)') '  параметров?   '
      read *,i
      if(i.eq.0) then
      ip=3
      hp=0.0
      np=1
      else
      write(*,'(a)') '  Введите индекс этого элемента, величину шага  '
      write(*,'(a)') '  по параметру, количество повторений цикла.    '
      read *,ip,hp,np
      endif
      endif
      write(*,'(a)') '  Введите количество краевых условий, заданных '
      write(*,'(a)') '  в начале отрезка.    '
      read *,ky
      write(*,'(a)') '  Вы будете распечатывать значения собственных '
      write(*,'(a)') '  решений?'
      read *,prmt(4)
      if(prmt(4).eq.0) then
      prmt(5)=0.0
      else
      write(*,'(a)') '  Введите значение шага печати.'
      read *,prmt(5)
      endif
   4  if((ipr1.eq.1).or.(ipr.eq.0)) then
c      write(*,'(a)') '  Введите имя файла данных:'
c      read(*,'(a)') fname
      if(lz.eq.0) fname='krit.dat'
      if(lz.eq.1) fname='yrv.dat'
      write(*,76) fname
  76  format(3x,'Имя файла данных:',a)
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
   20 format('   Сколько уравнений разветвления требуется вычислить? (их
     *yt более 10')
      read *,kkyrv
      if(kkyrv.gt.10) write(*,21)
   21 format(3x,'Исправить в описании osn.for размер массива kof')
      do i=1,kkyrv
      write(*,19) i
   19 format('  Сколько коэффициентов у ',i3,' уравнения разветвления')
      read *,kof(i)
      enddo
      kobr=1
      imax=1
  2   write(*,'(a)') '  Введите информацию о группе выбранных слагаемых'
      write(*,'(a)') '  Сначала индексы решения неоднородной задачи:'
      write(*,'(a)') '  верхний, а затем три нижних по порядку '
      read *,iaa(1,kobr),iaa(2,kobr),i,ij
      ijk=0
      do j=3,4
  3   if (i.eq.m) then
      iaa(j,kobr)=1
      else
      if (i.ne.m) then
      iaa(j,kobr)=2
      else
      print *,'  Вы неправильно ввели индекс.Введите снова'
      read *,i
      ijk=ijk+1
      if (ijk.gt.3) then
      print *,'  Вы невнимательны, постоянно ошибаетесь.'
      stop
      endif
      go to 3
      endif
      endif
      i=ij
      enddo
      write(*,'(a)') '  Вы набрали информацию о решении неоднородной '
      write(*,'(a)') '  задачи '
      write(*,17) (iaa(j,kobr),j=1,4)
  17  format(3x,'h=',i4,' j=',i4,' p=',i4,' q=',i4) 
      write(*,'(a)') '  Вы хотите исправить эту информацию,'
      write(*,'(a)') '  и ввести снова'
      read *,i
      if(i.eq.1) go to 2
      ii=1
  6   write(*,'(a)') '  Введите номер коэффициента, которому принадле-'
      write(*,'(a)') '  жит вычисляемое слагаемое '
      read *,inf(kobr,1,ii)
      do ij=1,2
      write(*,'(a)') '  Введите первый и второй индексы собственного '
      if(ij.eq.1) then
      jj=3
      write(*,'(a)') '  решения линейной задачи'
      else
      write(*,'(a)') '  решения сопряженной задачи'
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
      print *,'  Вы неправильно ввели ',ij,' индекс.Введите снова'
      read *,i
      ijk=ijk+1
      if (ijk.gt.3) then
      print *,'  Вы невнимательны, постоянно ошибаетесь.'
      stop
      endif
      go to 5
      endif
      endif
      if(ij.eq.1) then
      write(*,'(a)') '  Используется комплексно-сопряженное '
      write(*,'(a)') '  собственное решение линейной задачи?'
      read *,inf(kobr,2,ii)
      endif
      enddo
      write(*,'(a)') '  Выражение, зависящее от решения неоднородной'
      write(*,'(a)') '  задачи  U  и  входящее в выбранное слагаемое,'
      write(*,'(a)') '  определяет тип слагаемого. Укажите этот тип ,'
      write(*,'(a)') '  т.е. введите   '
      write(*,'(a)') '    1 - если используется U   '
      write(*,'(a)') '    2 - если используется U* (компл. сопряженное)'
      write(*,'(a)') '    3 - если используется J(U)   '
      write(*,'(a)') '    4 - если используется J(U*)   '
      read *,inf(kobr,7,ii)
      write(*,18)  ii
  18  format(3x,'Вы набрали информацию о ',i4,' слагаемом ')
      write(*,11) (inf(kobr,i,ii),i=1,7)
      write(*,'(a)') '  Вы хотите исправить эту информацию,'
      write(*,'(a)') '  и ввести снова'
      read *,i
      if(i.eq.1) go to 6
      write(*,'(a)') '  Вы будете вычислять еще одно из группы '
      write(*,'(a)') '   слагаемых?'
      read *,i
      if(i.eq.1) then
      ii=ii+1
      if(ii.gt.9) then
      write(*,'(a)') '  Вы выбрали очень большую группу слагаемых, раз-'
      write(*,'(a)') '  делите ее, начиная со следующего слагаемого.'
      i=1
      go to 7
      endif
      go to 6
      else
      if(imax.lt.ii) imax=ii
      iaa(5,kobr)=ii
      endif
      write(*,'(a)') '   Вы будете продолжать вычисления? '
      read *,i
  7   if(i.eq.1) then
      kobr=kobr+1
      if(kobr.gt.20) then
      write(*,'(a)') '  У вас очень много слагаемых. Требуется в файле'
      write(*,'(a)') '  OSN.FOR увеличить число 20 в описании '
      write(*,'(a)') '   массивов iaa и inf.'
      stop
      endif
      go to 2
      endif
      write(6,1) kkoef,imax,kobr,ipt,kkyrv
   1  format(5(i4))
      write(6,1) (kof(i),i=1,kkyrv)
      write(6,8) (xp(i),i=1,kp)
c если параметров больше 7,то они будут вводится более чем в одну строчку
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
      write(*,'(a)') '  Вы неправильно указали задачу, которую соби- '
      write(*,'(a)') '  раетесь решать. '
      stop
      END
