c  Выражение, зависящее от решения неоднородной задачи U и входящее
c  в слагаемое амплитудного коэффициента, имеет различный тип. Этот тип 
c  определяется параметром  nvar: 
c      nvar=1 - если используется U,
c      nvar=2 - если используется U* (компл. сопряженное),
c      nvar=3 - если используется J(U),
c      nvar=4 - если используется J(U*). 4.12.98   
      SUBROUTINE variant(u,kv,xp1,nvar,n)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      DIMENSION u[reference](1),xp1[reference](1)
c  nvar=1   u,xp1
      if(nvar.ne.1) then
c nvar=2,3,4
      if(nvar.eq.3) then
      call operko(x,u,u,xp1,2,n)
c  nvar=3 в u теперь J(u), соответствующее xp1
      else  
c nvar=2,4
      do 2 i=2,n,2
   2  u(i)=-u(i) 
      do 16 i=1,kv
   16 xp1(i)=-xp1(i) 
c  в u теперь u*, соответствующее xp1
      if(nvar.eq.4) then
      call operko(x,u,u,xp1,2,n)
c  nvar=4 в u теперь J(u*), соответствующее xp1
      endif
c  nvar=2 в u теперь u*, соответствующее xp1
      endif
      endif
      return
      END
