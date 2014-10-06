      subroutine analise(kkoef,cal,ca,nkoef)
      IMPLICIT REAL*8(a-b,d-h,O-Z)
      COMPLEX*16 ca[reference](1),cal
c  
      br=dreal(ca(1))
      bi=dimag(ca(1))
      cr=dreal(ca(2))
      ci=dimag(ca(2))
      dr=dreal(ca(3))
      di=dimag(ca(3))
      er=dreal(ca(3))
      ei=dimag(ca(3))
      fr=dreal(ca(3))
      fi=dimag(ca(3))
      ar=dreal(cal)
c---------------------------------------------------
      if(ar.gt.0) then
      write(6,6)  
      write(*,6)  
   6  format(3x,'Оcновной режим теряет устойчивость при дельта > 0. ')
c Real(a)>0
      else
c Real(a)<0
      write(6,7)  
      write(*,7)  
   7  format(3x,'Оcновной режим теряет устойчивость при дельта < 0. ')
      endif
c----------------------------------------------------
      r0=1 
c 2)===================
      r=br+dr
      r1=br-dr
      r2=r-cr-er+fr
      r3=r-cr-er-fr
      call fanalis(r,r1,r2,r3,r0,2)      
c 6)===================
      r=br
      r1=-(br-cr)
      r2=br-dr
      r3=br-er
      r0=1
      call fanalis(r,r1,r2,r3,r0,6)
      RETURN
      END
c==============================================================
      subroutine fanalis(r,r1,r2,r3,r0,i)
      IMPLICIT REAL*8(a-b,d-h,O-Z)
      if(r.lt.0) then
      write(6,1) i 
      write(*,1) i
      write(8,1) i	  
   1  format(3x,'Ответвление режима ',i2,' происходит в сверхкритиче',
     *'скую область.')
      else
      write(6,2) i
      write(*,2) i
      write(8,2) i
   2  format(3x,'Ответвление режима ',i2,' происходит в докритическую',
     *' область и он неустойчив.')
      endif
      if(r.lt.0) then
      if((r1.lt.0).and.(r2.gt.0).and.(r3.gt.0).and.(r0.gt.0)) then 
      write(6,3)
      write(8,3) 	  
      write(*,3) 
   3  format(3x,'Этот режим устойчив. ')
      else
      write(6,4) 
      write(8,4) 
      write(*,4) 
   4  format(3x,'Этот режим неустойчив. ')
      endif
      endif 
      return
      end
c============================================
