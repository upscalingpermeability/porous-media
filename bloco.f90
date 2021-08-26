      SUBROUTINE blocopoisson2(MAT,MX,MY)
!
        USE VARIAVEIS
      implicit none
!
      integer :: MX,MY,nt,nbt!,nf
      integer :: kd,kxp,kxm,kyp,kym
      real(8) :: xx,yy,xl,yl,hx,hy,rxy,ryx,hxy
      real(8) :: xkim,xkip,xkjm,xkjp,xk,q
      real(8) :: pleft,pright
      integer :: i,j,ip,im,jp,jm,k,kim,kip,kjm,kjp
      real(8), dimension(:,:), allocatable :: a
      real(8), dimension(:),   allocatable :: f,p
      REAL(8) :: MAT(MX*MY),AUX
!
!     dimensoes
!
      xl=1.d0
      yl=1.d0
!
!     numero de blocos
!
      nt=MX*MY
!
      hx=xl/MX
      hy=yl/MY
!
      hxy=hx*hy
      rxy=2.d0*hx/hy
      ryx=2.d0*hy/hx
!
!     fonte
!
!      nf=1
!
!     banda
!      
      nbt=2*MX+1
!
!     inicializando
!
      allocate(a(nt,nbt))
      allocate(f(nt))
      allocate(p(nt))
!
      a=0.d0
      f=0.d0
      p=0.d0
!
!...  le o campo de permeabilidade
!      
!      open(unit=15,file="campo.dat",status="unknown")
!
!      do i=1,nt
!      read(15,*) k,MAT(k)
!      end do      
!
      pleft =1.d0
      pright=0.d0
!
      xx=-hx/2.d0
      yy=-hy/2.d0
!
!     montagem do sistema 
!
!      write(*,*) "..... Montando o Sistema"
!
      do j=1,MY
!
      yy=yy+hy
!
      jp=j+1
      jm=j-1
!
      xx=0.d0
      do i=1,MX
!
      xx=xx+hx
!
      ip=i+1
      im=i-1
!
      k=i+MX*(j-1)
      kip=k+1 
      kim=k-1
      kjp=k+MX
      kjm=k-MX
!
!     condicoes
!
      if(i.eq. 1) kim = k
      if(i.eq.MX) kip = k
!
!     fluxo nulo
!
      if(j.eq. 1) kjm = k
      if(j.eq.MY) kjp = k
!
!     permeabilidades
!
      xk  =MAT(k)
      xkip=MAT(kip)
      xkim=MAT(kim)
      xkjp=MAT(kjp)
      xkjm=MAT(kjm)
!
!     media harmonica
!
      xkip=ryx*(xkip*xk)/(xkip+xk)
      xkim=ryx*(xkim*xk)/(xkim+xk)
      xkjp=rxy*(xkjp*xk)/(xkjp+xk)
      xkjm=rxy*(xkjm*xk)/(xkjm+xk)
!
!     estrutura em banda
!
      kd=MX+1
      kxp=kd+1
      kxm=kd-1
      kyp=kd+MX
      kym=1      
!
      if(i.eq. 1) kxm = kd
      if(i.eq.MX) kxp = kd
!
!     fluxo nulo
!
      if(j.eq. 1) kym = kd
      if(j.eq.MY) kyp = kd
!      
!     diagonal    
!
      a(k,kd) = a(k,kd)-(xkjp+xkjm+xkip+xkim)
!
!     vizinhos em x
!
      a(k,kxp) = a(k,kxp)+xkip
      a(k,kxm) = a(k,kxm)+xkim
!
      if(i.eq.1) then
      a(k,kxm) = a(k,kxm)-2.d0*xkim
      f(k)=f(k)-2.d0*xkim*pleft
      end if
!
      if(i.eq.MX) then
      a(k,kxp) = a(k,kxp)-2.d0*xkip
      f(k)=f(k)-2.d0*xkip*pright
      end if
!
!     vizinhos em y
!
      a(k,kyp) = a(k,kyp)+xkjp
      a(k,kym) = a(k,kym)+xkjm
!
!     termo de fonte
!
      f(k)=f(k)+hxy*q(xx,yy)
!
      end do ! i
!
      end do ! j
!
!      resolvendo o sistema
!    
!      write(*,*) "..... Resolvendo o Sistema"
!
      call solve(a,f,nt,nt,nbt)
      p=f    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTEGRAL NO BORDO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        U=0.D0
        AUX=2.D0/HX
        DO J=1,MY
           I=J*MX
           U=U+F(I)*AUX*MAT(I)*HY
        END DO
!      
!      imprimindo a solucao
!     
!      write(*,*) "..... Concluido"
!
!      do j=1,MY
!      do i=1,MX
!      k=i+MX*(j-1)
!      write(100,*) hx*(i-1)+hx/2.d0,hy*(j-1)+hy/2.d0,p(k)
!      end do
!!      
!      end do
!!
!      yy=0.d0
!      do j=1,MY
!      xx=0.d0
!      yy=yy+hy
!      do i=1,MX
!      k=i+MX*(j-1)
!      xx=xx+hx
!      write(200,*) xx-hx,yy-hy,p(k)
!      write(200,*) xx   ,yy-hy,p(k)
!      write(200,*) xx   ,yy   ,p(k)
!      write(200,*) xx-hx,yy   ,p(k)
!      write(200,*) xx-hx,yy-hy,p(k)
!      write(200,*)
!      write(200,*)
!      end do
!      end do
!
    end SUBROUTINE blocopoisson2
!
!----------------------------------------------------------------------
!
      function q(xx,yy)
      implicit none
!
      real(8) :: xx,yy,q,pi,dpi
!     
      pi =4.d0*datan(1.d0)
      dpi=2.d0*pi
! 
!      q=-8.d0*pi*pi*dsin(dpi*xx)*dcos(dpi*yy)
!      q= 2.d0*pi*pi*dsin(pi*xx)*dcos(pi*yy)
!
      q=0.d0
!
      end function
!
!******************************************************************************
!                                                 
!     Objetivo: resolver um sistema linear em banda
!
!----------------------------------------------------------------------
!
      subroutine solve(s,r,ns,idm,lb)                                        
      implicit real *8(a-h,o-z)                                          
      dimension s(idm,*),r(*)                                           
      ls=(lb-1)/2                                                        
      l1=ls+1                                                            
      l2=ls+2                                                            
      do 300 n=1,ns                                                      
      i=n                                                                
      do 290 l=2,l1                                                      
      i=i+1                                                              
      if(i.gt.ns) go to 300                                              
      ll=n-i+l1                                                          
      if(s(i,ll).eq.0.) go to 290                                        
      c=s(i,ll)/s(n,l1)                                                  
      j=ll-1                                                             
      do 270 k=l1,lb                                                     
      j=j+1                                                              
      if(s(n,k).eq.0.) go to 270                                         
      s(i,j)=s(i,j)-c*s(n,k)                                             
  270 continue                                                           
      r(i)=r(i)-c*r(n)                                                   
  290 continue                                                           
  300 continue                                                           
      r(ns)=r(ns)/s(ns,l1)                                               
      n=ns                                                               
  350 n=n-1                                                              
      if(n.le.0) go to 500                                               
      l=n                                                                
      do 400 k=l2,lb                                                     
      l=l+1                                                              
      if(l.gt.ns) go to 450                                              
      r(n)=r(n)-s(n,k)*r(l)                                              
  400 continue                                                           
  450 r(n)=r(n)/s(n,l1)                                                  
      go to 350                                                          
  500 return                                                             
      end             
