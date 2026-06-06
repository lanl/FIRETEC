      subroutine metryc()
      use metryic
      use gridsetup
      use msga
      use met2
      use io
      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k
      real :: f,sigma0
      real,external :: zcart,gdeform
      
      real,allocatable:: go(:,:)
      zcrdata=(/0.0 ,   2.0000,    4.000,     6.00,       8.000,
     +        10.000,  14.000 ,   18.000,    22.00 ,     26.00 ,
     +        30.00 ,  34.00  ,   40.00 ,    50.00 ,     70.00 ,
     +       100.00 , 130.00  ,  160.00 ,   200.00 ,    250.00 ,
     +       300.00 , 350.00  ,  450.00 ,   550.00 ,    750.00 ,
     +       950.00 ,1150.0   , 1400.0  ,  1700.0  ,   2000.0  ,
     +      2400.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.   ,
     +         0.,      0.    ,    0.   ,     0.   ,      0.    /)
      npoints=l
      allocate (go(1-ih:np+ih,1-ih:mp+ih))
 
! deformation of sigma coordinate is done using polynomial fit
! aa3*sigma**3+aa2*sigma**2+aa1*sigma=sigma0
! domain for both sigma sysyems is 0 <= sigma  <= zb,
!                                  0 <= sigma0 <= zb
! now aa1 is defined in gridlist (optional value and default is 0.1)
       aa1m2=aa1
!                          ! aa1 determines compression of sigma at surface
      !aa1=1.0              ! aa1 can vary from 0 to 1 (no stretching)
!     aa1=0.0
c     aa1=1./10. !used prior to 8/17/01
      !aa1=0.1 !/10.
      !
      !aa1=2./10.
                           ! if aa1=0 spline data will be used
      f=0.                 ! 0 <= f <= 1
                           ! f=0, pure cubic fit
                           ! f=1, pure quadratic fit
 
      aa2=f*(1-aa1)/zb           ! (don't change this scaling constraint)
      aa3=(1-aa2*zb-aa1)/zb**2   ! (don't change this scaling constraint)
      if (mpi_rank.eq.0) write(6,*) 'aa1 is ', aa1 
      if(aa1.eq.0.)then          ! use spline data instead of deformation
         do i=1,npoints          ! polynomial
            zdata(i)=z(i)*zb/z(npoints)
         enddo

c get spline coefficients for ginverse and dgdsigma
 
         call spline(zdata,zcrdata,npoints,99.e31,99.0e31,zcoeff)
      endif
 
      do k=1,l
        sigma0=gdeform(z(k),0)
        gmul(k)=(zb-sigma0)/gdeform(z(k),1)
        do j=1,mp
        do i=1,np
           gi(i,j,k)=zb/(zb-zs(i,j))/gdeform(z(k),1)
           h(i,j,k)=1./gi(i,j,k)
           go(i,j)=zb/(zb-zs(i,j))
        enddo
        enddo
      !     write(6,*) zb, gi(1,1,k), gdeform(z(k),1)
      enddo

      call updated(gi,gi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      call updated(h,h,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      call updated(go,go,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0)
      
      do 41 j=1,mp
      do 41 i=1+leftdedge,np-rightdedge                         !added d rrl
      c13(i,j)=.5*dxi*(1./go(i+1,j)-1./go(i-1,j))*go(i,j)
   41 continue

      if (leftdedge.eq.1) then                         !added d rrl
      do j=1,mp
         c13(1,j)=0.+ibcx*.5*dxi*(1./go(2,j)-1./go(-1,j))*go(1,j)
      end do
      end if
      if (rightdedge.eq.1) then                         !added d rrl
         do j=1,mp
            c13(np,j)=0.+ibcx*.5*dxi*(1./go(np+2,j)-1./
     .           go(np-1,j))*go(np+1,j)
         end do
      end if
  
      do 51 i=1,np
      do 51 j=1+botdedge,mp-topdedge                         !added d rrl
   51 c23(i,j)=.5*dyi*(1./go(i,j+j3)-1./go(i,j-j3))*go(i,j)

      if (botdedge.eq.1) then                         !added d rrl
      do  i=1,np
         c23(i,1)=0.+ibcy*.5*dyi*(1./go(i,1+j3)-1./go(i,-j3))*
     .        go(i,1)
      end do
      end if
      if (topdedge.eq.1) then                         !added d rrl
      do i=1,np
         c23(i,mp)=0.+ibcy*.5*dyi*(1./go(i,mp+1+j3)-1./go(i,mp-j3))*
     .        go(i,mp+1)
      end do
      end if

      call updated(c13,c13,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0)
      call updated(c23,c23,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0)
      deallocate (go)
      if (mpi_rank.eq.0) open(unit=77,file='gi.dat',
     +                    form='unformatted',status='unknown')
       call writeio(gi(1-ih,1-ih,1),77,1-ih,np+ih,1-ih,mp+ih,l,0)
      if (mpi_rank.eq.0) close(77)

      return
      end

      real function zcart(sigmax,i,j)
! sigma is sigma coordinate, zcart is cartesian vertical coordinate
      use metryic
      use met2

      Implicit None

      integer :: i,j
      real :: sigmax
      real,external :: gdeform

      zcart=gdeform(sigmax,0)*(zb-zs(i,j))/zb+zs(i,j)
      return
      end function zcart
      real function zcart2(sigmax,i,j)
! sigma is sigma coordinate, zcart is cartesian vertical coordinate
      use metryic
      use met2

      Implicit None

      integer :: i,j
      real :: sigmax
      real,external :: gdeform

      zcart2=gdeform(sigmax,0)*(zb-zsio(i,j))/zb+zsio(i,j)
      return
      end function zcart2
!*************************end function zcart*****************************!

      real function gdeform(sigmax,iflag)
! iflag=0, compute gdeform(sigma)=sigma0,
! iflag=1, return derivative of gdeform
! sigma0=aa3**sigma**3+aa2*sigma**2+aa1*sigma   ! cubic polynomial fit
      use met2

      Implicit None

      integer :: iflag
      real :: sigmax
      real :: answer

      if(aa1m2.eq.0)then           ! use spline
         call splint(zdata,zcrdata,zcoeff,npoints,sigmax,answer,iflag)
         gdeform=answer
      else                       ! use cubic polynomial
         if(iflag.eq.0)gdeform=aa3*sigmax**3+aa2*sigmax**2+aa1m2*sigmax
         if(iflag.eq.1)gdeform=3*aa3*sigmax**2+2*aa2*sigmax+aa1m2
      endif
      return
      end function gdeform
!*************************end function gdeform*****************************!

 
! subroutines spline and splint are from "Numerical Recipes". Splint
! is modified to give derivatives as well as interpolated values
 
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
 
     *u(i-1))/p
 11   continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
 12   continue
      return
      END

 
      SUBROUTINE splint(xa,ya,y2a,n,x,y,kderivative)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo,kderivative
      REAL a,b,h
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else 
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      if(kderivative.eq.0)then
         y=a*ya(klo)+b*ya(khi)+( (a**3-a)*y2a(klo) + (b**3-b)*y2a(khi) )
     +    *(h**2)/6.
      else
         y=(ya(khi)-ya(klo))/h
     +    -(  (3*a**2-1)*y2a(klo)-(3*b**2-1)*y2a(khi)  )*h/6
      endif
c234567
      return
      END
