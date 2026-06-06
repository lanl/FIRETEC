       !subroutine dragb(xv,il,iu,jl,ju,lls,nvp)
       subroutine dragb()
       use turba
       use updatedFields
       use gridsetup
       use msga
       use constants
       use xvo

       Implicit None

       !JAS 3/7/06 added explicit declarations to comply with implicit none
       !integer,intent(in) :: il,iu,jl,ju,lls,nvp

       !real xv(il:iu,jl:ju,lls,nvp)

       !JAS 3/7/06 added explicit declarations to comply with implicit none
       integer :: i,j,k
       real :: drag,sp,av,alphas,rpi
       rpi=1./3.14159
       do k=1,lfuel
       do j=1,mp
       do i=1,np
          sp=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
          if (idrag.eq.1) then
            alphas=rhof(i,j,k)/rhomicro(i,j,k)
            ! here rho was removed by FP (wrong)
            ! FP replace ss by sizescale 
            drag=-sqrtk(i,j,k)*xvb(i,j,k,6)/sb(i,j,k)
     &        -4.*alphas*rpi/sizescale(i,j,k)*sp*xvb(i,j,k,6)*.5
     &        +4.*alphas*rpi/sizescale(i,j,k)*sp*xvb(i,j,k,5)*.25*.5
          else if (idrag.eq.2) then
            av=2./sizescale(i,j,k)*rhof(i,j,k)/rhomicro(i,j,k)  !correct for cylinder if ss is radius
            !FP added 0.93 as in Deardorf
            drag=-0.93*sqrtk(i,j,k)*xvb(i,j,k,6)/sb(i,j,k)
c NB: there was a factor 4 in 2009 IJWF
     &        + 2.0*0.5*cd(i,j,k)*av*sp*xvb(i,j,k,5)
     &        - 2.0*0.5*cd(i,j,k)*av*sp*xvb(i,j,k,6)
          endif
          fkb(i,j,k)=fkb(i,j,k)+2.*drag*dt

c          drag=-0.93*sqrtk(i,j,k)*rkb
c     &         /(sb(i,j,k))*rho
c     &      -4.*alphas*rpi/ss*rho*sp*rkb*.5
c     &      +4.*alphas*rpi/ss*rho*sp*rka*.25*.5
c     &         -2.*alphas*rpi/ss*rho*sqrtk(i,j,k)*rkb  
c     &         +2.*alphas*rpi/ss*rho*sqrtk(i,j,k)*rka*.5
c 3/4     &   *2./(sb(i,j,k)+sa(i,j,k))*rho
c 3/4     .       -0.375*cd(i,j,k)*rho*rkb*sqrtk(i,j,k)*alphas/
c     .       sc(i,j,k)+
cc     .       0.375*cd(i,j,k)*rho*(sqrtk(i,j,k)*rka)*
c     .       0.375*cd(i,j,k)*rho*(sqrtk(i,j,k)*rka+sp*sp*sp)*
c     .       alphas/sb(i,j,k)
       enddo
       enddo
       enddo
       do k=lfuel+1,l
       do j=1,mp
       do i=1,np
          if (idrag.eq.1) then
            drag=-sqrtk(i,j,k)*xvb(i,j,k,6)/sb(i,j,k)
          else if (idrag.eq.2) then
            drag=-0.93*sqrtk(i,j,k)*xvb(i,j,k,6)/sb(i,j,k)
          endif
          fkb(i,j,k)=fkb(i,j,k)+2.*drag*dt
       enddo
       enddo
       enddo
   
       return
       end
