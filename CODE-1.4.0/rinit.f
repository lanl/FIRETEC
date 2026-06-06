      subroutine rinit(tcv)
      use gridsetup
      use xvo
      use xvi
      use xve
      use pres
      use relax
      use weights
      use constants
      use metryic
      use msga
      use fireteca
      use radiation !, only: irad,icallrad ! KOO  
      use ignite
      use fueldrag
      Implicit None

      !JAS 3/6/06 added explicit declarations to comply with implicit none
      real :: tcv
      integer :: i,j,k,iw,kv
      real,external :: zcart
  

       
      xvb=0.0
      xe=0.0
      xv=0.0
      ! time integration parameters
      tcv=0.
      wt(1)=0.3/real(nts) 
      wtf(1)=0.0
      do iw=2,nts
       wt(iw)=1.0/real(nts) 
       wtf(iw)=1.0
      enddo
      wt(nts+1)=0.7/real(nts) 
      wtf(nts+1)=1.0
      ! grid definition
      do i=1-ih,n+ih
       x(i)=(i-1)*dx-((n+1)/2.0-1)*dx
      enddo
      do j=1-ih,m+ih
       y(j)=(j-1)*dy-((m+1)/2.0-1)*dy
      enddo
      if(ibctopbot.eq.1) then
         do k=1,l+1
           zedge(k)=(k-1)*dz
         enddo 
         zb=zedge(l+1)
         do k=1,l
           z(k)=(k-1)*dz+0.5*dz
c           if(mpi_rank.eq.0)
c     .        write(6,*) 'k,zcart(z(k),1,1),zb,z(k)',k,zcart(z(k),1,1),zb,z(k)
          enddo
      else !ibctopbot.eq.0 (not used)
        do k=1,l
          z(k)=(k-1)*dz
        enddo
        zb=z(l)
      endif
      ! FP conditions for temporal update of xe in compress (ixevariation.eq.1)
      if(iwindfieldin.eq.1.or.uswitch.eq.1.or.vswitch.eq.1)
     + ixevariation=1
      ! FP redefinition of is,ie,js,je when equals to 0 (meaning the whole domain is considered
      ! for windfield stuff 
      if (iwindfieldin.eq.1.or.iwindfieldout.eq.1) then
          ! windfield in our out is done by default on the whole domain
          ! if ie,is,je,js have not been defined first
          if (ie.eq.0) ie=n 
          if (is.eq.0) is=1 
          if (je.eq.0) je=m 
          if (js.eq.0) js=1 
      endif

      !compute coordinate transformation related matrices
      call topo
      call metryc()

      !compute absorber at the sides and top
      call tinit()
      !compute rho, theta, pressure profiles for ambient variables xe
      ! ambient data at height zgroundref
      call setRhoThetaP()
      ! here, we define rhof,...cpsolid... when no restart or restart
      ! with new fuelfiles (irst.eq.2)
      if (irst.ge.1.or.iwindfieldin.eq.1) then !: definition of xvb
         call startFromFile()
         if (iwindfieldin.eq.1.and.irst.eq.0) itrestart=0
      endif
      if(irod.eq.1)then
         call ign_setup
         call rinitfire()
         call ign_cleanup
      endif
      if(irad.gt.0)  call frad_init() !require actualfueldepth
      ! computation of xe(1,2,3)
      call setVelocityProfile()
      ! here we set a large scale pressure gradient
      ! NB : xe(1,2) should be defined first!
      if (ilspgf.ge.1) call defineLargeScalePressureGradientForce() 
      ! computation of xe(5,6),sa, cd...
      if(iturb.ge.1) then
        call rinitturb()
        call init_fueldrag()
      endif
      ! definition of xvb:
      ! no restart:
      if (irst.eq.0.and.iwindfieldin.eq.0) then
        itrestart=0
        do k=1,l
          do j=1,mp
            do i=1,np
              xv(i,j,k,1)=xe(i,j,k,1)
              xv(i,j,k,2)=xe(i,j,k,2)
              xv(i,j,k,3)=xe(i,j,k,3)
              xv(i,j,k,4)=xe(i,j,k,4)
              xv(i,j,k,5)=xe(i,j,k,nv)
              do kv=1,nv 
                xvb(i,j,k,kv)=xe(i,j,k,kv)
              enddo
            enddo
          enddo
        enddo
        if (iperturb.ge.1) call setPerturbation() ! perturbation of xv and xvb
      endif ! irst.eq.0.andiwindfieldin.eq.0
      ! definition of higrad variable xv
      do k=1,l
        do j=1,mp
          do i=1,np
            xv(i,j,k,1)=xvb(i,j,k,1)
            xv(i,j,k,2)=xvb(i,j,k,2)
            xv(i,j,k,3)=xvb(i,j,k,3)
            xv(i,j,k,4)=xvb(i,j,k,4)
            xv(i,j,k,5)=xvb(i,j,k,nv)
          enddo
        enddo
      enddo

c updated of xvb and xv from xe
      do kv=1,nv
       call updated(xvb(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),
     &              np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
       if(kv.le.4)then
        call updated(xv(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),
     &              np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
       elseif(kv.eq.nv)then
        call updated(xv(1-ih,1-ih,1,5),xe(1-ih,1-ih,1,kv),
     &              np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
       endif 
      enddo
      
 

      if(islip.eq.1) then  ! not used for firetec
       do j=1,mp
        do i=1,np
         xv(i,j,1,1)=0.
         xv(i,j,1,2)=0.
         xvb(i,j,1,1)=0.
         xvb(i,j,1,2)=0.
         xe(i,j,1,1)=0.
         xe(i,j,1,2)=0.
        enddo
       enddo
       endif !if(islip.eq.1)
       return
      end ! end subroutine rinit




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc setVelocityProfile sets the inital values of xv(1,2,3) and
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine setVelocityProfile()
      use xve
      use xvo
      use gridsetup
      use msga
      use metryic
      use constants
      use turba
      Implicit none
      integer i,j,k,kv
      real,external :: zcart
      real::cosc,sinc,gam,zla,one  ! cosi and sin of ekman spiral for icorio.eq.2
      if (uswitch.eq.2.or.vswitch.eq.2) call empiricalProfileFromFuelData()
      if(icfmeflag.eq.1)then
        if(mpi_rank.eq.0)write(*,*)'going into readicfmew'
        call readicfmew(xe,1-ih,np+ih,1-ih,mp+ih,l,nv)
        call readicfmew(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)
        if(mpi_rank.eq.0)write(*,*)'exiting readicfmew'
      else ! icfmeflag.eq.0 
       if (icorio.eq.2) then
        ! special definition of xe (used for initial condition only,
        ! cause designed to be used with cyclic bc): xe used the ekman
        ! spiral as initial prof
        do k=1,l
         do j=1,mp
           do i=1,np
              zla=zcart(z(k),i,j) - zgroundref
              gam = 3.14/zab*zla
              cosc=1.0-exp(-gam)*cos(gam)  !u(z)/ug from ekman
              sinc=exp(-gam)*sin(gam)  !v(z)/ug from ekman
              one = sqrt(cosc*cosc+sinc*sinc)
              cosc = cosc/one
              sinc = sinc/one
              !NB uswitch.eq.1 is defined in xevariation.f
              xe(i,j,k,1)=0.0
              xe(i,j,k,2)=0.0
              if(uswitch==0)then
                  xe(i,j,k,1)=u0*cosc*xe(i,j,k,nv)
                  xe(i,j,k,2)=u0*sinc*xe(i,j,k,nv)
               else if (uswitch==2) then
                  xe(i,j,k,1)=u0*cosc*xe(i,j,k,nv)
     &         *uprofile(i,j,k)      !u0 represents speed at zu
                xe(i,j,k,2)= u0*sinc*xe(i,j,k,nv)
     &         *uprofile(i,j,k)      !u0 represents speed at zu
               endif
              if(vswitch==0) then
                  xe(i,j,k,1)=xe(i,j,k,1)-v0*sinc*xe(i,j,k,nv)
                  xe(i,j,k,2)=xe(i,j,k,2)+v0*cosc*xe(i,j,k,nv)
              else if (vswitch==2) then
                xe(i,j,k,1)= xe(i,j,k,1)-v0*sinc*xe(i,j,k,nv)
     &         *uprofile(i,j,k)      !v0 represents speed at zu
                xe(i,j,k,2)= xe(i,j,k,2)+v0*cosc*xe(i,j,k,nv)
     &         *uprofile(i,j,k)      !v0 represents speed at zu
              endif
          enddo
        enddo
       enddo
      else ! icorio.ne.2
        do k=1,l
         do j=1,mp
           do i=1,np
              !NB uswitch.eq.1 is defined in xevariation.f
              if(uswitch==0)then
                   xe(i,j,k,1)=u0*xe(i,j,k,nv)
               else if (uswitch==2) then
                   xe(i,j,k,1)=u0*xe(i,j,k,nv)
     &         *uprofile(i,j,k)      !u0 represents speed at zu
               endif
              if(vswitch==0) then
                  xe(i,j,k,2)=v0*xe(i,j,k,nv)
              else if (vswitch==2) then
                xe(i,j,k,2)= v0*xe(i,j,k,nv)
     &         *uprofile(i,j,k)      !v0 represents speed at zu
              endif
          enddo
        enddo
      enddo
      endif !icorio.ne.2
      endif !icfmeflag.eq.0
       !TODO : is potflow working or not?
       if (ipotflow.eq.1) call potflow(xe,1-ih,np+ih,1-ih,mp+ih,l,nv)
       do kv=1,3
       call updated(xe(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),
     &              np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
       enddo

       do k=1,l
         if (mpi_rank.eq.0)     write (6,*) 'initial wind profile',
     +         zcart(z(k),1,1)-zs(1,1), sqrt(xe(1,1,k,1)**2+
     +          xe(1,1,k,2)**2)/xe(i,j,k,nv)
       enddo 
       return
      end !subroutine setVelocityProfile
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c defineLargeScalePressureGradientForce
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine defineLargeScalePressureGradientForce()
      use xve
      use xvo
      use pres
      use gridsetup
      use msga
      use metryic
      use constants
      use lspgf
      Implicit none
      integer i,j,k,ierr
      real,external :: zcart
      real :: zla
      real :: gam
      real :: z1,z2,z12,zcoef
      real :: intsinthetal
      real ::dzcell 
      ! computation of geostrophic variables : NB : xe(:,:,l) should be
      ! constant
      rhoug=xe(1,1,l,1)
      rhovg=xe(1,1,l,2)
      ! quick local check that top xe is constant
      do j=1,mp
           do i=1,np
              if(rhoug.ne.xe(i,j,l,1).or.rhovg.ne.xe(i,j,l,2))  then
                  write(*,*) 'ERROR: top xe not compatible with lspgf'
                  stop
              endif
           enddo
      enddo 
      ! compute wind direction
      cosg = rhoug/sqrt(rhoug**2 + rhovg**2)
      sing = rhovg/sqrt(rhoug**2 + rhovg**2)
      ! lspgf time step
      tau = frqlspgf * dt
     
      ! definition of massFlow integrals
      if (abs(rhoug)*m >= abs(rhovg)*n) then ! massFlow computed along y axis (as if wind aligned with xaxis)
         iwindx = 1
         nMassFlux = m
         if (mpi_rank.eq.0) write (6,*) 'mass flux computed along yaxis
     +      over ', nMassFlux, ' cells'
      else  ! massFlux computed along x axis (as if wind aligned with y axis)
         iwindx = 0
         nMassFlux = n
         if (mpi_rank.eq.0) write (6,*) 'mass flux computed along xaxis
     +      over ' , nMassFlux, ' cells'
      endif
      allocate(massFluxTemp(nMassFlux))
      allocate(massFlux(nMassFlux))
      allocate(massFluxTarget(nMassFlux))

      massFluxTemp = 0.0
      massFlux = 0.0
      massFluxTarget = 0.0      
      intsinthetal = 0.0 ! integrated value of sintheta over the vertical 
      ! definition of sintheta (vertical profile of the lspgf extrapolated from ekmann spiral) 
      do k=1,l
       do j=1,mp
        do i=1,np
           zla=zcart(z(k),i,j) - zgroundref !FP took zgroundref as a reference for pressure gradient
c         ! LSPGF IS NOT COMPATIBLE WITH ROTATED GRAVITY BECAUSE W IS ) AT THE TOP
c         ! when gravity is rotated, zla vary with slopeangle
c         ! and slope azimuth:
           gam = 3.14/zab*zla
           sintheta(i,j,k)=0.5*exp(-gam)*sin(gam)/
     &        sqrt(1.0+exp(-2.0*gam)-2.0*exp(-gam)*cos(gam))
             dzcell = zcart(zedge(k+1),i,j) - zcart(zedge(k),i,j)
           if (izlspgf.eq.0) then
             intsinthetal = intsinthetal + sintheta(i,j,k) * dzcell/(n*m)
           endif
        enddo
       enddo
      enddo
      if (izlspgf.eq.0) then  ! mean of integral over z of sintheta*dz
        call mpi_allreduce(intsinthetal,intsintheta,1,mpi_real,
     +                  mpi_sum,mpi_comm_world,ierr)
      else !izlspgf.eq.1  : intsintheta=sintheta(zu)
          gam = 3.14/zab*(zu-zgroundref)
          intsintheta =0.5*exp(-gam)*sin(gam)/
     &        sqrt(1.0+exp(-2.0*gam)-2.0*exp(-gam)*cos(gam))
      endif

      ! definition of flspgf and sinthetaf
      if (irst.eq.0) then ! definition of flspgf (n*m) array
           flspgf = 8.0*fcor3  ! initial value 
           do k=1,l
             do j=1,mp
               do i=1,np
                  sinthetaf(i,j,k) = sintheta(i,j,k) * flspgf (i,j) 
               enddo
             enddo
           enddo
      else ! compute flspgf from sinthetaf in case of restart (from io.f)
        do j=1,mp
           do i=1,np
           flspgf(i,j) = sinthetaf(i,j,1)/sintheta(i,j,1)
           enddo
        enddo
      endif ! irst.eq.0
      if (mpi_rank.eq. 0) write(6,*) 'ilspgf: flspgfini',flspgf(1,1)
      

      ! test of value of zu (should be between bottom and top of the
      ! domain for izlspgf=1
      if (izlspgf.eq.1) then
         if (zu<=zcart(zedge(2),i,j)-zs(i,j).or.(zu>=zcart(zedge(l),i,j)-zs(i,j))) then
            write(6,*) 'izlspgf=1 incompatible with zu that should be
     +       above cell 1 and below cell l-1'
            stop
         endif
      endif
      ! computation of inital massFluxTemp on all sub domains (reduced with a mpi sum in array
      ! massFluxTarget later). Something similar is done in largeScalePGF to compute current value of
      ! massFlux (to compare it to the target define here).
      do k=1,l
       do j=1,mp
         do i=1,np
            z1=zcart(zedge(k),i,j)-zs(i,j)
            z2=zcart(zedge(k+1),i,j)-zs(i,j)
            if (izlspgf.eq.0) then ! massFlux on the whole domain
              call addToMassFluxTemp(i, j, k, 1, z2-z1)
            else ! massFlux at a given height (zu) 
              ! compute the contributon of cell (i,j,k)  at height zu
              ! z1 is cell center of cell k-1 (or domain bottom)
              if (k.eq.1) then
                z1=0.0
              else
                z1=zcart(z(k-1),i,j)-zs(i,j)
              endif
              z2=zcart(z(min(k+1,l)),i,j)-zs(i,j)
              if (z1 <= zu .and. zu <= z2) then
              ! cell i,j,k should contribute to massFlux
                z12 = zcart(z(k),i,j)-zs(i,j)
                ! zcoef is the weight of the cell to massFlux
                if (zu>=z12) then ! zu is between z12 and z2
                   zcoef = (z2 - zu) / (z2 - z12)
                else ! zu is between z1 and z12
                   zcoef = (zu - z1) / (z12 - z1)
                endif
                call addToMassFluxTemp(i, j, k, 1, zcoef)
               endif
            end if  ! izlspgf.eq.1


         enddo
       enddo
      enddo
      call mpi_allreduce(massFluxTemp,massFluxTarget,nMassFlux,mpi_real,
     +                  mpi_sum,mpi_comm_world,ierr)
       targMassFlux=0.0
       do j=1,nMassFlux
         targMassFlux = targMassFlux + massFluxTarget(j)
c          write(6,*) 'massflux initial:', massFluxTarget(j)
       enddo
        if (mpi_rank.eq.0)
     + write(6,*) 'total massflux initial:', targMassFlux
       
       return
      end !subroutine defineLargeScalePressureGradientForce

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c empiricalProfileFromFuelData computes a profile uprofile(k), normalized
c  at height zu.
c uprofile can be used to set initial profile in rinit.f for u and v with a value
c uO,v0 at height zu above ground
c based on empirical formulation of Kaimal and Finnigan 1994 and Raupach 1993
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine empiricalProfileFromFuelData()
      use metryic
      use gridsetup
      use fireteca
      use xve
      use turba
      use msga
      integer :: i,j,k,ift,ia,ja
       real,external :: zcart
       real :: c1,c2,c3,c4,LAIlow,LAIup,LAIlowt,LAIupt,LAIt,hmax,hmaxt,zk,ufueltopmax
         ! COMPUTATION OF LAIs and hmax
         hmax=0. !fuel max height (including actualfueldepth)
         LAIlow=0. ! LAI in cells k=1
         LAIup=0. !LAI in cells k>1
         if (ius.eq.0) ius=1
         if (iue.eq.0) iue=n
         if (jus.eq.0) jus=1
         if (jue.eq.0) jue=m
         do k=1,lfuel
          do j=1,mp
           ja = (mpos-1)*mp + j
           do i=1,np
            ia = (npos-1)*np + i
            if(sum(rhof(:,i,j,k)).gt.min_rhof.and.ia.ge.ius.and.ia.le.iue.and
     +            .ja.ge.jus.and.ja.le.jue) then
              if (k.eq.1) then 
                hmax=max(hmax,maxval(actualfueldepth(:,i,j,k)))
              else
                  hmax=max(hmax,zcart(zedge(k+1),i,j)-zs(i,j))
              endif
              if (k.eq.1) then
                LAIlow=LAIlow+sum(rhof(:,i,j,k))/(sum(rhomicro(:,i,j,k))/nfuel*sum(sizescale(:,i,j,k))/nfuel)
     +      *(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))
               else
                LAIup=LAIup+sum(rhof(:,i,j,k))/(sum(rhomicro(:,i,j,k))/nfuel*sum(sizescale(:,i,j,k))/nfuel)
     +      *(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))
              endif
           endif
           enddo
          enddo
         enddo
c           if (mpi_rank.eq.0) write (6,*) 'hmax',hmax,'LAI',LAI 
         call mpi_allreduce(hmax,hmaxt,1,mpi_real,mpi_max,
     +                   mpi_comm_world,ierror)
         call mpi_allreduce(LAIlow,LAIlowt,1,mpi_real,mpi_sum,
     +                   mpi_comm_world,ierror)
         call mpi_allreduce(LAIup,LAIupt,1,mpi_real,mpi_sum,
     +                   mpi_comm_world,ierror)
         LAIlowt=LAIlowt/real((iue-ius+1)*(jue-jus+1))
         LAIupt=LAIupt/real((iue-ius+1)*(jue-jus+1))
 
           if (mpi_rank.eq.0) write (6,*) 'hmaxt',hmaxt,'LAIlowt',LAIlowt,'LAIupt',LAIupt
     +       ,'in zone ',ius,iue,jus,jue
             if (LAIupt>0.1) then  !when a canopy is there: 
                LAIt=LAIupt+0.1*LAIlowt  
                ! so that the profile in the canopy is not too much affected by dense low vegetation
             else 
                LAIt=LAIupt+LAIlowt
             end if
          ! DEFINITION OF CONSTANT FROM LAIt
             c1=min(sqrt(0.003+0.15*LAIt), 0.3)      ! U*/uh
             c2=(1-exp(-sqrt(7.5*LAIt)))/sqrt(7.5*LAIt)      ! 1-d/h
             c3=c2*exp(-0.41/c1-log(2.)+0.5)          ! z0/h
             !c4 is coef within the canopy      
              !(Kaimal and Finnigan 1994 : range between 1.7 and 3.2 for LAI between 1 and 4) 
             if (LAIt.le.1) then
                   c4=1.7
             else if (LAIt.ge.4) then
                  c4=3.2
             else
                 c4=(3.2-1.7)/(4-1)*(LAIt-1)+1.7
             end if
           ! COMPUTATION OF UPROFILE at height zu (for normalization at height zu)
            if (zu.ge.hmaxt) then
               uprofilezu=c1/0.41*log((zu/hmaxt+c2-1)/c3)
               if (zu.lt.2.*hmaxt) then
                ufueltopmax=c1/0.41*log(c2/c3) !theoretical value at hmaxt         
                uprofilezu=uprofilezu+(1-ufueltopmax)*(2.-zu/hmaxt)**3
               end if
             else !under the canopy
              uprofilezu=exp(-c4*(1-zu/hmaxt))
              end if
         ! COMPUTATION OF PROFILE (NORMALIZED AT CANOPY HEIGHT AT THIS STAGE)
         do k=1,l
          do j=1,mp
           do i=1,np
            zk=zcart(z(k),i,j)-zs(i,j)
            if (zk.ge.hmaxt) then  !above fuel
             ! u(z)=c1/0.41*u(h)*ln((z/h+c2-1)/c3) where h= fueldepth
             ! this formula is valid for z>2h in the inertial sublayer
             uprofile(i,j,k)=c1/0.41*log((zk/hmaxt+c2-1)/c3)
             ! correction for the roughness sublayer (when actualfueldepth is greater than ztopcell/2)
             if (zk.lt.2.*hmaxt) then
                ufueltopmax=c1/0.41*log(c2/c3) !theoretical value at hmaxt         
                uprofile(i,j,k)=uprofile(i,j,k)+(1-ufueltopmax)*(2.-zk/hmaxt)**3
             end if
           else !under the canopy
              uprofile(i,j,k)=exp(-c4*(1-zk/hmaxt))
           end if
            ! correction top of first cell :
           if (k.eq.1.and.zcart(zedge(2),i,j)-zs(i,j).ge.hmaxt) then
              uprofile(i,j,1)=1.
           end if
           !RENORMALIZATION
            uprofile(i,j,k)=uprofile(i,j,k)/uprofilezu 
           !if (mpi_rank.eq.0) write (6,*) 'uprof',zk, uprofile(k)
         enddo
         enddo
         enddo
           !if (mpi_rank.eq.0) write (6,*) 'uprof at',zu,'=', uprofilezu

      end !subroutine empiricalProfileFromFuelData

c setPerturbation computes some perturbation to initialize relosved turbulence
c typically for cycli runs
c if iperturb=0 no perturb
c if iperturb=1 random perturbation of theta 
c if iperturb=2 pinwheel
c NB perturbation should be done on xvb and xv
      subroutine setPerturbation()
      use xvi
      use xvo
      use xve
      use gridsetup
      use msga
      use metryic
      use constants
      integer ::i,j,k,ia,ja,ipert,jpert,kpert
      real, external :: zcart
      real, allocatable :: randArray(:,:,:)
      real::rfluctuation ! fluctuation range
       if (mpi_rank.eq.0) write(6,*) 'iPerturbation=', iperturb
      if (iperturb.eq.1) then ! random perturbation on rhotheta
      allocate(randArray(n,m,l)) 
       rfluctuation=0.1
       if (mpi_rank.eq.0) write(6,*) 'perturbation magnitude=', rfluctuation
        ! that will be identical whatever the number of procs
       !call random_seed()
c RAND(0) does not use seed, (1 restart seed generator), other number is the seed itself
        do k=1,l
         do j=1,m
           do i=1,n
       !if (mpi_rank.eq.0) write(6,*) i,j,k
            randArray(i,j,k)=rand();  
            enddo
          enddo
        enddo
        do k=1,l
         do j=1,mp
           do i=1,np
             ia = (npos-1)*np + i
             ja = (mpos-1)*mp + j 
             xv(i,j,k,4)=xv(i,j,k,4) + rfluctuation*(randArray(ia,ja,k)-0.5)
             xvb(i,j,k,4)=xv(i,j,k,4)
             if ((mpi_rank==0).and.(k.eq.1).and.(i.le.2).and.(j.le.2)) 
     +            write(6,*) 'random theta for i=',i,' j=', j,' is '
     +            ,xv(i,j,k,4)/xv(i,j,k,5)
            enddo
          enddo
        enddo
       else if (iperturb.eq.2) then
         !*********Various Pinwheel initializations********************!
  
         !FP pinwheel
         kpert=10
         ipert=5
         jpert=5
         if ((mpi_rank==0)) then
           rfluctuation=0.05*sqrt(xe(ipert,jpert,kpert,1)**2+xe(ipert,jpert,kpert,2)**2)
           write(6,*) 'pinwheel on proc ', mpi_rank
           write(6,*) 'pinwheel i,j,k,height', ipert,jpert,kpert,zcart(z(kpert),1,1)
           write(6,*) 'pinwheel magnitude', zcart(z(10),1,1)
           xv(ipert,jpert,kpert,1)=xv(ipert,jpert,kpert,1)+rfluctuation
           xv(ipert+1,jpert,kpert,2)=xv(ipert+1,jpert,kpert,2)+rfluctuation
           xv(ipert+1,jpert+1,kpert,1)=xv(ipert+1,jpert+1,kpert,1)-rfluctuation
           xv(ipert,jpert+1,kpert,2)=xv(ipert,jpert+1,kpert,2)-rfluctuation
           xvb(ipert,jpert,kpert,1)=xvb(ipert,jpert,kpert,1)+rfluctuation
           xvb(ipert+1,jpert,kpert,2)=xvb(ipert+1,jpert,kpert,2)+rfluctuation
           xvb(ipert+1,jpert+1,kpert,1)=xvb(ipert+1,jpert+1,kpert,1)-rfluctuation
           xvb(ipert,jpert+1,kpert,2)=xvb(ipert,jpert+1,kpert,2)-rfluctuation
         end if !mpirank.eq.0
       else if (iperturb.eq.3) then
        !LANL pinwheel
        k=3
        do j=1,mp
        do i=1,np
          ia=(npos-1)*np+i
          ja=(mpos-1)*mp+j
          if((ia.eq.n/4.and.ja.eq.m/4).and.k.eq.3)then
            xv(i,j,k,1)=xv(i,j,k,1)-.001
            xv(i,j,k,2)=xv(i,j,k,2)+.001
            xv(i-1,j,k,1)=xv(i-1,j,k,1)-.001
            xv(i-1,j,k,2)=xv(i-1,j,k,2)-.001
            xv(i,j-1,k,1)=xv(i,j-1,k,1)+.001
            xv(i,j-1,k,2)=xv(i,j-1,k,2)+.001
            xv(i-1,j-1,k,1)=xv(i-1,j-1,k,1)+.001
            xv(i-1,j-1,k,2)=xv(i-1,j-1,k,2)-.001
          endif
         enddo
         enddo
       end if !iperturb.eq.3
      end   ! end setPerturbation

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  the routine setRhoThetaP computes rho,theta and pressure ambient values
c the reference height for pressground and tambient is zgroundref
c if itheta=0 atmosphere is neutral
c if itheta>0 a stable atmosphere is implemented according to parameters:
c   - zstabbot : bottom of stable layer
c   - zstabmiddle : above it, temperature growths rate is lower
c NB when the domain is rotated (slopeazimuth and slopeangle.ne.0), the definition of
c the ambient profiles includes this rotation    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setRhoThetaP()
      use xve
      use gridsetup
      use msga
      use pres
      use metryic
      use constants
      use fireteca
      Implicit none
      integer :: i,j,k,ia,ja
      real :: rhogasground,thetaground
      real :: atheta1,atheta2,btheta1,btheta2
      real :: zla,zstabbot, zstabmiddle, pstabbot,pstabmiddle
      real :: pretemp
      real :: rg0 
      real,external :: zcart, presComp
      !rhogasground=pressground/(tambient*rg)
      ! FP added initial computation of thermal properties
      call set_specifichumidity()
      call updateGasThermalProperties(specifichumidity)
      rg0=rg_over_cp_gas*cp_gas

      !rg0=specifichumidity*(gammav-1)*cvvapor+(1-specifichumidity)*287.04
      ! initialisation
      rhogasground=pressground/(tambient*rg0)
      !thetaground= tambient*(pressground*1.0e-5)**(-rg/cp)
      thetaground= tambient*(pressground*1.0e-5)**(-rg_over_cp_gas)
      if(mpi_rank.eq.0) then 
       write (6,*)  'RelativeHumidity(%),specifichumidity(kg/kg)='
     +,relativehumidity,specifichumidity
       write (6,*) 'compute potential temperature: itheta=',itheta
!       write (6,*) 'cp/cv,rg/cp,cp0,rg_over_prrcp0,rg0'
!     + ,cp_over_cv_gas,rg_over_cp_gas,cp_gas,rg_over_prrcp_gas,rg0
       write (6,*) 'thetaground=',thetaground
       write (6,*) 'rhogasground=',rhogasground
      endif
      if (itheta.ge.1) then
         if (itheta.eq.1) then
            zstabbot=400.0 ! above zstabbot:
            atheta1 = 8.0/50.0 !8k increase per 50m elevation
           zstabmiddle=450.0 !1000.0 ! above zstabmiddle:
           atheta2 = 6.0/1000. !6k increase per 1000m elevation
         end if  
         if (itheta.eq.2) then
            zstabbot=800.0 ! above zstabbot:
            atheta1 = 8.0/200.0 !8k increase per 200m elevation
           zstabmiddle=1000.0 ! above zstabmiddle:
           atheta2 = 3.0/1000. !3k increase per 1000m elevation
         end if  
         ! pressure at elevation zgroundref+zstabbot
         pstabbot=presComp(zstabbot, thetaground)
         btheta1 = thetaground -zstabbot *atheta1
         btheta2 = thetaground + atheta1 * (zstabmiddle-zstabbot)
     +                - zstabmiddle *atheta2
         ! computation of (at elevation zstabmiddle+zgroundref) 
         !by integration of pressure (analytic formulation)
         pretemp = (log(atheta1 * zstabmiddle + btheta1)
     +               - log(thetaground))/atheta1
         !FP09/2019 change 10e5 to 1e5)
!         pstabmiddle=(pstabbot**(rg/cp) - g/cp*(10e5)**(rg/cp)*
!     +                pretemp)**(cp/rg)
         pstabmiddle=(pstabbot**(rg_over_cp_gas)
     +                - g/cp_gas*(1.e5)**(rg_over_cp_gas)*
     +                pretemp)**(1./rg_over_cp_gas)
         if (mpi_rank.eq.0) then
          write(6,*) 'rate above ', zstabbot, 'm is ',atheta1,' K/m' 
          write(6,*) 'rate above ', zstabmiddle, 'm is ',atheta2,' K/m' 
         endif
      endif
 
      do k=1,l
        do j=1,mp
          do i=1,np
            zla=zcart(z(k),i,j) - zgroundref !FP took zgroundref as a reference for pressure
            ! when gravity is rotated, zla vary with slopeangle
            ! and slope azimuth:
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            zla=zla + sin(slopeangle*3.14159/180.0)*(
     &        cos(slopeazimuth*3.14159/180.)*ia*dx
     &        +sin(slopeazimuth*3.14159/180.)*ja*dy)
            ! NB: here xe(...4) is just theta not theta*rho
            if (itheta.eq.0.or.(itheta.ge.1.and.zla.le.zstabbot)) then
              ! atmosphere is stable, theta is constant
              xe(i,j,k,4)=thetaground  ! for thetaground of typically 290K
              pre(i,j,k)=presComp(zla, thetaground)
            elseif (itheta.ge.1.and.zla.le.zstabmiddle) then
              ! here atmosphere is very stable:
              ! theta=a*zla+b
              xe(i,j,k,4) = atheta1 * zla + btheta1
              ! integration of pressure (analytic formulation)
              pretemp = (log(xe(i,j,k,4))
     +          - log(thetaground))/atheta1
              !FP09/2019 change 10e5 to 1e5)
!             pre(i,j,k)=(pstabbot**(rg/cp) - g/cp*(10e5)**(rg/cp)*
!     +         pretemp)**(cp/rg)
              pre(i,j,k)=(pstabbot**(rg_over_cp_gas) 
     +          - g/cp_gas*(1.e5)**(rg_over_cp_gas)*
     +          pretemp)**(1./rg_over_cp_gas)
            elseif (itheta.ge.1.and.zla.gt.zstabmiddle) then
              ! theta=a*zla+b
              xe(i,j,k,4) = atheta2 * zla + btheta2
              ! integration of pressure (analytic formulation)
              pretemp = (log(xe(i,j,k,4))
     +          - log(atheta2 * zstabmiddle + btheta2))/atheta2
              !FP09/2019 change 10e5 to 1e5)
              pre(i,j,k)=(pstabmiddle**(rg_over_cp_gas) 
     +          - g/cp_gas*(1.e5)**(rg_over_cp_gas)*
     +          pretemp)**(1./rg_over_cp_gas)
            endif
            xe(i,j,k,nv)=pre(i,j,k)/(xe(i,j,k,4)*(pre(i,j,k)*1.0e-5)**rg_over_cp_gas*rg0)
            xe(i,j,k,4)=xe(i,j,k,4)*xe(i,j,k,nv)    !here xe(...4) is theta*rho
            if (irod.eq.1) then
              tempg(i,j,k)= xe(i,j,k,4)*(pre(i,j,k)*1.0e-5)**(rg_over_cp_gas) !copied for comparison
              tambientarray(i,j,k)=xe(i,j,k,4)/xe(i,j,k,nv)*(pre(i,j,k)*1.0e-5)**(rg_over_cp_gas)
            endif
            if(i.eq.1.and.j.eq.1.and.(mpi_rank.eq.0))
     +        write(6,*) 'z:',k,zla, xe(i,j,k,4)/xe(i,j,k,nv),
     +        xe(i,j,k,nv),pre(i,j,k)
          enddo
        enddo
      enddo  
      call updated(xe(1-ih,1-ih,1,4),xe(1-ih,1-ih,1,4),
     &  np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      call updated(xe(1-ih,1-ih,1,nv),xe(1-ih,1-ih,1,nv),
     &  np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      !call updated(tambientarray,tambientarray,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      !call updated(pre,pre,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      !call updated(tempg,tempg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      ! call rmaxmin1(xe(:,:,:,nv),'xe(nv)',1-ih,np+ih,1-ih,mp+ih,l)
       !call rmaxmin1(pre,'pre',1-ih,np+ih,1-ih,mp+ih,l)
       !call rmaxmin1(tempg,'tempg',1-ih,np+ih,1-ih,mp+ih,l)
      return

      end ! END SUBROUTINE setRhoThetaP


       !endif

! this function computes pressure as a function of elevation above zgroundref
! for a neutral atmosphere
      real function presComp(zla, thetaground)
        use metryic
        use constants
        Implicit None
        real zla, thetaground

        presComp=1e5*(-g/(cp_gas*thetaground)*zla  
     +                     +(pressground*1e-5)**(rg_over_cp_gas))**(1./rg_over_cp_gas)
        return
       end function presComp
