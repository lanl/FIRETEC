!-----------------------------------------------------------------------
! rinit initializes all arrays and variables for use throughout
! the simulation
!-----------------------------------------------------------------------
subroutine rinit
  use gridlist_variables, only : wp,nts,n,m,l,ih,dx,dy,dz,irst, &
    topofile,iemissions,ilspgf,iturb,islip,iperturb,ifbrand, &
    dt,nts,ixevariation,iheatsource,isensor
  use radiation_gridlist_variables, only : irad
  use gridsetup, only : np,mp,time,itrestart,dxi,dyi,dzi, &
    dti,dts
  use forcings, only : wt
  use metric_variables, only : x,y,z,zedge,zb
  use xvall, only : nv,xv,xe,iuvel,ivvel
  Implicit None

  ! Local Variables
  integer :: i,j,k,kv,iuv

  ! Executable Code 
  
  ! Time integration parameters
  wt(1)=0.3_wp/nts 
  wt(nts+1)=0.7_wp/nts
  
  ! Grid definition
  dts=dt/real(nts,kind=wp)        !dts is the small time step
  dxi=1._wp/dx 
  dyi=1._wp/dy 
  dzi=1._wp/dz 
  dti=1._wp/dt 
  do i=1-ih,n+ih
    x(i)=(i-1)*dx-((n+1)/2._wp-1)*dx
  enddo
  do j=1-ih,m+ih
    y(j)=(j-1)*dy-((m+1)/2._wp-1)*dy
  enddo
  do k=1,l+1
    zedge(k)=(k-1)*dz
  enddo 
  zb=zedge(l+1)
  do k=1,l
    z(k)=(k-1)*dz+0.5_wp*dz
  enddo
  
  ! Setup coordinate transformation related matrices
  if(topofile(1:1)/='*'.and.topofile/=' ') call topo
  call metric

  ! Setup absorber at the sides and top
  call tinit
  ! Initialize rho, theta, and pressure profiles for xe
  call setRhoThetaP
  
  ! Ignitialize fire and fuel parameters
  call rinitVeg
  if(irad.gt.0) call rinitrad ! checks irad internally after defineRadArray
  if(iemissions.gt.0) call rinitEmission(iemissions,xe,1,np,1,mp,l,nv)

  ! Initialize wind profiles
  call setVelocityProfile
  
  ! Ignitialize fuel parameters
  if(iturb.ge.1) call rinitSubTurb
  
  xv(1:np,1:mp,:,:)=xe
  if(iperturb.ge.1) call setPerturbation
  
  if(irst.eq.1) call irstreadio
  
  ! Set up a large scale pressure gradient
  if(ilspgf.ge.1) call rinitLargeScalePGF

  ! Initialize environmental temporal variations
  
  time=itrestart*dts*nts
  if(ixevariation.ge.1) call rinitXeVariation(ixevariation)
  
  do kv=1,nv
    call loop_update_flags(kv,iuv)
    call update(xv(:,:,:,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0,iuv)
  enddo

  if(islip.eq.1)then
    xv(:,:,1,iuvel)=0._wp
    xv(:,:,1,ivvel)=0._wp
  endif

  if(ifbrand.eq.1) call rinitFirebrands(irst)
  if(isensor.eq.1) call rinitSense
  if(iheatsource.ge.1) call rinitHeatSource()

end subroutine rinit

!-----------------------------------------------------------------------
! setRhoThetaP computes rho,theta and pressure ambient values
! the reference height for pressground and tambient is zgroundref
! if itheta=0 atmosphere is neutral
! if itheta>0 a stable atmosphere is implemented according to parameters
!   - zstabbot : bottom of stable layer
!   - zstabmiddle : above it, temperature growths rate is lower
!-----------------------------------------------------------------------
subroutine setRhoThetaP  !Theta stab is the lapse rate
  use gridlist_variables, only : zStabBot,zStabMid,thetaStabMid, &
    thetaStabTop,dx,dy,l,pressground,tambient,zgroundref,slopeangle, &
    slopeazimuth,wp
  use emission_general_variables, only : specificHumidity
  use emission_gridlist_variables, only : relativeHumidity
  use constants, only : pi,g,pref,Rgas
  use gridsetup, only : np,mp
  use xvall, only : xe,itheta,irho
  use thermo_variables, only : pre,tempg
  use forcings, only : xvLim
  use metric_variables, only : z
  use msga_variables, only : mpi_rank,npos,mpos
  use gridlist_variables, only : wp
  Implicit none

  ! Local Variables  
  integer :: i,j,k,ia,ja
  real(wp) :: rhogasground,thetaground
  real(wp) :: zla,pStabBot,pStabMid
  real(wp) :: yH2O,yO2,yN2
  real(wp) :: cpGas
  real(wp) :: mwGas 

  real(wp),external :: zcart
  
  ! Executable Code
  yH2O=specificHumidity
  yO2 =0.233_wp/(specificHumidity+1._wp)
  yN2 =1._wp-yH2O-yO2
  cpGas=yH2O*1996.0_wp+yO2*918.0_wp+yN2*1040.0_wp
  mwGas=1./(yH2O/18.015_wp+yO2/31.998_wp+yN2/28.014_wp)

  ! Set ground potential temperature and density
  rhogasground=pressground*mwGas/(tambient*Rgas)
  thetaground =Tambient*(pref/pressground)**(Rgas/cpGas/mwGas)
  if(mpi_rank.eq.0) then 
    print*,'relativeHumidity(%),specificHumidity(kg/kg)=', &
      relativeHumidity,specificHumidity
    print*,'thetaGround=',thetaGround
    print*,'rhoGasGround=',rhoGasGround
  endif

  ! Initialisation atmopsheric profile
  if(abs(thetaStabMid).lt.1.0e-8_wp) thetaStabMid=1.0e-8_wp
  if(abs(thetaStabTop).lt.1.0e-8_wp) thetaStabTop=1.0e-8_wp
  ! pressure at elevation zgroundref+zstabbot
  pStabBot=pref*((pressground/pref)**(Rgas/cpGas/mwGas)- &
    g/(cpGas*thetaground)*zStabBot)**(cpGas*mwGas/Rgas)
  ! computation of (at elevation zstabmiddle+zgroundref) 
  ! by integration of pressure (analytic formulation)
  pStabMid=pref*((pStabBot/pref)**(Rgas/mwGas/cpGas)- &
    g*log((thetaground+(zStabMid-zStabBot)*thetaStabMid)/thetaground)/ &
    (cpGas*thetaStabMid))**(cpGas*mwGas/Rgas)
 
  do k=1,l
    do j=1,mp
      ja=(mpos-1)*mp+j
      do i=1,np
        ia=(npos-1)*np+i
        zla=zcart(z(k),i,j)-zgroundref &  ! zgroundref as a reference for pressure
          +sin(slopeangle*pi/180.0_wp)* &     ! when gravity is rotated, zla vary with slopeangle and slope azimuth:
          (cos(slopeazimuth*pi/180.0_wp)*ia*dx &
          +sin(slopeazimuth*pi/180.0_wp)*ja*dy)
        if(zla.le.zStabBot)then
          ! atmosphere is stable, theta is constant
          xe(i,j,k,itheta)=thetaground  ! for thetaground of typically 300K
          pre(i,j,k)=pref*((pressground/pref)**(Rgas/mwGas/cpGas)- &
            g/(cpGas*thetaground)*zla)**(cpGas*mwGas/Rgas)
        elseif(zla.le.zStabMid)then
          ! here atmosphere is usually very stable (middle layer):
          xe(i,j,k,itheta)=thetaground+(zla-zStabBot)*thetaStabMid
          ! integration of pressure (analytic formulation)
          pre(i,j,k)=pref*((pStabBot/pref)**(Rgas/mwGas/cpGas)- &
            g*log(xe(i,j,k,itheta)/thetaground)/ &
            (cpGas*thetaStabMid))**(cpGas*mwGas/Rgas)
        else
          ! upper atmphsphere where instabilities are most common
          xe(i,j,k,itheta)=thetaground+(zStabMid-zStabBot)*thetaStabMid&
            +(zla-zStabMid)*thetaStabTop
          ! integration of pressure (analytic formulation)
          pre(i,j,k)=pref*((pStabMid/pref)**(Rgas/mwGas/cpGas)- &
            g*log(xe(i,j,k,itheta)/ &
            (thetaground+(zStabMid-zStabBot)*thetaStabMid))/ &
            (cpGas*thetaStabTop))**(cpGas*mwGas/Rgas)
        endif
        xe(i,j,k,irho)=pre(i,j,k)/(xe(i,j,k,itheta)* &
          (pre(i,j,k)/pref)**(Rgas/mwGas/cpGas)*Rgas/mwGas)
        tempg(i,j,k)=xe(i,j,k,itheta)* &
          (pre(i,j,k)/pref)**(Rgas/mwGas/cpGas)
        xe(i,j,k,itheta)=xe(i,j,k,itheta)*xe(i,j,k,irho)
        if(i.eq.1.and.j.eq.1.and.(mpi_rank.eq.0)) &
          print*,'z:',k,real(zla),real(xe(i,j,k,itheta)/xe(i,j,k,irho))&
            ,real(xe(i,j,k,irho)),real(pre(i,j,k))
      enddo
    enddo
  enddo 
  xvLim(1,itheta)=200._wp
  xvLim(2,itheta)=5000._wp
  xvLim(1,irho)=1.e-6_wp
  xvLim(2,irho)=10._wp 

end subroutine setRhoThetaP

!-----------------------------------------------------------------------
! setVelocityProfile sets the inital values of xe(1,2,3)
!-----------------------------------------------------------------------
subroutine setVelocityProfile()
  use gridlist_variables, only : ius,iue,jus,jue,n,m,l,zu, &
    zgroundref,zab,u0,v0,uswitch,vswitch,icorio,wp
  use vegetation_gridlist_variables, only : nfuel,rhoMicro
  use constants, only : pi
  use gridsetup, only : np,mp
  use fuel_variables, only : lfuel,min_rhoFuel,sizeScale,actualFuelDepth
  use msga_variables, only : npos,mpos,mpi_comm_world,ierror,mpi_sum, &
    mpi_rank,mpi_max
  use msga_variables, only : mpi_wp
  use metric_variables, only : zedge,zs,z
  use xvall, only : xe,iuvel,ivvel,irho,xvfuel,irhof
  use gridlist_variables, only : wp
  Implicit none

  ! Local Variables
  integer :: i,j,k,ia,ja
  real(wp) :: hmax,hmaxt,LAIlow,LAIlowt,LAIup,LAIupt,LAIt
  real(wp) :: c1,c2,c3,c4
  real(wp) :: uprofilezu,ufueltopmax
  real(wp) :: zk,zla,cosc,sinc,one,gam
  real(wp),external :: zcart
  real(wp),allocatable :: uprofile(:,:,:)

  ! Executable Code
  allocate(uprofile(np,mp,l))
  if(uswitch.eq.2.or.vswitch.eq.2)then ! Empirical profile of Kaimal and Finnigan 1994 and Raupach 1993
    ! COMPUTATION OF LAIs and hmax
    hmax=0.0_wp !fuel max height (including actualfueldepth)
    LAIlow=0.0_wp ! LAI in cells k=1
    LAIup=0.0_wp !LAI in cells k>1
    if(ius.eq.0) ius=1
    if(iue.eq.0) iue=n
    if(jus.eq.0) jus=1
    if(jue.eq.0) jue=m
    do k=1,lfuel
      do j=1,mp
        ja=(mpos-1)*mp+j
        do i=1,np
          ia=(npos-1)*np+i
          if(sum(xvfuel(:,i,j,k,irhof)).gt.min_rhoFuel.and.ia.ge.ius &
            .and.ia.le.iue.and.ja.ge.jus.and.ja.le.jue) then
            if(k.eq.1)then 
              LAIlow=LAIlow+sum(xvfuel(:,i,j,k,irhof))/ &
                (rhomicro*sum(sizeScale(:,i,j,k))/nfuel)* &
                (zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))
              hmax=max(hmax,maxval(actualFuelDepth(:,i,j)))
            else
              hmax=max(hmax,zcart(zedge(k+1),i,j)-zs(i,j))
              LAIup=LAIup+sum(xvfuel(:,i,j,k,irhof))/ &
                (rhomicro*sum(sizeScale(:,i,j,k))/nfuel)* &
                (zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))
            endif
          endif
        enddo
      enddo
    enddo
    call mpi_allreduce(hmax,hmaxt,1,mpi_wp,mpi_max, &
      mpi_comm_world,ierror)
    call mpi_allreduce(LAIlow,LAIlowt,1,mpi_wp,mpi_sum, &
      mpi_comm_world,ierror)
    call mpi_allreduce(LAIup,LAIupt,1,mpi_wp,mpi_sum, &
      mpi_comm_world,ierror)
    LAIlowt=LAIlowt/real((iue-ius+1)*(jue-jus+1),kind=wp)
    LAIupt=LAIupt/real((iue-ius+1)*(jue-jus+1),kind=wp)
 
    if (mpi_rank.eq.0) print*,'hmaxt',hmaxt,'LAIlowt',LAIlowt, &
      'LAIupt',LAIupt,'in zone ',ius,iue,jus,jue
    if(LAIupt>0.1_wp)then  !when a canopy is there: 
      LAIt=LAIupt+0.1_wp*LAIlowt ! so that the profile in the canopy is not too much affected by dense low vegetation
    else 
      LAIt=LAIupt+LAIlowt
    endif

    ! DEFINITION OF CONSTANT FROM LAIt
    c1=min(sqrt(0.003_wp+0.15_wp*LAIt),0.3_wp)           ! U*/uh
    c2=(1-exp(-sqrt(7.5_wp*LAIt)))/sqrt(7.5_wp*LAIt)  ! 1-d/h
    c3=c2*exp(-0.41_wp/c1-log(2.0_wp)+0.5_wp)             ! z0/h
     !(Kaimal and Finnigan 1994 : range between 1.7 and 3.2 for LAI between 1 and 4) 
    if(LAIt.le.1)then
      c4=1.7_wp
    elseif(LAIt.ge.4)then
      c4=3.2_wp
    else
      c4=(3.2_wp-1.7_wp)/(4-1)*(LAIt-1)+1.7_wp
    endif
    
    ! COMPUTATION OF UPROFILE at height zu (for normalization at height zu)
    if(zu.ge.hmaxt)then
      uprofilezu=c1/0.41_wp*log((zu/hmaxt+c2-1)/c3)
      if(zu.lt.2.0_wp*hmaxt)then
        ufueltopmax=c1/0.41_wp*log(c2/c3) !theoretical value at hmaxt         
        uprofilezu=uprofilezu+(1-ufueltopmax)*(2.0_wp-zu/hmaxt)**3
      endif
    else !under the canopy
      uprofilezu=exp(-c4*(1-zu/hmaxt))
    endif
    ! COMPUTATION OF PROFILE (NORMALIZED AT CANOPY HEIGHT AT THIS STAGE)
    do k=1,l
      do j=1,mp
        do i=1,np
          zk=zcart(z(k),i,j)-zs(i,j)
          if(zk.ge.hmaxt)then  !above fuel
            uprofile(i,j,k)=c1/0.41_wp*log((zk/hmaxt+c2-1)/c3) ! correction for the roughness sublayer (when actualfueldepth is greater than ztopcell/2)
            if(zk.lt.2.0_wp*hmaxt)then
              ufueltopmax=c1/0.41_wp*log(c2/c3) !theoretical value at hmaxt         
              uprofile(i,j,k)=uprofile(i,j,k) &
                +(1-ufueltopmax)*(2._wp-zk/hmaxt)**3
            endif
          else !under the canopy
            uprofile(i,j,k)=exp(-c4*(1-zk/hmaxt))
          endif
          ! correction top of first cell :
          zk=zcart(zedge(2),i,j)-zs(i,j)
          if((k.eq.1).and.(zk.ge.hmaxt)) uprofile(i,j,1)=1._wp
          !RENORMALIZATION
          uprofile(i,j,k)=uprofile(i,j,k)/uprofilezu 
        enddo
      enddo
    enddo
  endif

  if(icorio.eq.2)then
    ! special definition of xe (used for initial condition only,
    ! cause designed to be used with cyclic bc): xe used the ekman
    ! spiral as initial prof
    do k=1,l
      do j=1,mp
        do i=1,np
          zla=zcart(z(k),i,j)-zgroundref
          gam=pi/zab*zla
          cosc=1._wp-exp(-gam)*cos(gam)  !u(z)/ug from ekman
          sinc=exp(-gam)*sin(gam)  !v(z)/ug from ekman
          one=sqrt(cosc**2.0_wp+sinc**2.0_wp)
          cosc=cosc/one
          sinc=sinc/one
          ! NB uswitch.eq.1 is defined in xevariation.f
          if((uswitch.eq.0).or.(uswitch.eq.1))then
            xe(i,j,k,iuvel)=u0*cosc*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=u0*sinc*xe(i,j,k,irho)
          elseif(uswitch.eq.2)then
            xe(i,j,k,iuvel)=u0*cosc*xe(i,j,k,irho)*uprofile(i,j,k)  !u0 represents speed at zu
            xe(i,j,k,ivvel)=u0*sinc*xe(i,j,k,irho)*uprofile(i,j,k)  !u0 represents speed at zu
          endif
          if((vswitch.eq.0).or.(vswitch.eq.1)) then
            xe(i,j,k,iuvel)=xe(i,j,k,irho)-v0*sinc*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=xe(i,j,k,irho)+v0*cosc*xe(i,j,k,irho)
          elseif(vswitch.eq.2)then
            xe(i,j,k,iuvel)=xe(i,j,k,irho)-v0*sinc*xe(i,j,k,irho) &
              *uprofile(i,j,k)  !v0 represents speed at zu
            xe(i,j,k,ivvel)=xe(i,j,k,irho)+v0*cosc*xe(i,j,k,irho) &
              *uprofile(i,j,k)  !v0 represents speed at zu
          endif
        enddo
      enddo
    enddo
  else ! icorio.ne.2
    do k=1,l
      do j=1,mp
        do i=1,np
          if((uswitch.eq.0).or.(uswitch.eq.1))then
            xe(i,j,k,iuvel)=u0*xe(i,j,k,irho)
          elseif(uswitch.eq.2)then
            xe(i,j,k,iuvel)=u0*xe(i,j,k,irho)*uprofile(i,j,k)      !u0 represents speed at zu
          endif
          if((vswitch.eq.0).or.(vswitch.eq.1))then
            xe(i,j,k,ivvel)=v0*xe(i,j,k,irho)
          elseif(vswitch.eq.2)then
            xe(i,j,k,ivvel)=v0*xe(i,j,k,irho)*uprofile(i,j,k)      !v0 represents speed at zu
          endif
        enddo
      enddo
    enddo
  endif !icorio.ne.2
  
  do k=1,l
    if(mpi_rank.eq.0) print*,'initial wind profile', &
      zcart(z(k),1,1)-zs(1,1),xe(1,1,k,iuvel), &
      sqrt(xe(1,1,k,iuvel)**2+xe(1,1,k,ivvel)**2)/xe(1,1,k,irho)
  enddo 
  
  deallocate(uprofile)

end subroutine setVelocityProfile

!-----------------------------------------------------------------------
! setPerturbation computes some perturbation to initialize 
! resolved turbulence typically for cyclic runs
! if 1 random perturbation of theta
! if 2 or 3 pinwheel
! perturbation done on xv
!-----------------------------------------------------------------------
subroutine setPerturbation
  use gridlist_variables, only : n,m,l,iperturb,wp
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank,npos,mpos
  use metric_variables, only : z
  use xvall, only : xv,xe,itheta,irho,iuvel,ivvel
  Implicit None
  
  ! Local Variables
  integer :: i,j,k,ia,ja
  integer :: kpert,ipert,jpert
  real(wp) :: rndNum
  real(wp) :: rfluctuation
  real(wp),external :: zcart

  ! Executable Code
  if(mpi_rank.eq.0) print*,'iPerturbation=',iperturb
  if(iperturb.eq.1)then ! random perturbation on rhotheta
    do k=1,l
      do j=1,mp
        do i=1,np
          call random_number(rndNum)
          xv(i,j,k,itheta)=xv(i,j,k,itheta)+0.1_wp*(rndNum-0.5_wp)
          if(mpi_rank.eq.0.and.k.eq.1.and.i.le.2.and.j.le.2) &
            print*,'random theta for i=',i,' j=',j,' is ', &
            xv(i,j,k,itheta)/xv(i,j,k,irho)
        enddo
      enddo
    enddo
  elseif(iperturb.eq.2)then ! INRA Pinwheel initialization
    kpert=10
    ipert=5
    jpert=5
    if(mpi_rank.eq.0) then
      rfluctuation=0.05_wp*sqrt(xe(ipert,jpert,kpert,1)**2 &
        +xe(ipert,jpert,kpert,2)**2)
      write(6,*) 'pinwheel on proc ', mpi_rank
      write(6,*) 'pinwheel i,j,k,height', ipert,jpert,kpert, &
        zcart(z(kpert),1,1)
      write(6,*) 'pinwheel magnitude', zcart(z(10),1,1)
      xv(ipert,jpert,kpert,1)=xv(ipert,jpert,kpert,1)+rfluctuation
      xv(ipert+1,jpert,kpert,2)=xv(ipert+1,jpert,kpert,2)+rfluctuation
      xv(ipert+1,jpert+1,kpert,1)=xv(ipert+1,jpert+1,kpert,1) &
        -rfluctuation
      xv(ipert,jpert+1,kpert,2)=xv(ipert,jpert+1,kpert,2)-rfluctuation
    endif
  elseif(iperturb.eq.3)then ! LANL Pinwheel initialization
    k=3
    do j=1,mp
      do i=1,np
        ia=(npos-1)*np+i
        ja=(mpos-1)*mp+j
        if(ia.eq.n/4.and.ja.eq.m/4)then
          xv(i,j,k,iuvel)=xv(i,j,k,iuvel)-0.001_wp
          xv(i,j,k,ivvel)=xv(i,j,k,ivvel)+0.001_wp
          xv(i-1,j,k,iuvel)=xv(i-1,j,k,iuvel)-0.001_wp
          xv(i-1,j,k,ivvel)=xv(i-1,j,k,ivvel)-0.001_wp
          xv(i,j-1,k,iuvel)=xv(i,j-1,k,iuvel)+0.001_wp
          xv(i,j-1,k,ivvel)=xv(i,j-1,k,ivvel)+0.001_wp
          xv(i-1,j-1,k,iuvel)=xv(i-1,j-1,k,iuvel)+0.001_wp
          xv(i-1,j-1,k,ivvel)=xv(i-1,j-1,k,ivvel)-0.001_wp
        endif
      enddo
    enddo
  endif

end subroutine setPerturbation
