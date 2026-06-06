!----------------------------------------------------------------
! rinit initializes all arrays and variables for use throughout
! the simulation
!----------------------------------------------------------------
subroutine rinit
  use gridlist_variables, only : nts,n,m,l,ih,dx,dy,dz,irst, &
    iwindfieldin,iwindfieldout,uswitch,vswitch,topofile, &
    ie,is,je,js,ifire,ilspgf,iturb,islip,iperturb,ifbrand, &
    dts,nts,ntp
  use gridsetup, only : np,mp,time,itrestart,dxi,dyi,dzi, &
    dti,dt,dtp
  use forcings, only : wt,wtf
  use metric_variables, only : x,y,z,zedge,zb
  use xvall, only : ixevariation,nv,xv,xe,iuvel,ivvel
  Implicit None

  ! Local Variables
  integer :: i,j,k,kv

  ! Executable Code 
  
  ! Time integration parameters
  wt(1)=0.3/nts 
  wtf(1)=0.0
  wt(nts+1)=0.7/nts
  
  ! Grid definition
  dt=dts*real(nts)        !dt is the large time step
  dtp=dt/real(ntp)        !dt is the physics time step for FIRETEC
  dxi=1./dx 
  dyi=1./dy 
  dzi=1./dz 
  dti=1./dt 
  do i=1-ih,n+ih
    x(i)=(i-1)*dx-((n+1)/2.-1)*dx
  enddo
  do j=1-ih,m+ih
    y(j)=(j-1)*dy-((m+1)/2.-1)*dy
  enddo
  do k=1,l+1
    zedge(k)=(k-1)*dz
  enddo 
  zb=zedge(l+1)
  do k=1,l
    z(k)=(k-1)*dz+0.5*dz
  enddo
  
  ! Set xevariation and windfield defaults if not set throughout the domain
  if(iwindfieldin.eq.1.or.uswitch.eq.1.or.vswitch.eq.1) ixevariation=1
  ! FP redefinition of is,ie,js,je when equals to 0 (meaning the whole domain is 
  ! considered) for windfield stuff 
  if (iwindfieldin.eq.1.or.iwindfieldout.eq.1) then
    ! windfield in or out is done by default on the whole domain
    ! if ie,is,je,js have not been defined first
    if (ie.eq.0) ie=n 
    if (is.eq.0) is=1 
    if (je.eq.0) je=m 
    if (js.eq.0) js=1 
  endif

  ! Setup coordinate transformation related matrices
  if(topofile(1:1)/='*'.and.topofile/='') call topo
  call metric

  ! Setup absorber at the sides and top
  call tinit
  ! Initialize rho, theta, and pressure profiles for xe
  call setRhoThetaP
  
  ! Ignitialize fire and fuel parameters
  call rinitfire
  if(ifire.eq.1) call rinitrad ! checks irad internally after defineRadArray
  
  ! Initialize wind profiles
  call setVelocityProfile
  
  ! Set up a large scale pressure gradient
  if(ilspgf.ge.1) call rinitLargeScalePGF
  
  ! Ignitialize fuel parameters
  if(iturb.ge.1) call rinitturb
  
  if(irst.eq.0.and.iwindfieldin.eq.0)then
    xv(:,:,:,:)=xe(:,:,:,:)
    if(iperturb.ge.1) call setPerturbation
  else
    call startFromFile
    time=itrestart*dts*nts
  endif
  time=itrestart*dts*nts
  
  do kv=1,nv
    call update(xv(:,:,:,kv),xv(:,:,:,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
  enddo

  if(islip.eq.1)then
    xv(:,:,1,iuvel)=0.
    xv(:,:,1,ivvel)=0.
  endif

  if(ifbrand.eq.1) call rinitfbrand

end subroutine rinit

!----------------------------------------------------------------
! setRhoThetaP computes rho,theta and pressure ambient values
! the reference height for pressground and tambient is zgroundref
! if itheta=0 atmosphere is neutral
! if itheta>0 a stable atmosphere is implemented according to parameters:
!   - zstabbot : bottom of stable layer
!   - zstabmiddle : above it, temperature growths rate is lower
!----------------------------------------------------------------
subroutine setRhoThetaP
  use gridlist_variables, only : itheta,relativeHumidity,dx,dy,ifire,ih, &
    l,pressground,tambient,zgroundref,slopeangle,slopeazimuth
  use constants, only : pi,g
  use gridsetup, only : np,mp
  use xvall, only : xe,itemp,irho
  use thermo_variables, only : specificHumidity,rg_over_cp_gas,cp_gas, &
    pre,tempg
  use metric_variables, only : z
  use msga_variables, only : mpi_rank,npos,mpos
  Implicit none

  ! Local Variables  
  integer :: i,j,k,ia,ja
  real :: pw
  real :: rhogasground,thetaground
  real :: atheta1,atheta2,btheta1,btheta2
  real :: zla,zstabbot, zstabmiddle, pstabbot,pstabmiddle
  real :: pretemp
  real :: rg0 
  real,external :: zcart,presComp

  ! Executable Code
  ! Set specific humidity
  ! water pressure saturation in hPA (from Humidity Conversion Formulas, Vaisala)
  ! OK at 0.083% within the range 50-100°C
  pw=relativeHumidity*6.1164*10.**(7.5914*(tambient-273.16)/ &
    (tambient-32.44)) ! water vapor pressure (hPa)
  specificHumidity=0.62199*pw/(pressground-pw)
  call updateGasThermProps(0.233,specificHumidity)
  rg0=rg_over_cp_gas*cp_gas

  ! Initialisation
  rhogasground=pressground/(tambient*rg0)
  thetaground =tambient*(pressground*1.0e-5)**(-rg_over_cp_gas)
  if(mpi_rank.eq.0) then 
    print*,'RelativeHumidity(%),specifichumidity(kg/kg)=',relativeHumidity,specificHumidity
    print*,'compute potential temperature: itheta=',itheta
    print*,'thetaground=',thetaground
    print*,'rhogasground=',rhogasground
  endif
  if(itheta.ge.1)then
    if(itheta.eq.1)then
      zstabbot=400.0 ! above zstabbot:
      atheta1 = 8.0/50.0 !8k increase per 50m elevation
      zstabmiddle=450.0 !1000.0 ! above zstabmiddle:
      atheta2 = 6.0/1000. !6k increase per 1000m elevation
    elseif(itheta.eq.2)then
      zstabbot=800.0 ! above zstabbot:
      atheta1 = 8.0/200.0 !8k increase per 200m elevation
      zstabmiddle=1000.0 ! above zstabmiddle:
      atheta2 = 3.0/1000. !3k increase per 1000m elevation
    endif  
    ! pressure at elevation zgroundref+zstabbot
    pstabbot=1e5*(-g/(cp_gas*thetaground)*zstabbot+ &
      (pressground*1e-5)**rg_over_cp_gas)**(1./rg_over_cp_gas)
    btheta1 =thetaground-atheta1*zstabbot
    btheta2 =thetaground+atheta1*(zstabmiddle-zstabbot)- &
      zstabmiddle*atheta2
    ! computation of (at elevation zstabmiddle+zgroundref) 
    ! by integration of pressure (analytic formulation)
    pretemp =(log(atheta1*zstabmiddle+btheta1)-log(thetaground)) &
      /atheta1
    pstabmiddle=(pstabbot**(rg_over_cp_gas)- &
      g/cp_gas*(1.e5)**(rg_over_cp_gas)*pretemp)**(1./rg_over_cp_gas)
    if(mpi_rank.eq.0)then
      print*,'rate above ', zstabbot, 'm is ',atheta1,' K/m' 
      print*,'rate above ', zstabmiddle, 'm is ',atheta2,' K/m' 
    endif
  endif
 
  do k=1,l
    do j=1,mp
      ja=(mpos-1)*mp+j
      do i=1,np
        zla=zcart(z(k),i,j)-zgroundref !FP took zgroundref as a reference for pressure
        ia=(npos-1)*np+i
        zla=zla+sin(slopeangle*pi/180.)* &
          (cos(slopeazimuth*pi/180.)*ia*dx &
          +sin(slopeazimuth*pi/180.)*ja*dy)
        ! when gravity is rotated, zla vary with slopeangle and slope azimuth:
        if(itheta.eq.0.or.(itheta.ge.1.and.zla.le.zstabbot))then
          ! atmosphere is stable, theta is constant
          xe(i,j,k,itemp)=thetaground  ! for thetaground of typically 300K
          pre(i,j,k)=1e5*(-g/(cp_gas*thetaground)*zla &
            +(pressground*1e-5)**rg_over_cp_gas)**(1./rg_over_cp_gas)
        elseif(itheta.ge.1.and.zla.le.zstabmiddle)then
          ! here atmosphere is very stable:
          ! theta=a*zla+b
          xe(i,j,k,itemp)=atheta1*zla+btheta1
          ! integration of pressure (analytic formulation)
          pretemp=(log(xe(i,j,k,4))-log(thetaground))/atheta1
          pre(i,j,k)=(pstabbot**(rg_over_cp_gas)- &
            g/cp_gas*(1.e5)**(rg_over_cp_gas)* &
            pretemp)**(1./rg_over_cp_gas)
        elseif(itheta.ge.1.and.zla.gt.zstabmiddle)then
          ! theta=a*zla+b
          xe(i,j,k,itemp)=atheta2*zla+btheta2
          ! integration of pressure (analytic formulation)
          pretemp=(log(xe(i,j,k,itemp)) &
            -log(atheta2*zstabmiddle+btheta2))/atheta2
          pre(i,j,k)=(pstabmiddle**rg_over_cp_gas- &
            g/cp_gas*(1.e5)**rg_over_cp_gas* &
            pretemp)**(1./rg_over_cp_gas)
        endif
        xe(i,j,k,irho)=pre(i,j,k)/(xe(i,j,k,itemp)* &
          (pre(i,j,k)*1.e-5)**rg_over_cp_gas*rg0)
        xe(i,j,k,itemp)=xe(i,j,k,itemp)*xe(i,j,k,irho)
        if(ifire.eq.1) tempg(i,j,k)=xe(i,j,k,itemp)/xe(i,j,k,irho)* &
          (pre(i,j,k)*1.e-5)**rg_over_cp_gas
        if(i.eq.1.and.j.eq.1.and.(mpi_rank.eq.0)) &
          print*,'z:',k,zla, xe(i,j,k,itemp)/xe(i,j,k,irho), &
          xe(i,j,k,irho),pre(i,j,k)
      enddo
    enddo
  enddo  
  call update(xe(1-ih,1-ih,1,itemp),xe(1-ih,1-ih,1,itemp), &
    np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
  call update(xe(1-ih,1-ih,1,irho),xe(1-ih,1-ih,1,irho), &
    np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)

end subroutine setRhoThetaP

!----------------------------------------------------------------
! setVelocityProfile sets the inital values of xe(1,2,3)
!----------------------------------------------------------------
subroutine setVelocityProfile()
  use gridlist_variables, only : ius,iue,jus,jue,n,m,l,ih,nfuel, &
    rhomicro,zu,zgroundref,zab,u0,v0,uswitch,vswitch,icorio, &
    icfmeflag
  use constants, only : pi
  use gridsetup, only : np,mp
  use fuel_variables, only : lfuel,rhoFuel,min_rhoFuel,sizeScale,actualFuelDepth
  use msga_variables, only : npos,mpos,mpi_real,mpi_comm_world, &
    ierror,mpi_sum,mpi_rank,mpi_max
  use metric_variables, only : zedge,zs,z
  use xvall, only : xe,nv,iuvel,ivvel,irho
  Implicit none

  ! Local Variables
  integer :: i,j,k,kv,ia,ja
  real :: hmax,hmaxt,LAIlow,LAIlowt,LAIup,LAIupt,LAIt
  real :: c1,c2,c3,c4
  real :: uprofilezu,ufueltopmax
  real :: zk,zla,cosc,sinc,one,gam
  real,external :: zcart
  real,allocatable :: uprofile(:,:,:)

  ! Executable Code
  if(uswitch.eq.2.or.vswitch.eq.2)then ! Empirical profile of Kaimal and Finnigan 1994 and Raupach 1993
    ! COMPUTATION OF LAIs and hmax
    hmax=0. !fuel max height (including actualfueldepth)
    LAIlow=0. ! LAI in cells k=1
    LAIup=0. !LAI in cells k>1
    if(ius.eq.0) ius=1
    if(iue.eq.0) iue=n
    if(jus.eq.0) jus=1
    if(jue.eq.0) jue=m
    do k=1,lfuel
      do j=1,mp
        ja=(mpos-1)*mp+j
        do i=1,np
          ia=(npos-1)*np+i
          if(sum(rhoFuel(:,i,j,k)).gt.min_rhoFuel.and.ia.ge.ius.and.ia.le.iue.and. &
            ja.ge.jus.and.ja.le.jue) then
            if(k.eq.1)then 
              LAIlow=LAIlow+sum(rhoFuel(:,i,j,k))/(rhomicro*sum(sizeScale(:,i,j,k))/nfuel)* &
                (zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))
              hmax=max(hmax,maxval(actualFuelDepth(:,i,j)))
            else
              hmax=max(hmax,zcart(zedge(k+1),i,j)-zs(i,j))
              LAIup=LAIup+sum(rhoFuel(:,i,j,k))/(rhomicro*sum(sizeScale(:,i,j,k))/nfuel)* &
                (zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))
            endif
          endif
        enddo
      enddo
    enddo
    call mpi_allreduce(hmax,hmaxt,1,mpi_real,mpi_max,mpi_comm_world,ierror)
    call mpi_allreduce(LAIlow,LAIlowt,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(LAIup,LAIupt,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
    LAIlowt=LAIlowt/real((iue-ius+1)*(jue-jus+1))
    LAIupt=LAIupt/real((iue-ius+1)*(jue-jus+1))
 
    if (mpi_rank.eq.0) print*,'hmaxt',hmaxt,'LAIlowt',LAIlowt,'LAIupt',LAIupt, &
      'in zone ',ius,iue,jus,jue
    if(LAIupt>0.1)then  !when a canopy is there: 
      LAIt=LAIupt+0.1*LAIlowt ! so that the profile in the canopy is not too much affected by dense low vegetation
    else 
      LAIt=LAIupt+LAIlowt
    endif

    ! DEFINITION OF CONSTANT FROM LAIt
    c1=min(sqrt(0.003+0.15*LAIt),0.3)           ! U*/uh
    c2=(1-exp(-sqrt(7.5*LAIt)))/sqrt(7.5*LAIt)  ! 1-d/h
    c3=c2*exp(-0.41/c1-log(2.)+0.5)             ! z0/h
     !(Kaimal and Finnigan 1994 : range between 1.7 and 3.2 for LAI between 1 and 4) 
    if(LAIt.le.1)then
      c4=1.7
    elseif(LAIt.ge.4)then
      c4=3.2
    else
      c4=(3.2-1.7)/(4-1)*(LAIt-1)+1.7
    endif
    
    ! COMPUTATION OF UPROFILE at height zu (for normalization at height zu)
    if(zu.ge.hmaxt)then
      uprofilezu=c1/0.41*log((zu/hmaxt+c2-1)/c3)
      if(zu.lt.2.*hmaxt)then
        ufueltopmax=c1/0.41*log(c2/c3) !theoretical value at hmaxt         
        uprofilezu=uprofilezu+(1-ufueltopmax)*(2.-zu/hmaxt)**3
      endif
    else !under the canopy
      uprofilezu=exp(-c4*(1-zu/hmaxt))
    endif
    ! COMPUTATION OF PROFILE (NORMALIZED AT CANOPY HEIGHT AT THIS STAGE)
    allocate(uprofile(np,mp,l))
    do k=1,l
      do j=1,mp
        do i=1,np
          zk=zcart(z(k),i,j)-zs(i,j)
          if(zk.ge.hmaxt)then  !above fuel
            uprofile(i,j,k)=c1/0.41*log((zk/hmaxt+c2-1)/c3) ! correction for the roughness sublayer (when actualfueldepth is greater than ztopcell/2)
            if(zk.lt.2.*hmaxt)then
              ufueltopmax=c1/0.41*log(c2/c3) !theoretical value at hmaxt         
              uprofile(i,j,k)=uprofile(i,j,k)+(1-ufueltopmax)*(2.-zk/hmaxt)**3
            endif
          else !under the canopy
            uprofile(i,j,k)=exp(-c4*(1-zk/hmaxt))
          endif
          ! correction top of first cell :
          if(k.eq.1.and.zcart(zedge(2),i,j)-zs(i,j).ge.hmaxt) &
            uprofile(i,j,1)=1.
          !RENORMALIZATION
          uprofile(i,j,k)=uprofile(i,j,k)/uprofilezu 
        enddo
      enddo
    enddo
  endif

  if(icfmeflag.eq.1)then ! Read icfmew file for wind profile
    if(mpi_rank.eq.0)then 
      print*,'going into readicfmew'
      open (701,file='icfmeallprof.dat',form='unformatted',status='old')
    endif
      do kv=1,nv
        call readio(xe(:,:,:,kv),701,l)
      enddo
    if(mpi_rank.eq.0)then
      close (701)
      print*,'exiting readicfmew'
    endif
  else ! icfmeflag.eq.0 
    if(icorio.eq.2)then
      ! special definition of xe (used for initial condition only,
      ! cause designed to be used with cyclic bc): xe used the ekman
      ! spiral as initial prof
      do k=1,l
        do j=1,mp
          do i=1,np
            zla=zcart(z(k),i,j)-zgroundref
            gam=pi/zab*zla
            cosc=1.-exp(-gam)*cos(gam)  !u(z)/ug from ekman
            sinc=exp(-gam)*sin(gam)  !v(z)/ug from ekman
            one=sqrt(cosc**2.+sinc**2.)
            cosc=cosc/one
            sinc=sinc/one
            ! NB uswitch.eq.1 is defined in xevariation.f
            if(uswitch==0)then
              xe(i,j,k,iuvel)=u0*cosc*xe(i,j,k,irho)
              xe(i,j,k,ivvel)=u0*sinc*xe(i,j,k,irho)
            elseif(uswitch==2)then
              xe(i,j,k,iuvel)=u0*cosc*xe(i,j,k,irho)*uprofile(i,j,k)  !u0 represents speed at zu
              xe(i,j,k,ivvel)=u0*sinc*xe(i,j,k,irho)*uprofile(i,j,k)  !u0 represents speed at zu
            endif
            if(vswitch==0) then
              xe(i,j,k,iuvel)=xe(i,j,k,irho)-v0*sinc*xe(i,j,k,irho)
              xe(i,j,k,ivvel)=xe(i,j,k,irho)+v0*cosc*xe(i,j,k,irho)
            elseif(vswitch==2)then
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
            ! NB uswitch.eq.1 is defined in xevariation.f
            if(uswitch==0)then
              xe(i,j,k,iuvel)=u0*xe(i,j,k,irho)
            elseif(uswitch==2)then
              xe(i,j,k,iuvel)=u0*xe(i,j,k,irho)*uprofile(i,j,k)      !u0 represents speed at zu
            endif
            if(vswitch==0)then
              xe(i,j,k,ivvel)=v0*xe(i,j,k,irho)
            elseif(vswitch==2)then
              xe(i,j,k,ivvel)=v0*xe(i,j,k,irho)*uprofile(i,j,k)      !v0 represents speed at zu
            endif
          enddo
        enddo
      enddo
    endif !icorio.ne.2
  endif !icfmeflag.eq.0
  deallocate(uprofile)
       
  do kv=iuvel,ivvel
    call update(xe(:,:,:,kv),xe(:,:,:,kv), &
      np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
  enddo

  do k=1,l
    if(mpi_rank.eq.0) print*,'initial wind profile', &
      zcart(z(k),1,1)-zs(1,1),sqrt(xe(1,1,k,iuvel)**2+ &
      xe(1,1,k,ivvel)**2)/xe(i,j,k,irho)
  enddo 

end subroutine setVelocityProfile

!----------------------------------------------------------------
! setPerturbation computes some perturbation to initialize 
! resolved turbulence typically for cyclic runs
! if 1 random perturbation of theta
! if 2 or 3 pinwheel
! perturbation done on xv
!----------------------------------------------------------------
subroutine setPerturbation
  use gridlist_variables, only : n,m,l,iperturb
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank,npos,mpos
  use metric_variables, only : z
  use xvall, only : xv,xe,itemp,irho,iuvel,ivvel
  Implicit None
  
  ! Local Variables
  integer :: i,j,k,ia,ja
  integer :: kpert,ipert,jpert
  real :: rndNum
  real :: rfluctuation
  real,external :: zcart

  ! Executable Code
  if(mpi_rank.eq.0) print*,'iPerturbation=',iperturb
  if(iperturb.eq.1)then ! random perturbation on rhotheta
    do k=1,l
      do j=1,mp
        do i=1,np
          call random_number(rndNum)
          xv(i,j,k,itemp)=xv(i,j,k,itemp)+0.1*(rndNum-0.5)
          if(mpi_rank.eq.0.and.k.eq.1.and.i.le.2.and.j.le.2) &
            print*,'random theta for i=',i,' j=',j,' is ', &
            xv(i,j,k,itemp)/xv(i,j,k,irho)
        enddo
      enddo
    enddo
  elseif(iperturb.eq.2)then ! INRA Pinwheel initialization
    kpert=10
    ipert=5
    jpert=5
    if(mpi_rank.eq.0) then
      rfluctuation=0.05*sqrt(xe(ipert,jpert,kpert,1)**2+xe(ipert,jpert,kpert,2)**2)
      write(6,*) 'pinwheel on proc ', mpi_rank
      write(6,*) 'pinwheel i,j,k,height', ipert,jpert,kpert,zcart(z(kpert),1,1)
      write(6,*) 'pinwheel magnitude', zcart(z(10),1,1)
      xv(ipert,jpert,kpert,1)=xv(ipert,jpert,kpert,1)+rfluctuation
      xv(ipert+1,jpert,kpert,2)=xv(ipert+1,jpert,kpert,2)+rfluctuation
      xv(ipert+1,jpert+1,kpert,1)=xv(ipert+1,jpert+1,kpert,1)-rfluctuation
      xv(ipert,jpert+1,kpert,2)=xv(ipert,jpert+1,kpert,2)-rfluctuation
    endif
  elseif(iperturb.eq.3)then ! LANL Pinwheel initialization
    k=3
    do j=1,mp
      do i=1,np
        ia=(npos-1)*np+i
        ja=(mpos-1)*mp+j
        if(ia.eq.n/4.and.ja.eq.m/4)then
          xv(i,j,k,iuvel)=xv(i,j,k,iuvel)-.001
          xv(i,j,k,ivvel)=xv(i,j,k,ivvel)+.001
          xv(i-1,j,k,iuvel)=xv(i-1,j,k,iuvel)-.001
          xv(i-1,j,k,ivvel)=xv(i-1,j,k,ivvel)-.001
          xv(i,j-1,k,iuvel)=xv(i,j-1,k,iuvel)+.001
          xv(i,j-1,k,ivvel)=xv(i,j-1,k,ivvel)+.001
          xv(i-1,j-1,k,iuvel)=xv(i-1,j-1,k,iuvel)+.001
          xv(i-1,j-1,k,ivvel)=xv(i-1,j-1,k,ivvel)-.001
        endif
      enddo
    enddo
  endif

end subroutine setPerturbation
