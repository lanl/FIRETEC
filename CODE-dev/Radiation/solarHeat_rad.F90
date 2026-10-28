!-----------------------------------------------------------------------
! SolarHeating computes rates of vegetation's internal energy change due 
! to solar radiative heating
!-----------------------------------------------------------------------
subroutine solarHeating(force_sies,rhoFuel,rhoWater,atmWater, &
    ilh,iuh,jlh,juh,ku,il,iu,jl,ju,lls,nfuel)
  use gridlist_variables, only : wp,n,m,dx,dy,ibc_left,ibc_right, &
    ibc_top,ibc_bot,nprocx,nprocy
  use vegetation_gridlist_variables, only : rhoMicro
  use radiation_gridlist_variables, only : unixTimeStart,plotLong, &
    plotLat
  use gridsetup, only : np,mp,time
  use msga_variables, only : mpos,npos,numprocs,mpi_comm_world,ierror
  use radiation_montecarlo_variables, only : zCartEdgeAZGlobal
  use constants, only : tolerance,sigma
  use metric_variables, only : zedge
  use solar_variables, only : aod500,col_d_w,col_d_O3,cloud_frac, &
    leaf_albedo,leaf_emis
  use thermo_variables, only : pr,tempg
  use xvall, only : iH2O
  use msga_variables, only : mpi_wp
  use fuel_variables, only : sizeScale,lfuel,temps,kfindex
  use iso_fortran_env, only : int64
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,ku,il,iu,jl,ju,lls,nfuel
  real(wp),intent(in) :: rhoFuel(nfuel,il:iu,jl:ju,lls)
  real(wp),intent(in) :: rhoWater(nfuel,il:iu,jl:ju,lls)
  real(wp),intent(in) :: atmWater(ilh:iuh,jlh:juh,ku) ! Mass Fraction
  real(wp),intent(inout) :: force_sies(nfuel,il:iu,jl:ju,lls)

  integer :: i,j,k,ia,ja,it
  integer :: i1,j1,k1,kcell
  integer :: iproc,iprocx,jprocy
  integer :: nsize
  real(wp) :: DeepShadeLAI=10._wp
  real(wp) :: dir(3)
  real(wp) :: z1top,z2top,z3top,z4top,zlayer
  real(wp) :: xf,yf,zf
  real(wp) :: t,tx,ty,tz,tz1,tz2,delta
  real(wp) :: acoeftop,bcoeftop,ccoeftop
  real(wp) :: LAIlocal(np,mp,lfuel+nfuel)
  real(wp) :: LAIgather(np,mp,lfuel+nfuel,numprocs)
  real(wp) :: LAIall(n,m,lfuel+nfuel)
  real(wp) :: sun_azimuth,sun_elevation
  real(wp) :: sumShade(np,mp,lfuel+nfuel)
  real(wp) :: ghi
  real(wp) :: QswIN,QlwIN,QlwOUT
  real(wp) :: av,emis_g
  real(int64) :: unixT
  real(wp),external :: zcart

  ! Executable Code
  ! ----- Section 1 Compute and Communicate leaf area index
  do k=1,nfuel
    do j=1,mp
      do i=1,np
        ia=(npos-1)*np+i
        ja=(mpos-1)*mp+j
        zlayer=(sum(zCartEdgeAZGlobal(:,ia,ja,k)- &
          zCartEdgeAZGlobal(:,ia,ja,k-1)))/4._wp
        LAIlocal(i,j,k)=sum(rhoFuel(kfindex(i,j,k:nfuel),i,j,1)/ &
          sizeScale(kfindex(i,j,k:nfuel),i,j,1)*zlayer/rhoMicro)
      enddo
    enddo
  enddo
  LAIlocal(:,:,nfuel+1)=0._wp
  do k=nfuel+2,lfuel+nfuel
    do j=1,mp
      do i=1,np
        ia=(npos-1)*np+i
        ja=(mpos-1)*mp+j
        zlayer=(sum(zCartEdgeAZGlobal(:,ia,ja,k)- &
          zCartEdgeAZGlobal(:,ia,ja,k-1)))/4._wp
        LAIlocal(i,j,k)=sum(rhoFuel(:,i,j,k-nfuel)/rhoMicro &
          *zlayer/sizeScale(:,i,j,k-nfuel))
      enddo
    enddo
  enddo
  nsize=np*mp*(lfuel+nfuel)
  call mpi_allgather(LAIlocal,nsize,mpi_wp,LAIgather,nsize,mpi_wp, &
    mpi_comm_world,ierror)
  do iprocx=1,nprocx
    do jprocy=1,nprocy
      iproc=1+(iprocx-1)+(jprocy-1)*nprocx
      do k=1,lfuel+nfuel
        do j=1,mp
          do i=1,np
            ia=(iprocx-1)*np+i
            ja=(jprocy-1)*mp+j
            LAIall(ia,ja,k)=LAIgather(i,j,k,iproc)
          enddo
        enddo
      enddo
    enddo
  enddo

  ! ----- Section 2 Cast Shade on cells based on elevated fuels, sun
  ! ----- angle, and TODO:topography
  unixT=unixTimeStart+int(time,int64)
  call sun_position(plotLat,plotLong,unixT,sun_elevation,sun_azimuth)
  
  ! Ray direction
  dir(1)=cos(sun_elevation)*cos(sun_azimuth)
  dir(2)=cos(sun_elevation)*sin(sun_azimuth)
  dir(3)=sin(sun_elevation)
  sumShade(:,:,:)=0._wp
  do k=1,lfuel+nfuel
    if(k.lt.nfuel+1)then
      kcell=1
    else
      kcell=k-nfuel
    endif
    do j=1,mp
      do i=1,np
        if(any((rhoWater(:,i,j,kcell)+rhoFuel(:,i,j,kcell)) &
          .gt.tolerance))then
          i1=(npos-1)*np+i
          j1=(mpos-1)*mp+j
          k1=k
          xf=dx*(i1-0.5_wp)
          yf=dy*(j1-0.5_wp)
          zf=(sum(zCartEdgeAZGlobal(:,i1,j1,k)+ &
            zCartEdgeAZGlobal(:,i1,j1,k-1)))/8._wp
          do while(i1.ge.1.and.i1.le.n.and.j1.ge.1.and.j1.le.m.and. &
            k1.ge.1.and.k1.le.lfuel)

            if(dir(1).gt.0)then
              tx=(i1*dx-xf)/dir(1)
            else
              tx=(i1*dx-dx-xf)/dir(1)
            endif
            if(dir(2).gt.0)then
              ty=(j1*dy-yf)/dir(2)
            else
              ty=(j1*dy-dy-yf)/dir(2)
            endif
            tz=1.e6_wp
            tz1=0._wp
            tz2=0._wp
            z1top=zCartEdgeAZGlobal(1,i1,j1,k1)
            z2top=zCartEdgeAZGlobal(2,i1,j1,k1)
            z3top=zCartEdgeAZGlobal(3,i1,j1,k1)
            z4top=zCartEdgeAZGlobal(4,i1,j1,k1)
            acoeftop=(z1top+z4top-z2top-z3top)*dir(1)*dir(2)/(dx*dy)
            bcoeftop=(-dir(3)+dir(2)/dy*(z3top-z1top)+dir(1)/dx* &
              (z2top-z1top)+(z1top+z4top-z2top-z3top)* &
              (dir(2)/dy*(xf/dx-(i1-1))+dir(1)/dx*(yf/dy-(j1-1))))
            ccoeftop=(z1top-zf+(z3top-z1top)*(yf/dy-(j1-1))+ &
              (z2top-z1top)*(xf/dx-(i1-1))+(z1top+z4top-z2top-z3top)* &
              (xf/dx-(i1-1))*(yf/dy-(j1-1)))
            if(abs(acoeftop).lt.tolerance)then    !a.eq.0
              if(ccoeftop*bcoeftop.lt.0) tz=-ccoeftop/bcoeftop  ! b.eq.0 correct intersection with tz>0
            else   ! the code never goes here in case of flat topo
              delta=bcoeftop**2-4*acoeftop*ccoeftop
              if(delta.gt.0._wp)then     !two solutions
                tz1=-0.5_wp*(bcoeftop-sqrt(delta))/acoeftop
                tz2=-0.5_wp*(bcoeftop+sqrt(delta))/acoeftop
                if(tz1.gt.0._wp)then   ! choose the smallest positive root
                  tz=tz1 
                  if(tz2.gt.0._wp) tz=min(tz1,tz2)
                elseif(tz2.gt.0._wp)then
                  tz=tz2
                endif
              endif
            endif

            ! Intersection with cell face
            t=min(tx,ty,tz)
            sumShade(i,j,kcell)=sumShade(i,j,kcell)+ &
              t/(zcart(zedge(k1),i1,j1)+zcart(zedge(k1+1),i1,j1)) &
              /2._wp*LAIall(i1,j1,k1)/DeepShadeLAI
            if(abs(t-tx).lt.tolerance)then
              i1=i1+int(sign(1._wp,dir(1)))
            elseif(abs(t-ty).lt.tolerance)then
              j1=j1+int(sign(1._wp,dir(2)))
            else
              k1=k1+1
            endif
           
            xf=xf+dir(1)*t
            yf=yf+dir(2)*t
            zf=zf+dir(3)*t
            
            if(ibc_left.eq.0.and.i1.lt.1) then ! Cyclic boundary condition
              i1=i1+n
              xf=xf+n*dx
            endif
            if(ibc_right.eq.0.and.i1.gt.n) then ! Cyclic boundary condition
              i1=i1-n
              xf=xf-n*dx
            endif
            if(ibc_bot.eq.0.and.j1.lt.1) then ! Cyclic boundary condition
              j1=j1+m
              yf=yf+m*dy
            endif
            if(ibc_top.eq.0.and.j1.gt.m) then ! Cyclic boundary condition
              j1=j1-m
              yf=yf-m*dy
            endif
          enddo
        endif
      enddo
    enddo
  enddo
  do i=1,np
    do j=1,mp
      sumShade(i,j,1)=sumShade(i,j,1)/(nfuel+1) ! Correction for the multiple counting done in the loops above
    enddo
  enddo
  
  ! ----- Section 3 Calculate Solar Heating/Cooling
  call compute_ghi(ghi,sun_elevation,aod500,col_d_w,col_d_O3,cloud_frac)
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        ! Incoming shortwave radiation
        QswIN=(1._wp-leaf_albedo)*ghi* &
          (1._wp-min(sumShade(i,j,k),0.95_wp))/2._wp ! 5% of SW will always get through for ambient light, divide by 2 as only half
                                                     ! the surface area is assumed facing the sun
        
        ! Incoming longwave radiation
        emis_g=0.52_wp ! Emissivity of atmospheric gases
        if(iH2O.ne.1) &
          emis_g=emis_g+0.0065_wp*sqrt(atmWater(i,j,k)*pr(i,j,k)) ! Brunt equation for gas emissivity
        QlwIN=emis_g*sigma*tempg(i,j,k)**4
        do it=1,nfuel
          
          ! Outgoing longwave radiation
          QlwOUT=leaf_emis*sigma*temps(it,i,j,k)**4
          
          ! Add to forcing term
          av=2._wp*(rhoFuel(it,i,j,k)/rhoMicro)/sizeScale(it,i,j,k)
          force_sies(it,i,j,k)=force_sies(it,i,j,k) &
            +(QswIN+QlwIN-QlwOUT)*av &
            /(rhoFuel(it,i,j,k)+rhoWater(it,i,j,k))
        enddo
      enddo
    enddo
  enddo

end subroutine SolarHeating

!-----------------------------------------------------------------------
! sun_position calculates the zenith and azimuth (in radians) of the sun
! based on the geodetic position and UTC time.
!-----------------------------------------------------------------------
subroutine sun_position(lat_deg,long_deg,unixT,elevation,azimuth)
  use gridlist_variables, only : wp
  use constants, only : pi
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none

  ! Local Variables
  real(int64),intent(in) :: unixT    ! Unix Time (s)
  real(wp),intent(in)  :: lat_deg    ! Geodetic longitude (+E, degrees)
  real(wp),intent(in)  :: long_deg   ! Geodetic latitude (+N, degrees)
  real(wp),intent(out) :: elevation  ! Solar elevation angle (radians)
  real(wp),intent(out) :: azimuth    ! Solar azimuth from true north, increasing eastward (radians)

  real(wp) :: L0,M,C,lambda
  real(wp) :: yfac,EoT_min
  real(wp) :: Hr,Ha
  real(wp) :: phi,delta
  real(wp) :: e=0.016708634_wp ! Eccentricity of Earth's Orbit
  real(wp) :: eps=0.4091_wp    ! Mean obliquity (mean)

  ! Executable Code
  ! ---- Sun apparent longitude & obliquity
  L0 = modulo(REAL(unixT/5022440._wp-183.6045_wp,KIND=wp),2._wp*pi) ! Mean longitude (geom. not reduced, radians)
  M  = modulo(REAL(unixT/5022680._wp-182.2505_wp,KIND=wp),2._wp*pi) ! Mean anomaly (radians)
  C  = 0.03342_wp*sin(M)                  ! Equation of center (radians)
  lambda = L0+C                           ! True Longitude (radians)
  delta=asin(sin(eps)*sin(lambda))        ! Declination (radians)

  ! ---- Equation of time (minutes)
  yfac = tan(eps/2._wp)**2
  EoT_min = 720/pi*(yfac*sin(2*L0)-2*e*sin(M)+ &
    4*e*yfac*sin(M)*cos(2*L0)-yfac*yfac*sin(4*L0)/2-1.25_wp*e*e*sin(2*M))

  ! ---- Local true solar time and hour angle
  !   TST (min) = 60*UT + EoT + 4*lon
  Hr=MOD(REAL(MOD(unixT,REAL(86400._wp,KIND=int64)),KIND=wp)/60+ &
    EoT_min+4*long_deg,1440._wp)
  ! Handle wrapping: if Hr is negative after mod, add 1440
  if(Hr.lt.0._wp) Hr=Hr+1440._wp
  ! Hour angle: 0 at solar noon (720 min), negative before noon, positive after
  Ha=pi*(Hr/720._wp-1._wp)

  ! ---- Topocentric altitude and azimuth
  phi = lat_deg*pi/180._wp
  elevation = asin(sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(Ha))
  azimuth = atan2(sin(Ha),cos(Ha)*sin(phi)-tan(delta)*cos(phi))+pi

end subroutine sun_position

!-----------------------------------------------------------------------
! compute_ghi estimates the global horizontal irradiance (shortwave
! radiation flux [W/m^2]) at a given time and place.
! Notes:
!   - Clear-sky baseline: Haurwitz GHI_cs ≈ 1098*cosZ*exp(-0.059/cosZ).
!   - AOD/PW/O3 tweaks are simple Beer-type broadband factors (calibratable).
!   - Clouds via Kasten–Czeplak CMF = 1 - 0.75*C^3.
!-----------------------------------------------------------------------
subroutine compute_ghi(ghi,sea,aod500,pw,oz,cloud_frac)
  use gridlist_variables, only : wp
  use constants, only : pi
  implicit none

  ! Local Variables
  real(wp),intent(in) :: sea        ! solar elevation angle [radians; 0 = sun at horizon]
  real(wp),intent(in) :: aod500     ! Aerosol optical depth at 500-550 nm (>=0)
  real(wp),intent(in) :: pw         ! Precipitable water [m]
  real(wp),intent(in) :: oz         ! Total column ozone [m-atm ~ 0.3 cm = 300 DU]
  real(wp),intent(in) :: cloud_frac ! Total cloud fraction [0-1]  
  real(wp),intent(inout) :: ghi     ! GHI on a horizontal plane [W/m^2]
    
  real(wp) :: ghi_cs, m_air
  real(wp) :: t_aero, t_water, t_ozone, t_cloud
  
  ! Executable Code
  if(sea.le.0._wp)then
    ghi=0._wp
  else
    ! ---- Clear-sky baseline (Haurwitz) ----
    ghi_cs=1098._wp*sin(sea)*exp(-0.059_wp/sin(sea))

    ! ---- Airmass (Kasten–Young 1989) from zenith only ----
    ! Convert sza from radians to degrees for the formula
    m_air =1._wp/(sin(sea)+0.50572_wp* &
      ((6.07995_wp+sea*180._wp/pi)**(-1.6364_wp)))

    ! ---- Broadband optical tweaks (very compact, tunable) ----
    t_aero  =exp(-aod500*m_air)                 ! Treat AOD500 as effective broadband tau
    t_water =exp(-0.2986_wp*(pw*m_air)**0.3_wp) ! Empirical water vapor attenuation (kept modest)
    t_ozone =exp(-4._wp*oz*m_air)               ! Ozone absorption with effective k ≈ 0.04 cm^-1
    t_cloud =1._wp-0.75_wp*cloud_frac**3        ! Clouds: Kasten-Czeplak CMF
    
    ghi=ghi_cs*t_aero*t_water*t_ozone*t_cloud
  endif

end subroutine compute_ghi
