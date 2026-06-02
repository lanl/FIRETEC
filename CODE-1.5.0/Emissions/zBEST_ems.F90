!-----------------------------------------------------------------------
! ZBEST (Zonal-Based Emissions Source Term Model) calculates forcing 
! terms for source emissions with a physics-based estimation where 
! combustion properties are averaged throughout a cell and combined 
! to do a 1-Dimensional emission evolution
!-----------------------------------------------------------------------
subroutine ZBEST(forcing,xv,frhof,temps,frhow,rhof,il,iu,jl,ju,lls, &
    nvp,nfuel,dt)
  use gridlist_variables, only : rhoMicro,dx,dy,prec
  use xvall, only : iM0,iM1,iuvel,ivvel,iwvel,iO2,irho,iH2O
  use emission_general_variables, only : spEmit
  use fuel_variables, only : rnfuel,rno,hf,sizeScale
  use constants, only : pi,Na,Rgas,kB
  use thermo_variables, only : pr
  use metric_variables, only : zedge
  use zcart_function, only : zcart
  use computePsi_function, only : computePsi
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nfuel,nvp
  real(prec),intent(inout) :: forcing(il:iu,jl:ju,lls,nvp)
  real(prec),intent(in)    :: frhof(nfuel,il:iu,jl:ju,lls)  ! Consumption rate of fuel
  real(prec),intent(in)    :: temps(nfuel,il:iu,jl:ju,lls)  ! Consumption rate of fuel
  real(prec),intent(in)    :: frhow(nfuel,il:iu,jl:ju,lls)  ! Consumption rate of moisture
  real(prec),intent(in)    :: rhof(nfuel,il:iu,jl:ju,lls)   ! Bulk fuel density
  real(prec),intent(in)    :: xv(il:iu,jl:ju,lls,nvp)
  real(prec),intent(in)    :: dt

  integer :: ift,i,j,k
  real(prec) :: TempI=600.      ! Incipient Zone Temperature
  real(prec) :: TempF=1200.     ! Flame Zone Temperature
  real(prec) :: TempR=1800.     ! Reaction Zone Temperature
  real(prec) :: SP=0.56         ! Fuel's Sooting Potential
  real(prec) :: Zst=0.135       ! Stoichiometric mixture fraction (C_33 H_48 O_19)
  real(prec) :: AO2=7.98E-1     ! pre-exponential constant (kg K^0.5/Pa m2 s)
  real(prec) :: EO2 = -1.77E8   ! activation energy (J/kmole)
  real(prec) :: AOH = 1.89E-3   ! pre-exponentail constant (kg K^0.5/Pa m2 s)
  real(prec) :: mC  = 12.01     ! Molecular weight of carbon (kg)
  real(prec) :: xtol = 0.111    ! 0.132 mass fraction   # mole fraction of toluene surrogates in the tar
  real(prec) :: xnaph= 0.342    ! 0.457 mass fraction   # mole fraction of naphthalene surrogates in the tar
  real(prec) :: rhoV = 0.45     ! Density of volatiles (kg/m3)
  real(prec) :: eps = 2.2*0.001 ! Van der Waals enhancement factor*Collision Efficiency
  real(prec) :: rhost= 1850.0   ! Particle density (kg/m3)
  real(prec) :: Diff,Deq,Q,dstar,nu,VolC,ff
  real(prec) :: yH2,yCH4,yCO2,yCO,SP_mtar,SP_xtar,Ntar
  real(prec) :: PO2r,POHr,Oxidation
  real(prec) :: speed,gam,krxn,df,psif!,scalar_beta
  !real(prec) :: shapeFactor
  real(prec) :: rtF,rtR
  real(prec) :: N,M
  real(prec) :: EF

  ! Executable Code
  ! Flame Characteristic constants
  Diff=2./3.*(kB*TempR/pi)**1.5/(101325.*(350.e-12)**2.*sqrt(3.*mC/Na))    ! Average diffusivity of species in reaction zone (m2/s)
  nu  =(1.7616E-5*(TempR/273.16)**1.5*(273.16+110.4)/(TempR+110.4))/rhoV ! Kinematic viscosity (m2/s)
  
  yH2     = (1.-SP)*0.018 ! Mass fraction of H2
  yCH4    = (1.-SP)*0.085 ! Mass fraction of CH4
  yCO2    = (1.-SP)*0.298 ! Mass fraction of CO2
  yCO     = 1.-SP-yH2-yCH4-yCO2 ! Mass fraction of CO
  SP_mtar = 128.*xnaph+92.*xtol+94.*(1.-xnaph-xtol) 
  SP_xtar = SP/SP_mtar/(SP/SP_mtar+yH2/2.+yCH4/16.+yCO2/44.+yCO/28.) ! Mole fraction of precursors
  
  do ift=1,nfuel
    do i=il,iu
      do j=jl,ju
        do k=1,lls
          ! Do we bother with this?
          if(frhof(ift,i,j,k).le.0) cycle   
         
          ! Oxygen Consumption
          forcing(i,j,k,iO2)=forcing(i,j,k,iO2) &
              -frhof(ift,i,j,k)*rno/rnfuel
        
          !*************************************************************
          ! Flame characteristics
          !*************************************************************
          VolC=dx*dy*(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)) ! Cell volume (m3)
          ff  =frhof(ift,i,j,k)*dt*VolC ! Mass of raw fuel consumed (kg)
        
          ! Total flame length (Heskestad 1983 Luminous heights of turbulent diffusion flames)
          Deq  =sqrt(4.*ff/(pi*rhoMicro*sizeScale(ift,i,j,k))) ! Diameter equivalence of the fuel bed (m)
          Q    =ff*hf/dt                                       ! Energy release rate for burning fuels (J/s)
          dstar=0.01476*Q**0.4-1.02*Deq                        ! Total length of the flame (m)
        
          ! Reaction zone thickness (Bilger 1976 Reaction zone thickness and formation of nitric oxide in turbulent diffusion flames)
          speed=sqrt(xv(i,j,k,iuvel)**2.+xv(i,j,k,ivvel)**2. &
            +xv(i,j,k,iwvel)**2.)/xv(i,j,k,irho)
          gam =nu/2.*(speed/dstar)**2.              ! Kolmogoroff strain rate (1/s)
          psif=computePsi(temps(ift,i,j,k))
          krxn=frhof(ift,i,j,k)/(dt*psif*rhof(ift,i,j,k))! Oxidation rate constant
          df  =min(dstar,sqrt(Diff/gam)*(gam/krxn)**(1./3.))    ! reaction zone thickness (m)
          !scalar_beta=sqrt(2.*Diff*(((0.1392/(rnfuel*0.455)-0.1667 &
          !  /(rno*0.032))-(0.0/(rnfuel*0.455)-0.1936/(rno*0.032))) &
          !  /df)**2./gam)
          !df  =min(dstar,df*(rno*rnfuel*scalar_beta)**(-1./3.)) ! corrected reaction zone thickness (m)
    
          ! Zone residence times
          rtF =(dstar-df)/speed ! Flame residence time (s)
          rtR =df/speed         ! Reaction zone residence time (s)
          
          !*************************************************************
          ! Soot Inception
          !*************************************************************
          ! Rates of precursor reacting
          Ntar=pr(i,j,k)*SP_xtar/(kB*TempI) ! Concentration of precursors (#/m3)
          N   =Ntar*(xnaph+xtol)/4. !!! ASSUMPTION !!! naphthalene and toluene precursors nucleate while phenol crack (#/m3)
          M   =N*4.*SP_mtar/Na ! kg/m3
    
          !*************************************************************
          ! Flame zone mechanisms
          !*************************************************************
          ! Coagulation
          N = (N**(-5./6.)+eps*5./6.*(6./rhost)**(2./3.)* &
            sqrt(8.*kB*TempF)*(M/pi)**(1./6.)*rtF)**(-6./5.)
        
          !*************************************************************
          ! Reaction zone mechanisms
          !*************************************************************
          ! Oxidation Kinetics
          PO2r=pr(i,j,k)*0.1775*(xv(i,j,k,iO2)/0.233*(1.-Zst/2.)+Zst/2.)
          POHr=pr(i,j,k)*2.7e-4*(xv(i,j,k,iO2)/0.233*(1.-Zst/2.)+Zst/2.)
          Oxidation=(AO2*PO2r*exp(EO2/Rgas/TempR)+AOH*POHr)/sqrt(TempR)
          M = (max(0.,M**(1./3.)-Oxidation/3.*(pi*N)**(1./3.)* &
            (6./rhost)**(2./3.)*rtR))**3.
          EF=SP*M/(Ntar*SP_mtar/Na)

          if(EF.gt.0)then
            forcing(i,j,k,iM1)=forcing(i,j,k,iM1)+EF*frhof(ift,i,j,k)
            forcing(i,j,k,iM0)=forcing(i,j,k,iM0)+ &
              EF*frhof(ift,i,j,k)*N/M
          endif
        enddo
      enddo
    enddo
  enddo
  
  if(ANY(spEmit.eq."H2O"))then
    do ift=1,nfuel
      do i=il,iu
        do j=jl,ju
          do k=1,lls
            forcing(i,j,k,iH2O)=forcing(i,j,k,iH2O) &
              +frhof(ift,i,j,k)*0.56+frhow(ift,i,j,k)
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine ZBEST
