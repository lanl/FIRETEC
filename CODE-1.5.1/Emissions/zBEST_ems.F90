!-----------------------------------------------------------------------
! ZBEST (Zonal-Based Emissions Source Term Model) calculates forcing 
! terms for source emissions with a physics-based estimation where 
! combustion properties are averaged throughout a cell and combined 
! to do a 1-Dimensional emission evolution
!-----------------------------------------------------------------------
subroutine ZBEST(forcing,xv,frhof,temps,rhof,ilh,iuh,jlh,juh,ku, &
    il,iu,jl,ju,lls,nvp,nfuel,dt)
  use gridlist_variables, only : dx,dy,wp
  use vegetation_gridlist_variables, only : rhoMicro
  use xvall, only : iM0,iM1,iuvel,ivvel,iwvel,iO2,irho,iH2O
  use emission_general_variables, only : spEmit
  use fuel_variables, only : rnfuel,rno,hf,sizeScale
  use constants, only : pi,Na,Rgas,kB
  use thermo_variables, only : pr
  use metric_variables, only : zedge
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,ku,il,iu,jl,ju,lls,nfuel,nvp
  real(wp),intent(inout) :: forcing(ilh:iuh,jlh:juh,ku,nvp)
  real(wp),intent(in)    :: frhof(nfuel,il:iu,jl:ju,lls)  ! Consumption rate of fuel
  real(wp),intent(in)    :: temps(nfuel,il:iu,jl:ju,lls)  ! Consumption rate of fuel
  real(wp),intent(in)    :: rhof(nfuel,il:iu,jl:ju,lls)   ! Bulk fuel density
  real(wp),intent(in)    :: xv(ilh:iuh,jlh:juh,ku,nvp)
  real(wp),intent(in)    :: dt

  integer :: ift,i,j,k
  real(wp) :: TempI=600.0_wp      ! Incipient Zone Temperature
  real(wp) :: TempF=1200.0_wp     ! Flame Zone Temperature
  real(wp) :: TempR=1800.0_wp     ! Reaction Zone Temperature
  real(wp) :: SP=0.56_wp         ! Fuel's Sooting Potential
  real(wp) :: Zst=0.135_wp       ! Stoichiometric mixture fraction (C_33 H_48 O_19)
  real(wp) :: AO2=7.98E-1_wp     ! pre-exponential constant (kg K^0.5/Pa m2 s)
  real(wp) :: EO2 = -1.77E8_wp   ! activation energy (J/kmole)
  real(wp) :: AOH = 1.89E-3_wp   ! pre-exponentail constant (kg K^0.5/Pa m2 s)
  real(wp) :: mC  = 12.01_wp     ! Molecular weight of carbon (kg)
  real(wp) :: xtol = 0.111_wp    ! 0.132 mass fraction   # mole fraction of toluene surrogates in the tar
  real(wp) :: xnaph= 0.342_wp    ! 0.457 mass fraction   # mole fraction of naphthalene surrogates in the tar
  real(wp) :: rhoV = 0.45_wp     ! Density of volatiles (kg/m3)
  real(wp) :: eps = 2.2_wp*0.001_wp ! Van der Waals enhancement factor*Collision Efficiency
  real(wp) :: rhost= 1850.0_wp   ! Particle density (kg/m3)
  real(wp) :: Diff,Deq,Q,dstar,nu,VolC,ff
  real(wp) :: yH2,yCH4,yCO2,yCO,SP_mtar,SP_xtar,Ntar
  real(wp) :: PO2r,POHr,Oxidation
  real(wp) :: speed,gam,krxn,df,psif!,scalar_beta
  real(wp) :: rtF,rtR
  real(wp) :: N,M
  real(wp) :: EF
  real(wp),external :: zcart

  ! Executable Code
  ! Flame Characteristic constants
  Diff=2._wp/3._wp*(kB*TempR/pi)**1.5_wp/ &
    (101325._wp*(350.e-12_wp)**2._wp*sqrt(3._wp*mC/Na))    ! Average diffusivity of species in reaction zone (m2/s)
  nu  =(1.7616E-5_wp*(TempR/273.16_wp)**1.5_wp*(273.16_wp+110.4_wp)/ &
    (TempR+110.4_wp))/rhoV ! Kinematic viscosity (m2/s)
  
  yH2     = (1.0_wp-SP)*0.018_wp ! Mass fraction of H2
  yCH4    = (1.0_wp-SP)*0.085_wp ! Mass fraction of CH4
  yCO2    = (1.0_wp-SP)*0.298_wp ! Mass fraction of CO2
  yCO     = 1.0_wp-SP-yH2-yCH4-yCO2 ! Mass fraction of CO
  SP_mtar = 128._wp*xnaph+92._wp*xtol+94._wp*(1._wp-xnaph-xtol) 
  SP_xtar = SP/SP_mtar/ &
    (SP/SP_mtar+yH2/2._wp+yCH4/16._wp+yCO2/44._wp+yCO/28._wp) ! Mole fraction of precursors
  
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        do ift=1,nfuel
          ! Do we bother with this?
          if(frhof(ift,i,j,k).le.0._wp) cycle   
         
          ! Oxygen Consumption
          forcing(i,j,k,iO2)=forcing(i,j,k,iO2) &
              -frhof(ift,i,j,k)*rno/rnfuel
        
          !*************************************************************
          ! Flame characteristics
          !*************************************************************
          VolC=dx*dy*(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)) ! Cell volume (m3)
          ff  =frhof(ift,i,j,k)*dt*VolC ! Mass of raw fuel consumed (kg)
        
          ! Total flame length (Heskestad 1983 Luminous heights of turbulent diffusion flames)
          Deq  =sqrt(4._wp*ff/(pi*rhoMicro*sizeScale(ift,i,j,k))) ! Diameter equivalence of the fuel bed (m)
          Q    =ff*hf/dt                                       ! Energy release rate for burning fuels (J/s)
          dstar=0.01476_wp*Q**0.4_wp-1.02_wp*Deq                        ! Total length of the flame (m)
        
          ! Reaction zone thickness (Bilger 1976 Reaction zone thickness and formation of nitric oxide in turbulent diffusion flames)
          speed=sqrt(xv(i,j,k,iuvel)**2.0_wp+xv(i,j,k,ivvel)**2.0_wp &
            +xv(i,j,k,iwvel)**2.0_wp)/xv(i,j,k,irho)
          gam =nu/2.0_wp*(speed/dstar)**2.0_wp              ! Kolmogoroff strain rate (1/s)
          call computePsif(psif,temps(ift,i,j,k))
          krxn=frhof(ift,i,j,k)/(dt*psif*rhof(ift,i,j,k))! Oxidation rate constant
          df  =min(dstar,sqrt(Diff/gam)*(gam/krxn)**(1.0_wp/3.0_wp))    ! reaction zone thickness (m)
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
          N   =Ntar*(xnaph+xtol)/4.0_wp !!! ASSUMPTION !!! naphthalene and toluene precursors nucleate while phenol crack (#/m3)
          M   =N*4.0_wp*SP_mtar/Na ! kg/m3
    
          !*************************************************************
          ! Flame zone mechanisms
          !*************************************************************
          ! Coagulation
          N = (N**(-5._wp/6._wp)+eps*5._wp/6._wp*(6._wp/rhost)**(2._wp/3._wp)* &
            sqrt(8._wp*kB*TempF)*(M/pi)**(1._wp/6._wp)*rtF)**(-6._wp/5._wp)
        
          !*************************************************************
          ! Reaction zone mechanisms
          !*************************************************************
          ! Oxidation Kinetics
          PO2r=pr(i,j,k)*0.1775_wp*(xv(i,j,k,iO2)/0.233_wp* &
            (1._wp-Zst/2._wp)+Zst/2._wp)
          POHr=pr(i,j,k)*2.7e-4_wp*(xv(i,j,k,iO2)/0.233_wp* &
            (1._wp-Zst/2._wp)+Zst/2._wp)
          Oxidation=(AO2*PO2r*exp(EO2/Rgas/TempR)+AOH*POHr)/sqrt(TempR)
          M = (max(0._wp,M**(1._wp/3._wp)- &
            Oxidation/3._wp*(pi*N)**(1._wp/3._wp)* &
            (6._wp/rhost)**(2._wp/3._wp)*rtR))**3
          EF=SP*M/(Ntar*SP_mtar/Na)

          if(EF.gt.0.0_wp)then
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
              +frhof(ift,i,j,k)*0.56_wp
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine ZBEST
