c234567**********************************************
      subroutine con
c this subroutine sets the values of the constants that
c are going to be used throughout HIGRAD and FIRETEC
      use gridsetup
      use constants
      use nonlocal
      use turba, only :diffcst,kbcratio
!      use turba, only :rturbprandtl
      Implicit None

      r13=0.33333            ! fraction 1/3
      rnfuel=0.4552      ! This is the stoichiometric coef. for wood in 
c                          simple wood model  (Drysdale pp.179) 
      rno=0.5448         ! This is the stoichiometric coef. for wood in 
c                          simple wood model  (Drysdale pp.179)
c      cfhydro=.50            ! this was the value before correction of tke boundary condition
c     cfhydro and cfchar are equivalant to CF *rhoref*sx**2*.09
      cfhydro=.30*3.        ! Coef. of mass reaction rates.  (found through 
                             !      numerical experiments with grass)
                             ! cfhydro of .8 seems a little low with code of
                             ! 5/24/01, flat ground, no wind, rhof=1 kg/m^3 
                             ! and 5% moisture
      cfchar=.015*6.      ! Coef. of mass reaction rates.  (found through 
      rhohydrothresh=.4
      rhoref=1.          ! reference density for use in reaction rates (kg/m**3)
      sx=0.005            ! length scales for mixing reactions          (m)
      fueldepth=1.6     ! value for tall grass
      rkmin=1.e-06       ! minimum turbulence level (to prevent neg. values)
      rkmax=80.          ! maximum turbulence level (to diagnose blow-up)
      c1a=1.             ! turbulence coef.
      tc=300.            ! tcold for FIRETEC bouy. turbulence creation (K)
      hf=8913.48e3         ! heat of reaction for simple wood  (J/Kg of products)
      thetag=.75        ! fraction of energy that is deposited in the gas phase
                        ! for the computation of thetasolid
      ss=0.0005            ! default size scale for the solids (small fuels)     (m)
      tfire=1000.         !initial average temp (k) of fire

c  atmospheric parameters
      ! moved to gridlist
      !if (ilocation.eq.0) then 
      !  tambient=290.      ! Ambient temperature                         (K)
      !  pressground=1.0e5  ! Pa
      !  zgroundref=0.   ! reference elevation for tambient and pressground
      !else if (ilocation.eq.1) then
      !  tambient=300.      ! Ambient temperature                         (K)
      !  pressground=0.86e5  ! Pa
      !  zgroundref=0.   ! reference elevation for tambient and pressground
      !endif
c radiation parameters
      rke=2.             ! radiation emmision coefficient              
      sigma=5.678e-8    ! Stefan-Boltzmann constant               (W/(m^2K^4))
      rsourceht=200.     ! maximum height the the radiation is tracked (m)
      ep=1.e-10          ! epsilon (min threshold for various things)

c material parameters
      gamma=1.4          ! cp/cv for ambient air at ~ 350 K
      gammo=1.38         ! cp/cv for oxygen at ~ 350 K
      !  cpwood and rhomicro now defined in gridlist
      !if (iinra.eq.0) then
      !   cpwood=2500.         ! cp wood (estimated from Incorpera, Dewitt) (J/KG/K)
      !   rhomicrovalue=500.
      !else if (iinra.eq.1) then
      !   cpwood=1800.         ! cp jack pine needles(estimated n Albini and stocks) (J/KG/K)
      !   rhomicrovalue=700.
      !end if
      cpwater=4200.        ! cp for water (KJ/Kg/K)
      hwevap=2257.e3       ! energy needed for water evap at 373 K (Moran)(J/Kg)
      gammav=1.7        ! cp/cv for water vapor
      twvap=373.         ! temp. of water evap. for amb. cond (sea level) (K)
      tcrit=600.         ! temp. at which wood begins to pyrolyze
      tramp=500.        ! the tail of the ramp straight line (if extrapolated to psi=0)
      psijoint=.30        ! height of the ramp joint
      tfstep=310.        !tail of the ramp curve (gives marginal fire with no 
      tstep=(373-tambient)*(tfstep-tambient)/(tcrit-tambient)  
                    ! (above) for 2 m res. step used in water ramp functions           (K)
      cvvapor=1160       ! cv for water vapor  (rrl 9/9/99)          (J/Kg/K)
c                           estimated based on dI/dT at 373 K (Moran Shapiro)
      cvoxygen=670.      ! cv for oxygen (rrl 9/9/99)                (J/Kg/K)
c                           estimated for ideal gas (Moran, Shapiro)  (J/Kg/K)

      rv=461.        ! gas constant of moist air                     (J/Kg-K)
      t00=273.16     ! freezing temperature                          (K)
      ee0=611.       ! saturated vapor pressure                      (Nt/m2)
      hlat=2.5e6     ! latent heat of condensation                   (J/Kg)
      g=9.8          ! gravity                                       (m/s^2)
      rg=287.04      ! gas constant for dry air                      (J/Kg-K)
      prndt=1.       ! Prandtal number
      !bv=sqrt(st*g)  ! Bruint-Vas frequency
      
c  nonlocal constants
      
      rnhc=0.68
      
      rng=0.404                 ! This is the stoichiometric coef. for combustable gas
                                ! in the nonlocal model
      rnonl=0.596               ! This is the stoichiometric coef. for oxygen
                                ! in the nonlocal model
      cg=10.0                  ! Taken from Rod's thesis pg 60 (0.088)

      hfgas= 8789.245e3         ! This is the heat of combustion of our generic
                                ! gas in the nonlocal model (J/kg) products.
                                
      hfsolid=-200.e3           ! this is the heat of reaction for k3 (J/kg)
                                ! I am following the firetec convention that
                                ! positive values of hf are exothermic
      
      cv=717        ! specific heat of air at constant volume        (J/Kg/K)
      cp=rg+cv      ! specific heat of air at constant pressure     (J/Kg-K)
      cap=rg/cp
      pi=acos(-1.)   ! pi
      prrcp=1.e5**(rg/cp)
      fcor2=fcr0*cos(pi/180.*ang)
      fcor3=fcr0*sin(pi/180.*ang)
      
      !TODO : check turbulent prandlt number with JLD
      ! in Mell et al 2009, pr=0.5, now defined in gridlist
      !rturbprandtl=2.0 !0.66 !1.43  !inverse of the turbulent prandtl number
      diffcst=0.09
      kbcratio=0.2  ! kb/kc
      !bv=sqrt(st*g)
      return
      end
      !******************************************************************!

      !*****************************************************!
      subroutine set_fueldepth()
      use gridsetup
      use constants
      use turba
      use msga
      use metryic
      use fireteca, only:min_rhof
      Implicit None


      integer :: i,j,k
      integer:: lfueltmp
      real :: max_z,tmp_z
      real,external :: zcart

      min_rhof = 1.0e-6  
      lfueltmp=2
      max_z = 2.5     !default minimum fueldepth = 2.5 meters 
      !max_z = 29.5     !default minimum fueldepth = 2.5 meters 
c      max_z = 3.8     !default minimum fueldepth = 2.5 meters 
        do k=1,l
         do j=1,mp
          do i=1,np
            if(rhof(i,j,k).gt.min_rhof)then
             tmp_z = zcart(zedge(k+1),i,j)-zs(i,j)
             if(tmp_z.gt.max_z)then
              max_z = tmp_z
             endif
             if (k.gt.lfueltmp)  lfueltmp=k
            endif
          enddo
         enddo
             !write(6,*) 'k,zcart(zedge(k+1))',k,tmp_z
        enddo
        call mpi_allreduce(max_z,fueldepth,1,mpi_real,mpi_max,
     +                   mpi_comm_world,ierror)
        call mpi_allreduce(lfueltmp,lfuel,1,mpi_integer,mpi_max,
     +                   mpi_comm_world,ierror)

      end subroutine set_fueldepth
      !*****************************************************!

