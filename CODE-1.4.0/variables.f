c2345678***************************************************
      module gridsetup
      Implicit None

      integer::j3=1,irst,nt,nts,ibctopbot
c              !j3-1 is 3-d and 0 is 2-d
c                  !irst is
c                          !nt is the number of time steps 
c                             !nts is the number of subtime steps
      integer::isplit,nprocx,nprocy,n,m,l,lfuel,nproc
c              !isplit is
c                    !number of processors in the x direction
c                           !number of processors in the y direction
c                                  !number of real cells in x dir.
c                                    !number of real cells in y dir.
c                                      !number of real cells in z dir.
c                                      ! number of cells for definition of fuel 
c                                        !number of proc.
      integer::nv
c             !number of variables in the xvb array
c                !number of 
      real::dx,dy,dz,aa1,dts,time=0, restarttime
c                   !the subcycle time step
c                      ! time in seconds
      real::zab,zabt,tow
c          !vertical minimum distance from ground of the absorber  (m)
c          !vertical minimum distance from ground of the absorber for rho and theta (m)
c              !relaxation time constant                           (s)
      integer::ntp,ibcx,ibcy,ih,irlx,irly,ibclatopen,ibctopopen
c             !number of subcycles for FIRETEC (for reactions~.5nts)
c                 !flags for cyclic bc's (1 is turned on)
c
c                           !number of ghost points subdomains
c                              !relaxation bc's in x
c                                   !relaxation bc's in y direction
c ! if ibclatopen is 1,east,north and south  are open for u or v (exit only)
c ! if ibclatopen is 2,east,north and south  are open for u or v (exit and come in)
c ! if ibctopopen is 1 top is open (w is relaxed to 0)
      integer::ibxo,ibyo
c             !flax for relaxation bc's for MPDATA
c                  !flag for relaxation bc's for MPDATA (y direction)
      integer::nr,iab
      integer::irod,iturb,irhovapor,inonlocal,isubgridgas
c                  !switch for FIRETEC
c                       !switch for turbulence model
c                     irad : switch for radiation heat transfer <- moved to radiation.mod. KOO
c                                  !switch for transport of H20 vapor
c										!switch for nonlocal chemistry
c											!switch for subgrid gas chemistry
      integer:: isootmodel=0  ! 0 is rrl' original soot model, 1 is inra's soot model
      integer::idrag,isa !drag
      integer::ixevariation
      ! tell is xe is variating with time
      integer::itheta,iperturb
              ! selection of potential temperature profile
                ! selection of perturbation of initial state (theta, pinwheel...)
      integer::icorio,ilspgf,frqlspgf,izlspgf,islip
c                   !flag for calculation of coriolis (1 from top, 2 from xe)
c		  ! flag for large scale pressure gradient force (1 ajusted homogeneously along x and y, 
c                                  2 ajustment vary along the perpendicular to wind to limit streaks, 3 smoothing of the mode 2)
c                   !frq for updating large scale pressure gradient force (ilspgf>=1)
c                 ! izlspgf = 0 compute mass flow over the vertical, izlspgf=1 : compute at height zu (for target wind speed)
c                          !islip =1 is no slip bcs at z=0
      integer :: itrestart,ittot
      integer::ipotflow ! call for potflow
      real::slopeangle  ! slope in degree for rotated gravity
      real::slopeazimuth  ! slope azimuth in degree ifor rotated gravity (0: the slope is onx, 90: the slope is on y)
      integer::uswitch,vswitch
      real :: uramp,vramp,uramptime,vramptime
      real::u0,v0 !,st
      real::zu,uprofilezu ! height of reference wind speed for uswich=2  
      integer::ius=0! extension of reference zone in the domain for uswitch=2 (default is 1)
      integer::iue=0! extension of reference zone in the domain for uswitch=2 (default is n)
      integer::jus=0 !! extension of reference zone in the domain for uswitch=2 (default is 1
       integer::jue=0 ! extension of reference zone in the domain for uswitch=2 (default  is m)
       integer::windspeedupfactor=1 ! ratio between wind time step versus fire time step
      integer::iord,isor,nonos,idiv,nfct,nonosold
      integer::impdataold=1 ! FP2019: decide if we want the old mpdata
      integer::ifuel,ifuelinra,ivegread,kmax,kf77max
      integer::fuelinranumber 
c          !pdf selector
c             !fuel type specifier
                 !read fuels data from trees.f (1) or not (0)
c                   !maximum number o partitions in a cell for Voff
c                        !switch for fortran 77
      integer::nx,ny,nz
c             !first real x cell in radiation super grid
c                !first real y cell in radiation super grid
c                   !first real z cell in radiation super grid
      integer::nnx,mny,Lnz
c             !last real x cell in radiation super grid
c                 !last real y cell in radiation super grid
c                     !last real z cell in radiation super grid
      integer::np,mp,nm,nml,ml,npmp,npmpl,mpl
c             !number of cells per processor in the x direction
c                !number of cells per processor in the y direction
c                   !number of real cells in horizontal plane
c                      !number of total cells
c                          !number of cells in plane perp to x direction
c                             !number of cells in horiz. plane per processor
c                                  !total number of cells per processor
c                                        !number of cells in plane perp.
c                                         to x per processor

      real::dt,dtp,dxi,dyi,dzi,dti
c          !large time step (for the main code)
c             !time scale for some of the physics in FIRETEC
c                 !1/dx                                            (1/m)
c                     !1/dy                                        (1/m)
c                         !1/dz                                    (1/m)
c                             !1/dt                                (1/s)
   
      real::gc1s,gc2s,gc3s
c          !gc1s is dt/dx for the MOA subcycles      (s/m)
c               !gc2s is dt/dy for the MOA subcycles      (s/m)
c                    !gc3s is dt/dz for the MOA subcycles      (s/m)


      real::gh1s,gh2s,gh3s
c          !gh1s is .5*dt/dx for the MOA subcycles                 (s/m)
c               !gh2s is .5*dt/dx for the MOA subcycles            (s/m)
c                    !gh3s is .5*dt/dx for the MOA subcycles       (s/m)

      real::gc1,gc2,gc3
c          !gc1 is dt/dx                                           (s/m)
c              !gc2 is dt/dy                                       (s/m)
c                  !gc3 is dt/dz                                   (s/m)

      real::gh1,gh2,gh3
c          !gc1 is .5*dt/dx                                        (s/m) 
c              !gc2 is .5*dt/dy                                    (s/m)
c                  !gc3 is .5*dt/dz                                (s/m)

      integer::frqoutput
c             !frequency of outputs                             (1/timesteps)
      integer::frqoutputsub ! frequency of output on a subdomain (1/timesteps)      
      integer::frqfilstr !frequency of filtering                             (1/timesteps)
      integer:: issub,jssub,nisub,njsub,nksub 
      ! start index of the subdomain, length of subdomain in number of cell
      integer::icfmeflag
               !icfmeflag switches icfme functionality on if 1, off if 0 
               !JAS 1/20/06

      integer:: iunstable
c               !switch to turn on heat source at the ground
      integer:: ioextra ! entails to have extra io not required for restart...
      integer::iwindfieldout,iwindfieldin,itwindfield,itinterp,
     &         ibcells,jbcells,is,ie,js,je
c             !run a windfield      JMC 8/15/6
c                         !timestep to start files
c                                    !delta timestep for windfield interpolation
c             !number of x boundary cells
c                     !number of y boundary cells
c                             !start a fire JMC 8/15/6
c         coordinates of the candidate zone for data extraction
      
      ! Legacy flags to be cleaned up
      integer::idirt
        
      ! personnal flags
      integer::ifp                 ! francois 's personal flag
 
      ! BELOW THIS LINE:
      ! personnal variable parameter set
      ! parameters associated to iheatsource 
      integer::iheatsource 
       ! 0 no heat source
       ! 1 heat source that does not burn the fuel
       ! 2 heat source that does burn the fuel for k>=2
      real::hsros ! spread rate of heat source (along the axis
      real::hsint ! intensity of heat source (along the axis) kW/m
      integer::ihsmass ! mass source flag associated to heat source
      integer::lHeatSource ! height of heat source in number of cells
      end module gridsetup 
!2345678***************************************************
      module firebranda                                !KOO

      Implicit None
      real, allocatable:: brand(:,:),fb_imm(:,:)
!                          (nbmax,FBDIM)   
      real, allocatable:: fb_ign_time(:,:,:)
!                         starting time [sec] of ignition by firebrands
      integer, allocatable:: ibrand(:,:),ifb_imm(:,:)
!                              (nbmax,IFBDIM)
      integer, allocatable:: nbrand(:,:,:),nfb_imm(:,:,:),nbproc(:)
!                                  (n,m,l)                      mpi
      integer, parameter::FBDIM=20, IFBDIM=6
! FBDIM: 1=x,  2=y,  3=z,  4=V_x, 5=V_y, 6=V_z, 
!        7=time,  8=radi,  9=tck, 10=temp 
!       11=xo,12=yo,13=zo,14=density,15=densityo,16=anglea
!       17=timeo,18=radio,19=tcko,20=tempo
!       21=....
! IFBDIM: 1=globalID, 2=b_status, 3=b_shape, 4=mpi_rank, 5=x-cycle,
! 6=y-cycle
!       b_status: 0=burnout
!                 1=in flight, 2=hit upper b'ry 
!                 3=landingnded on no fuel, 4=hit domain b'ry 
!                 5=landed on fuel (effective)
!                 6=landed on burning site (ineffective)
!                 11-18=be sent to other processor 
!                (temporary -> to be 1 after recieved) 
!                   
!       b_shape:  0=dsk_ub, 1=dsk_dh/dt, 2=dsk_dr/dt
!                 3=cyl_ub, 4=cyl_dh/dt, 5=cyl_dr/dt
!                 6=sph_ub, 7=sph_dr/dt  8=massless particles 
!  expended to 3D- z cordinate for streathcing 
       real, allocatable::  zposition(:,:,:),cl_hgt(:,:,:), dzk(:,:,:)
!                                    (ix,iy,iz), 
       integer ifbrand,nbmax,nb_imm,nb_crt,nb_ign
       integer nb_per_cell,fb_start,lau_frq
       integer lau_low_limit,lau_up_limit,ishape
       integer ivar_fbsize,iunit1,iunit2,iland_fuel
!! ivar_fbsize: 1=various size, 0=fixed size
       real Cd_dn,Cd_cn,Cd_sp,tck,radi,aka,lrrat_d,lrrat_c
       real temphot,fden_limit,rho_fb,v_in_rat
       real tck_limit_d,rad_limit_c,rad_limit_sp
       real :: burnout   ! added in fbrandlist 072414 KOO
       logical w_ter,tx_out,nfb_io 

       end module firebranda
      
c2345678***************************************************
c module bc contains boundary condition arrays for turbulence
      module bc
      Implicit None
      real, allocatable::tauw(:,:),uzs(:,:),vzs(:,:)
c    !for storing momentum fluxes for bottom boundary (turbulence)
c              !x component of the shear at z=0
c                       !y component of the shear at z=0

      end module bc
c2345678***************************************************
c module workavg contains variables that are computed in MOA (inner loop)
      module workavg
      Implicit None
      real,allocatable::  uavg(:,:,:),vavg(:,:,:),
     .                    wavg(:,:,:),oavg(:,:,:)
c                       ! uavg, vavg, wavg and oavg are cell centered average
c                         velocities computed in the inner loop   (m/s)
      real,allocatable:: f1avg(:,:,:),f2avg(:,:,:),
     .                   f3avg(:,:,:)
c                      ! f1avg, f2avg, and f3avg are cell centered forcing 
c                       funcions that are used in inner loop and in turbulence

      end module workavg

c2345678***************************************************
c module forcinner contains forcing functions for the inner loop
      module forcinner
      Implicit None
      real, allocatable::f1(:,:,:),f2(:,:,:),f3(:,:,:)
c                       !f1, f2, and f3 are edge centered forces
      !real, allocatable::rg_over_cp(:,:,:) ! for higrad
      real, allocatable::rg_over_prrcp(:,:,:) ! for higrad
      real, allocatable::cp_over_cv(:,:,:) ! for higrad
      end module forcinner

c2345678***************************************************
c module metryic contains arrays for the coordinate transformations
      module metryic
      Implicit None
      real,allocatable:: c13(:,:),
     .     c23(:,:),
     .      gi(:,:,:),
c          !1/(determinant of metric tensor)**.5`
     .      gmul(:),h(:, :,:)
      real,allocatable:: zs(:,:),zsio(:,:),x(:),
     .                                y(:),z(:),zedge(:)
c     .                               zcoords(:,:,:),zcoordsf(:,:,:) 
      real :: zb
      end module metryic

c2345678***************************************************
      module met2
      Implicit None
        real::aa1m2,aa2,aa3,zdata(100),zcrdata(100),zcoeff(100)
        integer::npoints
      end module met2
c2345678***************************************************
      module relax
      Implicit None
c      real,allocatable:: tau(:,:,:)
c      real,allocatable:: relx(:),rely(:)
      real,allocatable:: relaxxv(:,:,:,:) !FP
c		relaxation coef for i,j,k,inv inv=1 to 4 (u,v,w and others)
      end module relax

c***********************************************************
      module updatedfields
      Implicit None
      ! the following arrays are defined for k=0 to l and has their halocell updated
      ! they are defined in fieldUpdate for use in turb package 
      real,allocatable:: u(:,:,:),v(:,:,:),w(:,:,:) 
      real,allocatable:: theta(:,:,:),tke_a(:,:,:),tke_b(:,:,:) 
      real,allocatable:: rtke_abc(:,:,:),ox(:,:,:),vap(:,:,:)
      end module updatedfields
c**********************************************************

c2345678***************************************************
      module turba
      Implicit None
      ! the following arrays are defined for k=1 to l 
      ! they are defined in fieldUpdate for use in turb package 
      real,allocatable:: K_axy(:,:,:),K_b(:,:,:),K_az(:,:,:)  ! rho*vt
      real,allocatable:: sqrtG_Kxy(:,:,:),sqrtG_Kz(:,:,:)  !sqrtG=1/gi
      real,allocatable:: sqrtG_KG33(:,:,:) ! (J*K*Jt)33
      !stresrij arrays:
      real,allocatable::AB11c(:,:,:),AB12c(:,:,:),AB13c(:,:,:),
     .  AB21c(:,:,:),AB22c(:,:,:),AB23c(:,:,:),
     .  AB31c(:,:,:),AB32c(:,:,:),AB33c(:,:,:)
      real, allocatable::Ci(:,:,:),Cj(:,:,:),Ck(:,:,:)

      ! source terms
      real,allocatable:: fka(:, :,:)
      real,allocatable:: fkb(:, :,:)
      real,allocatable:: fib(:, :,:)
      real,allocatable:: sqrtk(:, :,:) 
      real,allocatable:: rkc(:, :,:) 
      real,allocatable:: rhof(:,:,:,:)
      real,allocatable:: sizescale(:,:,:,:) ! ss as an array FP
      real,allocatable:: actualfueldepth(:,:,:,:)
      real,allocatable:: rhofinitial(:,:,:,:)
      real,allocatable:: cd(:,:,:)
      real,allocatable:: sa(:,:,:),saxy(:,:,:),saz(:,:,:)
      real,allocatable:: sb(:,:,:),rhomicro(:,:,:,:) !FP added a
           !dimension to rhomicro for multifuel 
      real::sc
      integer::kfuelmax ! max height of the fuel
      integer:: isoturb
      integer:: ilapdo
      real::diffcst ! diffusion cst (=0.09)
      !real::rturbprandtl ! inverse of the turbulent prandtl number (usually between 1 and 3)
      !FP09/2019 hard coded an reasonable value for the inverse prandtl
      !number, in case it is not defined in gridlist to avoid
      !rturbprandtl=0  by default
      real::rturbprandtl=2 ! inverse of the turbulent prandtl number (usually between 1 and 3)
      real::kbcratio
      end module turba

c2345678***************************************************
      module turbb
      Implicit None
      real,allocatable::  d13(:,:,:),
     .                    d23(:,:,:)
      real,allocatable::  d13a(:,:,:), 
     .                    d12a(:,:,:),
     .                    d23a(:,:,:), 
     .                    d11a(:,:,:),
     .                    d22a(:,:,:),   
     .                    d33a(:,:,:)
      real,allocatable::  d13b(:,:,:), 
     .                    d12b(:,:,:),
     .                    d23b(:,:,:), 
     .                    d11b(:,:,:),
     .                    d22b(:,:,:),   
     .                    d33b(:,:,:)
      real,allocatable::  d13t(:,:,:), 
     .                    d12t(:,:,:),
     .                    d23t(:,:,:), 
     .                    d11t(:,:,:),
     .                    d22t(:,:,:),   
     .                    d33t(:,:,:)
      end module turbb

c2345678***************************************************
      module fireteca
      Implicit None
      real,allocatable:: convht(:,:,:,:)   !wss
      real,allocatable:: tambientarray(:,:,:)
      real,allocatable:: temps(:,:,:,:),tempg(:,:,:)
      real,allocatable::  foxb(:,:,:) ! diffusion of oxygen
      real,allocatable::  firad(:,:,:)
      real,allocatable::  eastFlux(:,:,:) ! radiative flux to east cell face
      integer,allocatable:: ifirestart(:,:,:)
      real,allocatable:: frhosiesrad(:,:,:,:)
      real,allocatable:: rmoist(:,:,:,:),
     .       rhos(:,:,:,:),
     .   rhowater(:,:,:,:),
     .    cpsolid(:,:,:,:),
     .       sies(:,:,:,:),
     .       psiwmax(:,:,:,:),
     .  frhovaporb(:,:,:)
      real :: min_rhof 
      real::xhsmin,xhsmax,yhsmin,yhsmax !iheatsource extension (in m)
      integer:: nfuel ! number of fuel types in simulation
      end module fireteca

c2345678***************************************************
      module xvi
      Implicit None
      real,allocatable:: xv(:,:,:,:)
      end module xvi

c2345678***************************************************
      module xvo
      Implicit None
      real,allocatable:: xvb(:,:,:,:)
      real,allocatable:: xvfuel(:,:,:,:,:)
      real,allocatable:: tkewght(:)
      end module xvo

c2345678***************************************************
      module xve
      Implicit None
      real,allocatable:: xe(:,:,:,:)     ! environmental field
      real,allocatable:: uprofile(:,:,:)     ! uprofile from LAI (FP)
        ! normalized at height zu
      end module xve

c2345678***************************************************
      module pres
      Implicit None
      real,allocatable:: pr(:,:,:),pre(:, :,:) 
      end module pres

c**************************************************************
      module lspgf
      real, allocatable:: massFluxTemp(:) ! mass flux in the direction of geostrophic wind (ug,vg), temporary on local domain
      real, allocatable:: massFlux(:) ! mass flux in the direction of geostrophic wind (ug,vg), after mpi_allreduce
      real, allocatable:: massFluxTarget(:) ! mass flux in the direction of geostrophic wind (ug,vg), at initialization (target)
      real :: totMassFlux, targMassFlux

      real :: rhoug,rhovg ! geostrophic wind
      real :: cosg,sing ! geostrophic wind cos and sin
      integer :: iwindx 
      ! if iwindx==1 mass flux is computed for each j (as if wind aligned with x axis) : n*vg<m*ug
      ! if iwindx==0 mass flux is computed for each i (as if wind aligned with y axis) : n*vg>m*ug
      integer :: nMassFlux ! size of massFluxArray (either n or m depending of iwindx)
      real :: tau ! time between two updates of flspgf
      real :: intsintheta ! vertical integral of sintheta
      real,allocatable:: sintheta(:,:,:) ! evaluation of v/sqrt(u2+v2)
      real,allocatable:: sinthetaf(:,:,:) ! evaluation of f*v/sqrt(u2+v2)
      real,allocatable:: flspgf(:,:) ! evaluation of f
      ! based on ekman theory to impose a pressure gradient (ilspgf)
      end module lspgf

c2345678***************************************************
      module xvbin    !JMC  10/21/6  This module is for reading 
                      !windfields
      Implicit None
      integer :: itabsold,itabsnew
      real,allocatable:: xvbdataold(:,:,:,:)
      real,allocatable:: xvbdatanew(:,:,:,:)
c      real,allocatable:: prold(:,:,:)
c      real,allocatable:: prnew(:,:,:)
      character(len=257)::fxvbdataname=''        !JMC
      character(len=257)::xvbdataname  !='xvbdata/xvbdata'        !JMC
      character(len=257)::windfieldstartfile !FIXME MJH changed (larger
!filename 1/26/18
   
      contains
      !---------------------------------------------------
      real function lin_interp(xi,xl,xh,yl,yh)
      Implicit None
      real :: xi,xl,xh,yl,yh

                lin_interp=
     &               (xi-xl)/real(xh-xl)
     &              *(yh-yl)
     &              +yl

      end function lin_interp
      !---------------------------------------------------
      end module xvbin

c2345678***************************************************

      module advo
      Implicit None
      real,allocatable:: u1(:,:,:),u2(:,:,:), u3(:,:,:) ! face contravariant velocities*gc1s/gi
      real,allocatable:: fd1(:,:,:),fd2(:,:,:), fd3(:,:,:) ! donorcell fluxes
      real,allocatable:: v1(:,:,:) ! antidiffusive advective velocities(iord.eq.2)
      real,allocatable:: v2(:,:,:) 
      real,allocatable:: v3(:,:,:) 
      real,allocatable:: pmx(:,:,:)! max and min xv (iord.eq.2)
      real,allocatable:: pmn(:,:,:)
      real,allocatable:: cp(:,:,:) ! non oscillary coeff (iord.eq.2)
      real,allocatable:: cn(:,:,:)
    
      end module advo

c2345678***************************************************
      module weights
      Implicit None
      real,allocatable:: wt(:),wtf(:)
      end module weights

c2345678***************************************************
      module constants
      Implicit None
      real :: ffparam  ! reaction rate constant
      real ::cfhydro,cfchar,rhohydrothresh
      real ::rv,t00,ee0,hlat,g !,rg
      real ::prndt,cp_air,cv_air,rg_air,prrcp_air,pi !cap,pi,prrcp !,bv
      real :: cv_vapor,cp_vapor
      real ::cp_over_cv_gas,cp_gas,rg_over_cp_gas,rg_over_prrcp_gas
      real ::rnfuel,rhoref,sx,rho,rno
      real ::tramp,tfstep,aramp,bramp
      real ::psijoint,tjoint,tjoint2
      real ::fueldepth
      real ::rkmin,rkmax,rhomin,rhomax
      real ::c1a,tc,hf
      real ::thetag,ss,rke,sigma
      real :: tambient
      real :: relativehumidity=0.0 !in % of partial pressure of water
                  !vapor/equilibrium pressure
      real :: specifichumidity !(computed from relativehumidity)
      ! mass of water vapor to air parcel's total (including water vapor)
      ! FP09/2019 : in his preliminary implementation, we assume that
      ! specific humidity is constant over the vertical, given the RH
      ! near the ground
      real :: pressground,zgroundref 
      real ::rsourceht
      real ::ep
      real ::rhomicrovalue
      real ::tfire,cpwood,cpwater,hwevap,gammav
      real ::twvap,tcrit,tstep
      real ::cvvapor,cvoxygen
      real :: waterProductToWoodMassRatio
      real ::thermcondair=33.8e-3 
      real ::ang=55. ,fcr0=1.4584e-4, fcor2,fcor3
      real :: r13

      end module constants

c2345678***************************************************
      module restarta
      Implicit None
c
c  setup read write data files
c
      character(len=257)::restartfile
c  topography file name
      character(len=60)::topofile
c
c  number of timesteps between writes to plot files (frpplot) and
c  restart files (frqrestart), number of writes to go to a plot
c  or restart file before opening a new one (nplotwrites and
c  nrestartwrites), and base file names for plot files (outname)
c  and restart files (restartname)
      character(len=257)::outname='comp.out'        !wss
      character(len=257)::outnamesub='compsub.out'  !FIXME MJH 1/26/18
c increased file length
      ! name of subdomain file for extraction on a subdomain
      end module restarta

c2345678***************************************************
      module nonlocal
      Implicit None
      real ,allocatable :: fg(:,:,:)
      real ,allocatable :: fhc(:,:,:)
      real ,allocatable :: fhcb(:,:,:)
      real ,allocatable :: psig(:,:,:)
      real :: rnhc=0.68      !these constants are later initialized in subroutine con()
      real :: rng=0.404
      real :: rnonl=0.596
      real :: hfgas=8789245.0
      real :: hfsolid=-200000.0
      real :: cg=20.0
      real :: gamma,gammo
      integer :: ifuelcount=0
      ! this arrays are not array anymore in normal version
      real,allocatable::  frho(:,:,:)
      real,allocatable::  fox(:,:,:)
      real,allocatable::frhof(:,:,:,:)
      real,allocatable:: frhosies(:,:,:,:),
     .   frhowater(:,:,:,:),
     .   frhovapor(:,:,:),fi(:,:,:,:),
     .       psiw(:,:,:,:),psif(:,:,:,:)
      real,allocatable:: cf(:,:,:,:),thetasolid(:,:,:,:)
      real,allocatable:: ff(:,:,:,:),fw(:,:,:,:),
     .      qflux(:, :,:)
      end module nonlocal
c2345678***************************************************
      module filesub
      !JMC this module contains a subroutine that will append
      !the time stamp to the end of a file name.
      implicit none
      contains
!----------------------------------------------------------
        subroutine namefile(itn,fnamein,fnameout)
        implicit none
        integer  itn,indx
        character(len=257):: fnamein,fnameout !FIXME MJH 1/26/18
        character(len=10):: end

        write (end,1001) itn
 1001   format(i0)

!        if (itn.lt.10) then
!          write (end,1001) itn
! 1001     format(i1.1)
!        else if (itn.lt.100) then
!          write (end,1101) itn
! 1101     format(i2.2)
!        else if (itn.lt.1000) then
!          write (end,1201) itn
! 1201     format(i3.3)
!        else if (itn.lt.10000) then
!          write (end,1301) itn
! 1301     format(i4.4)
!        else if (itn.lt.100000) then
!          write (end,1401) itn
! 1401     format(i5.5)
!        else
!          write (end,1501) itn
! 1501     format(i6.6)
!        endif
        indx=index(fnamein,' ')
        fnameout=fnamein(1:indx-1)//'.'//end
        end subroutine namefile
!----------------------------------------------------------
      end module filesub


!*****************************  MODULE fueldrag   ********************************!
      module fueldrag

      use fireteca, only:nfuel,min_rhof
      implicit none

      save

      integer, allocatable :: fuelinds(:)
      real, allocatable :: zonehts(:),zonedzs(:),zonerhos(:),iftwght(:)
      real, allocatable :: zoneus(:),zonevs(:),deltaus(:),deltavs(:)
      real :: ftop,uftop,vftop,zonerhotot
      real :: min_dz = 1e-10

      contains
!********************************* subroutine init_fueldrag ********************!
      subroutine init_fueldrag()

      Implicit None


      allocate(zonehts(nfuel))
      allocate(zonedzs(nfuel))
      allocate(zonerhos(nfuel))
      allocate(iftwght(nfuel))
      allocate(fuelinds(nfuel))
      allocate(zoneus(nfuel))
      allocate(zonevs(nfuel))
      allocate(deltaus(nfuel))
      allocate(deltavs(nfuel))

      end subroutine init_fueldrag
!********************************* end subroutine init_fueldrag ********************!
!********************************* subroutine calcZoneHts ********************!
      subroutine calcZoneHts(zbot,ztop,frhodg,fdepthdg)

      use fireteca, only:nfuel,min_rhof
      Implicit None

      integer :: i,j
      real :: ht,lastht
      real, intent(in) :: frhodg(:),fdepthdg(:),zbot,ztop

      fuelinds = 0
      zonedzs=1E-8
      zonehts=1E-8     
      lastht = zbot
      do j=1,nfuel
        ht = ztop
        do i=1,nfuel
          if(frhodg(i).gt.min_rhof/nfuel)then
            if((fdepthdg(i).le.ht).and.(fdepthdg(i).ge.lastht)
     &        .and.(count(fuelinds.eq.i).eq.0))then
              fuelinds(j) = i ! KOO081814
              ht = fdepthdg(i)
            endif !fdepthdg...
          endif !frhodg
        enddo !i
        lastht = ht
      enddo !j

      do j=1,nfuel
        if(fuelinds(j).ne.0)then
          iftwght(j) = (ztop-zbot)/(fdepthdg(fuelinds(j))+1E-8)
          zonehts(j) = fdepthdg(fuelinds(j));
          if(j.eq.1)then
            zonedzs(j) = zonehts(j)
          else
            zonedzs(j) = (zonehts(j)-zonehts(j-1));
          endif
          ftop = zonehts(j)
        endif
      enddo !j

      end subroutine calcZoneHts
!********************************* end subroutine calcZoneHts ********************!
!********************************* subroutine calcZoneRhos ********************!
      subroutine calcZoneRhos(zbot,ztop,frhodg,fdepthdg)

      Implicit None

      integer :: i,j
      real, intent(in) :: frhodg(:),fdepthdg(:),zbot,ztop

      zonerhos = 0.0
      do j=1,nfuel
       if(zonedzs(j).gt.min_dz)then
       do i=j,nfuel
        if(fuelinds(i).ne.0)then
         zonerhos(j) = zonerhos(j)+frhodg(fuelinds(i))*iftwght(fuelinds(i))
        endif
       enddo !i
       endif
      enddo !j

      zonerhotot = 0.0
      do j=1,nfuel
       if(zonedzs(j).gt.min_dz)then
        zonerhotot = zonerhotot+zonerhos(j)
       endif
      enddo
      end subroutine calcZoneRhos
!********************************* end subroutine calcZoneRhos ********************!
!********************************* subroutine calcZoneVels ********************!
      subroutine calcZoneVels(zbot,ztop,uftop,vftop)

      Implicit None

      integer :: i,j
      real, intent(in) :: zbot,ztop,uftop,vftop
      do j=nfuel,1,-1
       if(zonerhos(j).gt.min_rhof/nfuel)then
        if(zonedzs(j).gt.min_dz)then
         deltaus(j) = (uftop/ftop)*(zonerhos(j)/zonedzs(j))* 
     &                                (ftop/zonerhotot)*zonedzs(j)
         deltavs(j) = (vftop/ftop)*(zonerhos(j)/zonedzs(j))* 
     &                                (ftop/zonerhotot)*zonedzs(j)
        else !this fueltype and the previous have identical fueldepth
             !so use previous fueltype's zonedzs, and add previous
             !fueltype's zonerho to recalc both deltaus
         deltaus(j) = (uftop/ftop)*((zonerhos(j-1)+zonerhos(j))/zonedzs(j-1))* 
     &                                (ftop/zonerhotot)*zonedzs(j-1)
         deltavs(j) = (vftop/ftop)*((zonerhos(j-1)+zonerhos(j))/zonedzs(j-1))* 
     &                                (ftop/zonerhotot)*zonedzs(j-1)
         deltaus(j-1)=deltaus(j)
         deltavs(j-1)=deltavs(j)
        endif
       else
        deltaus(j) = 0.0
        deltavs(j) = 0.0
       endif
      enddo !j  

      do j=1,nfuel
       if(zonerhos(j).gt.min_rhof/nfuel)then
        if(j.eq.1)then
         zoneus(j) = deltaus(j)/2.0
         zonevs(j) = deltavs(j)/2.0
        else
         if(zonedzs(j).gt.1E-8)then
          zoneus(j) = zoneus(j-1)+deltaus(j-1)/2.0 + deltaus(j)/2.0
          zonevs(j) = zonevs(j-1)+deltavs(j-1)/2.0 + deltavs(j)/2.0
         else !this fueltype and the previous have identical fueldepth
             !so use previous fueltype's zoneus
          zoneus(j) = zoneus(j-1)
          zonevs(j) = zonevs(j-1)
         endif !zonedzs.gt.0.0
        endif !j=1 else ...
       else
        zoneus(j) = 0.0
        zonevs(j) = 0.0
       endif !zonerhos.gt.min_rhof
      enddo !j

      end subroutine calcZoneVels
!********************************* end subroutine calcZoneVels ********************!
!*********************************  subroutine calcXvfuels ********************!
      subroutine calcXvfuels(xvfu,xvfv,tkew,frhodg,fdepthdg)

      Implicit None

      integer :: i,j
      real, intent(in) :: frhodg(:),fdepthdg(:)
      real :: xvfu(:),xvfv(:),tkew(:)
      real :: velsq_sum

      !write(6,*) "The size of the ", &
      !           "argument vector xvfu, in calcXvfuels is:",size(xvfu)
      !write(6,*)
      xvfu = 0.0
      xvfv = 0.0
      tkew = 0.0 
      velsq_sum = 0.0 
      !find a total energy for this cell
      do j=1,nfuel
       velsq_sum = velsq_sum+(zoneus(j)**2+zonevs(j)**2)*zonedzs(j)
       if(zonehts(j).EQ.0.0) zonehts(j)=1e-8
      enddo !j
       
       velsq_sum=velsq_sum/(ftop+1e-8) 
       if(velsq_sum.eq.0) velsq_sum=1e-8
!       print*,'velsq_sum',velsq_sum,ftop    
      !obtain the hieght weighted velocity components for each fueltype
      do j=1,nfuel
!       print*,j,'zonehts',zonehts(j)
       if(zonehts(j).gt.0.0)then
        do i=1,nfuel
         if((fdepthdg(i).ge.zonehts(j)).and.(frhodg(i).gt.min_rhof/nfuel))then
          !xvfu(i) = xvfu(i)+zoneus(j)*zonedzs(j)/zonehts(j)   !JAS 7/30/07
          !xvfv(i) = xvfv(i)+zonevs(j)*zonedzs(j)/zonehts(j)   !just testing
          xvfu(i) = xvfu(i)+zoneus(j)*zonedzs(j)/fdepthdg(i)
          xvfv(i) = xvfv(i)+zonevs(j)*zonedzs(j)/fdepthdg(i)
        tkew(i) = tkew(i)+((zoneus(j)**2+zonevs(j)**2)/velsq_sum)*
      !&                              zonedzs(j)/zonehts(j) 
     &                              zonedzs(j)/fdepthdg(i) 
         endif
        enddo !i
       endif !(zonehts.gt.0.0)
      enddo !j

      end subroutine calcXvfuels
!********************************* end subroutine calcXvfuels ********************!

!*********************************  subroutine cleanup_fueldrag ********************!
      subroutine cleanup_fueldrag()

      Implicit None

      deallocate(zonehts)
      deallocate(zonedzs)
      deallocate(zonerhos)
      deallocate(iftwght)
      deallocate(fuelinds)
      deallocate(zoneus)
      deallocate(zonevs)
      deallocate(deltaus)
      deallocate(deltavs)

      end subroutine cleanup_fueldrag
!********************************* end subroutine cleanup_fueldrag ********************!

      end module fueldrag
!********************************* end module fueldrag ********************************!


