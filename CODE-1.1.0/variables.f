c2345678***************************************************
      module gridsetup
      Implicit None

      integer::j3=1,irst,nrst,nt,nts,ibctopbot
c              !j3-1 is 3-d and 0 is 2-d
c                  !irst is
c                      !nrst is
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
      integer::nv,nv2d
c             !number of variables in the xvb array
c                !number of 
      real::dx,dy,dz,dts,time=0, restarttime
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
      integer::idrag,isa,iinra,ilocation,ist=0
      !drag
      ! choice of sa
      ! inra flag
      !location flag
      !ist : choice of turbulence diffusion for o2 and theta (based on sb,kb or sa, sqrtk)
      integer::ixevariation
      ! tell is xe is variating with time
      integer::itheta,iperturb
              ! selection of potential temperature profile
                ! selection of perturbation of initial state (theta, pinwheel...)
      integer::icorio,ilspgf,frqlspgf,islip
c                   !flag for calculation of coriolis (1 from top, 2 from xe)
c		  ! flag for large scale pressure gradient force (1 fixed, 2 ajusted)
c                   !frq for updating large scale pressure gradient force (ilspgf=2)
c                          !islip =1 is no slip bcs at z=0
      integer :: itrestart,ittot
      integer::ipotflow ! call for potflow
      real::slopeangle  ! slope in degree for rotated gravity
      real::slopeazimuth  ! slope azimuth in degree ifor rotated gravity (0: the slope is onx, 90: the slope is on y)
      integer::uswitch,vswitch
      real :: uramp,vramp,uramptime,vramptime
      real::u0,v0,st
      real::zu,uprofilezu ! height of reference wind speed for uswich=2  
      integer::ius=0! extension of reference zone in the domain for uswitch=2 (default is 1)
      integer::iue=0! extension of reference zone in the domain for uswitch=2 (default is n)
      integer::jus=0 !! extension of reference zone in the domain for uswitch=2 (default is 1
       integer::jue=0 ! extension of reference zone in the domain for uswitch=2 (default  is m)
       integer::windspeedupfactor=1 ! ratio between wind time step versus fire time step
      integer::iord,isor,nonos,idiv,nfct,nonosold
      integer::ifuel,ifuelinra,idirt,ivegread,kmax,kf77max
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

       ! personnal flags
       integer:: irrl               !rod's personal flag
       integer ::ijww               !judy's personal flag
       integer ::ijmc                ! jesse's personal flag
       integer::ijas               !jeremy's personal flag
       integer::ieko              !eunmo's personal flag
       integer::ifp                 ! francois 's personal flag
       integer::ijld                 ! jean-luc 's personal flag
       integer::iekl             ! etienne's personal flag

 
        ! BELOW THIS LINE:
        ! personnal variable parameter set
        integer::fuelnumber 
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
      
c2345678***************************************************
c module bc contains boundary condition arrays for turbulence
      module bc
      Implicit None
      real, allocatable::
     .tauw(:,:),uzs(:,:),vzs(:,:)
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
        real::aa1,aa2,aa3,zdata(100),zcrdata(100),zcoeff(100)
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
      real,allocatable:: rhof(:,:,:)
      real,allocatable:: sizescale(:,:,:) ! ss as an array FP
      real,allocatable:: actualfueldepth(:,:,:)
      real,allocatable:: rhofinitial(:,:,:)
      real,allocatable::  cd(:, :,:)
      real,allocatable::  sa(:,:,:),saxy(:, :,:),saz(:, :,:)
      real,allocatable::  sb(:, :,:),
     .   rhomicro(:, :,:)
      real::     sc
      integer::kfuelmax ! max height of the fuel
      integer:: isoturb
      integer:: ilapdo
      real::diffcst ! diffusion cst (=0.09)
      real::rturbprandtl ! inverse of the turbulent prandtl number (usually between 1 and 3)
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
      real,allocatable:: convht(:,:,:)   !wss
      real,allocatable:: tambientarray(:,:,:)
      real,allocatable:: temps(:, :,:),tempg(:,:,:)
      real,allocatable::  foxb(:, :,:)
      real,allocatable::  firad(:, :,:)
      real,allocatable::  eastFlux(:, :,:) ! radiative flux to east cell face
      integer,allocatable:: ifirestart(:, :,:)
      real,allocatable:: frhosiesrad(:,:,:)
      real,allocatable:: rmoist(:, :,:),
     .       rhos(:,:,:),rhodirt(:,:,:),
     .   rhowater(:, :,:),
     .    cpsolid(:, :,:),
     .       sies(:, :,:),
     .       psiwmax(:, :,:),
     .  frhovaporb(:,:,:)
      real :: min_rhof 
      real::xhsmin,xhsmax,yhsmin,yhsmax !iheatsource extension (in m)
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
      real,allocatable:: xvfuel(:,:,:,:)
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
      real,allocatable:: sintheta(:,:,:) ! evaluation of v/sqrt(u2+v2)
      real:: xe1tot ! total initial amount of mass flux through west boundary
      real:: xe2tot ! total initial amount of mass flux through north boundary
      real,allocatable::xe1j(:),xe2i(:)
      ! based on ekman theory to impose a pressure gradient (ilspgf)
      end module pres

c2345678***************************************************
      module xvbin    !JMC  10/21/6  This module is for reading 
                      !windfields
      Implicit None
      integer :: itabsold,itabsnew
      real,allocatable:: xvbdataold(:,:,:,:)
      real,allocatable:: xvbdatanew(:,:,:,:)
c      real,allocatable:: prold(:,:,:)
c      real,allocatable:: prnew(:,:,:)
      character(len=60)::fxvbdataname=''        !JMC
      character(len=60)::xvbdataname  !='xvbdata/xvbdata'        !JMC
      character(len=60)::windfieldstartfile
   
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
      real,allocatable:: uab(:,:,:),vab(:,:,:),
     .                   oab(:,:,:)
      real,allocatable:: ua(:,:,:),va(:,:,:),
     .                   oa(:,:,:)
      end module advo

c2345678***************************************************
      module weights
      Implicit None
      real,allocatable:: wt(:),wtf(:)
      end module weights

c2345678***************************************************
      module constants
      Implicit None
      real ::cfhydro,cfchar,rhohydrothresh
      real ::rv,t00,ee0,hlat,g,rg
      real ::prndt,bv,cp,cv,cap,pi,prrcp
      real ::rnfuel,rhoref,sx,rho,rno
      real ::tramp,tfstep,aramp,bramp
      real ::psijoint,tjoint,tjoint2
      real ::fueldepth
      real ::rkmin,rkmax,rhomin,rhomax
      real ::c1a,tc,hf
      real ::thetag,ss,rke,sigma
      real :: tambient, pressground,zgroundref 
      real ::rsourceht,gamma,gammo
      real ::ep
      real ::rhomicrovalue
      real ::tfire,cpwood,cpwater,hwevap,gammav,cpdirt
      real ::twvap,tcrit,tstep
      real ::cvvapor,cvoxygen
      real ::thermcondair=33.8e-3 
      real ::ang=40. ,fcr0=1.4584e-4, fcor2,fcor3
      real :: r13

      end module constants

c2345678***************************************************
      module restarta
      Implicit None
c
c  setup read write data files
c
      character(len=60)::restartfile
c  topography file name
      character(len=60)::topofile
c
c  number of timesteps between writes to plot files (frpplot) and
c  restart files (frqrestart), number of writes to go to a plot
c  or restart file before opening a new one (nplotwrites and
c  nrestartwrites), and base file names for plot files (outname)
c  and restart files (restartname)
      character(len=60)::outname='comp.out'        !wss
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
      integer :: ifuelcount=0
      ! this arrays are not array anymore in normal version
      real,allocatable::  frho(:, :,:)
      real,allocatable::  fox(:, :,:)
      real,allocatable::frhof(:, :,:)
      real,allocatable:: frhosies(:,:,:),
     .   frhowater(:,:,:),
     .   frhovapor(:,:,:),fi(:,:,:),
     .       psiw(:, :,:),psif(:,:,:)
      real,allocatable:: cf(:,:,:),thetasolid(:,:,:)
      real,allocatable:: ff(:, :,:),fw(:, :,:),
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
        character(len=60):: fnamein,fnameout
        character(len=10):: end
        if (itn.lt.10) then
          write (end,1001) itn
 1001     format(i1.1)
        else if (itn.lt.100) then
          write (end,1101) itn
 1101     format(i2.2)
        else if (itn.lt.1000) then
          write (end,1201) itn
 1201     format(i3.3)
        else if (itn.lt.10000) then
          write (end,1301) itn
 1301     format(i4.4)
        else if (itn.lt.100000) then
          write (end,1401) itn
 1401     format(i5.5)
        else
          write (end,1501) itn
 1501     format(i6.6)
        endif
        indx=index(fnamein,' ')
        fnameout=fnamein(1:indx-1)//'.'//end
        end subroutine namefile
!----------------------------------------------------------
      end module filesub
