      module treedatinfo

c  max # trees allowed in the sampled/measured (excel) data
        parameter (ntr=400000)
c
c  type/source of data
c    itdsource=0 -> Forest Service Rocky Mountain Research Station data(20 X 50)
c    itdsource=1 -> LANL XL06 data 
c    itdsource=2 -> Phil Dennisons Calabasas fire data 
c    itdsource=3 -> General California chaparral (100 X 100)
c    itdsource=4 -> Las Vegas/Gallinas - douglas fir (100 X 100)
c    itdsource=5 -> Las Vegas/Gallinas - aspen (100 X 100)
c    itdsource=6 -> Angel Fire tree data from Sandia (25 X 10)
c    itdsource=7 -> Canada tree data ( 100m X  100m)
        parameter (itdsource=7)
         real, allocatable:: xloc(:) 
        real, allocatable::yloc(:) 
        real, allocatable::dbh(:)
        real, allocatable:: height(:)
        real, allocatable:: htlc(:)
        real, allocatable:: cr1(:)
        integer, allocatable:: treetype(:)
        integer, allocatable:: health(:)

                                                        
c  arrays and variables for reading sampled tree data
c        real, allocatable ::  xloc(:),yloc(:),dbh(:),height(:),
c     +       htlc(:),cr1(:),cr2(:)
c        integer, allocatable :: treetype(:)

c        character dummy*1 ! name of tree data file

c
      end module treedatinfo

      module higradinfo
c
c  x,y,z size for HIGRAD grid and x,y,z grid resolution (meters)
c  NOTE  ******  DO NOT CHANGE NXH, NYH, NZH HERE WITHOUT CHANGING
c                N, M, L IN SUBROUTINES METRYC AND ZCART
        parameter (nxh=750,nyh=1000,nzh=61,dxh=2.0,dyh=2.0,dzh=15.0)
c
c  x,y grid resolution for small resolution cells used to resolve
c  a tree within a HIGRAD cell
        parameter (dxs=.2,dys=.2)
c
c  max # trees allowed in the HIGRAD grid (once fully populated)
        parameter (ntrh=400000)
c
c  area of HIGRAD grid (min -> max) to populate with trees, so can
c  have segments of the grid with no trees
        parameter (itminxh=1,itmaxxh=750,itminyh=1,itmaxyh=1000) 
c
c  area of HIGRAD grid (min -> max) to be treated as a grass area
c  (this can overlap with tree areas, or be just grass) - for no
c  grass, set igminxh=0
        parameter (igminxh=1,igmaxxh=750,igminyh=1,igmaxyh=1000) 
c
c  area of HIGRAD grid (min -> max) to have the possibility of litter
c  (this can overlap with tree or grass) 
        parameter (ilminxh=1,ilmaxxh=750,ilminyh=1,ilmaxyh=1000) 
c
c  consider northern/southern exposure effects (more trees on northern
c  facing slopes, fewer on southern) - 0=>dont consider, 1=>do consider  
c  (NOTE: this will be ignored if topography is not being used)
c        parameter (treeslopefactor=0)
c
c  if considering north/south effect, percentage more trees on north face
c  and less on south
c        parameter (slopeeffect=20)
c
c  if considering north/south effect, what is the orientation of the HIGRAD
c  grid - 0=>north is in the U direction; 1=>north is in the V direction
        parameter (gridorientation=1)
c
c  arrays and variables used in populating the higrad grid with trees
        integer itreeh(ntrh),itreehtmp(ntrh),inewtr
c        integer itreeh2(ntrh)
c        real xloch(ntrh),yloch(ntrh)



c
c  variable for topography file name and array for topography data
c        character (len=40) :: topofile      ! name of topography file
c        real zs(nxh,nyh)                    ! topography data
c
      end module higradinfo

      module thininfo
c
        use higradinfo
c
c  thinning method
c    ithinmethod=0 => no thinning
c    ithinmethod=1 => thin trees based on tree diameter
c    ithinmethod=2 => thin trees to leave uneven sized round-ish patches
c    ithinmethod=3 => thin based on diameter and to leave uneven patches
c    ithinmethod=4 => thin to remove uneven sized eliptical shaped patches
c    ithinmethod=5 => thin to leave uneven sized eliptical shaped patches
c    ithinmethod=7 => thin according to Angelfire data functions
c    ithinmethod=8 => remove trees in linear strips
        parameter (ithinmethod=0)
c
c  area of HIGRAD grid (min -> max) to thin (the area must have trees 
c  populated on it) - ignored if no thinning is being done or a mask that
c  indicates areas to be thinned is being used
        parameter(ithinminxh=1,ithinmaxxh= 750,ithinminyh=1,
     +                                            ithinmaxyh= 1000) 
c
c  use a mask file (an unformatted 2D file of integers, size (nxh,nyh),
c  where 0 => dont thin in this HIGRAD cell, 1 => do thin in this cell) 
c  to indicate areas to be thinned - ithinmask=0 => no mask file, 
c  ithinmask=1 => use mask file (will be prompted for name)
c        parameter(ithinmask=0)
c
c  use a mask file (same as above) to determine a select area to compute
c  canopy bulk density for - ithinmaskbd=0 => no mask file,
c  ithinmaskbd=1 => use mask file (will be prompted for name)
c        parameter(ithinmaskbd=0)
c
c  for thinning based on tree diameter, diameter under which trees 
c  will be removed (dbh less than this value)
c        parameter (dbhthin=28.0)
c
c  for patch thinning, use circles of radius "patchradius" meters - for
c  eliptical patches, make them "patchradius" meters in x and "patchradiusy"
c  meters in y (NOTE:  For ithinmethod=2,3,or 5, you are specifying the size
c  of the patches to be left; for ithinmethod=4, you are specifying the size
c  of the patches to be removed.)
c        parameter (patchradius=35.0,patchradiusy=12.0)
c       parameter (patchradius=12.0,patchradiusy=3.0)
c
c  for patch thinning, estimate 1 patch every patchperiod meters
c        parameter (patchperiod=50.0)
c
c  for patch thinning, use computer generated seed for randomly placing 
c  patches or specify the seed, which will allow you to duplicate a 
c  specific random placing - ipatchseed=0 => use computer generated 
c  seed; ipatchseed=1 => use integer seed values specified in array
c  ipseed 
c        parameter (ipatchseed=0)
c        integer, dimension(4) :: ipseed = (/2,4,6,8/) c compiler error - needs 12 not 4
c          integer,dimension(12)::ipseed=(/2,4,6,8,10,12,14,16
c     &,18,20,22,24/)

c
      end module thininfo
 
      module constants
  
c
c  fraction of the tree canopy shape that is concave downward - this
c  really equates to how high are the branches off of the ground (in
c  general, use .8 for ponderosa and Douglas fir, .5 for pinon/juniper, 
c  .99 for coast live oak and California chaparral, .5 for aspen)
        parameter (clfactor=.8)
c
c  average tree canopy density (fine fuel) - this is the value that
c  the Forest Service typically gives us - compute maximum density
c  from the average - ????????  should this vary for healthy/unhealthy ???
c  (in general, use .4 for ponderosa, .7 for pinon/juniper, 1.3 for
c  California chaparral, .5 for Douglas fir, .2 for aspen)

         parameter (rhoavg=0.449,rhomax=(3.0/2.0)*rhoavg) ! this is conifer - based on white spruce ! 0.4 for fire run 
         parameter (rhodeadavg=0.449,rhodeadmax=(3.0/2.0)*rhoavg) ! this is dead conifer - based on 80% foliage loss in dead trees - 60% of mass is foliage . (0.449*0.6)*0.2 + (0.449*0.4)
         parameter(rhoavgpotr=0.104,rhomaxpotr=3.0/2.0*rhoavgpotr) ! Aspen
         parameter (rhoavgburn=0.04,rhomaxburn=(3.0/2.0)*rhoavgburn) 
c         parameter (rhoavg=0.5,rhomax=(3.0/2.0)*rhoavg) ! Judy's original
c   was 0.5 GM Nov 27


c  small depth of fuel upwind of the ignition
c        parameter (burnedsfcdepth=0.001)

c  actual open grassland grass density (within the fuelbed) (kg/m3)
        parameter (actualrhograss=2.13) ! 0.32/0.15=2.13 
c! 3.5 tons per ha == 0.32 kg/m2 /0.5 m ==> 0.635 kg/m3 !was 0.635 - using grass as the only surface fuel for wind
c     grass was 0.635 kg/m3 and 0.5m high originally
c  actual open grassland grass height (meters)
        parameter (actualgrassht=0.15)! was 0.5 - using grass as the only surface fuel for the wind run
c was 0.7 - changed Nov 27 GM to have 0.25 kg/m^2          

c
c  estimated maximum litter density (kg/m^3) (within the fuelbed)
        parameter (actualrholitter=11.0) ! 0.33 kg/m2 / 0.03 m = 11 ! from ft mac
c        parameter (actualrholitter=22.3) ! 0.33 + 0.34 kg/m2 / 0.03 m = 22.3 ! from ft mac (0.33) ! we add 0.34 kg/m2 since 70% of trees are dead and they've lost 80% of their folliage

c
c  estimated maximum litter height (m)
        parameter (actuallitterht=0.03)

c  estimated maximum fwd density (kg/m^3) (within the fuelbed)
        parameter (actualrhofwd=1.67) ! 0.05 kg/m2 / 0.03 m = 1.67 ! from ft mac
c
c  estimated maximum fwd height (m)
        parameter (actualfwdht=0.03)

c
c  nominal fuel moisture for canopy and grass (1.0 -> 100%)
       parameter (FFMC=96) !from ft mac wildfire-mnp-report.pdf
       parameter (fmoistcan=0.852,fmoistpotr= 2.) ! based on white spruce 0.852 ! guess potr=200
       parameter (fmoistlitter=((147.2*(101-FFMC)/(59.5+FFMC))/100.0)
     +  )
       parameter (fmoistseed=fmoistcan,fmoistmulch=fmoistlitter,
     +   fmoistfwd=fmoistlitter,fmoistgrass=fmoistlitter,
     +     fmoistdeadcan=fmoistlitter) ! These will probably need to be adjusted


c  nominal fuel size scale (solids) - meters
         parameter (fsizescaleconifer=0.0006,fsizescaleburned=0.0016, 
     +        fsizescalepotr=0.001,  ! based on 1 cm branch 
     +        fsizescalelitter=0.00045, !calculated from Ft Mcmurray using litter, sc1 fwd and sc2 fwd
     +        fsizescalegrass=0.0005) ! from FIRETEC default
c         parameter(fsizescale=0.0004)
c  was 0.0004 GM Nov 27
c  data for headers in output files that are generated (treesfueldata.dat) -
c  convertto = 'big_endian' or 'little_endian', ireal = 4 or 8 (real numbers
c  single or double precision)
        character (len=50) :: convertto='little_endian'
        parameter (ireal=4)

c
      end module constants

      module rhoinfo
        use higradinfo
        save

        type :: rho_contrib
          integer :: itdtreenum
          real :: rho_amount
          type (rho_contrib), pointer :: rnext
        end type rho_contrib

        type :: rho_pointer
          type (rho_contrib), pointer :: rptr
        end type rho_pointer

        type (rho_pointer), dimension (nxh,nyh,nzh) :: rhead
        type (rho_pointer), dimension (nxh,nyh,nzh) :: rtail
        type (rho_contrib), pointer :: rho_info
        type (rho_contrib), pointer :: rtempptr

        logical :: rhofirst =.true.

      end module rhoinfo

      program trees
c
c  program to populate a grid with trees from the sampled/measured data
c
c  Judy Winterkamp  9/02 
c
      use treedatinfo
      use higradinfo 
      use thininfo
      use constants
      use rhoinfo
c
c  arrays for HIGRAD data

c
c  array for average height of fuel within a cell and initialize
        real,allocatable,dimension(:,:,:) :: actualfueldepth
        real,allocatable:: ztopcellout(:,:,:)
        real,allocatable:: zbottomcellout(:,:,:) 
        real,allocatable:: rhof(:,:,:)         ! fuel density for each HIGRAD cell - conifer initially
       real,allocatable:: rhofdead(:,:,:)         ! fuel density for each HIGRAD cell - conifer initially
        real,allocatable:: rhofdecid(:,:,:)     ! fuel density of deciduous
        real,allocatable:: rhofgrass(:,:,:)    ! array for saving grass contributions
        real,allocatable:: rhoflitter(:,:,:)   ! array for saving litter contrib
        real,allocatable:: rhoffwd(:,:,:)
        real,allocatable:: sizescale(:,:,:)   ! fuel size scale for each HIGRAD cell
        real,allocatable:: moist(:,:,:)        ! fuel moisture or each HIGRAD cell
      real,allocatable :: x(:,:,:),y(:,:,:)
      real,allocatable :: zs(:,:)
c      real z(nzh),zedge(nzh+1)  ! for metryc
      real,allocatable :: z(:),zedge(:)
      real,allocatable :: ignline(:),temp(:)
      real,allocatable:: gi(:,:,:)
      real gmul(nzh)      ! for metryc
      real, allocatable ::  c13(:,:),c23(:,:)      ! for metryc
      character (len=50) :: dataname      ! data type name for output file  

        real,allocatable::xloch(:),yloch(:) 
C               xloch(ntrh),yloch(ntrh)
        integer ntreedat,ntreedat2 !,itreenum2(ntr)
        integer, allocatable::  itreenum(:)
        real xmax,ymax 


         integer,allocatable ::  genie(:,:)
c       real effrhomax
c  initialize
      allocate(xloc(ntr),yloc(ntr),dbh(ntr),height(ntr))
      allocate(cr1(ntr),treetype(ntr),health(ntr),itreenum(ntr))
      allocate(htlc(ntr),rhof(nxh,nyh,nzh))
      allocate(rhofdead(nxh,nyh,nzh))
      allocate(rhofgrass(nxh,nyh,nzh))
      allocate(rhoflitter(nxh,nyh,nzh))
      allocate(rhoffwd(nxh,nyh,nzh))
      allocate(rhofdecid(nxh,nyh,nzh))
      allocate(sizescale(nxh,nyh,nzh))
      allocate(moist(nxh,nyh,nzh))
      allocate(x(nxh,nyh,nzh),y(nxh,nyh,nzh))
      allocate(actualfueldepth(nxh,nyh,nzh))
      allocate(ztopcellout(nxh,nyh,nzh))
      allocate(zbottomcellout(nxh,nyh,nzh))
      allocate(gi(nxh,nyh,nzh))
      allocate(zs(nxh,nyh),genie(nxh,nyh))
      allocate(z(nzh),zedge(nzh+1))
      allocate(c13(nxh,nyh),c23(nxh,nyh))
      allocate(xloch(ntrh),yloch(ntrh))
      allocate(ignline(nxh))
      allocate(temp(2))


      rhof=0.0
      rhofdead=0.0
      rhofdecid=0.0
      rhofgrass=0.0
      rhoflitter=0.0
      rhoffwd=0.0
c      rhofmoist=0.00001
c      rhofseedling=0.0
c      rhofmoss=0.0
c       rhofmulch=0.0
c       mulch_fuel=0.0
c       mulch_depth=0.0
c      sizescale=0.0
        zs=0.0


c    The mask file was written in R so (1,1) corresponds to the top left corner NOT cartesian coordinates.
      
           open (3,file='Mixedwood_genie_swind.txt',
     +        form='formatted',status='old')
         ios=1
          iy=nyh
         do while (ios.ge.0.and.iy.ge.1)
        read (3,*,iostat=ios) genie(:,iy)
         iy=iy-1
         end do
           close (3)


c Genie	Vegetation
c 0	Roads
c 1	grassy corridor
c 2	river
c 3	shadow (conifer)
c 4	short grass
c 5	spruce
c 6	yellow aspen
c 7	green aspen
c 8	roof
c 9	roof/shoulder
c 10	border/transition - I called it grass - could be dead aspen?? grass under is reasonable.

c      fsizescaleconifer,fsizescaleburned,fsizescalepotr,fsizescalelitter,fsizescalegrass

      sizescale=fsizescalelitter

      do i=1,nxh
        do j=1,nyh
           do k=1,nzh ! conifer size scale below 13.25 m - aspen size scale above 
            if(k.le.10.and.k.gt.1)then  
              sizescale(i,j,k)=fsizescaleconifer 
            else if(k.gt.10)then 
              sizescale(i,j,k)=fsizescalepotr
            end if
            if(genie(i,j).eq.4.and.k.gt.1)then
              sizescale(i,j,k)=fsizescalegrass
            end if
           end do
         end do
      end do

      actualfueldepth=0.0
      ztopcellout=0.0
      zbottomcellout=0.0
      moist=0.00001
      call rhoinfo_init
c
c  call routine to read tree data

c      if (itdsource.eq.7) then
c        call read_canada
c      endif

      open (1,file='Mixedwood_stems_swind_domain_input.txt'
     +  ,form='formatted',status='old')
      ntreedat=1
      do while (ios.ge.0)
        read (1,*,iostat=ios) xloc(ntreedat),
     +       yloc(ntreedat),height(ntreedat),cr1(ntreedat),
     +       htlc(ntreedat),treetype(ntreedat),health(ntreedat)
c1100   format(i6,f5.2,f5.2,a3,a6,f5.1,f5.1,a8,a3,f6.1,f4.1,f5.1,
c    +         f5.1,f4.1,f4.1,a5,a1)
ccc gm        if (ntreedat.le.100)
ccc gm     &         height(ntreedat),cr1(ntreedat),htlc(ntreedat)
        if (ios.lt.0) then
          ntreedat=ntreedat-1     ! end-of-file
        else if (.not.(xloc(ntreedat).eq.0.and.
     +                 yloc(ntreedat).eq.0)) then
          ntreedat=ntreedat+1     ! not a blank line in excel file
        endif
      enddo
      close (1)
      

      htavg=0.0
      do i=1,ntreedat
        itreenum(i)=i
c        htlc(i)=0.1
        htavg=htavg+height(i)
ccc gm        if (height(i).gt.2.0) print *,'height ',i,height(i)
      enddo
       
      xmax=maxval(xloc(1:ntreedat))
      ymax=maxval(yloc(1:ntreedat))
       print*,'xmax =',xmax,'ymax =',ymax
      hmin=minval(height(1:ntreedat))
      hmax=maxval(height(1:ntreedat))
      htavg=htavg/float(ntreedat)

        
       print *,'End read_canada - ntreedat = ',ntreedat
ccc gm      print *,'height min,max,avg ',hmin,hmax,htavg
c


      print *,'After reading tree data  - ntreedat = ',ntreedat

c       print*,sizescale(:,100,1)
c
c  check parameter settings for errors
      if (itminxh.lt.1.or.itmaxxh.gt.nxh.or.itminxh.gt.itmaxxh.or.
     +    itminyh.lt.1.or.itmaxyh.gt.nyh.or.itminyh.gt.itmaxyh.or.
c    +    igminxh.lt.0.or.igmaxxh.gt.nxh.or.igminxh.gt.igmaxxh.or.
c    +    igminyh.lt.1.or.igmaxyh.gt.nyh.or.igminyh.gt.igmaxyh.or.
     +    dxh.le.dxs.or.dyh.le.dys.or.int(dxh/dxs)*dxs.ne.dxh.or.
     +    int(dyh/dys)*dys.ne.dyh) then
        print *,' Error in parameter settings - stop and fix'
        stop
      endif
c
c  read topography data
c      write (*,1200) 
c 1200 format(1x,'Topography file (2D binary) or flat: ',$)
c      read (*,*) topofile
c      if (topofile.eq.'flat') then


        zs=0.0

      ios=1
c      topofile=""
      open (1,file="Mixedwood_dem_swind.txt",form='formatted'
     +   ,status='old')
       iy=nyh
      do while (ios.ge.0.and.iy.ge.1)
        read (1,*,iostat=ios) zs(:,iy)
          iy=iy-1     ! not a blank line in excel file
      enddo
      close (1)
c      endif

      ios=1
      open (1,file="ignite.dat",form='formatted'
     +   ,status='old')
      do i=1,750
        ignline(i) = 0.0
      enddo
      do i=1,7
        read(1,*)
      enddo
      do while (ios.ge.0)
        read (1,*,iostat=ios) temp(:)
        ignline(temp(1)) = MAX(ignline(temp(1)),temp(2))
      enddo
      close (1)

      open (1,file="mcmurray_topo_swind_domain.bin",form='unformatted',
     + status='unknown')
        write (1) zs
      close (1)

              
c
c  populate the HIGRAD grid with trees or read old tree locations
c  and info from an old run of this program

        do i=1,ntreedat
          xloch(i)=xloc(i)
          yloch(i)=yloc(i)
          itreeh(i)=itreenum(i)
        enddo
        inewtr=ntreedat

c       inewtr=ntreedat

      print *,'After populating - inewtr = ',inewtr

       print*,'xmax=',xmax,'ymax=',ymax
c
c  generate file for plotting tree locations with xmgrace
      open (1,file='treesloc.dat',form='formatted',status='unknown')
      do i=1,inewtr
        write (1,*) xloch(i),',',yloch(i)
      enddo
      close (1)
c
c  initialize
      do iz=1,nzh
        do iy=1,nyh
          do ix=1,nxh
            x(ix,iy,iz)=(ix-1)*dxh-((nxh+1)/2-1)*dxh ! (ix-1)*2.-((160+1)/2-1)*2.
          enddo
        enddo
      enddo
       
      do iz=1,nzh
        do iy=1,nyh
          ycoord=(iy-1)*dyh-((nyh+1)/2-1)*dyh
          do ix=1,nxh
            y(ix,iy,iz)=ycoord
          enddo
        enddo
      enddo
      do iz=1,nzh
        z(iz)=(iz-1)*dzh+0.5*dzh
      enddo
      do iz=1,nzh+1
        zedge(iz)=(iz-1)*dzh
      enddo
      zb=zedge(nzh+1)
      call metryc(x,y,z,zs,gi,gmul,c13,c23,dxh,dyh,dzh,0.0,zb,1)

c      do iz=1,nzh
c       print*,zcart(zedge(iz),1,1,zs,zb)
c      end do

      treehtavg=0.0
      htlcavg=0.0
      navg=0
      rmcanopysum=0.0
      thinmcanopysum=0.0
c
c  for each tree...

      do it=1,inewtr
c
c  initialize tree characteristics - tree height (treeht), height live
c  crown (htlctr), height from base of live crown (hb), crown length (cl),
c  and crown radius (cr)  NOTE:  The tree numbers are searched to find a 
c  tree number match - it is not assumed that tree #1 is in the first 
c  array position, tree #2 in second position, etc.
c 
c  note:  in an email from John Bailey of the Forest Service he says -
c  ""crown radius" is the tree's longest dimension (the longest living 
c  branch), which is CR1, and 90 degrees clockwise to that dimension is 
c  CR2.  Crowns are rarely circular, typically elliptical....hence, two 
c  measurements are needed for volume."  For now, we average cr1 and cr2.
        if (itreeh(it).eq.itreenum(itreeh(it))) then
          itdindex=itreeh(it)
        else
          itdindex=1
          do while (itreeh(it).ne.itreenum(itdindex))
            itdindex=itdindex+1
            if (itdindex.gt.ntreedat) then
              print *,' Error locating tree number ',itreeh(it)
              stop
            endif
          enddo
        endif
        treeht=height(itdindex)
        htlctr=htlc(itdindex)
        hb=(treeht-htlctr)*(1.0-clfactor)
        cl=(treeht-htlctr)*clfactor
        cr=(cr1(itdindex))
        ispecies=treetype(itdindex)
c
c  sum tree height and htlc for averages of these quantities and save
c  max tree height
        treehtavg=treehtavg+treeht
        htlcavg=htlcavg+htlctr
        navg=navg+1 
        if (it.eq.1) then
          treehtmax=treeht
        else
          treehtmax=max(treehtmax,treeht)
        endif
c       
c  determine which HIGRAD cells the tree is in, and the minimum and
c  maximum x and y cells that contain part of the tree (note that
c  "corners" of this x/y area may not contain any part of the tree)
        ixcell=nint(xloch(it)/dxh)
        if (ixcell.le.0) ixcell=1
        iycell=nint(yloch(it)/dyh)
        if (iycell.le.0) iycell=1
        nminx=ixcell-int(cr/dxh+1)
        if (nminx.lt.1) nminx=1
        nmaxx=ixcell+int(cr/dxh+1)
        if (nmaxx.gt.nxh) nmaxx=nxh
        nminy=iycell-int(cr/dyh+1)
        if (nminy.lt.1) nminy=1
        nmaxy=iycell+int(cr/dyh+1)
        if (nmaxy.gt.nyh) nmaxy=nyh
        do iz=1,nzh
          ztopcell=zcart(zedge(iz+1),ixcell,iycell,zs,zb)-
     +                                      zs(ixcell,iycell)
          zbotcell=zcart(zedge(iz),ixcell,iycell,zs,zb)-
     +                                      zs(ixcell,iycell)
          if (htlctr.le.ztopcell.and.
     +                         htlctr.ge.zbotcell) nminz=iz
          if (treeht.le.ztopcell.and.
     +                         treeht.ge.zbotcell) nmaxz=iz
ccc gm          if (it.eq.1)
ccc gm    +        print *,' T0',zbotcell,ztopcell,htlctr,treeht,nmaxz,nminz
        enddo
c
c  save minimum and maximum vertical cells for later use
        if (it.eq.1) then
          icminzh=nminz
          icmaxzh=nmaxz
        else
          icminzh=min(nminz,icminzh)
          icmaxzh=max(nmaxz,icmaxzh)
        endif
ccc gm        if (it.le.5) print *,' T1',it,treeht,htlctr,hb,cl,cr,
ccc gm     +               ixcell,iycell,nminx,nmaxx,nminy,nmaxy,
ccc gm     +               nminz,nmaxz,icminzh,icmaxzh
c

c Here we need to determine the value of moisture - is it deciduous or conifer?



c  for all cells that contain part of the tree...
        do iz=nminz,nmaxz
          do iy=nminy,nmaxy
            do ix=nminx,nmaxx
            if(iy.le.ignline(ix))then

            else if(treetype(it).eq.1)then  ! treetype ==1 is POTR

              rmcanopy=0.0
              ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
              zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
c  divide the HIGRAD cell into smaller areas, and sum over all the smaller
c  areas to get the mass of the canopy in the HIGRAD cell 
              ixloop=nint(dxh/dxs)
              iyloop=nint(dyh/dys)
ccc gm              if (it.eq.1.and.ix.eq.nminx+2.and.iy.eq.nminy+2.and
ccc gm     +            .iz.eq.nminz+2) print *,' T2',ztopcell,zbotcell,
ccc gm     +             ixloop,iyloop
              do iyy=1,iyloop
                ylocs=iy*dyh+(real(iyy)-.5)*dys
                do ixx=1,ixloop
                  xlocs=ix*dxh+(real(ixx)-.5)*dxs
                  rprime2=(xlocs-xloch(it))**2+(ylocs-yloch(it))**2
                  ftop=((-1.0*cl)/(cr*cr))*rprime2+treeht
                  fbot=(hb/(cr*cr))*rprime2+htlctr
                  zmid=.5*(min(ftop,ztopcell)+max(fbot,zbotcell))
                  rhoscale=((zmid+(cl/(cr*cr)*rprime2)-htlctr)/
     +                      (treeht-htlctr))*rhomaxpotr
c  vcell = volume of canopy in small cell, rmcell= real mass of canopy
c  in small cell, rmcanopy = real mass of canopy in HIGRAD cell
                  vcell=dxs*dys*max(0.0,(min(ztopcell,ftop)-
     +                                   max(zbotcell,fbot)))
                  rmcell=vcell*rhoscale
                  rmcanopy=rmcanopy+rmcell
                enddo
              enddo
c
c  compute fuel density for this cell, and then add it to the density of
c  other trees in the HIGRAD cell, limiting it by rhomax
              rhoftree=rmcanopy/(dxh*dyh*(ztopcell-zbotcell)) !mass/volume
              rhofdecid(ix,iy,iz)=rhofdecid(ix,iy,iz)+rhoftree
c              if(genie(ix,iy).ne.12)then
              rhofdecid(ix,iy,iz)=rhofdecid(ix,iy,iz)+rhoftree
c              else if(genie(ix,iy).eq.12)then
c              rhofdecid(ix,iy,iz)=rhofdecid(ix,iy,iz)+(rhoftree*0.9) ! if the area is burned there is 90% deciduous fuel
c              end if

c   !!!!!! Let's come back to this... do we need to limit this value?? !!!!!!!!!!!!!!!!!!!
              if (rhofdecid(ix,iy,iz).gt.rhomaxpotr) then
                rhoftree=rhoftree-(rhofdecid(ix,iy,iz)-rhomaxpotr)
                rhofdecid(ix,iy,iz)=rhomaxpotr
              endif


          else if(treetype(it).eq.2.and.health(it).eq.1) then   ! treetype==2 is conifer and health==1 is alive


C             if(genie(ix,iy).eq.4)then  
C              effrhomax=rhomaxburn
C             else
C              effrhofmax=rhomax
C             endif


              rmcanopy=0.0
              ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
              zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
c  divide the HIGRAD cell into smaller areas, and sum over all the smaller
c  areas to get the mass of the canopy in the HIGRAD cell 
              ixloop=nint(dxh/dxs)
              iyloop=nint(dyh/dys)
ccc gm              if (it.eq.1.and.ix.eq.nminx+2.and.iy.eq.nminy+2.and
ccc gm     +            .iz.eq.nminz+2) print *,' T2',ztopcell,zbotcell,
ccc gm     +             ixloop,iyloop
              do iyy=1,iyloop
                ylocs=iy*dyh+(real(iyy)-.5)*dys
                do ixx=1,ixloop
                  xlocs=ix*dxh+(real(ixx)-.5)*dxs
                  rprime2=(xlocs-xloch(it))**2+(ylocs-yloch(it))**2
                  ftop=((-1.0*cl)/(cr*cr))*rprime2+treeht
                  fbot=(hb/(cr*cr))*rprime2+htlctr
                  zmid=.5*(min(ftop,ztopcell)+max(fbot,zbotcell))
                  rhoscale=((zmid+(cl/(cr*cr)*rprime2)-htlctr)/
     +                      (treeht-htlctr))*rhomax
c  vcell = volume of canopy in small cell, rmcell= real mass of canopy
c  in small cell, rmcanopy = real mass of canopy in HIGRAD cell
                  vcell=dxs*dys*max(0.0,(min(ztopcell,ftop)-
     +                                   max(zbotcell,fbot)))
                  rmcell=vcell*rhoscale
                  rmcanopy=rmcanopy+rmcell
                enddo
              enddo
c
c  compute fuel density for this cell, and then add it to the density of
c  other trees in the HIGRAD cell, limiting it by rhomax
              rhoftree=rmcanopy/(dxh*dyh*(ztopcell-zbotcell)) !mass/volume

c              if(genie(ix,iy).ne.12)then
              rhof(ix,iy,iz)=rhof(ix,iy,iz)+rhoftree
c              else if(genie(ix,iy).eq.12)then
c              rhof(ix,iy,iz)=rhof(ix,iy,iz)+(rhoftree*0.1) ! if the area is burned there is 10% conifer fuel
c              end if
              
c              if (genie(ix,iy).eq.3)then
                 if (rhof(ix,iy,iz).gt.rhomax) then
                   rhoftree=rhoftree-(rhof(ix,iy,iz)-rhomax)
                   rhof(ix,iy,iz)=rhomax
                 endif

          else if(treetype(it).eq.2.and.health(it).eq.2) then   ! treetype==2 is conifer and health==2 is dead


C             if(genie(ix,iy).eq.4)then  
C              effrhomax=rhomaxburn
C             else
C              effrhofmax=rhomax
C             endif


              rmcanopy=0.0
              ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
              zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
c  divide the HIGRAD cell into smaller areas, and sum over all the smaller
c  areas to get the mass of the canopy in the HIGRAD cell 
              ixloop=nint(dxh/dxs)
              iyloop=nint(dyh/dys)
ccc gm              if (it.eq.1.and.ix.eq.nminx+2.and.iy.eq.nminy+2.and
ccc gm     +            .iz.eq.nminz+2) print *,' T2',ztopcell,zbotcell,
ccc gm     +             ixloop,iyloop
              do iyy=1,iyloop
                ylocs=iy*dyh+(real(iyy)-.5)*dys
                do ixx=1,ixloop
                  xlocs=ix*dxh+(real(ixx)-.5)*dxs
                  rprime2=(xlocs-xloch(it))**2+(ylocs-yloch(it))**2
                  ftop=((-1.0*cl)/(cr*cr))*rprime2+treeht
                  fbot=(hb/(cr*cr))*rprime2+htlctr
                  zmid=.5*(min(ftop,ztopcell)+max(fbot,zbotcell))
                  rhoscale=((zmid+(cl/(cr*cr)*rprime2)-htlctr)/
     +                      (treeht-htlctr))*rhodeadmax
c  vcell = volume of canopy in small cell, rmcell= real mass of canopy
c  in small cell, rmcanopy = real mass of canopy in HIGRAD cell
                  vcell=dxs*dys*max(0.0,(min(ztopcell,ftop)-
     +                                   max(zbotcell,fbot)))
                  rmcell=vcell*rhoscale
                  rmcanopy=rmcanopy+rmcell
                enddo
              enddo
c
c  compute fuel density for this cell, and then add it to the density of
c  other trees in the HIGRAD cell, limiting it by rhomax
              rhoftree=rmcanopy/(dxh*dyh*(ztopcell-zbotcell)) !mass/volume

c              if(genie(ix,iy).ne.12)then
              rhofdead(ix,iy,iz)=rhofdead(ix,iy,iz)+rhoftree
c              else if(genie(ix,iy).eq.12)then
c              rhof(ix,iy,iz)=rhof(ix,iy,iz)+(rhoftree*0.1) ! if the area is burned there is 10% conifer fuel
c              end if
              

                if (rhofdead(ix,iy,iz).gt.rhodeadmax) then
                   rhoftree=rhoftree-(rhofdead(ix,iy,iz)-rhodeadmax)
                   rhofdead(ix,iy,iz)=rhodeadmax
                end if
      end if


c
c  save information about this trees contribution to the density for this
c  HIGRAD cell (note that rhoftree may be 0.0 if this cell does not really
c  contain part of the tree - as noted above, "corners" of the x/y tree
c  area may actually be empty)
              if (rhoftree.ne.0.0) 
     +          call rhoinfo_add(ix,iy,iz,itreenum(itdindex),rhoftree)
c
c  sum canopy mass for average later
              rmcanopysum=rmcanopysum+rmcanopy
          if(rmcanopy.ne.rmcanopy) print*,rmcanopy,ix,iy,iz,it
            enddo
          enddo
        enddo

       enddo    
      deallocate(xloc,yloc,dbh,height,xloch,yloch)
      deallocate(cr1,treetype,health,htlc)
 
c  compute canopy bulk density from the average tree height and htlc
      treehtavg=treehtavg/real(navg)
      htlcavg=htlcavg/real(navg)
      avgcanht=treehtavg-htlcavg
      volcanopy=avgcanht*((itmaxxh-itminxh)*dxh)*
     +                            ((itmaxyh-itminyh)*dyh)
      bdcanopy=rmcanopysum/volcanopy
c      print*,'rmcanopysum',rmcanopysum
c      print*,'volcanopy',volcanopy
      print *,'Canopy Bulk Density dead spruce = ',bdcanopy
      if (ithinmask.eq.1.or.ithinmaskbd.eq.1) then
        bdthincanopy=thinmcanopysum/(avgcanht*(nmaskcells*dxh*dyh))
        print *,'Masked Thinning Area Canopy Bulk Density = ',
     +                                                bdthincanopy
      endif
      print *,'Tree Height Max = ',treehtmax



c here we add grass, litter etc - they are mass weighted at the end


c     if (igminxh.eq.0) goto 20
        do iy=igminyh,igmaxyh
          do ix=igminxh,igmaxxh

c  izg locates the top cell of the grass at a particular x and y on the 
c  HIGRAD grid
             
              izg=1
              zgtopcell=zcart(zedge(izg+1),ix,iy,zs,zb)-zs(ix,iy)
              do while (zgtopcell.lt.actualgrassht)
                izg=izg+1
                zgtopcell=zcart(zedge(izg+1),ix,iy,zs,zb)-zs(ix,iy)
              enddo

c  loop from the ground up to the top of the grass and fill each cell with 
c  the appropriate density of grass depending on the amount of the cell that
c  is full and the shading
              do iz=1,izg
                actualfuelbottom= 0.     ! inserting 0 for ground
                fueldistributionslope=0. ! pos # -> load higher at top (kg/m3/m)


c                Genie raster
c# 0 roads/parkinglots
c# 1 grassy corridor
c# 2 river
c# 3 shadow on veg - CONIFER
c# 4 short grass/meadow... but could be forest in meadow
c# 5 dead spruce - could be meadow or dead spruce
c# 6 Yellow aspen
c# 7 green aspen
c# 8 roof
c# 9 roof/shoulder
c# 10 border/transition - could be dead tree/shaded veg/ - I called this grass... looks like could be dead aspen? also seen along
c# the edge of river/slump areas... grass is reasonable for these areas.                 

c ! putting litter and fwd in under trees only 

                if(genie(ix,iy).eq.3.or.genie(ix,iy).eq.5.or.
     +            genie(ix,iy).eq.6.or.genie(ix,iy).eq.7)then
                actualgroundloadlitter=actualrholitter ! fuel load @ ground if linear etrap
     +                       -.5*actuallitterht*fueldistributionslope
                ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
                zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
                rholitter=(min(actuallitterht,ztopcell)
     +                    -max(actualfuelbottom,zbotcell))/
     +                    (ztopcell-zbotcell)*
     +                    (actualgroundloadlitter +
     +                    .5*fueldistributionslope*
     +                    (min(actuallitterht,ztopcell)+
     +                    max(actualfuelbottom,zbotcell)))

                rhoflitter(ix,iy,iz)=rhoflitter(ix,iy,iz)+rholitter
     
                actualgroundloadfwd=actualrhofwd ! fuel load @ ground if linear etrap
     +                       -.5*actualfwdht*fueldistributionslope
                ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
                zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
                rhofwd=(min(actualfwdht,ztopcell)
     +                    -max(actualfuelbottom,zbotcell))/
     +                    (ztopcell-zbotcell)*
     +                    (actualgroundloadfwd +
     +                    .5*fueldistributionslope*
     +                    (min(actualfwdht,ztopcell)+
     +                    max(actualfuelbottom,zbotcell)))

                rhoffwd(ix,iy,iz)=rhoffwd(ix,iy,iz)+rhofwd

 


              end if
c ! put grass everywhere except non-fuel areas

              if(genie(ix,iy).eq.1.or.genie(ix,iy).eq.15.or.
     +          (genie(ix,iy).ge.3.and.
     +          genie(ix,iy).le.7))then  

                actualgroundloadgrass=actualrhograss ! fuel load @ ground if linear etrap
     +                   -.5*actualgrassht*fueldistributionslope 
                ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
                zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
                rhograss=(min(actualgrassht,ztopcell)-
     +                       max(actualfuelbottom,zbotcell))/
     +                       (ztopcell-zbotcell)*
     +                       (actualgroundloadgrass +
     +                       .5*fueldistributionslope*
     +                       (min(actualgrassht,ztopcell)+
     +                       max(actualfuelbottom,zbotcell)))

                rhofgrass(ix,iy,iz)=rhofgrass(ix,iy,iz)+rhograss

            end if
                ztopcellout(ix,iy,iz)=ztopcell
                zbottomcellout(ix,iy,iz)=zbottomcell

                if(rhoflitter(ix,iy,iz).ne.0.or.rhoffwd(ix,iy,iz)
     +            .ne.0.or.rhofgrass(ix,iy,iz).ne.0.or.rhof(ix,iy,iz)
     +            .ne.0.or.rhofdead(ix,iy,iz).ne.0)then
                      actualfueldepth(ix,iy,iz)=
     +                    ((min(actuallitterht,ztopcell)-zbotcell)
     +                         *rhoflitter(ix,iy,iz)
     +                    +(min(actualfwdht,ztopcell)-zbotcell)
     +                        *rhoffwd(ix,iy,iz)
     +                    +(min(actualgrassht,ztopcell)
     +                        -zbotcell)*rhofgrass(ix,iy,iz)
     +                    +rhofdead(ix,iy,iz)*(ztopcell-zbotcell)
     +                    +rhof(ix,iy,iz)*(ztopcell-zbotcell)
     +                    +rhofdecid(ix,iy,iz)*(ztopcell-zbotcell))
     +                    /(rhoflitter(ix,iy,iz)+rhoffwd(ix,iy,iz)
     +               +rhofgrass(ix,iy,iz)+
     +                rhof(ix,iy,iz)+rhofdead(ix,iy,iz)+
     +                rhofdecid(ix,iy,iz))
                rhof(ix,iy,iz)=rhof(ix,iy,iz)+rhoflitter(ix,iy,iz)+
     +          rhoffwd(ix,iy,iz)+rhofgrass(ix,iy,iz)
                end if

c              if(actualfueldepth(ix,iy,iz).le.0)print*,ix,iy,iz

              enddo

          enddo
        enddo
                        print*,'sfc fuels added'

 20   continue
 
 

      do iz=1,nzh
        do iy=1,nyh
          do ix=1,nxh
            if (rhof(ix,iy,iz).gt.0.0) then
              rhoftree=rhof(ix,iy,iz)-(rhofgrass(ix,iy,iz)+
     +                                 rhoflitter(ix,iy,iz)+
     +                                 rhoffwd(ix,iy,iz))
c          if(rhoftree.lt.0)print*,ix,iy,iz
              moist(ix,iy,iz)=((rhoftree*fmoistcan)+
     +                         (rhofdead(ix,iy,iz)*fmoistdeadcan)+
     +                       (rhofgrass(ix,iy,iz)*fmoistgrass)+
     +                       (rhoflitter(ix,iy,iz)*fmoistlitter)+
     +                       (rhoffwd(ix,iy,iz)*fmoistfwd)+
     +                       (rhofdecid(ix,iy,iz)*fmoistpotr))/
     +                       (rhoftree+rhofdead(ix,iy,iz)+
     +                          rhofgrass(ix,iy,iz)+
     +                          rhoflitter(ix,iy,iz)+ 
     +                          rhoffwd(ix,iy,iz) +rhofdecid(ix,iy,iz))

            endif
          enddo
        enddo
      enddo
      k=20

c      print*, fmoistcan,fmoistdeadcan,fmoistgrass,fmoistlitter,fmoistfwd
c     + ,fmoistpotr
      open(3,file='rhofgrass.dat',form='formatted',status='unknown')
       do J=1,nyh
      WRITE(3,"(1000(F8.5,','))")rhofgrass(:,J,1)
       end do
      close(3)
c      print*,rhofgrass(:,:,1)
      open(3,file='rhoflitter.dat',form='formatted',status='unknown')
       do J=1,nyh
      WRITE(3,"(1000(F8.5,','))")rhoflitter(:,J,1)
       end do
      close(3)
      open(3,file='rhoffwd.dat',form='formatted',status='unknown')
       do J=1,nyh
      WRITE(3,"(1000(F8.5,','))")rhoffwd(:,J,1)
       end do
      close(3)

      open(3,file='checkmoist.dat',form='formatted',status='unknown')
       do J=1,nyh
      WRITE(3,"(1000(F8.5,','))")moist(:,J,k)
       end do
      close(3)

      open(3,file='checkdead.dat',form='formatted',status='unknown')
       do J=1,nyh
      WRITE(3,"(1000(F8.5,','))")rhofdead(:,J,k)
       end do
      close(3)

      open(3,file='checkrhof.dat',form='formatted',status='unknown')
       do J=1,nyh
      WRITE(3,"(1000(F8.5,','))")rhof(:,J,k)
       end do
      close(3)

      open(3,file='checkpotr.dat',form='formatted',status='unknown')
       do J=1,nyh
      WRITE(3,"(1000(F8.5,','))")rhofdecid(:,J,k)
       end do
      close(3)

      open(3,file='checkafd.dat',form='formatted',status='unknown')
       do J=1,nyh
      WRITE(3,"(1000(F8.5,','))")actualfueldepth(:,J,1)
       end do
      close(3)
c      print*,actualfueldepth(:,50,1)
C      open(3,file='checktop.dat',form='formatted',status='unknown')
c       do J=1,nyh
c      WRITE(3,"(1000(F8.5,','))")ztopcellout(:,J,1)
c       end do
c      close(3)
C      open(3,file='checkbottom.dat',form='formatted',status='unknown')
C       do J=1,nyh
C      WRITE(3,"(1000(F8.5,','))")zbottomcellout(:,J,1)
C       end do
C      close(3)
Cc
cC       print*, moist(:,:,1)
c
c
c        do ix=1,nxh
c          do iy=1,nyh
c            do iz=1,nzh
c        if(abs(moist(ix,iy,iz)).lt.huge(moist(ix,iy,iz)))print*,ix,iy,iz
c            enddo
c          enddo
c        enddo

        deallocate(rhofgrass,rhoflitter,rhoffwd)
c     Here we are adding the deciduous leafless trees into rhof.. The deciduous fuel is kept seperate since 
c     it has extremely high moisture content.

      do iz=1,nzh
        do iy=1,nyh
          do ix=1,nxh
           rhof(ix,iy,iz)=rhof(ix,iy,iz)+rhofdecid(ix,iy,iz)+
     +      rhofdead(ix,iy,iz)
c     +      + rhofdead(ix,iy,iz)
          enddo
        enddo
      enddo
        deallocate(rhofdecid,rhofdead)

c         print*,"grass,litter,fwd,seedling,can:",
c     +  fmoistgr,fmoistlitter,fmoistfwd,fmoistseed,
c     +  fmoistcan

      
c
c  write out fuels data in "old" FIRETEC format
      open (1,file='treesrhof.dat',form='unformatted',status='unknown')
      write (1) rhof
      close (1)
      open (1,file='treesfueldepth.dat',form='unformatted'
     &      ,status='unknown')
      write (1) actualfueldepth
      close (1)

      open (1,file='treesss.dat',form='unformatted',status='unknown')
      write (1) sizescale
      close (1)
      open (1,file='treesmoist.dat',form='unformatted',status='unknown')
      write (1) moist
      close (1)

c      open(3,file='checkdepth.dat',form='formatted',status='unknown')
c       do J=1,nyh
c      WRITE(3,"(1000(F6.4,','))")((actualfueldepth(:,J,1)))
c       end do
c      close(3)

c      open(3,file='checkmoisture.dat',form='formatted',status='unknown')
c       do J=1,nyh
c      WRITE(3,"(1000(F8.5,','))")moist(:,J,1)
c      end do
c      close(3)



      ! open(3,file='rhofmulch.dat',form='formatted',status='unknown')
       ! do J=1,nyh
      ! WRITE(3,"(1000(F8.5,','))")rhofmulch(:,J,1)
      ! end do
      ! close(3)

c     open(3,file='rhofcheck.dat',form='formatted',status='unknown')
c      do J=1,nyh
c     WRITE(3,"(1000(F8.5,','))")rhof(:,J,1)
c     end do
c     close(3)



c       print*,'rhof',rhof(100:140,100,1)
c       print*,'rhofseedling',rhofseedling(100:140,100,1)
c              print*,'rhofgrass',rhofgrass(100:140,100,1)
c              print*,'rhoflitter',rhoflitter(100:140,100,1)           
c              print*,'rhoffwd',rhoffwd(100:140,100,1)      

c  print summary
      print *,' '
      sum=0
      v1=rhof(1,1,2)
      v2=rhof(1,1,2)
      do k=2,nzh
        do j=1,nyh
          do i=1,nxh
            sum=sum+rhof(i,j,k)
            if (rhof(i,j,k).lt.v1) v1=rhof(i,j,k)
            if (rhof(i,j,k).gt.v2) v2=rhof(i,j,k)
          enddo
        enddo
      enddo
      print *,'Canopy rhof min,max,avg = ',v1,v2,
     &                             sum/float(nxh*nyh*(nzh-1))
      sum=0
      v1=rhof(1,1,1)
      v2=rhof(1,1,1)
      do j=1,nyh
        do i=1,nxh
          sum=sum+rhof(i,j,1)
          if (rhof(i,j,1).lt.v1) v1=rhof(i,j,1)
          if (rhof(i,j,1).gt.v2) v2=rhof(i,j,1)
        enddo
      enddo
      print *,'Ground rhof min,max,avg = ',v1,v2,sum/float(nxh*nyh)
      sum=0
      v1=moist(1,1,2)
      v2=moist(1,1,2)
      do k=2,nzh
        do j=1,nyh
          do i=1,nxh
            sum=sum+moist(i,j,k)
            if (moist(i,j,k).lt.v1) v1=moist(i,j,k)
            if (moist(i,j,k).gt.v2) v2=moist(i,j,k)
          enddo
        enddo
      enddo
      print *,'Canopy moist min,max,avg = ',v1,v2,
     &                             sum/float(nxh*nyh*(nzh-1))
      sum=0
      v1=moist(1,1,1)
      v2=moist(1,1,1)
      do j=1,nyh
        do i=1,nxh
          sum=sum+moist(i,j,1)
          if (moist(i,j,1).lt.v1) v1=moist(i,j,1)
          if (moist(i,j,1).gt.v2) v2=moist(i,j,1)
        enddo
      enddo
      print *,'Ground moist min,max,avg = ',v1,v2,sum/float(nxh*nyh)
      sum=0
      v1=actualfueldepth(1,1,1)
      v2=actualfueldepth(1,1,1)
      do j=1,nyh
        do i=1,nxh
          sum=sum+actualfueldepth(i,j,1)
          if (actualfueldepth(i,j,1).lt.v1) v1=actualfueldepth(i,j,1)
          if (actualfueldepth(i,j,1).gt.v2) v2=actualfueldepth(i,j,1)
        enddo
      enddo
      print *,'Ground actualfueldepth min,max,avg = ',v1,v2,
     +         sum/float(nxh*nyh)
      sum=0
      v1=actualfueldepth(1,1,2)
      v2=actualfueldepth(1,1,2)
      do k=1,nzh
        do j=1,nyh
          do i=1,nxh
            sum=sum+actualfueldepth(i,j,k)
            if (actualfueldepth(i,j,k).lt.v1) v1=actualfueldepth(i,j,k)
            if (actualfueldepth(i,j,k).gt.v2) v2=actualfueldepth(i,j,k)
          enddo
        enddo
      enddo
      print *,'Total actualfueldepth min,max,avg = ',v1,v2,
     &                             sum/float(nxh*nyh*nzh)
      sum=0
      v1=sizescale(1,1,2)
      v2=sizescale(1,1,2)
      do k=1,nzh
        do j=1,nyh
          do i=1,nxh
            sum=sum+sizescale(i,j,k)
            if (sizescale(i,j,k).lt.v1) v1=sizescale(i,j,k)
            if (sizescale(i,j,k).gt.v2) v2=sizescale(i,j,k)
          enddo
        enddo
      enddo
      print *,'Total sizescale min,max,avg = ',v1,v2,
     &                             sum/float(nxh*nyh*nzh)

      print *,' '

        deallocate(sizescale,rhof,moist,x,y,actualfueldepth,z,
     +    zedge,genie,zs,gi,c13,c23,ztopcellout,zbottomcellout)
c
      end

      subroutine read_canada
c
      use higradinfo
      use treedatinfo

          
c      treefile='mixedwood_matrix_70pct_trees.txt'
c      open (1,file='Mixedwood_stems_swind_domain.csv'
c     +  ,form='formatted',status='old')
c      ntreedat=1
c      do while (ios.ge.0)
c        read (1,*,iostat=ios) xloc(ntreedat),
c     +       yloc(ntreedat),height(ntreedat),cr1(ntreedat),
c     +       htlc(ntreedat),treetype(ntreedat),health(ntreedat)
cc1100   format(i6,f5.2,f5.2,a3,a6,f5.1,f5.1,a8,a3,f6.1,f4.1,f5.1,
c    +         f5.1,f4.1,f4.1,a5,a1)
ccc gm        if (ntreedat.le.100)
ccc gm     &    print *,'rec ',ntreedat,xloc(ntreedat),yloc(ntreedat),
ccc gm     &         height(ntreedat),cr1(ntreedat),htlc(ntreedat)
c        if (ios.lt.0) then
c          ntreedat=ntreedat-1     ! end-of-file
c        else if (.not.(xloc(ntreedat).eq.0.and.
c     +                 yloc(ntreedat).eq.0)) then
c          ntreedat=ntreedat+1     ! not a blank line in excel file
c        endif
c      enddo
c      close (1)
c      
cc      print*,treetype(10000:10005)
cc      print*,height(10000:10005)
cc      print*,xloc(10000:10005)
cc      print*,htlc(10000:10005)
cc      print*,cr1(10000:10005)
c
c      htavg=0.0
c      do i=1,ntreedat
c        itreenum(i)=i
cc        htlc(i)=0.1
c        htavg=htavg+height(i)
cccc gm        if (height(i).gt.2.0) print *,'height ',i,height(i)
c      enddo
c       
c      xmax=maxval(xloc(1:ntreedat))
c      ymax=maxval(yloc(1:ntreedat))
cc       print*,'xmax =',xmax,'ymax =',ymax
c      hmin=minval(height(1:ntreedat))
c      hmax=maxval(height(1:ntreedat))
c      htavg=htavg/float(ntreedat)
c
c        
c       print *,'End read_canada - ntreedat = ',ntreedat
cccc gm      print *,'height min,max,avg ',hmin,hmax,htavg
ccc
      return
      end




      subroutine rhoinfo_init
c
      use rhoinfo
c
      do k=1,nzh
        do j=1,nyh
          do i=1,nxh
            nullify(rhead(i,j,k)%rptr)
            nullify(rtail(i,j,k)%rptr)
          enddo
        enddo
      enddo
      rhofirst=.false.
c
      return
      end
c
      subroutine rhoinfo_add(ixh,iyh,izh,itreenum,rhoamount)
c
      use rhoinfo
c
      if (.not.associated(rhead(ixh,iyh,izh)%rptr)) then
        allocate(rho_info,stat=istat)
        rhead(ixh,iyh,izh)%rptr => rho_info
        rtail(ixh,iyh,izh)%rptr => rho_info
        rho_info%itdtreenum = itreenum
        rho_info%rho_amount = rhoamount
        nullify(rho_info%rnext)
      else
        allocate(rho_info,stat=istat)
        rtempptr => rtail(ixh,iyh,izh)%rptr
        rtempptr%rnext => rho_info
        rtail(ixh,iyh,izh)%rptr => rho_info
        rho_info%itdtreenum = itreenum
        rho_info%rho_amount = rhoamount
        nullify(rho_info%rnext)
      endif
c
      return
      end
c
      subroutine rhoinfo_write(lun)
c
      use rhoinfo
c
      do k=1,nzh
        do j=1,nyh
          do i=1,nxh
            rtempptr => rhead(i,j,k)%rptr
            if (associated(rtempptr)) then
              icount=0
              do while (associated(rtempptr))
                icount=icount+1
                rtempptr => rtempptr%rnext
              enddo
              write (1) i,j,k,icount
c             print *,i,j,k,icount
              rtempptr => rhead(i,j,k)%rptr
              do while (associated(rtempptr))
                write (1) rtempptr%itdtreenum,rtempptr%rho_amount
c               print *,rtempptr%itdtreenum,rtempptr%rho_amount
                rtempptr => rtempptr%rnext
              enddo
            endif
          enddo
        enddo
      enddo
c
c  write a record to indicate end-of-data (all 0's)
      write (1) 0,0,0,0
c
      return
      end

      subroutine metryc(x,y,z,zs,gi,gmul,c13,c23,dx,dy,dz,dt,zb,j3)
c  ***************
c  NOTE:  n, m, l must be the same as nxh,nyh,nzh in subroutine trees
c  ***************
      parameter(n=750,m=1000,l=61)
      dimension x(n),y(m),z(l)
      dimension go(n,m)
      real c13(n,m),c23(n,m)
      real gmul(l),zs(n,m),gi(n,m,l)
      common/cyclbc/ ibcx,ibcy,irlx,irly
      common/met2/ aa1,aa2,aa3,zdata(100),zcrdata(100),
     +             zcoeff(100),npoints
      data zcrdata/0.0 ,    2.0000,    4.000,     6.00,       8.000,
     +            10.000,  14.000 ,   18.000,    22.00 ,     26.00 ,
     +            30.00 ,  34.00  ,   40.00 ,    50.00 ,     70.00 ,
     +           100.00 , 130.00  ,  160.00 ,   200.00 ,    250.00 ,
     +           300.00 , 350.00  ,  450.00 ,   550.00 ,    750.00 ,
     +           950.00 , 1150.0   , 1400.0  , 1700.0  ,   2000.0  ,
     +          2400.,69*0./
      dxi=1./dx
      dyi=1./dy 
      npoints=31
c      allocate(gi(n,m,l))
 
! deformation of sigma coordinate is done using polynomial fit
! aa3*sigma**3+aa2*sigma**2+aa1*sigma=sigma0
! domain for both sigma sysyems is 0 <= sigma  <= zb,
!                                  0 <= sigma0 <= zb

!                          ! aa1 determines compression of sigma at surface
c     aa1=1.0              ! aa1 can vary from 0 to 1 (no stretching)
!     aa1=0.0
      aa1=1./10.
                           ! if aa1=0 spline data will be used
      f=0.0                 ! 0 <= f <= 1
                           ! f=0, pure cubic fit
                           ! f=1, pure quadratic fit
 
      aa2=f*(1-aa1)/zb           ! (don't change this scaling constraint)
      aa3=(1-aa2*zb-aa1)/zb**2   ! (don't change this scaling constraint)
 
      if(aa1.eq.0.)then          ! use spline data instead of deformation
         do i=1,npoints          ! polynomial
            zdata(i)=z(i)*zb/z(npoints)
         enddo

c get spline coefficients for ginverse and dgdsigma
 
         call spline(zdata,zcrdata,npoints,99.e31,99.0e31,zcoeff)
      endif
 
      do 2 k=1,l
      sigma0=gdeform(z(k),0)
      gmul(k)=(zb-sigma0)/gdeform(z(k),1)
      do 2 j=1,m
      do 2 i=1,n
    2 gi(i,j,k)=zb/(zb-zs(i,j))/gdeform(z(k),1)
 
      do 3 j=1,m
      do 3 i=1,n
    3 go(i,j)=zb/(zb-zs(i,j))
 
      do 40 j=1,m
      do 41 i=2,n-1
   41 c13(i,j)=.5*dxi*(1./go(i+1,j)-1./go(i-1,j))*go(i,j)
      c13(1,j)=0.+ibcx*.5*dxi*(1./go(2,j)-1./go(n-1,j))*go(1,j)
   40 c13(n,j)=c13(1,j)
      do 50 i=1,n
      do 51 j=1+j3,m-j3
   51 c23(i,j)=.5*dyi*(1./go(i,j+j3)-1./go(i,j-j3))*go(i,j)
      c23(i,1)=0.+ibcy*.5*dyi*(1./go(i,1+j3)-1./go(i,m-j3))*
     $          go(i,1)
   50 c23(i,m)=c23(i,1)

c      deallocate(gi)
 
      return
      end

 
      function zcart(sigma,i,j,zs,zb)
! sigma is sigma coordinate, zcart is cartesian vertical coordinate
c  ***************
c  NOTE:  n, m, l must be the same as nxh,nyh,nzh in subroutine trees
c  ***************
      parameter(n=750,m=1000,l=61)
       real zs(n,m)
      common/met2/ aa1,aa2,aa3,zdata(100),zcrdata(100),
     +             zcoeff(100),npoints
      zcart=gdeform(sigma,0)*(zb-zs(i,j))/zb+zs(i,j)
      return 
      end 

 
      function gdeform(sigma,iflag)
! iflag=0, compute gdeform(sigma)=sigma0,
! iflag=1, return derivative of gdeform
! sigma0=aa3**sigma**3+aa2*sigma**2+aa1*sigma   ! cubic polynomial fit
      common/met2/ aa1,aa2,aa3,zdata(100),zcrdata(100),
     +             zcoeff(100),npoints
      if(aa1.eq.0)then           ! use spline
         call splint(zdata,zcrdata,zcoeff,npoints,sigma,answer,iflag)
         gdeform=answer
      else                       ! use cubic polynomial
         if(iflag.eq.0)gdeform=aa3*sigma**3+aa2*sigma**2+aa1*sigma
         if(iflag.eq.1)gdeform=3*aa3*sigma**2+2*aa2*sigma+aa1
      endif
      return
      end

 
      SUBROUTINE splint(xa,ya,y2a,n,x,y,kderivative)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else 
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
c      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      if(kderivative.eq.0)then
         y=a*ya(klo)+b*ya(khi)+( (a**3-a)*y2a(klo) + (b**3-b)*y2a(khi) )
     +    *(h**2)/6.
      else
         y=(ya(khi)-ya(klo))/h
     +    -(  (3*a**2-1)*y2a(klo)-(3*b**2-1)*y2a(khi)  )*h/6
      endif
      return
      END

 
! subroutines spline and splint are from "Numerical Recipes". Splint
! is modified to give derivatives as well as interpolated values
 
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
 
     *u(i-1))/p
 11   continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
 12   continue
      return
      END
