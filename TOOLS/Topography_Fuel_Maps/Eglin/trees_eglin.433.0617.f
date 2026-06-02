      module treedatinfo
c
c  x,y size for tree data (meters) (actual sampled/measured data from 
c  the Forest Service, LANL, etc)
        parameter (nxtd=60,nytd=106)
c
c  max # trees allowed in the sampled/measured (excel) data
        parameter (ntr=120000)
c
c  type/source of data
c    itdsource=0 -> Forest Service Rocky Mountain Research Station data
c    itdsource=1 -> LANL XL06 data 
c    itdsource=2 -> Phil Dennison's Calabasas fire data 
c    itdsource=3 -> General California chaparral 
c    itdsource=4 -> Las Vegas/Gallinas - douglas fir 
c    itdsource=5 -> Las Vegas/Gallinas - aspen
c    itdsource=6 -> eglin 0107U001 data 6/15/12
        parameter (itdsource=6)
c
c  rotate tree data (0 => no, 1 => yes - so, if 20m X 50m data,
c  0 => 50 aligns with Y axis, 1=> 50 aligns with X axis)
        parameter (itdrotate=0)
c
c  itdloc=0 => spread FS tree data in random locations within any nxtd X
c  nytd section of the HIGRAD grid; itdloc=1 => spread FS tree data in
c  uniform nxtd X nytd patterns; itdloc=2 => read tree locations and 
c  information from a previous run of this program from file oldtreefile;
c  itdloc=3 => do not populate, just use trees as read in
c    and
c  use computer generated seed for randomly placing trees or specify the
c  seed, which will allow you to duplicate a specific random placing (only 
c  used when itdloc=0) - itdranseed=0 => use computer generated
c  seed; itdranseed=1 => use integer seed values specified in array
c  itdseed 
        parameter (itdloc=0,itdranseed=1)
        integer, dimension(4) :: itdseed = (/8,6,4,2/)
c
c  arrays and variables for reading sampled tree data
        integer itreenum(ntr),ntreedat
        real xloc(ntr),yloc(ntr),dbh(ntr),dmr(ntr),height(ntr),
     +       htlc(ntr),htdc(ntr),htdt(ntr),cr1(ntr),cr2(ntr)
        character species(ntr)*8,live(ntr)*1,spp(ntr)*3,status(ntr)*6,
     +            bole_scar(ntr)*8,bss(ntr)*3,keens(ntr)*5,
     +            position(ntr)*1
        character treefile*40,dummy*1 ! name of tree data file
        character oldtreefile*60      ! name of file with old tree info

      end module treedatinfo
 
      module higradinfo
c
c  x,y,z size for HIGRAD grid and x,y,z grid resolution (meters)
c  NOTE  ******  DO NOT CHANGE NXH, NYH, NZH HERE WITHOUT CHANGING
c                N, M, L IN SUBROUTINES METRYC AND ZCART
        parameter (nxh=800,nyh=150,nzh=41,dxh=2.0,dyh=2.0,dzh=15.0)
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
        parameter (itminxh=1,itmaxxh=800,itminyh=1,itmaxyh=150) 
c
c  area of HIGRAD grid (min -> max) to be treated as a grass area
c  (this can overlap with tree areas, or be just grass) - for no
c  grass, set igminxh=0
        parameter (igminxh=1,igmaxxh=800,igminyh=1,igmaxyh=150) 
c
c  area of HIGRAD grid (min -> max) to have the possibility of litter
c  (this can overlap with tree or grass) 
        parameter (ilminxh=1,ilmaxxh=800,ilminyh=1,ilmaxyh=150) 
c
c  consider northern/southern exposure effects (more trees on northern
c  facing slopes, fewer on southern) - 0=>don't consider, 1=>do consider  
c  (NOTE: this will be ignored if topography is not being used)
        parameter (treeslopefactor=0)
c
c  if considering north/south effect, percentage more trees on north face
c  and less on south
        parameter (slopeeffect=0)
c
c  if considering north/south effect, what is the orientation of the HIGRAD
c  grid - 0=>north is in the U direction; 1=>north is in the V direction
        parameter (gridorientation=1)
c
c  arrays and variables used in populating the higrad grid with trees
        integer itreeh(ntrh),itreehtmp(ntrh),inewtr 
        real xloch(ntrh),yloch(ntrh)
        real xlochtmp(ntrh),ylochtmp(ntrh)
c
c  arrays for HIGRAD data
        real rhof(nxh,nyh,nzh)         ! fuel density for each HIGRAD cell
c
c  SPECIAL for eglin
        real rhofllpine(nxh,nyh,nzh)
        real rhoftoak(nxh,nyh,nzh)
        real rhofpersim(nxh,nyh,nzh)

        real rhofgrass(nxh,nyh,nzh)    ! array for saving grass contributions
        real rhoflitter(nxh,nyh,nzh)   ! array for saving litter contributions
        real sizescale(nxh,nyh,nzh)    ! fuel size scale for each HIGRAD cell
        real moist(nxh,nyh,nzh)        ! fuel moisture or each HIGRAD cell
c
c  add for eglin
        integer open(nxh,nyh)
c
c  array for average height of fuel within a cell and initialize
        real, dimension (nxh,nyh,nzh) :: actualfueldepth = -1.0
c
c  variable for topography file name and array for topography data
        character (len=40) :: topofile      ! name of topography file
        real zs(nxh,nyh)                    ! topography data
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
        parameter (ithinmethod=0)
c
c  area of HIGRAD grid (min -> max) to thin (the area must have trees 
c  populated on it) - ignored if no thinning is being done or a mask that
c  indicates areas to be thinned is being used
        parameter(ithinminxh=1,ithinmaxxh=800,ithinminyh=1,
     +                                            ithinmaxyh=150) 
c
c  use a mask file (an unformatted 2D file of integers, size (nxh,nyh),
c  where 0 => don't thin in this HIGRAD cell, 1 => do thin in this cell) 
c  to indicate areas to be thinned - ithinmask=0 => no mask file, 
c  ithinmask=1 => use mask file (will be prompted for name)
        parameter(ithinmask=0)
c
c  use a mask file (same as above) to determine a select area to compute
c  canopy bulk density for - ithinmaskbd=0 => no mask file,
c  ithinmaskbd=1 => use mask file (will be prompted for name)
        parameter(ithinmaskbd=0)
c
c  for thinning based on tree diameter, diameter under which trees 
c  will be removed (dbh less than this value)
        parameter (dbhthin=28.0)
c
c  for patch thinning, use circles of radius "patchradius" meters - for
c  eliptical patches, make them "patchradius" meters in x and "patchradiusy"
c  meters in y (NOTE:  For ithinmethod=2,3,or 5, you are specifying the size
c  of the patches to be left; for ithinmethod=4, you are specifying the size
c  of the patches to be removed.)
        parameter (patchradius=35.0,patchradiusy=12.0)
c       parameter (patchradius=12.0,patchradiusy=3.0)
c
c  for patch thinning, estimate 1 patch every patchperiod meters
        parameter (patchperiod=50.0)
c
c  for patch thinning, use computer generated seed for randomly placing 
c  patches or specify the seed, which will allow you to duplicate a 
c  specific random placing - ipatchseed=0 => use computer generated 
c  seed; ipatchseed=1 => use integer seed values specified in array
c  ipseed 
        parameter (ipatchseed=0)
        integer, dimension(4) :: ipseed = (/2,4,6,8/)
c
c  arrays and variables used in thinning
        real patchx(ntrh),patchy(ntrh)
        integer thinmask(nxh,nyh),nmaskcells
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
c  for eglin, will use rhoavg=0.3 for Longleaf pine and Persimmon, and
c             rhoavg=0.4 for Turkey oak  (**  SET BELOW  **)
        parameter (rhoavg=0.4,rhomax=(3.0/2.0)*rhoavg)
c
c  actual open grassland grass density (within the fuelbed) (kg/m3)
c  for eglin...  3.5 tons/acre = 0.7865 kg/m2
c                0.7865 kg/m2 / .5 m (fuel bed depth) = 1.573 kg/m3
c                   (***  SET BELOW  ***)
c       parameter (actualrhograss=1.573)
c
c  actual open grassland grass height (meters)   (in the clearings)
c  for eglin, this is .5m everywhere  (***  SET BELOW  ***)
        parameter (actualgrassht=0.5)
c
c  exponential factor reducing grass density due to canopy shading
c  for eglin...  surface fuels are spread uniformly below
c       parameter (grassconstant=0.0)
c
c  estimated maximum litter density (kg/m^3) (within the fuelbed) 
c  for eglin...  no litter beyond the "grass" above
c       parameter (actualrholitter=1e-6)
c
c  estimated maximum litter height (m)   (under the trees)
c       parameter (actuallitterht=0.0001)
c
c  exponential factor enhancing litter density due to canopy shading
c       parameter (rlitterconstant=0.0)
c
c  nominal fuel moisture for canopy and grass (1.0 -> 100%)
c  for eglin, will use 133% for Longleaf pine, 170% for Persimmon, and
c             200% for Turkey oak, 8% for grass and litter  (** SET BELOW **)
        parameter (fmoistcan=1.0,fmoistgr=.08,fmoistlitter=.08)
c
c  nominal fuel size scale (solids) - meters
c  for eglin, will use 0.0005 for Longleaf pine and Persimmon, and
c            0.0002 for Turkey oak (from Kevin Hiers, 
c                          "/Q. laevis/  SAV 96.7 (7.1 SD) cm2/cm3", so
c            97 cm2/cm3 = 9700 m2/m3 and SAV=2/SS so 2/9700=0.0002)  ( ** SET BELOW **)
        parameter (fsizescale=0.0005)
c
c  SPECIAL for eglin - data type
c    ieglintype=0 => overstory only
c    ieglintype=1 => overstory and well-managed midstory (later called "open")
c    ieglintype=2 => overstory and unmanaged midstory (later called "midstory")
        parameter(ieglintype=2)
c
c  SPECIAL for eglin - indicate whether it is the "wet" (0) or "dry" (1) season -
c  in the "dry" season, all midstory DIOVIR are killed and the moisture for
c  all QUELAE is set to 15%
        parameter (ieglindry=0)
c
c  data for headers in output files that are generated (treesfueldata.dat) -
c  convertto = 'big_endian' or 'little_endian', ireal = 4 or 8 (real numbers
c  single or double precision)
        character (len=50) :: convertto='big_endian'
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
      real x(nxh,nyh,nzh),y(nxh,nyh,nzh)  ! for metryc
      real gi(nxh,nyh,nzh),gmul(nzh)      ! for metryc
      real c13(nxh,nyh),c23(nxh,nyh)      ! for metryc
      real zedge(nzh+1)             ! bottom of cell before stretching
      character (len=50) :: dataname      ! data type name for output file 
c
c  JLW temp
      real road(nxh,nyh)
c
c  initialize
      rhof=0.0
c
c  SPECIAL for eglin
      rhofllpine=0.0
      rhoftoak=0.0
      rhofpersim=0.0

      rhofgrass=0.0
      rhoflitter=0.0
      sizescale=fsizescale
      actualfueldepth=0.0
      moist=0.00001
      call rhoinfo_init
c
c  call routine to read tree data
      if (itdsource.eq.0) then
        call read_fsdata
      else if (itdsource.eq.1) then
        call read_lanldata
      else if (itdsource.eq.2) then
        call read_caladata
      else if (itdsource.eq.3) then
        call make_chapdata
      else if (itdsource.eq.4) then
        call make_dougfir
      else if (itdsource.eq.5) then
        call make_aspen
      else if (itdsource.eq.6) then
        call read_eglin
      endif
      print *,'After reading tree data  - ntreedat = ',ntreedat
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
      write (*,1200) 
 1200 format(1x,'Topography file (2D binary) or flat: ',$)
      read (*,*) topofile
      if (topofile.eq.'flat') then
        zs=1.0
      else
        open (1,file=topofile,form='unformatted',status='old')
        read (1) zs
        close (1)
      endif
c
c  populate the HIGRAD grid with trees or read old tree locations
c  and info from an old run of this program
      if (itdloc.le.1) then 
        call populate
      else if (itdloc.eq.2) then
        call readoldtrees
      else if (itdloc.eq.3) then
        do i=1,ntreedat
          xloch(i)=xloc(i)
          yloch(i)=yloc(i)
          itreeh(i)=itreenum(i)
        enddo
        inewtr=ntreedat
      endif
      print *,'After populating - inewtr = ',inewtr
c
c  "thin" trees if requested
      if (ithinmethod.ne.0.or.ithinmaskbd.eq.1) then
        call thintrees
        if (ithinmethod.ne.0)print *,'After thinning - inewtr = ',inewtr
      endif
c
c  generate file for plotting tree locations with xmgrace
c     open (1,file='treesloc.dat',form='formatted',status='unknown')
c     do i=1,inewtr
c       write (1,*) xloch(i),',',yloch(i)
c     enddo
c     close (1)
c
c  initialize
      do iz=1,nzh
        do iy=1,nyh
          do ix=1,nxh
            x(ix,iy,iz)=(ix-1)*dxh-((nxh+1)/2-1)*dxh
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
      do iz=1,nzh+1
        zedge(iz)=(iz-1)*dzh
      enddo
      zb=zedge(nzh+1)
      call metryc(x,y,zedge,zs,gi,gmul,c13,c23,dxh,dyh,dzh,dt,zb,1)
c
      call canopy(x,y,zedge,zb,.true.,bdcanopy,treehtmax)

      print *,'Canopy Bulk Density = ',bdcanopy
      if (ithinmask.eq.1.or.ithinmaskbd.eq.1) then
        bdthincanopy=thinmcanopysum/(avgcanht*(nmaskcells*dxh*dyh))
        print *,'Masked Thinning Area Canopy Bulk Density = ',
     +                                                bdthincanopy
      endif
      print *,'Tree Height Max = ',treehtmax
      overstoryht=treehtmax
c
c  SPECIAL eglin midstory section...
      if (ieglintype.gt.0) then
        ntreedat=0
c
c  create midstory shrub data if needed
        if (ieglintype.eq.1) then
          call eglin_midstory_wellmanaged
        else if (ieglintype.eq.2) then
          call eglin_midstory_unmanaged
        endif

        do i=1,ntreedat
          xloch(i)=xloc(i)
          yloch(i)=yloc(i)
          itreeh(i)=itreenum(i)
        enddo
        inewtr=ntreedat
c
c  loop to process each shrub as done above for the overstory
        call canopy(x,y,zedge,zb,.false.,bdcanopy,treehtmax)

        print *,'Canopy Midstory Bulk Density = ',bdcanopy
        print *,'Midstory Tree Height Max = ',treehtmax
        if (overstoryht.gt.treehtmax) treehtmax=overstoryht
      endif
c
c  end SPECIAL eglin midstory section...

c
c  SPECIAL for eglin case - deal with "thinning" by rhomax (which is
c  0.4 for llpine and persim and 0.7 for toak) for each tree type here
c  NOTE:  At this point, rhof is assumed to be the sum of all tree types,
c  and not limited by any "rhomax"
      do iz=1,nzh
        do iy=1,nyh
          do ix=1,nxh
            if (rhof(ix,iy,iz).gt.0.0) then

              rhomaxcell=(rhofllpine(ix,iy,iz)*(3.0/2.0)*0.3+
     &                    rhoftoak(ix,iy,iz)*(3.0/2.0)*0.4+
     &                    rhofpersim(ix,iy,iz)*(3.0/2.0)*0.3)/
     &                                           rhof(ix,iy,iz)
c             rhomaxcell=(3.0/2.0)*(rhofllpine(ix,iy,iz)*0.4+
c    &                              rhoftoak(ix,iy,iz)*0.7+
c    &                              rhofpersim(ix,iy,iz)*0.4)/
c    &                                           rhof(ix,iy,iz)

              rllpine=rhofllpine(ix,iy,iz)*min(1.0,
     &                              rhomaxcell/rhof(ix,iy,iz))
              rtoak=rhoftoak(ix,iy,iz)*min(1.0,
     &                              rhomaxcell/rhof(ix,iy,iz))
              rpersim=rhofpersim(ix,iy,iz)*min(1.0,
     &                              rhomaxcell/rhof(ix,iy,iz))

c             if (iz.eq.2.and.rllpine.ne.0.0.and.rtoak.ne.0.0.
c    &            and.rpersim.ne.0.0) 
              if (iz.eq.2.and.rhof(ix,iy,iz)-
     &             (rllpine+rtoak+rpersim).gt.0.1) then
                 print *,rhof(ix,iy,iz),rllpine+rtoak+rpersim,
     &                             rhofllpine(ix,iy,iz),rllpine,
     &                             rhoftoak(ix,iy,iz),rtoak,
     &                             rhofpersim(ix,iy,iz),rpersim,
     &                             rhomaxcell,
     &                         rhof(ix,iy,iz)-(rllpine+rtoak+rpersim)
                 print *,' '
               endif    
              rhof(ix,iy,iz)=rllpine+rtoak+rpersim
              rhofllpine(ix,iy,iz)=rllpine
              rhoftoak(ix,iy,iz)=rtoak
              rhofpersim(ix,iy,iz)=rpersim
            endif
          enddo
        enddo
      enddo
crrl            if (rhof(ix,iy,iz).gt.0.0) then
crrl              rllpine=(min(rhofllpine(ix,iy,iz),(3.0/2.0)*0.4)*
crrl     +                                  rhofllpine(ix,iy,iz))/
crrl     +             (rhofllpine(ix,iy,iz)+rhoftoak(ix,iy,iz)+
crrl     +                                    rhofpersim(ix,iy,iz))
crrl              rtoak=(min(rhoftoak(ix,iy,iz),(3.0/2.0)*0.7)*
crrl     +                                    rhoftoak(ix,iy,iz))/
crrl     +             (rhofllpine(ix,iy,iz)+rhoftoak(ix,iy,iz)+
crrl  +                                    rhofpersim(ix,iy,iz))
crrl        if (rtoak.gt.1.05) print *,'RTOAK ',ix,iy,iz,rtoak,
crrl     &                  rhoftoak(ix,iy,iz),rhofllpine(ix,iy,iz),
crrl     &                  rhofpersim(ix,iy,iz)
crrl              rpersim=(min(rhofpersim(ix,iy,iz),(3.0/2.0)*0.4)*
crrl     +                                 rhofpersim(ix,iy,iz))/
crrl     +             (rhofllpine(ix,iy,iz)+rhoftoak(ix,iy,iz)+
crrl     +                                    rhofpersim(ix,iy,iz))

c
c  in order to calculate the load for grass and litter on the ground and
c  a mass weighted height for the grass and litter, we identify areas where
c  there would be grass or litter and modify the amount based on canopy cover
c     if (igminxh.eq.0) goto 20
c       do iy=igminyh,igmaxyh
c         do ix=igminxh,igmaxxh
c  determine the amount of canopy directly above a given x and y location
c           rhocolumn=0.0
c           
c           if (itminxh.le.ix.and.itmaxxh.ge.ix.and.    ! trees and grass
c    +          itminyh.le.iy.and.itmaxyh.ge.iy) then
c             do ic=icminzh,icmaxzh
c               ztopcell=zcart(zedge(ic+1),ix,iy,zs,zb)-zs(ix,iy)
c               zbotcell=zcart(zedge(ic),ix,iy,zs,zb)-zs(ix,iy)
c               rhocolumn=rhocolumn+
c    +                 (rhof(ix,iy,ic)*(ztopcell-zbotcell))
c             enddo
c           else
c             rhocolumn=0.0
c           endif
c  rhocolumn will be used to scale the absence of grass and the abundance 
c  of litter - determine if there is grass possible in the cell
c           if (igminxh.le.ix.and.igmaxxh.ge.ix.and.    !identify areas with grass
c    +          igminyh.le.iy.and.igmaxyh.ge.iy) then
c  izg locates the top cell of the grass at a particular x and y on the 
c  HIGRAD grid
c             izg=1
c             zgtopcell=zcart(zedge(izg+1),ix,iy,zs,zb)-zs(ix,iy)
c             do while (zgtopcell.lt.actualgrassht)
c               izg=izg+1
c               zgtopcell=zcart(zedge(izg+1),ix,iy,zs,zb)-zs(ix,iy)
c             enddo
c  in the absence of canopy the grass load in a cell (k<=izg) will be 
c  calculated using
c   (actualrhograss kg/m3*(max(grass height(m),top of cell)-bottom of cell))/dz
c  shade factor is the fraction of the grass that remains under the
c  computed rhocolumn

c             shadefactor=exp(-1.0*grassconstant*
c    +                    (rhocolumn/max(rhomax*avgcanht,rhocolumn)))
c
c  loop from the ground up to the top of the grass and fill each cell with 
c  the appropriate density of grass depending on the amount of the cell that
c  is full and the shading
c             do iz=1,izg
c               actualfuelbottom= 0.     ! inserting 0 for ground
c               fueldistributionslope=0. ! pos # -> load higher at top (kg/m3/m)
c               actualgroundload=actualrhograss ! fuel load @ ground if linear etrap
c    +                   -.5*actualgrassht*fueldistributionslope 
c               rhomicro=500.            ! kg/m3
c               ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
c               zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
c               rhograss=shadefactor*(min(actualgrassht,ztopcell)-
c    +                       max(actualfuelbottom,zbotcell))/
c    +                       (ztopcell-zbotcell)*
c    +                       (actualgroundload +
c    +                       .5*fueldistributionslope*
c    +                       (min(actualgrassht,ztopcell)+
c    +                       max(actualfuelbottom,zbotcell)))
c               actualfueldepth(ix,iy,iz)=
c    +                       ((min(actualgrassht,ztopcell)
c    +                       -zbotcell)*rhograss
c    +                       +rhof(ix,iy,iz)
c    +                       *min(treehtmax,ztopcell)-zbotcell)
c    +                       /(rhograss+rhof(ix,iy,iz))
c               rhof(ix,iy,iz)=rhof(ix,iy,iz)+rhograss
c               rhofgrass(ix,iy,iz)=rhofgrass(ix,iy,iz)+rhograss
c             enddo
c           endif
c
c  now do a similar process for litter except we are working under the 
c  assumption that grass load diminishes under trees and litter increases -
c  check to see if we are in a litter area 
c           if (ilminxh.le.ix.and.ilmaxxh.ge.ix.and.    !identify areas with litter
c    +          ilminyh.le.iy.and.ilmaxyh.ge.iy) then
c  izl locates the top cell of the litter at a particular x and y on the higrad grid
c             izl=1
c             zltopcell=zcart(zedge(izl+1),ix,iy,zs,zb)-zs(ix,iy)
c             do while (zltopcell.lt.actuallitterht)
c               izl=izl+1
c               zltopcell=zcart(zedge(izl+1),ix,iy,zs,zb)-zs(ix,iy)
c             enddo
c  coverfactor is the fraction of the grass that remains under the
c  computed rhocolumn

c             coverfactor=1.-exp(-1.0*rlitterconstant*
c    +                 (rhocolumn/max(rhomax*avgcanht,rhocolumn)))
c
c  loop from the ground up to the top of the litter and fill each cell with 
c  the appropriate density of litter depending on the amount of the cell that
c  is full and the shading
c             do iz=1,izl
c               actualfuelbottom= 0.     ! inserting 0 for ground
c               fueldistributionslope=0. ! pos # -> load higher at top (kg/m3/m)
c               actualgroundload=actualrholitter ! fuel load @ ground if linear etrap
c    +                       -.5*actuallitterht*fueldistributionslope
c               rhomicro=500.            ! kg/m3
c               ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
c               zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
c               rholitter=coverfactor
c    +                    *(min(actuallitterht,ztopcell)
c    +                    -max(actualfuelbottom,zbotcell))/
c    +                    (ztopcell-zbotcell)*
c    +                    (actualgroundload +
c    +                    .5*fueldistributionslope*
c    +                    (min(actuallitterht,ztopcell)+
c    +                    max(actualfuelbottom,zbotcell)))
c               actualfueldepth(ix,iy,iz)=
c    +                    ((min(actuallitterht,ztopcell)-zbotcell)
c    +                    *rholitter
c    +                    +rhof(ix,iy,iz)*
c    +                    min(actualfueldepth(ix,iy,iz),ztopcell)
c    +                    -zbotcell)
c    +                    /(rholitter+rhof(ix,iy,iz))
c               rhof(ix,iy,iz)=rhof(ix,iy,iz)+rholitter
c               rhoflitter(ix,iy,iz)=rhoflitter(ix,iy,iz)+rholitter
c             enddo
c           endif
c
c         enddo
c       enddo
c20   continue

c
c  SPECIAL for eglin...
c
c 
c  spread the grass load evenly (no litter) - from the XCEL eglin data, the
c  3.5 tons/acre = 0.7865 kg/m2 
      ztopcell=zcart(zedge(2),1,1,zs,zb)-zs(1,1)
      zbotcell=zcart(zedge(1),1,1,zs,zb)-zs(1,1)
      hgt=ztopcell-zbotcell
      rhograss=0.7865/hgt

      do iy=1,nyh
        do ix=1,nxh
          rhoftree=rhof(ix,iy,1)
          actualfueldepth(ix,iy,1)=(
     +                    ((min(actualgrassht,ztopcell)-zbotcell)*
     +                    rhograss)+
     +                    ((min(treehtmax,ztopcell)-zbotcell)*
     +                    rhoftree))/
     +                    (rhograss+rhoftree)
          rhof(ix,iy,1)=rhof(ix,iy,1)+rhograss
          rhofgrass(ix,iy,1)=rhograss
          rhoflitter(ix,iy,1)=0.0
        enddo
      enddo
      print *,'Spread ',rhograss,rholitter,hgt

c
c  call routine to write out detailed information about which trees
c  contributed how much to the density in each HIGRAD cell, etc, for
c  later analysis (after the FIRETEC run)
c     call writetrees    ! this won't work with eglin mid-story stuff!
c
c  special - blank out road for Phil Dennison's Calabasas run
c     open (1,file='roads.sun.dat',form='unformatted',
c    +      status='unknown')
c     read (1) road
c     close (1)
c     do iy=1,nyh
c       do ix=1,nxh
c         if (road(ix,iy).eq.1.0) then
c           do iz=1,nzh
c             rhof(ix,iy,iz)=0.0
c           enddo
c         endif
c       enddo
c     enddo
c
c  eglin - correct way to do moist
      moist=0.00001
      do iz=1,nzh
        do iy=1,nyh
          do ix=1,nxh
            if (rhof(ix,iy,iz).gt.0.0) then
              if (ieglindry.eq.0) then
                moist(ix,iy,iz)=((rhofllpine(ix,iy,iz)*1.33)+
     +                          (rhoftoak(ix,iy,iz)*2.0)+
     +                          (rhofpersim(ix,iy,iz)*1.7)+
     +                          (rhofgrass(ix,iy,iz)*0.08)+
     +                          (rhoflitter(ix,iy,iz)*0.08))/
     +               (rhofllpine(ix,iy,iz)+rhoftoak(ix,iy,iz)+
     +                rhofpersim(ix,iy,iz)+rhofgrass(ix,iy,iz)+
     +                rhoflitter(ix,iy,iz))
              elseif (ieglindry.eq.1) then
                moist(ix,iy,iz)=((rhofllpine(ix,iy,iz)*1.33)+
     +                          (rhoftoak(ix,iy,iz)*0.15)+
     +                          (rhofpersim(ix,iy,iz)*1.7)+
     +                          (rhofgrass(ix,iy,iz)*0.08)+
     +                          (rhoflitter(ix,iy,iz)*0.08))/
     +               (rhofllpine(ix,iy,iz)+rhoftoak(ix,iy,iz)+
     +                rhofpersim(ix,iy,iz)+rhofgrass(ix,iy,iz)+
     +                rhoflitter(ix,iy,iz))
              else
                print *,'Error - need to set ieglindry'
                stop
              endif
            endif
          enddo
        enddo
      enddo
c
c  eglin - set sizescale
      sizescale=fsizescale
      do iz=1,nzh
        do iy=1,nyh
          do ix=1,nxh
            if (rhofllpine(ix,iy,iz)+rhoftoak(ix,iy,iz)+
     +                 rhofpersim(ix,iy,iz).gt.0.0) then
              sizescale(ix,iy,iz)=((rhofllpine(ix,iy,iz)*0.0005)+
     +                             (rhoftoak(ix,iy,iz)*0.0002)+
     +                             (rhofpersim(ix,iy,iz)*0.0005))/
     +                (rhofllpine(ix,iy,iz)+rhoftoak(ix,iy,iz)+
     +                 rhofpersim(ix,iy,iz))
            endif
          enddo
        enddo
      enddo
c
c  write out fuels data in "old" FIRETEC format
      open (1,file='treesrhof.dat',form='unformatted',status='unknown')
      write (1) rhof
      close (1)
      open (1,file='treesfueldepth.dat',form='unformatted'
     &      ,status='unknown')
      write (1) actualfueldepth
      close (1)

c  for eglin, sizescale set above
c     sizescale=fsizescale

      open (1,file='treesss.dat',form='unformatted',status='unknown')
      write (1) sizescale
      close (1)
      open (1,file='treesmoist.dat',form='unformatted',status='unknown')
      write (1) moist
      close (1)
c
c  print summary
      print *,' '
      canload=0.0
      sum=0
      v1=rhof(1,1,2)
      v2=rhof(1,1,2)
      do k=2,nzh
        do j=1,nyh
          do i=1,nxh
            sum=sum+rhof(i,j,k)
            hgt=zcart(zedge(k+1),i,j,zs,zb)-zcart(zedge(k),i,j,zs,zb)
            canload=canload+(rhof(i,j,k)*hgt)
            if (rhof(i,j,k).lt.v1) v1=rhof(i,j,k)
            if (rhof(i,j,k).gt.v2) v2=rhof(i,j,k)
          enddo
        enddo
      enddo
      print *,'Canopy fuel load (kg/m2) = ',canload/(nxh*nyh)
      print *,'Canopy rhof min,max,avg = ',v1,v2,
     &                             sum/float(nxh*nyh*(nzh-1))
      grload=0.0
      sum=0
      v1=rhof(1,1,1)
      v2=rhof(1,1,1)
      do j=1,nyh
        do i=1,nxh
          sum=sum+rhof(i,j,1)
          hgt=zcart(zedge(2),i,j,zs,zb)-zcart(zedge(1),i,j,zs,zb)
          grload=grload+(rhof(i,j,1)*hgt)
          if (rhof(i,j,1).lt.v1) v1=rhof(i,j,1)
          if (rhof(i,j,1).gt.v2) v2=rhof(i,j,1)
        enddo
      enddo
      print *,'Ground fuel load (kg/m2) = ',grload/(nxh*nyh)
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
c
      end

      subroutine canopy(x,y,zedge,zb,overstory,bdcanopy,treehtmax)
c
c  loop to process each tree canopy

      use treedatinfo
      use higradinfo
      use thininfo
      use constants

      real x(nxh,nyh,nzh),y(nxh,nyh,nzh)  ! for metryc
      real zedge(nzh+1)                   ! bottom of cell befo
      logical overstory
c
      print *,'In subroutine canopy - inewtr = ',inewtr
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
        cr=(cr1(itdindex)+cr2(itdindex))/2.0
c
c  SPECIAL for eglin case
        if (species(itdindex).eq.'QUELAE') then
          eglinmax=(3.0/2.0)*0.4
        else
          eglinmax=(3.0/2.0)*0.3
        endif
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
          if (it.eq.1)
     +        print *,' T0',zbotcell,ztopcell,htlctr,treeht
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
        if (it.le.5) print *,' T1',it,treeht,htlctr,hb,cl,cr,
     +               ixcell,iycell,nminx,nmaxx,nminy,nmaxy,
     +               nminz,nmaxz,icminzh,icmaxzh
c
c  for all cells that contain part of the tree...
        do iz=nminz,nmaxz
          do iy=nminy,nmaxy
            do ix=nminx,nmaxx
              rmcanopy=0.0
              ztopcell=zcart(zedge(iz+1),ix,iy,zs,zb)-zs(ix,iy)
              zbotcell=zcart(zedge(iz),ix,iy,zs,zb)-zs(ix,iy)
c  divide the HIGRAD cell into smaller areas, and sum over all the smaller
c  areas to get the mass of the canopy in the HIGRAD cell 
              ixloop=nint(dxh/dxs)
              iyloop=nint(dyh/dys)
              if (it.eq.1.and.ix.eq.nminx+2.and.iy.eq.nminy+2.and
     +            .iz.eq.nminz+2) print *,' T2',ztopcell,zbotcell,
     +             ixloop,iyloop
              do iyy=1,iyloop
                ylocs=iy*dyh+(real(iyy)-.5)*dys
                do ixx=1,ixloop
                  xlocs=ix*dxh+(real(ixx)-.5)*dxs
                  rprime2=(xlocs-xloch(it))**2+(ylocs-yloch(it))**2
                  ftop=((-1.0*cl)/(cr*cr))*rprime2+treeht
                  fbot=(hb/(cr*cr))*rprime2+htlctr
                  zmid=.5*(min(ftop,ztopcell)+max(fbot,zbotcell))
c
c  SPECIAL for eglin
                  rhoscale=((zmid+(cl/(cr*cr)*rprime2)-htlctr)/
     +                      (treeht-htlctr))*eglinmax
c                 rhoscale=((zmid+(cl/(cr*cr)*rprime2)-htlctr)/
c    +                      (treeht-htlctr))*rhomax

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
c
c  SPECIAL for eglin case...
              rhof(ix,iy,iz)=rhof(ix,iy,iz)+rhoftree
              if (.not.(species(itdindex).eq.'PINPAL'.or.
     &                  species(itdindex).eq.'QUELAE'.or.
     &                  species(itdindex).eq.'DIOVIR')) then
                print *,'Unknown species ',species(itdindex),itdindex
                species(itdindex)='QUELAE'
              endif
              if (species(itdindex).eq.'PINPAL') then
                rhofllpine(ix,iy,iz)=rhofllpine(ix,iy,iz)+rhoftree
              else if (species(itdindex).eq.'QUELAE') then
                rhoftoak(ix,iy,iz)=rhoftoak(ix,iy,iz)+rhoftree
              else if (species(itdindex).eq.'DIOVIR') then
                rhofpersim(ix,iy,iz)=rhofpersim(ix,iy,iz)+rhoftree
              else
                print *,'Unknown species ',species(itdindex),itdindex
                stop
              endif
c
c  save information about this tree's contribution to the density for this
c  HIGRAD cell (note that rhoftree may be 0.0 if this cell does not really
c  contain part of the tree - as noted above, "corners" of the x/y tree
c  area may actually be empty)
c             if (rhoftree.ne.0.0) 
c    +          call rhoinfo_add(ix,iy,iz,itreenum(itdindex),rhoftree)
c
c  sum canopy mass for average later
              rmcanopysum=rmcanopysum+rmcanopy
c             print *,rmcanopysum,it,itdindex,treeht,htlctr,hb,cl,cr
              if (thinmask(ix,iy).eq.1) thinmcanopysum=thinmcanopysum+
     +                                                          rmcanopy
            enddo
          enddo
        enddo

c
c  SPECIAL eglin midstory section...
c
c  if these trees are the overstory...
        if (overstory) then
c
c  mark a 4m radius around each tree - open=0 => tree within 4m
          if (it.eq.1) open=1
          rad=4.0

          irxcell=nint(xloch(it)/dxh)
          if (irxcell.le.0) irxcell=1
          irycell=nint(yloch(it)/dyh)
          if (irycell.le.0) irycell=1
          ilow=irxcell-int(rad/dxh)
          if (ilow.lt.1) ilow=1
          ihigh=irxcell+int(rad/dxh)
          if (ihigh.gt.nxh) ihigh=nxh
          jlow=irycell-int(rad/dyh)
          if (jlow.lt.1) jlow=1
          jhigh=irycell+int(rad/dyh)
          if (jhigh.gt.nyh) jhigh=nyh
          do j=jlow,jhigh
            do i=ilow,ihigh
              cx=i*dxh
              cy=j*dyh
              if (((cx-xloch(it))**2+
     &             (cy-yloch(it))**2).le.rad**2) then
                open(i,j)=0
              endif
            enddo
          enddo
c         if (it.le.100) print *,'OPEN',it,ilow,ihigh,jlow,jhigh,
c    &          ((open(i,j),i=ilow,ihigh),j=jlow,jhigh),
c    &          xloch(it),yloch(it),rad
        endif

      enddo     ! end of "for each tree..." loop
c
c  compute canopy bulk density from the average tree height and htlc
      treehtavg=treehtavg/real(navg)
      print *,'Tree height average = ',treehtavg
      htlcavg=htlcavg/real(navg)
      avgcanht=treehtavg-htlcavg
      volcanopy=avgcanht*((itmaxxh-itminxh)*dxh)*
     +                            ((itmaxyh-itminyh)*dyh)
      bdcanopy=rmcanopysum/volcanopy
      if (ithinmask.eq.1.or.ithinmaskbd.eq.1) then
        bdthincanopy=thinmcanopysum/(avgcanht*(nmaskcells*dxh*dyh))
      endif

      return
      end

      subroutine read_fsdata
c
      use treedatinfo
c
c  read excel spreadsheet data (that has been converted to text)
c  that contains the actual Forest Service tree data from the Rocky
c  Mountain Research Station
c
c     write (*,1000) 
c1000 format(1x,'Tree data file (Excel converted to text): ',$)
c     read (*,*) treefile
c
      treefile='FStreedata.txt'
      open (1,file=treefile,form='formatted',status='old')
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      read (1,'(a)') dummy
      ntreedat=1
      do while (ios.ge.0)
        read (1,1100,iostat=ios) itreenum(ntreedat),xloc(ntreedat),
     +       yloc(ntreedat),spp(ntreedat),status(ntreedat),
     +       dbh(ntreedat),dmr(ntreedat),bole_scar(ntreedat),
     +       bss(ntreedat),height(ntreedat),htlc(ntreedat),
     +       htdc(ntreedat),htdt(ntreedat),cr1(ntreedat),
     +       cr2(ntreedat),keens(ntreedat),position(ntreedat) 
 1100   format(i6,f5.2,f5.2,a3,a6,f5.1,f5.1,a8,a3,f6.1,f4.1,f5.1,
     +         f5.1,f4.1,f4.1,a5,a1)
        print *,'rec ',ntreedat,height(ntreedat),xloc(ntreedat),
     +                 yloc(ntreedat)
        if (ios.lt.0) then
          ntreedat=ntreedat-1     ! end-of-file
        else if (.not.(xloc(ntreedat).eq.0.and.
     +                 yloc(ntreedat).eq.0)) then
          ntreedat=ntreedat+1     ! not a blank line in excel file
        endif
      enddo
      close (1)
c
      do i=1,ntreedat
        species(i)='Pipo'
        live(i)='L'
      enddo
c
      return
      end

      subroutine read_lanldata
c
      use treedatinfo
c
      character shrubfile*40,spec*4,transect*1,livedead*1
      character*8 str1,str2,str3,str4,str5
c
c  read excel spreadsheet data (that has been converted to text) that
c  contains the measured LANL plot XL06 tree data
c
c     write (*,1000)
c1000 format(1x,'XL06 tree data file (Excel -> text): ',$)
c     read (*,*) treefile
c
      treefile='xl06.trees.txt'
      open (1,file=treefile,form='formatted',status='old')
      read (1,'(a)') dummy
      ntreedat=0
      do while (ios.ge.0)
        read (1,1100,iostat=ios) isubplot,iquad,spec,livedead,str1,
     +               str2,str3,str4,str5
 1100   format(8x,8x,i1,5x,i1,7x,a4,4x,a1,7x,a8,a8,a8,a8,8x,a8)
c       read (1,1100,iostat=ios) isubplot,iquad,spec,livedead,xdbh,
c    +               xdrc,iheight,ihtlc,xcrwid
c1100   format(8x,8x,i1,5x,i1,7x,a4,4x,a1,7x,f8.1,f8.1,i8,i8,8x,f8.1)
        print *,'LANL ',isubplot,iquad,spec,livedead,
     +                  cklanlnum(str1),cklanlnum(str2),
     +                  cklanlnum(str3),cklanlnum(str4),
     +                  cklanlnum(str5)
        if (ios.lt.0.or.livedead.eq.'d'.or.livedead.eq.'D') then
          continue     ! end-of-file or dead tree (so ignore)
        else if (.not.(isubplot.eq.0.and.iquad.eq.0)) then
          ntreedat=ntreedat+1     ! not a blank line in excel file
          itreenum(ntreedat)=ntreedat
          species(ntreedat)=spec
          live(ntreedat)=livedead
          if (ntreedat.eq.1) then
            if (itdranseed.eq.0) then
              call random_seed()
            else
              call random_seed(put=itdseed)
            endif
          endif
          call random_number(randomnumber)
          xloc(ntreedat)=randomnumber*nxtd
          call random_number(randomnumber)
          yloc(ntreedat)=randomnumber*nytd
          print *,'LOCS ',xloc(ntreedat),yloc(ntreedat)
c convert inches -> meters (1 inch = .0254 meters) 
          dbh(ntreedat)=cklanlnum(str1)*.0254
c convert feet -> meters (1 foot = .3048 meters) 
          height(ntreedat)=cklanlnum(str3)*.3048
          htlc(ntreedat)=cklanlnum(str4)*.3048
          cr1(ntreedat)=(cklanlnum(str5)*.3048)/2
          if (cr1(ntreedat).eq.0) cr1(ntreedat)=.4
          cr2(ntreedat)=cr1(ntreedat)
        endif
      enddo
      close (1)
      print *,'# LANL trees = ',ntreedat
c
c  read excel spreadsheet data (that has been converted to text) that
c  contains the measured LANL plot XL06 shrub data
c
c     write (*,1200)
c1200 format(1x,'XL06 large shrub data file (Excel -> text): ',$)
c     read (*,*) shrubfile
c
      shrubfile='xl06.largeshrub.txt'
      open (1,file=shrubfile,form='formatted',status='old')
      read (1,'(a)') dummy
      ios=0
      do while (ios.ge.0)
        read (1,1300,iostat=ios) isubplot,iquad,transect,spec,
     +               livedead,str1,str2,str3,str5,str5
 1300   format(8x,8x,7x,i1,7x,i1,a1,7x,8x,a4,4x,a1,7x,a8,a8,
     +         a8,a8,a8)
        print *,'LANL SHRUB ',isubplot,iquad,transect,spec,livedead,
     +                  str1,str2,str3,str4,str5,ntreedat
        if (ios.lt.0.or.livedead.eq.'d'.or.livedead.eq.'D') then
          print *,' largeshrub ',ios,livedead
          continue     ! end-of-file or dead tree (so ignore)
        else if (.not.(isubplot.eq.0.and.iquad.eq.0)) then
         do irepeat=1,5
          ntreedat=ntreedat+1     ! not a blank line in excel file
          print *,' largeshrub 2',ntreedat
          itreenum(ntreedat)=ntreedat
          species(ntreedat)=spec
          live(ntreedat)=livedead
          if (ntreedat.eq.1) then
            if (itdranseed.eq.0) then
              call random_seed()
            else
              call random_seed(put=itdseed)
            endif
          endif
          call random_number(randomnumber)
          xloc(ntreedat)=randomnumber*nxtd
          call random_number(randomnumber)
          yloc(ntreedat)=randomnumber*nytd
          print *,'LOCS ',xloc(ntreedat),yloc(ntreedat)
c convert inches -> meters (1 inch = .0254 meters) 
          dbh(ntreedat)=cklanlnum(str1)*.0254
c convert feet -> meters (1 foot = .3048 meters) 
          height(ntreedat)=cklanlnum(str3)*.3048
          htlc(ntreedat)=cklanlnum(str4)*.3048
          cr1(ntreedat)=(cklanlnum(str5)*.3048)/2
          cr2(ntreedat)=cr1(ntreedat)
         enddo
        endif
      enddo
      close (1)
      print *,'# LANL trees and large shrubs = ',ntreedat
c
      return
      end

      function cklanlnum(string)
c
c  routine to convert a character string from the LANL xl06 data
c  to a real number
c  
      character*8 string
      real cklanlnum,xnum
c
      if (index(string,'.').ne.0) then
        read (string,'(f8.1)') cklanlnum
      else
        read (string,'(i8)') num
        cklanlnum=float(num)
      endif
c 
      return
      end

      subroutine read_caladata
c
      use treedatinfo
c
c  read excel spreadsheet data (that has been converted to text)
c  that contains the actual tree data collected by Phil Dennison 
c  for the area burned in the first hour of the Calabasas fire
c
c     write (*,1000) 
c1000 format(1x,'Tree data file (Excel converted to text): ',$)
c     read (*,*) treefile
c
c     treefile='dennison.txt'
      treefile='judy.txt'
      open (1,file=treefile,form='formatted',status='old')
      print *,' opened dennison.txt'
      read (1,'(a)') dummy
      print *,'Dum',dummy
      ntreedat=1
      do while (ios.ge.0)
        read (1,1100,iostat=ios) itreenum(ntreedat),xloc(ntreedat),
     +       yloc(ntreedat),height(ntreedat),cr1(ntreedat),
     +       cr2(ntreedat) 
c1100   format(24x,i8,40x,f11.3,f14.3,8x,f10.1,f8.3,f8.3)
 1100   format(i7,f10.3,f12.3,f8.1,f6.3,f6.3)
        print *,'rec ',ntreedat,height(ntreedat),xloc(ntreedat),
     +                 yloc(ntreedat),itreenum(ntreedat),ios
        if (ios.lt.0) then
          ntreedat=ntreedat-1     ! end-of-file
        else if (.not.(xloc(ntreedat).eq.0.and.
     +                 yloc(ntreedat).eq.0)) then
          ntreedat=ntreedat+1     ! not a blank line in excel file
        endif
      enddo
      close (1)
c
      do i=1,ntreedat
c  xloc and yloc read in UTM - convert to meters from lower left
c  corner at UTM 344370, 3779757 
        xloc(i)=xloc(i)-344370.0
        yloc(i)=yloc(i)-3779757.0
        htlc(i)=1.5
        if (htlc(i).ge.height(i)) htlc(i)=.1
        species(i)='Quag'
        live(i)='L'
      enddo
c
      return
      end

      subroutine make_chapdata
c
      use treedatinfo
c
c  # fuel types, dead and live fuel moistures
      parameter (nf=8,fdm=0.1,flm=1.2)
c
      dimension ftyp(nf),dead1(nf),dead10(nf),herb(nf)
      dimension ratiod(nf),ratiol(nf),savdead1(nf)
      dimension savherb(nf),fdepth(nf),exmoist(nf)
c
c  these fuel loads come from program prefire.f
c  ftyp values - 4 => Anderson chaparral
c                16 => NPS ceanothus
c                17 => NPS young chamise
c                18 => sagebrush/buckwheat
      data ftyp/3.,4.,9.,14.,15.,16.,17.,18./
c  dead 1-hr fuel load (tons/acre)
      data dead1/3.01,5.01,2.92,3.,2.,2.25,1.3,5.5/
c  dead 10-hr fuel load (tons/acre)
      data dead10/0.0,4.01,0.41,4.5,3.,4.8,1.0,0.8/
c  live herb. fuel load (tons/acre)
      data herb/0.,0.,0.,1.45,0.5,3.0,2.0,2.5/
c  ratio of total dead fuel
      data ratiod/3.01,11.03,3.48,8.55,6.0,8.85,3.3,6.4/ 
c  ratio of total live fuel
      data ratiol/0.0,5.0,0.0,6.45,2.5,5.8,4.0,3.25/
c  1-hr dead SAV (ft2/ft3)
      data savdead1/1500.,2000.,2500.,350.,640.,500.,640.,640./
c  live herb. SAV (ft2/ft3)
      data savherb/0.,190.,0.,1500.,2200.,1500.,2200.,1500./ 
c  fuel bed depth (feet)
      data fdepth/2.5,6.,0.2,3.,3.,6.,4.,3./
c  moist. of extinction (%) 
      data exmoist/25.,20.,25.,15.,13.,15.,20.,25./ 
c
c  from prefire.f, rhof is calculated as:
c     rhof=dead1(m)+(dead10(m)*0.5)+herb(m)    ! m is ftype
c     rhof=rhof*0.224167                       ! tons/acre to kg/m2
c     rhof=rhof/(fdepth(m)*0.3049)             ! kg/m2 to kg/m3
c     rhof=rhof/dz                             ! cell depth dz=5.
c
c  ASSUMPTIONS (primary reference a document found on the web "Chaparral
c  in Southern and Central Coastal California in the Mid-1990s: Area,
c  Ownership, Condition, and Change" by Jeremy S. Fried, Charles L.
c  Bolsinger, and Debby Beardsley of the USDA Forest Service Pacific
c  Northwest Research Station, Bulletin PNW-RB-240) - 100m x 100m area, 
c  60% shrub cover, 50% chamise and 50% ceanothus, shrub heights of 3.4-
c  6.5 ft. (average 5.0 or 5.7), and a crown diameter of 2.5 meter (page
c  I found on the web said ceanothus crassifolius was 5-10'H x 6-8'W).
c  So...  using ((100**2) * %cover) / area of bush = (10000*.6)/1.22718=
c  4889.25 (I used 4900).  Also, kg/(m**2) * 1/(%cover) => kg/(m**2) in
c  covered area (A) so average rhof=(A)/actualfuelheight.  Thus, we have
c  4900 shrubs of varying heights in the 100mx100m area, and then populate
c  them
      if (nxtd.ne.100.or.nytd.ne.100) 
     +  print *,' ERROR - for chaparral, nxtd and nytd must be 100m' 
      ntreedat=0
      ht=3.4
      do i=1,10
        do j=1,100
          ntreedat=ntreedat+1
          itreenum(ntreedat)=ntreedat
          species(ntreedat)='CHAP'
          live(ntreedat)='L'
          if (ntreedat.eq.1) then
            if (itdranseed.eq.0) then
              call random_seed()
            else
              call random_seed(put=itdseed)
            endif
          endif
          call random_number(randomnumber)
          xloc(ntreedat)=randomnumber*nxtd
          call random_number(randomnumber)
          yloc(ntreedat)=randomnumber*nytd
c convert feet -> meters (1 foot = .3048 meters)
          dbh(ntreedat)=.4*.3048
          height(ntreedat)=ht*.3048
          htlc(ntreedat)=1.0*.3048
          cr1(ntreedat)=1.25
          cr2(ntreedat)=1.25
        enddo
        ht=ht+.15
        print *,ntreedat,ht
      enddo
c
      do i=1,3900
        ntreedat=ntreedat+1
        itreenum(ntreedat)=ntreedat
        species(ntreedat)='CHAP'
        live(ntreedat)='L'
        call random_number(randomnumber)
        xloc(ntreedat)=randomnumber*nxtd
        call random_number(randomnumber)
        yloc(ntreedat)=randomnumber*nytd
c convert feet -> meters (1 foot = .3048 meters)
        call random_number(randomnumber)
        dbh(ntreedat)=.4*.3048
        ht=5.0+randomnumber
        height(ntreedat)=ht*.3048
        htlc(ntreedat)=1.0*.3048
        cr1(ntreedat)=1.25
        cr2(ntreedat)=1.25
        if (i.le.30)  print *,ntreedat,ht
      enddo
      treefile='chap'
c
      print *,'# chaparral trees = ',ntreedat
c
      end

      subroutine make_dougfir
c
      use treedatinfo
c
c  ASSUMPTIONS - 100m x 100m area,`primarily Douglas Fir, canopy closure
c  of 59%, tree heights of 40'-80', height to live crown ~15'-25', crown  
c  "spread" of ~ 15'-25', diameter at breast height (dbh) ~30" or 2.5'.
c  Note:  For this first pass, we will use 60% cover because the trees
c         are being randomly placed, so there will be overlap of the tree
c         canopies, so using a higher closure number will help to make
c         it more realistic (hopefully!).
c  So...  The area is 100m x 100m = 10,000 m^2.  60% canopy closure means
c  that 6000m^2 is covered by trees.  Tree radius is 7.5'-12.5' (use 10'), 
c  so 3.048m, so the area of a tree is pi*r^2 = 29.186m^2,  Thus, the number 
c  of trees in 6000m^2 is 6000m^2 / 29.186m^2 = ~205.58 trees.  Now, you must
c  seriously adjust this for tree shading, north/south effect, etc.

      if (nxtd.ne.100.or.nytd.ne.100)  print *,
     +    'ERROR - for Doug fir fuels, nxtd and nytd must be 100m' 
      ntreedat=0
      ht=40.0
      do i=1,10
        do j=1,56               !increase because not uping for north/south
cJLW    do j=1,20
          ntreedat=ntreedat+1
          itreenum(ntreedat)=ntreedat
          species(ntreedat)='FIR'
          live(ntreedat)='L'
          if (ntreedat.eq.1) then
            if (itdranseed.eq.0) then
              call random_seed()
            else
              call random_seed(put=itdseed)
            endif
          endif
          call random_number(randomnumber)
          xloc(ntreedat)=randomnumber*nxtd
          call random_number(randomnumber)
          yloc(ntreedat)=randomnumber*nytd
c convert feet -> meters (1 foot = .3048 meters)
          if (mod(ntreedat,2).eq.0) then
            dbh(ntreedat)=2.0*.3048
            height(ntreedat)=ht*.3048
            htlc(ntreedat)=7.5*.3048
            cr1(ntreedat)=7.5*.3048
            cr2(ntreedat)=7.5*.3048
          else
            dbh(ntreedat)=2.5*.3048
            height(ntreedat)=ht*.3048
            htlc(ntreedat)=12.5*.3048
            cr1(ntreedat)=12.5*.3048
            cr2(ntreedat)=12.5*.3048
          endif
        enddo
        ht=ht+4.0
        print *,ntreedat,ht
      enddo
      treefile='dougfir'
c
      print *,'# Douglas fir (Las Vegas/Gallinas area) trees = ',ntreedat
c
      end

      subroutine make_aspen
c
      use treedatinfo
c
c  ASSUMPTIONS - 100m x 100m area,`primarily Douglas Fir, canopy closure
c  of 59%, tree heights of 40'-80', height to live crown ~6'-8', crown
c  "spread" of ~ 15'-25', diameter at breast height (dbh) ~30" or 2.5'.
c  Note:  For this first pass, we will use 60% cover because the trees
c         are being randomly placed, so there will be overlap of the tree
c         canopies, so using a higher closure number will help to make
c         it more realistic (hopefully!).
c  So...  The area is 100m x 100m = 10,000 m^2.  60% canopy closure means
c  that 6000m^2 is covered by trees.  Tree radius is 7.5'-12.5' (use 10'),
c  so 3.048m, so the area of a tree is pi*r^2 = 29.186m^2,  Thus, the number
c  of trees in 6000m^2 is 6000m^2 / 29.186m^2 = ~205.58 trees (I used 200).

      if (nxtd.ne.100.or.nytd.ne.100)  print *,
     +    'ERROR - for aspen fuels, nxtd and nytd must be 100m'
      ntreedat=0
      ht=40.0
      do i=1,10
        do j=1,20
          ntreedat=ntreedat+1
          itreenum(ntreedat)=ntreedat
          species(ntreedat)='FIR'
          live(ntreedat)='L'
          if (ntreedat.eq.1) then
            if (itdranseed.eq.0) then
              call random_seed()
            else
              call random_seed(put=itdseed)
            endif
          endif
          call random_number(randomnumber)
          xloc(ntreedat)=randomnumber*nxtd
          call random_number(randomnumber)
          yloc(ntreedat)=randomnumber*nytd
c convert feet -> meters (1 foot = .3048 meters)
          if (mod(ntreedat,2).eq.0) then
            dbh(ntreedat)=2.0*.3048
            height(ntreedat)=ht*.3048
            htlc(ntreedat)=6.0*.3048
            cr1(ntreedat)=7.5*.3048
            cr2(ntreedat)=7.5*.3048
          else
            dbh(ntreedat)=2.5*.3048
            height(ntreedat)=ht*.3048
            htlc(ntreedat)=8.0*.3048
            cr1(ntreedat)=12.5*.3048
            cr2(ntreedat)=12.5*.3048
          endif
        enddo
        ht=ht+4.0
        print *,ntreedat,ht
      enddo
      treefile='aspen'
c
      print *,'# aspen (Las Vegas/Gallinas area) trees = ',ntreedat
c
      end

      subroutine read_eglin
c
      use treedatinfo

      integer, dimension(4) :: imseed = (/3,7,5,9/)
c
c  Read overstory "tree" files from Eglin Air Force folks - 106.2m x 61m area
c
c     write (*,1000)
c1000 format(1x,'Tree data file: ',$)
c     read (*,*) treefile
      treefile='eglin0107U001.csv'

      open (1,file=treefile,form='formatted',status='old')
      ntreedat=1
      read (1,*,iostat=ios) iskip
      read (1,*,iostat=ios) iskip
      read (1,*,iostat=ios) iskip
      read (1,*,iostat=ios) iskip
      read (1,*,iostat=ios) iskip
      do while (ios.ge.0)
        read (1,*,iostat=ios) itreenum(ntreedat),xloc(ntreedat),
     +       yloc(ntreedat),species(ntreedat),dummy,dummy,
     +       dbh(ntreedat),height(ntreedat),htlc(ntreedat),
     +       cr1(ntreedat)
c       print *,xloc(ntreedat),yloc(ntreedat),
c    +       species(ntreedat),dbh(ntreedat),height(ntreedat)
        if (ios.lt.0) then
          ntreedat=ntreedat-1     ! end-of-file

c***  SPECIAL to fix eglin data for trees with 0 crown radius  ***********
        else if (cr1(ntreedat).eq.0.0) then
          print *,'Tree w/ crown radius=0 at ',ntreedat

        else if (.not.(xloc(ntreedat).eq.0.and.
     +                 yloc(ntreedat).eq.0)) then
          ntreedat=ntreedat+1     ! not a blank line in excel file
        endif
      enddo

      if (ntreedat.eq.0) ntreedat=1
      close (1)
c
c  randomly place trees within the 60m x 106m area (nxtd x nytd)
      call random_seed(put=imseed)
      do i=1,ntreedat
        itreenum(i)=i    ! cause skipped data with cr1=0.0
        call random_number(randomnumber)
        xloc(i)=randomnumber*nxtd
        call random_number(randomnumber)
        yloc(i)=randomnumber*nytd
        cr2(i)=cr1(i)
        print *,i,xloc(i),yloc(i),species(i),dbh(i),height(i),
     +          cr1(i),nxtd,nytd
      enddo
      print *,'read_eglin - # trees = ',ntreedat

      return
      end

      subroutine eglin_midstory_wellmanaged
c
c  NOTE:  This routine will over-write the overstory data!!
c
c  Create the "well managed" eglin midstory (from the Xcel file 
c                                Eglin 0107U001 data_midstory.xlsx) -
c    165 PINPAL (longleaf pine) per acre
c    223 DIOVIR (common persimmon) per acre  (unless "dry" (ieglindry=1), 
c        then 0 DIOVIR - all are killed)
c    18 QUELAE (turkey oak) per acre
c  nxh x nyh grid = nint(nxh*dxh*nyh*dyh m2 * 0.000247105) = # acres
c  So, first place (# acres)*165 PINPAL in gaps >4m from a tree 
c  ("Longleaf pine regenerates in uneven-aged clumps in canopy gaps.")
c  Then place (# acres)*223 DIOVIR (if not "dry") and (# acres)*18 
c  QUELAE randomly.
c
      use treedatinfo
      use higradinfo
      use constants

      integer, dimension(4) :: imseed = (/3,7,5,9/)
      
      cvtm2toacres=2.47105e-04
      npinpal=165
      ndiovir=223
      nquelae=18

      ntreedat=0

      iopen=0
      do j=1,nyh
        do i=1,nxh
          if (open(i,j).eq.1) iopen=iopen+1
        enddo
      enddo
      print *,'Open cells for midstory placement = ',iopen

      call random_seed(put=imseed)
      nxt=nxh*dxh
      nyt=nyh*dxh
c
c  place PINPAL midstory in "open" locations (canopy gaps) - 5 different
c  types of PINPAL bushes (from Xcel "midstory" spreadsheet)
      nbush=nint(nxh*dxh*nyh*dyh*cvtm2toacres*npinpal/5.0)
      print *,'# PINPAL = ',nbush*5,nint(nxh*dxh*nyh*dyh*cvtm2toacres)
      do i=1,5
       do itt=1,nbush
        ic=0
        do while (ic.le.1000)
          call random_number(randomnumber)
          xl=randomnumber*nxt
          ixcell=nint(xl/dxh)
          call random_number(randomnumber)
          yl=randomnumber*nyt
          iycell=nint(yl/dyh)
          if (xl.ge.1.0.and.xl.le.float(nxt).and.
     &        yl.ge.1.0.and.yl.le.float(nyt).and.
     &        open(ixcell,iycell).eq.1) then
c           print *,'Found open location ',ixcell,iycell,ic 
            exit
          else
            ic=ic+1
            if (ic.gt.1000) then
              print *,'Can not find open location for PINPAL'
              stop
            endif
          endif
        enddo
        ntreedat=ntreedat+1
        itreenum(ntreedat)=ntreedat
        xloc(ntreedat)=xl
        yloc(ntreedat)=yl
        species(ntreedat)='PINPAL'
        if (i.eq.1) then
          dbh(ntreedat)=2.0
          height(ntreedat)=1.5
          htlc(ntreedat)=1.0
          cr1(ntreedat)=0.5
          cr2(ntreedat)=0.5
        else if (i.eq.2.or.i.eq.3) then
          dbh(ntreedat)=4.0
          height(ntreedat)=2.5
          htlc(ntreedat)=1.5
          cr1(ntreedat)=1.5
          cr2(ntreedat)=1.5
        else
          dbh(ntreedat)=6.0
          height(ntreedat)=3.5
          htlc(ntreedat)=2.0
          cr1(ntreedat)=2.0
          cr2(ntreedat)=2.0
        endif
       enddo
      enddo
c
c  if not the "dry" season" case, place DIOVIR midstory in random 
c  locations - 2 different types of DIOVIR bushes (from Xcel "midstory"
c  spreadsheet)
      if (ieglindry.eq.0) then
        nbush=nint(nxh*dxh*nyh*dyh*cvtm2toacres*ndiovir/2.0)
        print *,'# DIOVIR = ',nbush*2
        do i=1,2
          do itt=1,nbush
            ntreedat=ntreedat+1
            call random_number(randomnumber)
            xloc(ntreedat)=randomnumber*nxt
            call random_number(randomnumber)
            yloc(ntreedat)=randomnumber*nyt
            itreenum(ntreedat)=ntreedat
            species(ntreedat)='DIOVIR'
            if (i.eq.1) then
              dbh(ntreedat)=1.5
              height(ntreedat)=1.5
              htlc(ntreedat)=1.0
              cr1(ntreedat)=1.0
              cr2(ntreedat)=1.0
            else
              dbh(ntreedat)=2.0
              height(ntreedat)=2.0
              htlc(ntreedat)=1.0
              cr1(ntreedat)=1.0
              cr2(ntreedat)=1.0
            endif
          enddo
        enddo
      else
        print *,'# DIOVIR = ',0
      endif
c
c  place QUELAE midstory in random locations (from Xcel "midstory"
c  spreadsheet)
      nbush=nint(nxh*dxh*nyh*dyh*cvtm2toacres*nquelae)
      print *,'# QUELAE = ',nbush
      do i=1,nbush
        ntreedat=ntreedat+1
        call random_number(randomnumber)
        xloc(ntreedat)=randomnumber*nxt
        call random_number(randomnumber)
        yloc(ntreedat)=randomnumber*nyt
        itreenum(ntreedat)=ntreedat
        species(ntreedat)='QUELAE'
        dbh(ntreedat)=2.0
        height(ntreedat)=2.0
        htlc(ntreedat)=1.0
        cr1(ntreedat)=2.0
        cr2(ntreedat)=2.0
      enddo

      print *,'# eglin midstory (well-managed) shrubs = ',ntreedat

      return
      end

      subroutine eglin_midstory_unmanaged
c
c  NOTE:  This routine will over-write the overstory data!!
c
c  Create the "unmanaged" eglin midstory (from the Xcel file 
c               Eglin 0107U001 data_midstory_Modified_07_18_12.xlsx) -
c    165 PINPAL (longleaf pine) per acre - 3 sizes
c    223 DIOVIR (common persimmon) per acre  (unless "dry" (ieglindry=1), 
c        then 0 DIOVIR - all are killed) - 2 sizes
c    398 QUELAE (turkey oak) per acre - 3 sizes
c  nxh x nyh grid = nint(nxh*dxh*nyh*dyh m2 * 0.000247105) = # acres
c  So, first place (# acres)*165 PINPAL in gaps >4m from a tree 
c  ("Longleaf pine regenerates in uneven-aged clumps in canopy gaps.")
c  Then place (# acres)*223 DIOVIR (if not "dry") and (# acres)*398 
c  QUELAE randomly.
c  
c  NOTE: Data specifies "density" of different sizes of each species, 
c        thus the nbush below.
c
      use treedatinfo
      use higradinfo
      use constants

      integer, dimension(4) :: imseed = (/3,7,5,9/)
      
      cvtm2toacres=2.47105e-04
      nacres=nint(nxh*dxh*nyh*dyh*cvtm2toacres)
      npinpal=165
      ndiovir=223
      nquelae=398

      ntreedat=0

      iopen=0
      do j=1,nyh
        do i=1,nxh
          if (open(i,j).eq.1) iopen=iopen+1
        enddo
      enddo
      print *,'Open cells for midstory placement = ',iopen

      call random_seed(put=imseed)
      nxt=nxh*dxh
      nyt=nyh*dxh
c
c  place PINPAL midstory in "open" locations (canopy gaps) - 3 different
c  types of PINPAL bushes (from Xcel "midstory" spreadsheet)
      print *,'# PINPAL = ',npinpal*nacres
      do i=1,3
       if (i.eq.1) then
         nbush=100*nacres
       else if (i.eq.2) then
         nbush=40*nacres
       else if (i.eq.3) then
         nbush=25*nacres
       endif
       do itt=1,nbush
        ic=0
        do while (ic.le.1000)
          call random_number(randomnumber)
          xl=randomnumber*nxt
          ixcell=nint(xl/dxh)
          call random_number(randomnumber)
          yl=randomnumber*nyt
          iycell=nint(yl/dyh)
          if (xl.ge.1.0.and.xl.le.float(nxt).and.
     &        yl.ge.1.0.and.yl.le.float(nyt).and.
     &        open(ixcell,iycell).eq.1) then
c           print *,'Found open location ',ixcell,iycell,ic 
            exit
          else
            ic=ic+1
            if (ic.gt.1000) then
              print *,'Can not find open location for PINPAL'
              stop
            endif
          endif
        enddo
        ntreedat=ntreedat+1
        itreenum(ntreedat)=ntreedat
        xloc(ntreedat)=xl
        yloc(ntreedat)=yl
        species(ntreedat)='PINPAL'
        if (i.eq.1) then
          dbh(ntreedat)=2.0
          height(ntreedat)=1.5
          htlc(ntreedat)=1.0
          cr1(ntreedat)=0.5
          cr2(ntreedat)=0.5
        else if (i.eq.2) then
          dbh(ntreedat)=4.0
          height(ntreedat)=2.5
          htlc(ntreedat)=1.5
          cr1(ntreedat)=1.5
          cr2(ntreedat)=1.5
        else
          dbh(ntreedat)=6.0
          height(ntreedat)=3.5
          htlc(ntreedat)=2.0
          cr1(ntreedat)=2.0
          cr2(ntreedat)=2.0
        endif
       enddo
      enddo
c
c  if not the "dry" season" case, place DIOVIR midstory in random 
c  locations - 2 different types of DIOVIR bushes (from Xcel "midstory"
c  spreadsheet)
      if (ieglindry.eq.0) then
        print *,'# DIOVIR = ',ndiovir*nacres
        do i=1,2
          if (i.eq.1) then
            nbush=200*nacres
          else if (i.eq.2) then
            nbush=23*nacres
          endif
          do itt=1,nbush
            ntreedat=ntreedat+1
            call random_number(randomnumber)
            xloc(ntreedat)=randomnumber*nxt
            call random_number(randomnumber)
            yloc(ntreedat)=randomnumber*nyt
            itreenum(ntreedat)=ntreedat
            species(ntreedat)='DIOVIR'
            if (i.eq.1) then
              dbh(ntreedat)=1.5
              height(ntreedat)=1.5
              htlc(ntreedat)=1.0
              cr1(ntreedat)=1.0
              cr2(ntreedat)=1.0
            else
              dbh(ntreedat)=2.0
              height(ntreedat)=2.0
              htlc(ntreedat)=1.0
              cr1(ntreedat)=1.0
              cr2(ntreedat)=1.0
            endif
          enddo
        enddo
      else
        print *,'# DIOVIR = ',0
      endif
c
c  place QUELAE midstory in random locations - 3 different types
c  of QUELAE bushes (from Xcel "midstory" spreadsheet)
      print *,'# QUELAE = ',nquelae*nacres
      do i=1,3
       if (i.eq.1) then
         nbush=304*nacres
       else if (i.eq.2) then
         nbush=65*nacres
       else if (i.eq.3) then
         nbush=29*nacres
       endif
       do itt=1,nbush
         ntreedat=ntreedat+1
         call random_number(randomnumber)
         xloc(ntreedat)=randomnumber*nxt
         call random_number(randomnumber)
         yloc(ntreedat)=randomnumber*nyt
         itreenum(ntreedat)=ntreedat
         species(ntreedat)='QUELAE'
         if (i.eq.1) then
           dbh(ntreedat)=2.0
           height(ntreedat)=2.0
           htlc(ntreedat)=1.0
           cr1(ntreedat)=2.0
           cr2(ntreedat)=2.0
         else if (i.eq.2) then
           dbh(ntreedat)=6.0
           height(ntreedat)=3.5
           htlc(ntreedat)=1.0
           cr1(ntreedat)=1.0
           cr2(ntreedat)=1.0
         else
           dbh(ntreedat)=15.0
           height(ntreedat)=4.5
           htlc(ntreedat)=2.0
           cr1(ntreedat)=2.0
           cr2(ntreedat)=2.0
         endif
       enddo
      enddo

      print *,'# eglin midstory shrubs = ',ntreedat

      return
      end

      subroutine populate
c
      use treedatinfo
      use higradinfo
c
      real nsslopefactor,slopeeffectval
c
c  indicate if rotating tree data was requested
      if (itdrotate.eq.1) then
        nxt=nytd
        nyt=nxtd
      else
        nxt=nxtd
        nyt=nytd
      endif
c
c  populate an area the size of the HIGRAD grid (or less if the x,y 
c  min and max indicate that areas of the grid are to be left empty)
c  with replicas of the smaller measured/sampled tree "grid" (really
c  tree locations, etc), either exactly replicating or spreading
c  the trees randomly within each area the size of the measured tree
c  "grid"
      xmeterh=(itmaxxh-itminxh+1)*dxh
      nxloop=int(xmeterh/nxt)
      if (mod(xmeterh,real(nxt)).ne.0) nxloop=nxloop+1
      ymeterh=(itmaxyh-itminyh+1)*dyh
      nyloop=int(ymeterh/nyt)
      if (mod(ymeterh,real(nyt)).ne.0) nyloop=nyloop+1
      ixadd=(itminxh-1)*dxh
      ixmaxadd=itmaxxh*dxh
      iyadd=(itminyh-1)*dyh
      iymaxadd=itmaxyh*dyh
      inewtr=0

      if (itdloc.eq.0) then 
        if (itdranseed.eq.0) then
          call random_seed()
        else
          call random_seed(put=itdseed)
        endif
      endif
      print *,'before loop ',xmeterh,nxloop,ymeterh,nyloop,ixadd,
     +                       ixmaxadd,iyadd,iymaxadd
      do j=1,nyloop
        do i=1,nxloop
          do it=1,ntreedat
            inewtr=inewtr+1
            if (inewtr.gt.ntrh) then
              print *,' Error - too many trees in HIGRAD grid'
              stop
            endif
            if (itdloc.eq.0) then
              call random_number(randomnumber)
              xl=randomnumber*nxt+ixadd
              call random_number(randomnumber)
              yl=randomnumber*nyt+iyadd
            else if (itdrotate.eq.0) then
              xl=xloc(it)+ixadd
              yl=yloc(it)+iyadd
            else
              xl=yloc(it)+ixadd
              yl=xloc(it)+iyadd
            endif
            if (xl.lt.ixmaxadd.and.yl.lt.iymaxadd) then
              xloch(inewtr)=xl
              yloch(inewtr)=yl
              itreeh(inewtr)=itreenum(it)
            else
              inewtr=inewtr-1
            endif
          enddo
          ixadd=ixadd+nxt
        enddo
        ixadd=(itminxh-1)*dxh
        iyadd=iyadd+nyt
      enddo
c
c  if requested, consider northern/southern exposure and reduce the
c  number of trees on the south face by nsslopefactor
      if (treeslopefactor.gt.0.and.topofile.ne.'flat') then
        if (gridorientation.eq.0) then
          northx=1
          northy=0
        else
          northx=0
          northy=1
        endif
        slopeeffectval=slopeeffect/100.0
        icurrent=0
        do it=1,inewtr
          ix=nint(xloch(it)/dxh)
          iy=nint(yloch(it)/dyh)
c  leave the tree if it is right on the edge of the HIGRAD grid
          if (.not.(ix.eq.1.or.ix.eq.nxh.or.
     +              iy.eq.1.or.iy.eq.nyh)) then 
            xslope=(zs(ix+1,iy)-zs(ix-1,iy))/(2.0*dxh)
            yslope=(zs(ix,iy+1)-zs(ix,iy-1))/(2.0*dyh)
            value=1.0/sqrt(xslope**2+yslope**2)
            nsslopefactor= -slopeeffectval*((northx*xslope*value)+
     +                                   (northy*yslope*value))
         if (it.eq.48) print *,'JLW48 ',-slopeeffectval,
     +   northx,xslope,value,northy,yslope,value,nsslopefactor
            call random_number(randomnumber)
            retainpercent=(1.+nsslopefactor)/(1.0+slopeeffectval)
            if (randomnumber.ge.retainpercent) itreeh(it)=-999
            if (it.le.100) print *,it,xloch(it),yloch(it),xslope,
     +                  yslope,nsslopefactor,retainpercent,
     +                  randomnumber,slopeeffectval,value,ix,iy,
     +                  zs(ix+1,iy),zs(ix-1,iy),zs(ix,iy+1),
     +                  zs(ix,iy-1)
          endif
          if (itreeh(it).ne.-999) then
            icurrent=icurrent+1
            xloch(icurrent)=xloch(it)
            yloch(icurrent)=yloch(it)
            itreeh(icurrent)=itreeh(it)
          endif
        enddo 
        print *,'Northern/southern exposure code - ',inewtr,
     +          ' trees to ',icurrent
        inewtr=icurrent
      endif
c
      return
      end

      subroutine readoldtrees
c
      use treedatinfo
      use higradinfo
c
c  temp arrays and variables for reading tree data
        integer itreenumtmp(ntr),ntreedattmp,nxtdtmp,nytdtmp
        real dbhtmp(ntr),dmrtmp(ntr),heighttmp(ntr),htlctmp(ntr),
     +       htdctmp(ntr),htdttmp(ntr),cr1tmp(ntr),cr2tmp(ntr)
        character speciestmp(ntr)*4,livetmp(ntr)*1,spptmp(ntr)*3,
     +            statustmp(ntr)*6,bole_scartmp(ntr)*8,
     +            bsstmp(ntr)*3,keenstmp(ntr)*5,positiontmp(ntr)*1
        character treefiletmp*40
        character answer*5

c
c  read part of a tree file generated by a previous run of this
c  program to populate the HIGRAD grid, making sure that the Forest
c  Service file is the same
      write (*,1000) 
 1000 format(1x,'Old tree data file: ',$)
      read (*,*) oldtreefile
c
      open (unit=1,file=oldtreefile,form='unformatted',status='old')
c
c  read sampled tree data information, making sure that it matches
c  the data that has already been read this run of this program 
      idiff=0
      read (1) nxtdtmp,nytdtmp,treefiletmp,ntreedattmp
      do i=1,ntreedattmp
        read (1) itreenumtmp(i),spptmp(i),statustmp(i),dbhtmp(i),
     +            dmrtmp(i),bole_scartmp(i),bsstmp(i),heighttmp(i),
     +            htlctmp(i),htdctmp(i),htdttmp(i),cr1tmp(i),cr2tmp(i),
     +            keenstmp(i),positiontmp(i),speciestmp(i),livetmp(i)
        if (itreenumtmp(i).ne.itreenum(i)) idiff=idiff+1
        if (spptmp(i).ne.spp(i)) idiff=idiff+1
        if (statustmp(i).ne.status(i)) idiff=idiff+1
        if (dbhtmp(i).ne.dbh(i)) idiff=idiff+1
        if (dmrtmp(i).ne.dmr(i)) idiff=idiff+1
        if (bole_scartmp(i).ne.bole_scar(i)) idiff=idiff+1
        if (bsstmp(i).ne.bss(i)) idiff=idiff+1
        if (heighttmp(i).ne.height(i)) idiff=idiff+1
        if (htlctmp(i).ne.htlc(i)) idiff=idiff+1
        if (htdctmp(i).ne.htdc(i)) idiff=idiff+1
        if (htdttmp(i).ne.htdt(i)) idiff=idiff+1
        if (cr1tmp(i).ne.cr1(i)) idiff=idiff+1
        if (cr2tmp(i).ne.cr2(i)) idiff=idiff+1
        if (keenstmp(i).ne.keens(i)) idiff=idiff+1
        if (positiontmp(i).ne.position(i)) idiff=idiff+1
        if (speciestmp(i).ne.species(i)) idiff=idiff+1
        if (livetmp(i).ne.live(i)) idiff=idiff+1
      enddo
      if (nxtdtmp.ne.nxtd.or.nytdtmp.ne.nytd.or.treefiletmp.ne.
     +    treefile.or.ntreedattmp.ne.ntreedat.or.idiff.ne.0) then
        print *,' Differences in sampled tree data header info - ',
     +          'new/old: ',treefiletmp,treefile,nxtdtmp,nxtd,
     +          nytdtmp,nytd,ntreedattmp,ntreedat
        print *,' Number of differences in Forest Service data = ',
     +          idiff
        write (*,1100)
 1100   format(1x,'Continue? (Y or N): ',$)
        read (*,*) answer
        if (answer(1:1).eq.'n'.or.answer(1:1).eq.'N') stop
      endif
c
c  read HIGRAD grid info and then info about trees as they were populated
c  into the HIGRAD grid
      read (1) nxhtmp,nyhtmp,nzhtmp,dxhtmp,dyhtmp,dzhtmp,inewtr
      if (nxhtmp.ne.nxh.or.nyhtmp.ne.nyh.or.nzhtmp.ne.nzh.or.
     +    dxhtmp.ne.dxh.or.dyhtmp.ne.dyh.or.dzhtmp.ne.dzh) then
        print *,' Differences in HIGRAD grid sizes - new/old: ',
     +    nxhtmp,nxh,nyhtmp,nyh,nzhtmp,nzh,dxhtmp,dxh,dyhtmp,
     +    dyh,dzhtmp,dzh
        stop
      endif
      do i=1,inewtr
        read (1) xloch(i),yloch(i),itreeh(i)
      enddo
c
      print *,' Tree locations read from file ',oldtreefile
c
      return
      end

      subroutine thintrees
c
      use treedatinfo
      use higradinfo
      use thininfo
c
      character (len=40) :: thinmaskfile 
c
c  if a mask showing which areas are to be thinned is being used...
      if (ithinmask.eq.1.or.ithinmaskbd.eq.1) then
        write (*,'(a28,$)') 'Name of thinning mask file: '
        read (*,*) thinmaskfile
        open (1,file=thinmaskfile,form='unformatted',status='old')
        read (1) thinmask
        close (1)
        nmaskcells=0
        do j=1,nyh
          do i=1,nxh
            if (thinmask(i,j).eq.1) nmaskcells=nmaskcells+1
          enddo
        enddo
        if (ithinmethod.eq.0) return
        iminxh=itminxh
        imaxxh=itmaxxh
        iminyh=itminyh
        imaxyh=itmaxyh
c
c  otherwise, ensure that the area to be thinned has trees
      else 
        thinmask=0
        if (ithinminxh.lt.itminxh.or.ithinmaxxh.gt.itmaxxh.or.
     +      ithinminyh.lt.itminyh.or.ithinmaxyh.gt.itmaxyh) 
     +         print *,' Possible error - area to be thinned is ',
     +          'larger than area with trees'
        iminxh=ithinminxh
        imaxxh=ithinmaxxh
        iminyh=ithinminyh
        imaxyh=ithinmaxyh
      endif
c
      xlochtmp=xloch
      ylochtmp=yloch
      itreehtmp=itreeh
      ioldnewtr=inewtr
      inewtr=0
c
c  initialize for thinning to patches
      if (ithinmethod.ge.2) then
        ipatchesx=int(real(imaxxh-iminxh)*dxh/patchperiod)+1
        ipatchesy=int(real(imaxyh-iminyh)*dyh/patchperiod)+1
        numpatches=ipatchesx*ipatchesy
        if (numpatches.gt.ntrh) then
          write (*,*) 'Too many patches ',numbatches,
     +                ' - ntrh= ',ntrh
          stop
        endif
        if (ipatchseed.eq.0) then
          call random_seed()
        else
          call random_seed(put=ipseed)
        endif
        do ipatch=1,numpatches
          call random_number(randomnumber)
          patchx(ipatch)=real(imaxxh-iminxh)*dxh*
     +                        randomnumber+real(iminxh)*dxh
          call random_number(randomnumber)
          patchy(ipatch)=real(imaxyh-iminyh)*dyh*
     +                        randomnumber+real(iminyh)*dyh
        enddo
      endif
c
c  for each tree...
      do it=1,ioldnewtr
        icuttree=0
c
c  search Forest Service data tree numbers to find a tree number match
        if (itreehtmp(it).eq.itreenum(itreehtmp(it))) then
          itdindex=itreehtmp(it)
        else
          itdindex=1
          do while (itreehtmp(it).ne.itreenum(itdindex))
            itdindex=itdindex+1
            if (itdindex.gt.ntreedat) then
              print *,' Error locating tree number ',itreehtmp(it)
              stop
            endif
          enddo
        endif
c
c  if tree is within area to be thinned...
        ix=nint(xloch(it)/dxh)
        iy=nint(yloch(it)/dyh)
        if ((ithinmask.eq.0.and.ix.ge.iminxh.and.ix.le.imaxxh.and.
     +                            iy.ge.iminyh.and.iy.le.imaxyh).or.
     +      (ithinmask.eq.1.and.thinmask(ix,iy).eq.1)) then
c
c  thin by diameter, if requested
          if ((ithinmethod.eq.1.or.ithinmethod.eq.3.).and.
     +                         dbh(itdindex).lt.dbhthin) icuttree=1
c
c  thin to uneven patches, if requested
          if ((ithinmethod.eq.2.or.ithinmethod.eq.3).and.
     +                                      icuttree.eq.0) then
            icuttree=1
            do ipatch=1,numpatches
              distfrompatch=sqrt((patchx(ipatch)-xloch(it))**2
     +                         + (patchy(ipatch)-yloch(it))**2)
              if (distfrompatch.le.patchradius) icuttree=0
            enddo
          endif
          if (ithinmethod.eq.4.and.icuttree.eq.0) then
            do ipatch=1,numpatches
              distfrompatch=((patchx(ipatch)-xloch(it))/patchradius)**2
     +                   + ((patchy(ipatch)-yloch(it))/patchradiusy)**2
              if (distfrompatch.le.1) icuttree=1
            enddo
          endif
 
        endif
c
c  do not keep tree information if it has been "thinned"
        if (icuttree.eq.0) then
          inewtr=inewtr+1
          xloch(inewtr)=xlochtmp(it)
          yloch(inewtr)=ylochtmp(it)
          itreeh(inewtr)=itreehtmp(it)
        endif
c
      enddo
c
      return
      end

      subroutine writetrees
c
      use treedatinfo
      use higradinfo
      use thininfo
      use rhoinfo
c
      character newtreefile*120, date*8, tmp*5
      logical yesno
c
c  create name for the file using the name of the Forest Service data,
c  HIGRAD grid size, current date, and thinning used
      newtreefile=treefile
      ilen=len_trim(newtreefile)+1
      if (nxh.lt.10) then
        write (tmp,'(i1.1)') nxh
      else if (nxh.lt.100) then
        write (tmp,'(i2.2)') nxh
      else
        write (tmp,'(i3.3)') nxh
      endif
      newtreefile(ilen:)='_'//tmp
      ilen=len_trim(newtreefile)+1
      if (nyh.lt.10) then
        write (tmp,'(i1.1)') nyh
      else if (nyh.lt.100) then
        write (tmp,'(i2.2)') nyh
      else
        write (tmp,'(i3.3)') nyh
      endif
      newtreefile(ilen:)='X'//tmp
      ilen=len_trim(newtreefile)+1
      if (nzh.lt.10) then
        write (tmp,'(i1.1)') nzh
      else if (nzh.lt.100) then
        write (tmp,'(i2.2)') nzh
      else
        write (tmp,'(i3.3)') nzh
      endif
      newtreefile(ilen:)='X'//tmp
      ilen=len_trim(newtreefile)+1
      write (tmp,'(i1.1)') ithinmethod
      newtreefile(ilen:)='_'//tmp
      ilen=len_trim(newtreefile)+1
      call date_and_time(date)
      newtreefile(ilen:)='_'//date(5:8)//date(3:4)
      ilen=len_trim(newtreefile)
c
c  open the file if appropriate
 10   inquire (file=newtreefile,exist=yesno)
      if (yesno) then
        write (*,1000) newtreefile(1:ilen)
 1000   format(1x,'File ',a,' already exists.')
        write (*,1100)
 1100   format(1x,'Overwrite? (Y or N): ',$)
        read (*,*) tmp
        if (tmp(1:1).eq.'n'.or.tmp(1:1).eq.'N') then
          write (*,1200)
 1200     format(1x,'New file name: ',$)
          read (*,*) newtreefile
          goto 10
        endif
      endif
      open (1,file=newtreefile,form='unformatted',status='unknown')
c
c  write sampled tree data information (left out xloc and yloc because
c  it has no use)
      write (1) nxtd,nytd,treefile,ntreedat
      do i=1,ntreedat
        write (1) itreenum(i),spp(i),status(i),dbh(i),dmr(i),
     +            bole_scar(i),bss(i),height(i),htlc(i),htdc(i),
     +            htdt(i),cr1(i),cr2(i),keens(i),position(i),
     +            species(i),live(i)
      enddo
c
c  write HIGRAD grid info and then info about trees as they were populated
c  into the HIGRAD grid
      write (1) nxh,nyh,nzh,dxh,dyh,dzh,inewtr,dxs,dys,itminxh,itmaxxh,
     +          itminyh,itmaxyh,igminxh,igmaxxh,igminyh,igmaxyh,ilminxh,
     +          ilmaxxh,ilminyh,ilmaxyh
      do i=1,inewtr
        write (1) xloch(i),yloch(i),itreeh(i)
      enddo
c
c  write data used for FIRETEC run
      write (1) rhof
      write (1) actualfueldepth
      write (1) sizescale
      write (1) moist
c
c  write out info about which trees contributed to the density in each
c  HIGRAD cell
      call rhoinfo_write(1)
c
      close(1)
c
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
      parameter(n=800,m=150,l=41)
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
 
      return
      end

 
      function zcart(sigma,i,j,zs,zb)
! sigma is sigma coordinate, zcart is cartesian vertical coordinate
c  ***************
c  NOTE:  n, m, l must be the same as nxh,nyh,nzh in subroutine trees
c  ***************
      parameter(n=800,m=150,l=41)
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
      if (h.eq.0.) pause 'bad xa input in splint'
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
