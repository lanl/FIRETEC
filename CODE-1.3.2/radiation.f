! radiation was reorganized by fp (07 and 08/2012), including modification of convergence scheme, that
! was dependant on procs number, and flux accuracy that was dependent on fire total energy (and size)
! FRAZAZ stuff using Modest was removed by FP (WRONG FOR TWO REASONS: THEORETICAL AND IMPLEMENTATION) (06/2012)
! The secondargument in updated2 (ae) is removed for all updated2 calls KOO 
!
      module radiation
   
      use gridsetup 
      use msga
      use turba, only:rhomicro  
      use metryic
      use fireteca, only:min_rhof
      Implicit None

      real,allocatable::
     .     E(:,:,:),
     .     Ef(:,:,:),     !wss
     .     Es(:,:,:),     !wss
     .     EAZ(:,:,:), !FPAZ E define on the "All Zone" grid for MonteCarlo, on np*mp domain
     .     E2AZ(:,:,:), !FPAZ E define on the "All Zone" grid for MonteCarlo, on n*m domain
     .     xflux(:,:),
     .     yflux(:,:),
     .     zflux(:,:),
     .     keff(:,:,:),
     .     keffx(:,:,:),
     .     keffy(:,:,:),
     .     keffz(:,:,:),
     .     diffterm(:,:,:),
     .     hx(:,:,:),
     .     hy(:,:,:),
     .     hz(:,:,:),
     .     Ex(:,:,:),
     .     Ey(:,:,:),
     .     Ez(:,:,:),
     .     sinkgas(:,:,:),
     .     sinksol(:,:,:),
     .     gi0(:,:,:),
     .     gmul0(:),
     .     zcoordsf(:,:,:) ! for mpi_fluxes
      integer::irad ! from gridsetup, 0=no rad, 1=diffu, 2=MC
      integer::icallrad ! from firetec, frq to call rad. diffu:10, MC:1
      integer::irandseed
      integer::iseed,iseed0
      integer::iradeastflux ! flag for computation of flux to east boundary of i,j,k cell
      integer::nn,mm,LL
      real:: spval
      integer,allocatable:: lix(:),liy(:),liz(:)
      real:: cst=1.0
      real :: ckeff,gimin,gimax
      integer :: nuemannbcs
      
      !Common variable between Monte Carlo and Diffusion approximation methods
      real :: rstep
      real :: csurr,camb,crad,pathcoef,ssoot,rmaxsootcon    !,rhowoodmicro
      real :: facp,facm,frac,wdth,absemis0,papvmax,papvgasmax
      logical, allocatable :: gasonly(:,:)
      real, allocatable :: sourcesol(:,:,:),
     &                     sourcegas(:,:,:),
     &                     srcgas(:,:,:),
     &                     srcsol(:,:,:),
     &                     papv(:,:,:),
     &                     papvgas(:,:,:),
     &                     papvtot(:,:,:),
     &                     rnetsol(:,:,:),
     &                     rnetgas(:,:,:), 
     &                     papvtotAZ(:,:,:),
     &                     papvtot2AZ(:,:,:),
     &                      volumeAZ(:,:,:),
     &                     t4barsolid(:,:,:),
     &                     t4bar(:,:,:),
     &                     absemis(:,:,:),
     &                     aemit(:,:,:),
     &                      zcartedgeAZ1(:,:,:),
     &                      zcartedgeAZ2(:,:,:),
     &                      zcartedgeAZ3(:,:,:),
     &                      zcartedgeAZ4(:,:,:),
     &                     zcartedge2AZ(:,:,:,:) 
       ! on small domains : these are elevation of top of the edge between cells
       !zcartedgeAZ1(i,j,k) : (i-1,j), (i,j-1), (i-1,j-1) and (i,j) (loser left corner in plan view
       !zcartedgeAZ2(i,j,k) : (i,j), (i+1,j-1), (i+1,j-1) and (i,j)  (lower right corner in plan view)
       !zcartedgeAZ3(i,j,k) : (i-1,j+1), (i,j), (i-1,j) and (i,j+1) (upper left corner in plan view)
       !zcartedgeAZ4(i,j,k) : (i,j), (i+1,j), (i,j) and (i+1,j+1)   (upper right corner in plan view)
       ! on the whole domain , zcartedge2AZ(l,i,j,k) with l between 1 and 4 are zcartedgeAZ1(i,j,k),zcartedgeAZ2(i,j,k)...
      integer,allocatable :: nzonebelowcell1(:,:) 
! number of fuel zones which end below the top of cell k=1 (if all of the fuels in cell have actual fuel height > top of cell 1 then this value will be 0
      REAL*8, allocatable ::  percabsAZold(:,:,:),
     &                         percabsAZnew(:,:,:)

      INTEGER, allocatable :: nphotonAZ(:,:,:)
           ! number of photons absorbed in a cell on entire domain by a single proc (n*m)
      INTEGER, allocatable :: nphotonAZgather(:,:,:)
           ! number of photons absorbed in a cell on entire domain by all procs (n*m)
      INTEGER, allocatable :: nphotemisAZ(:,:,:)
                            !number photons emitted in cell (entire domain) (n*m)
      INTEGER, allocatable :: nphotemisAZsplit(:,:,:)
      INTEGER, allocatable :: nphotemisAZtemp(:,:,:)
                            !number photons emitted in cell (entire domain), split per proc
      INTEGER, allocatable :: nphotonEastFlux(:,:,:) ! to compute flux to east face of a cell
      INTEGER, allocatable :: nphotontempEastFlux(:,:,:) ! to compute flux to east face of a cell
      INTEGER nphotwrong   ! photon miscalculated due to topography
      real*8 ::Etot !,Emin
      real::Emin,papvthresh,papvthresh2 ! minimal energy for a photon
      integer::nphotbatch
      integer ::lmc ! equivalent of l for montecarlo
      INTEGER, allocatable :: nperproc(:) !number photons per processor
      integer, allocatable ::  kfuel(:,:) !highest k value with fuel
      integer, allocatable :: seed(:)
      integer :: seedsize

      integer :: llowlim,lhighlim   !*** consider using highlim to cap the MC domain instead of l
      integer :: lphase ! number of fuel family for MonteCarlo
! these to real entail to avoid very flat cells for cell 0 or 1 all zone in Montecarlo
      real :: minfueldepth ! cell 0 all zone can not be lower than 1/10 of the high the real cell in MonteCarlo
      real :: maxfueldepth1,maxfueldepth2,maxfueldepth3,maxfueldepth4
! if actual>maxfuldepth, cell 0 and 1 all zone are a mixture of gas and solid
      integer :: neumannbcs
 
      save

      contains
!!!--------------------------subroutine frad_init  -----------------------------------!!!
      subroutine frad_init()
      use metryic

      use turba !FPAZ to have actualfueldepth
      use constants
      Implicit None

      real, allocatable ::  zcartedgetmp2AZ(:,:,:)
      integer :: i,j,k,kk,cnt,cnt_rate,cnt_max
      real ::rzonesabove,zincrement1,zincrement2,zincrement3,zincrement4
      real, external :: gdeform
      real,external :: zcart
      
      IF(ibctopbot .EQ. 0) THEN
       WRITE(6,*) 'Value ibctopbot - STOP'
       STOP
      ENDIF
 
      !specify type of boundary conditions to use with the diffusion approximation radiation scheme
      neumannbcs = 0
      !Determine lsize for MC od Diffusion method
      if(irad.eq.1)then
       llowlim = nz
       lhighlim = LL
      elseif(irad.eq.2)then
       llowlim = 1
       lhighlim = l   !**** consider using this for the upper limit in the MC calc
      endif

      !Allocate MC/Diffusion common arrays
      allocate(t4barsolid(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); t4barsolid=0.0
      allocate(t4bar(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); t4bar=0.0
      allocate(sourcesol(1-ih:np+ih,1-ih:mp+ih,lhighlim)); sourcesol=0.0
      allocate(sourcegas(1-ih:np+ih,1-ih:mp+ih,lhighlim)); sourcegas=0.0
      allocate(absemis(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); absemis=0.0
      allocate(aemit(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); aemit=0.0

      

      !Allocate diffusion method radiation arrays
      if(irad.eq.1)then 
       allocate(rnetsol(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); rnetsol=0.0
       allocate(rnetgas(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); rnetgas=0.0
       allocate(kfuel(1-ih:np+ih,1-ih:mp+ih)); kfuel = 0
       allocate(E(1-ih:np+ih,1-ih:mp+ih,lhighlim));  Ef=0.
       allocate(Ef(1-ih:np+ih,1-ih:mp+ih,lhighlim)); Ef=0.
       allocate(Es(1-ih:np+ih,1-ih:mp+ih,lhighlim)); Es=0.
       allocate(xflux(1-ih:n+ih,1-ih:m+ih));     xflux=0.
       allocate(yflux(1-ih:n+ih,1-ih:m+ih));     yflux=0.
       allocate(zflux(1-ih:n+ih,1-ih:m+ih));     zflux=0.
       allocate(lix(1-ih:np+ih));     lix=0.
       allocate(liy(1-ih:mp+ih));     liy=0.
       allocate(liz(lhighlim));     liz=0.
       allocate(diffterm(1-ih:np+ih,1-ih:mp+ih,lhighlim));  diffterm=0.
       allocate(keff(1-ih:np+ih,1-ih:mp+ih,lhighlim));     keff=0.
       allocate(keffx(1-ih:np+ih,1-ih:mp+ih,lhighlim));     keffx=0.
       allocate(keffy(1-ih:np+ih,1-ih:mp+ih,lhighlim));     keffy=0.
       allocate(keffz(1-ih:np+ih,1-ih:mp+ih,lhighlim));     keffz=0.
       allocate(hx(1-ih:np+ih,1-ih:mp+ih,lhighlim));     hx=0.
       allocate(hy(1-ih:np+ih,1-ih:mp+ih,lhighlim));     hy=0.
       allocate(hz(1-ih:np+ih,1-ih:mp+ih,lhighlim+1));     hz=0.
       allocate(Ex(1-ih:np+ih,1-ih:mp+ih,lhighlim))
       allocate(Ey(1-ih:np+ih,1-ih:mp+ih,lhighlim))
       allocate(Ez(1-ih:np+ih,1-ih:mp+ih,lhighlim+1))
       allocate(sinksol(1-ih:np+ih,1-ih:mp+ih,lhighlim)); sinksol=0.
       allocate(sinkgas(1-ih:np+ih,1-ih:mp+ih,lhighlim)); sinkgas=0.
       allocate(gi0(1-ih:np+ih,1-ih:mp+ih,lhighlim));  gi0 = 0.0
       allocate(gmul0(lhighlim));                       gmul0=0.0
       allocate(papv(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); papv = 0.0
       allocate(papvgas(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); papvgas = 0.0
       allocate(papvtot(1-ih:np+ih,1-ih:mp+ih,llowlim:lhighlim)); papvtot = 0.0
       allocate(srcsol(1-ih:np+ih,1-ih:mp+ih,lhighlim)); srcsol=0.0
       allocate(srcgas(1-ih:np+ih,1-ih:mp+ih,lhighlim)); srcgas=0.0
       allocate(zcoordsf(1-ih:np+ih,1-ih:mp+ih,L+1));     zcoordsf=0.
      endif !irad.eq.1

      !Allocate Monte Carlo  method radiation arrays
      if(irad.eq.2)then 
       lphase=1  ! (lphase= number of solid phase in the first layer) FPAZ
       ! if lphase=1, on the all array finishing by AZ (for All Zones),
       ! arrayAZ(:,:,1)=value of array in the solid phase (below actualfueldepth) : solid+gas
       ! arrayAZ(:,:,2)=value of array in the gas phase (above actualfueldepth) : gas only
       ! for k=3 to l+lphase, arrayAZ(:,:,k)=array(:,:,k-1) : solid+gas
       
       allocate(papv(np,mp,lhighlim)); papv = 0.0  !projected area per unit volume of the solid
       allocate(papvgas(np,mp,lhighlim)); papvgas = 0.0 !projected area per unit volume of the gas
       allocate(papvtot(np,mp,lhighlim)); papvtot = 0.0 !total gas and solid emitting of abosrbing area within a cell or zone
       allocate(papvtotAZ(np,mp,llowlim:lhighlim+lphase)); papvtotAZ = 0.0 !FPAZ
       allocate(papvtot2AZ(n,m,llowlim:lhighlim+lphase)); papvtot2AZ = 0.0 !FPAZ
       allocate(volumeAZ(np,mp,llowlim:lhighlim+lphase)); volumeAZ = 0.0 !FPAZ volume of cell
       ! heights of the edges on the All Zone grid
       allocate(zcartedgeAZ1(np,mp,llowlim-1:lhighlim+lphase)); zcartedgeAZ1 = 0.0 !FPAZ
       allocate(zcartedgeAZ2(np,mp,llowlim-1:lhighlim+lphase)); zcartedgeAZ2 = 0.0 !FPAZ
       allocate(zcartedgeAZ3(np,mp,llowlim-1:lhighlim+lphase)); zcartedgeAZ3 = 0.0 !FPAZ
       allocate(zcartedgeAZ4(np,mp,llowlim-1:lhighlim+lphase)); zcartedgeAZ4 = 0.0 !FPAZ
       allocate(zcartedge2AZ(4,n,m,llowlim-1:lhighlim+lphase)); zcartedge2AZ = 0.0 !FPAZ
       allocate(zcartedgetmp2AZ(n,m,llowlim-1:lhighlim+lphase)); zcartedgetmp2AZ= 0.0 !FPAZ
       allocate(nphotemisAZ(n,m,l+lphase))
       allocate(nphotemisAZsplit(n,m,l+lphase))
       allocate(nphotemisAZtemp(n,m,l+lphase))
       allocate(nphotonAZgather(n,m,l+lphase))
       allocate(nphotonAZ(n,m,l+lphase))
       if (iradeastflux.eq.1) then
           allocate(nphotonEastFlux(n,m,l+lphase))
           allocate(nphotontempEastFlux(n,m,l+lphase))
           nphotonEastFlux=0
           nphotontempEastFlux=0
       endif
       allocate(EAZ(np,mp,l+lphase))
       allocate(E2AZ(n,m,l+lphase))
       allocate(percabsAZold(n,m,l+lphase))
       allocate(percabsAZnew(n,m,l+lphase))
       allocate(nperproc(0:nproc-1))
       allocate(seed(numprocs))
       allocate(nzonebelowcell1(np,mp))
       !TODO : get a value for Emin
       Emin=1000.0
       papvthresh = 2.0*min_rhof/(ss*rhomicrovalue)
       papvthresh2 = 2.0*0.01/(ss*rhomicrovalue)
      endif !irad.eq.2

      !Set constants
      pathcoef=1.0   !coefficient for pathlength?
      ! now defined in gridlist
      !if (iinra.eq.0) then
      !   crad=20.0      !coefficient for papv calculations in radiation 
      !else
      !   crad=50.0      !coefficient for papv calculations in radiation 
      !endif
      csurr=.8       !surrounding weighting coeffient for gaussian subgrid temperature distribution
      camb=.2        !local weighting coeffient for gaussian subgrid temperature distribution
      facp=.5*(3.+2.449489743)
      facm=.5*(3.-2.449489743)
      rmaxsootcon=.2 !maximum soot concentration
      ssoot=1.e-5    !soot sizescale?
      absemis0=0.3    !


      !Setup one-time calculated arrays
      ckeff=1000.0

      
      if(irad.eq.1)then
       do k=1,LL
       do j=1,mp
       do i=1,np
         kk=k-nz
         if(kk.lt.nz)kk=nz
         if(kk.gt.l)kk=l
         gmul0(k)=gmul(kk)
         gi0(i,j,k)=gi(i,j,kk)
      enddo
      enddo
      enddo
       call updated(gi0,gi0,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1,0)
       !the next loop is taken from rinit because used here only i(FP)
       do k=1,l+1
       do j=1,mp
        do i=1,np
         zcoordsf(i,j,k)=zcart(zedge(k),i,j)-zs(i,j)
        enddo
        enddo
       enddo
       ! FP added l+1
       call updated(zcoordsf,zcoordsf,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1,0)
       call stretch_keff()
       write(6,*) 'finished irad=1 stuff'

      elseif(irad.eq.2)then   
      !initialize the random number generator for Monte Carlo
      seedsize=numprocs
       if(mpi_rank.eq.0)then
        call system_clock(cnt,cnt_rate,cnt_max)
!        write(6,*) 'cnt=',cnt,'mpi_rank =',mpi_rank
       endif
       call MPI_Bcast(cnt,1,mpi_integer,0,mpi_comm_world,ierror)
       if(irandseed.eq.1)then
        iseed0=cnt
        iseed=iseed0
       else
        iseed0=iseed
       endif
       if(mpi_rank.eq.0) write(6,*) 'iseed=',iseed
             
! FPAZ topcell
      ! FPAZ  define zcartedge from top to cell k=1 (real mesh) to the top 
       do k=1,l
        do j=1,mp
         do i=1,np
          zcartedgeAZ1(i,j,k+lphase)=0.25*(zcart(zedge(k+1),i-1,j-1)+  
     +    zcart(zedge(k+1),i-1,j)+zcart(zedge(k+1),i,j-1)+zcart(zedge(k+1),i,j))  !top of k cell, x-, y-
          zcartedgeAZ2(i,j,k+lphase)=0.25*(zcart(zedge(k+1),i+1,j-1)+  
     +    zcart(zedge(k+1),i+1,j)+zcart(zedge(k+1),i,j-1)+zcart(zedge(k+1),i,j))  !top of k cell, x+, y-
          zcartedgeAZ3(i,j,k+lphase)=0.25*(zcart(zedge(k+1),i-1,j+1)+  
     +    zcart(zedge(k+1),i-1,j)+zcart(zedge(k+1),i,j+1)+zcart(zedge(k+1),i,j))  !top of k cell, x-, y+
          zcartedgeAZ4(i,j,k+lphase)=0.25*(zcart(zedge(k+1),i+1,j+1)+  
     +    zcart(zedge(k+1),i+1,j)+zcart(zedge(k+1),i,j+1)+zcart(zedge(k+1),i,j))  !top of k cell, x+, y+

           enddo !i
        enddo !j
       enddo !k
          ! ground level           
        do j=1,mp
         do i=1,np
          zcartedgeAZ1(i,j,0)=0.25*(zcart(zedge(1),i-1,j-1)+  
     +    zcart(zedge(1),i-1,j)+zcart(zedge(1),i,j-1)+zcart(zedge(1),i,j))  !top of k cell, x-, y-
          zcartedgeAZ2(i,j,0)=0.25*(zcart(zedge(1),i+1,j-1)+  
     +    zcart(zedge(1),i+1,j)+zcart(zedge(1),i,j-1)+zcart(zedge(1),i,j))  !top of k cell, x+, y-
          zcartedgeAZ3(i,j,0)=0.25*(zcart(zedge(1),i-1,j+1)+  
     +    zcart(zedge(1),i-1,j)+zcart(zedge(1),i,j+1)+zcart(zedge(1),i,j))  !top of k cell, x-, y+
          zcartedgeAZ4(i,j,0)=0.25*(zcart(zedge(1),i+1,j+1)+  
     +    zcart(zedge(1),i+1,j)+zcart(zedge(1),i,j+1)+zcart(zedge(1),i,j))  !top of k cell, x+, y+
       nzonebelowcell1(i,j)=0
       enddo !i
       enddo
! for individual fuel zones in cell k=1  ***needs modification and loops for multiple fuel types****
       do k=1,lphase
          do j=1,mp
             do i=1,np
                minfueldepth=0.1*(zcartedgeAZ1(i,j,lphase+1)-zcartedgeAZ1(i,j,0))
                zcartedgeAZ1(i,j,k)= zcartedgeAZ1(i,j,0)+max(actualfueldepth(i,j,1),minfueldepth) 
                minfueldepth=0.1*(zcartedgeAZ2(i,j,lphase+1)-zcartedgeAZ2(i,j,0))
                zcartedgeAZ2(i,j,k)= zcartedgeAZ2(i,j,0)+max(actualfueldepth(i,j,1),minfueldepth) 
                minfueldepth=0.1*(zcartedgeAZ3(i,j,lphase+1)-zcartedgeAZ3(i,j,0))
                zcartedgeAZ3(i,j,k)= zcartedgeAZ3(i,j,0)+max(actualfueldepth(i,j,1),minfueldepth) 
                minfueldepth=0.1*(zcartedgeAZ4(i,j,lphase+1)-zcartedgeAZ4(i,j,0))
                zcartedgeAZ4(i,j,k)= zcartedgeAZ4(i,j,0)+max(actualfueldepth(i,j,1),minfueldepth) 

             enddo !i
          enddo !j
       enddo !k
!correction for actualfueldepth
         ! check that zcartedgeAZ(i,j,1) is relevant due to actualfueldeth
       do k=1,lphase
          do j=1,mp
             do i=1,np
         if (zcartedgeAZ1(i,j,k).le.zcartedgeAZ1(i,j,0).or.
     +     zcartedgeAZ2(i,j,k).le.zcartedgeAZ2(i,j,0).or.
     +     zcartedgeAZ3(i,j,k).le.zcartedgeAZ3(i,j,0).or.
     +     zcartedgeAZ4(i,j,k).le.zcartedgeAZ4(i,j,0)) then
           write(6,*) 'negative or nul value of actualfueldepth',actualfueldepth(i,j,1)
           stop
         end if
                maxfueldepth1=0.9*(zcartedgeAZ1(i,j,lphase+1)-zcartedgeAZ1(i,j,0))
                maxfueldepth2=0.9*(zcartedgeAZ2(i,j,lphase+1)-zcartedgeAZ2(i,j,0))
                maxfueldepth3=0.9*(zcartedgeAZ3(i,j,lphase+1)-zcartedgeAZ3(i,j,0))
                maxfueldepth4=0.9*(zcartedgeAZ4(i,j,lphase+1)-zcartedgeAZ4(i,j,0))
         if (zcartedgeAZ1(i,j,k).lt.zcartedgeAZ1(i,j,0)+maxfueldepth1.and.
     +     zcartedgeAZ2(i,j,k).lt.zcartedgeAZ2(i,j,0)+maxfueldepth2.and.
     +     zcartedgeAZ3(i,j,k).lt.zcartedgeAZ3(i,j,0)+maxfueldepth3.and.
     +     zcartedgeAZ4(i,j,k).lt.zcartedgeAZ4(i,j,0)+maxfueldepth4) then
           nzonebelowcell1(i,j)=1+nzonebelowcell1(i,j)
          endif
        enddo
        enddo
       enddo
! check the value of zcartedgeAZ
!     we create artificially two zones that will have both solid and gas
!   indices should be checked here
          do j=1,mp
         do i=1,np
!         do k=nzonebelowcell1(i,j)+1,lphase+1
!         zcartedgeAZ(i,j,k)=zcartedgeAZ(i,j,nzonebelowcell1(i,j)+1)+
!     +   (zcartedgeAZ(i,j,lphase+1)-zcartedgeAZ(i,j,nzonebelowcell1(i,j)+1))*
!     +    (k-nzonebelowcell1(i,j)+1)/(lphase-nzonebelowcell1(i,j))
!         if (nzonebelowcell1(i,j).eq.0) then
          rzonesabove=1/float(lphase-nzonebelowcell1(i,j)+1)
          zincrement1=rzonesabove*(zcartedgeAZ1(i,j,lphase+1)
     &               -zcartedgeAZ1(i,j,nzonebelowcell1(i,j)))
          zincrement2=rzonesabove*(zcartedgeAZ2(i,j,lphase+1)
     &               -zcartedgeAZ2(i,j,nzonebelowcell1(i,j)))
          zincrement3=rzonesabove*(zcartedgeAZ3(i,j,lphase+1)
     &               -zcartedgeAZ3(i,j,nzonebelowcell1(i,j)))
          zincrement4=rzonesabove*(zcartedgeAZ4(i,j,lphase+1)
     &               -zcartedgeAZ4(i,j,nzonebelowcell1(i,j)))
          do k=nzonebelowcell1(i,j)+1,lphase
         
             zcartedgeAZ1(i,j,k)=zcartedgeAZ1(i,j,k-1)+zincrement1
             zcartedgeAZ2(i,j,k)=zcartedgeAZ2(i,j,k-1)+zincrement2
             zcartedgeAZ3(i,j,k)=zcartedgeAZ3(i,j,k-1)+zincrement3
             zcartedgeAZ4(i,j,k)=zcartedgeAZ4(i,j,k-1)+zincrement4
            enddo !k
         enddo !i
        enddo !j
! FPAZ  volume
       do k=1,l+lphase
        do j=1,mp
         do i=1,np
         volumeAZ(i,j,k)=dx*dy*.25*(zcartedgeAZ1(i,j,k)+zcartedgeAZ2(i,j,k)+
     +   zcartedgeAZ3(i,j,k)+zcartedgeAZ4(i,j,k)-
     +   zcartedgeAZ1(i,j,k-1)-zcartedgeAZ2(i,j,k-1)-
     +   zcartedgeAZ3(i,j,k-1)-zcartedgeAZ4(i,j,k-1))
         enddo !i
        enddo !j
       enddo !k
! this section to build zcartedge2AZ (n,m,l+lphase) can probably be improved
       
       call allgather3d0(zcartedgeAZ1,zcartedgetmp2AZ,l+lphase)
       do k=0,l+lphase
        do j=1,m
         do i=1,n
          zcartedge2AZ(1,i,j,k)=zcartedgetmp2AZ(i,j,k)
         enddo !i
        enddo !j
       enddo !k
       call allgather3d0(zcartedgeAZ2,zcartedgetmp2AZ,l+lphase)
       do k=0,l+lphase
        do j=1,m
         do i=1,n
          zcartedge2AZ(2,i,j,k)=zcartedgetmp2AZ(i,j,k)
         enddo !i
        enddo !j
       enddo !k
       call allgather3d0(zcartedgeAZ3,zcartedgetmp2AZ,l+lphase)
       do k=0,l+lphase
        do j=1,m
         do i=1,n
          zcartedge2AZ(3,i,j,k)=zcartedgetmp2AZ(i,j,k)
         enddo !i
        enddo !j
       enddo !k
       call allgather3d0(zcartedgeAZ4,zcartedgetmp2AZ,l+lphase)
       do k=0,l+lphase
        do j=1,m
         do i=1,n
          zcartedge2AZ(4,i,j,k)=zcartedgetmp2AZ(i,j,k)
         enddo !i
        enddo !j
       enddo !k
!       k=1
!        do j=1,m
!         do i=1,n-1
!          if (zcartedge2AZ(1,i+1,j,k).ne.zcartedge2AZ(2,i,j,k)) write (6,*) mpi_rank,'wrong zcartedge',i,j
!          if (zcartedge2AZ(3,i+1,j,k).ne.zcartedge2AZ(4,i,j,k)) write (6,*) mpi_rank,'wrong zcartedge',i,j
!         enddo !i
!        enddo !j
!       stop
      endif  !if irad.eq.1 or 2

    
      if(irad.eq.2)then
      deallocate(zcartedgetmp2AZ)
      deallocate(zcartedgeAZ1)
      deallocate(zcartedgeAZ2)
      deallocate(zcartedgeAZ3)
      deallocate(zcartedgeAZ4)
      endif

      end subroutine frad_init
!!!-----------------------end subroutine frad_init  -----------------------------------!!!

!!!----------------------- subroutine papv_defs -----------------------------------!!!
      subroutine papv_defs()
    
      use metryic
      use fireteca
      use turba 
      use xvo
      use constants
 
      Implicit None

      real :: dV,papvtmp,papvsolreal !FPAZ :papvsol real in the bottom layer
      integer :: i,j,k,kreal,ia,ja !,msgcnt
      real ::pr, fvsoot !sootvolumefraction
      !calculate current tempgas  (moved from source_defs because required by new soot model
        do k=1,l
         do j=1,mp
          do i=1,np
              pr=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
              tempg(i,j,k)=xvb(i,j,k,4)/xvb(i,j,k,nv)*(pr*1.e-5)**(rg/cp)
          enddo
         enddo
        enddo
      call updated(temps,temps,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      call updated(tempg,tempg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      call updated(tambientarray,tambientarray,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      papv = 0.0
      papvgas = 0.0
      papvtot = 0.0
      papvtotAZ=0.0
      papvtot2AZ=0.0
      papvmax=-99e9
      papvgasmax=-99e9
!      if(mpi_rank .EQ. 0) WRITE(6,62) MAXVAL(sizescale) !KOO ss->sssolid
! 62   FORMAT('Caution: the max value of sizescale ss for radiation on
!     +        processor 0 is',1X,1pe17.10)
      if(irad.eq.1)then
       do k=lhighlim,llowlim,-1
       !adjust k for  real domain if using diffusion method (kreal indexes real domain arrays)
        kreal = k-nz+1
       do j=1,mp
        do i=1,np
         if(rhof(i,j,kreal).eq.0.0)kfuel(i,j)=max(k-1,nz) 
        enddo
       enddo
      enddo
      endif
      do k=llowlim,lhighlim
       !adjust k for  real domain if using diffusion method (kreal indexes real domain arrays)
       if(irad.eq.1)then
        kreal = k-nz+1
       else
        kreal = k
       endif
       do j=1,mp
        do i=1,np
         !set absemis and aemit
         ! TODO : check with JLD how to improve this
         absemis(i,j,k)=absemis0
         aemit(i,j,k)=(rmaxsootcon*
     &                 (sqrt(max((.21-xvb(i,j,kreal,7)/xvb(i,j,kreal,nv)),0.0)/.21)))
     &                 *xvb(i,j,kreal,nv)
         if(rhof(i,j,kreal).gt.min_rhof) ! FP09/2019 replaced the
            !incorrect index k here by kreal
     &      papv(i,j,k) = (2.0/sizescale(i,j,kreal))*(rhof(i,j,kreal)/(4.0*rhomicro(i,j,kreal)))
         papvgas(i,j,k) = crad*aemit(i,j,k)*absemis(i,j,k)*0.25
         if (isootmodel.eq.1) then
           !HERE FP DEFINE A NEW MODEL FOR ABSORPTION COEFFICIENT:
           !soot volume fraction is defined as correlated to mixture fraction, which is a linear fonction
           ! of tempg and O2 deficit
           fvsoot=0.899e-6*(max(3.05e-5*(tempg(i,j,kreal)-tambientarray(i,j,kreal))
     +             +0.05*(1-xvb(i,j,kreal,7)/(0.21*xvb(i,j,kreal,nv)))
     +             -1.48e-4,0.0))**0.671 
           papvgas(i,j,k)=266.0*7.0*fvsoot*tempg(i,j,kreal)  
         endif
         if (irad.eq.1) then
           if(papv(i,j,k).gt.papvmax)
     +      papvmax = papv(i,j,k)
           if(papvgas(i,j,k).gt.papvgasmax)
     +      papvgasmax = papvgas(i,j,k)
         endif
         ia=(npos-1)*np+i
         ja=(mpos-1)*mp+j
         papvtot(i,j,k) = papv(i,j,k)+papvgas(i,j,k)
        enddo  !do i
       enddo  !do j
      enddo  !do k
      !TODO : temporary output: to remove
!      if (mpi_rank.eq.0) open(unit=21,file='radfile',form='unformatted',
!     +                         status='unknown')
!       call writeiorad(papv,21,1-ih,np+ih,1-ih,mp+ih,l,0) 
!       call writeiorad(papvgas,21,1-ih,np+ih,1-ih,mp+ih,l,0) 
!       if (mpi_rank.eq.0) close(21) 

      if (irad.eq.1) then
       papvtmp = 0.0
       call mpi_allreduce(papvmax,papvtmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
       papvmax = papvtmp
       papvtmp = 0.0
       call mpi_allreduce(papvgasmax,papvtmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
       papvgasmax = papvtmp
      else if(irad.eq.2)then
       do k=2,l  !papvtotAZ
        do j=1,mp
         do i=1,np 
           papvtotAZ(i,j,k+lphase)=papv(i,j,k)+papvgas(i,j,k)
         enddo  !do i
        enddo  !do j
       enddo  !do k
       do j=1,mp
        do i=1,np
         ia=(npos-1)*np+i
         ja=(mpos-1)*mp+j
         dV=volumeAZ(i,j,2)+volumeAZ(i,j,1)
         if (nzonebelowcell1(i,j).eq.0) then
 	   papvtotAZ(i,j,1)=papvgas(i,j,1)+papv(i,j,1) !FPAZ first layer=gas+solid
 	   papvtotAZ(i,j,2)=papvgas(i,j,1)+papv(i,j,1) !FPAZ second layer=gas+solid
         else
 	   papvsolreal=papv(i,j,1)*(volumeAZ(i,j,2)+volumeAZ(i,j,1))
     +       /volumeAZ(i,j,1)   !papvsol*DVsol=papv*DV
 	   !FP multi: be careful if actualfueldepth=0
 	   papvtotAZ(i,j,1)=papvgas(i,j,1)+papvsolreal !FPAZ first layer=gas+solid
 	   papvtotAZ(i,j,2)=papvgas(i,j,1) !FPAZ second layer=gas only
         end if
        enddo  !do i
       enddo  !do j
       call allgather3d(papvtotAZ,papvtot2AZ,l+lphase)
      endif
 
      end subroutine papv_defs
!!!----------------------- end subroutine papv_defs -----------------------------------!!!

!!!----------------------- subroutine source defs -----------------------------------!!!
      subroutine source_defs()
    
      use metryic
      use fireteca 
      use turba 
      use xvo
      use constants
      implicit None


      integer :: i,j,k,kreal,ia,ja
      real :: localterm,nbrterm,yterm!,pr
      real*8:: tmpval
      rstep=1.0


      sourcesol = 0.0
      sourcegas = 0.0
      !do k=llowlim,lhighlim
      do k=llowlim,lhighlim-1   !rrl  9_26_08  (avoid reaching out the top of the grid)
      !adjust kreal if using diffusion method (kreal indexes real domain arrays)
       if(irad.eq.1)then
        kreal = k-nz+1
       else
        kreal = k
       endif
       do j=1,mp
        ja=(mpos-1)*mp+j
        do i=1,np
         ia=(npos-1)*np+i
         !define sourcesol 
         if((temps(i,j,kreal).gt.(tambientarray(i,j,kreal)+rstep)))then
          !nbrterm accounts for neighbor influence on t4barsolid
          if(k.gt.llowlim)then
           frac=(1.0/6.0)*j3+.25*(1-j3)    !sets avg weight for 2-d or 3-d
           nbrterm = csurr*frac*abs(temps(i,j,kreal)-temps(i,j,kreal-1))
          else
           frac=(1.0/5.0)*j3+.25*(1-j3)    !sets avg weight for 2-d or 3-d
           nbrterm = 0.0
          endif
          !FP09/2019 replaced incorrect tempg by temps here
           yterm=abs(temps(i,j,kreal)-temps(i,j-j3,kreal))
     &         +abs(temps(i,j,kreal)-temps(i,j+j3,kreal))
           nbrterm = nbrterm+csurr*frac*
     &        (abs(temps(i,j,kreal)-temps(i-1,j,kreal))
     &        +abs(temps(i,j,kreal)-temps(i+1,j,kreal))
     &        +abs(temps(i,j,kreal)-temps(i,j,kreal+1)) 
     &        +yterm )
          !localterm accounts for local influence on t4barsolid
          localterm = camb*(temps(i,j,kreal)-tambientarray(i,j,kreal))
          wdth=localterm+nbrterm
          t4barsolid(i,j,k)=sqrt((temps(i,j,kreal)**2+facp*wdth**2)
     &                          *(temps(i,j,kreal)**2+facm*wdth**2))
          sourcesol(i,j,k)=4.0*papv(i,j,k)*pathcoef*
     &                     (t4barsolid(i,j,k)-tambientarray(i,j,kreal)**2)
     &                    *(t4barsolid(i,j,k)+tambientarray(i,j,kreal)**2)*
     &                     5.67051e-8
         endif  !if(temps.gt.tambientarray)
 
         !define sourcegas
         if((tempg(i,j,kreal).gt.(tambientarray(i,j,kreal)+rstep)))then
          !nbrterm accounts for neighbor influence on t4barsolid
          if(k.gt.llowlim)then
           frac=(1.0/6.0)*j3+.25*(1-j3)    !sets avg weight for 2-d or 3-d
           nbrterm = csurr*frac*abs(tempg(i,j,kreal)-tempg(i,j,kreal-1))
          else
           frac=(1.0/5.0)*j3+.25*(1-j3)    !sets avg weight for 2-d or 3-d
           nbrterm = 0.0
          endif
          yterm=abs(tempg(i,j,kreal)-tempg(i,j-j3,kreal))
     &         +abs(tempg(i,j,kreal)-tempg(i,j+j3,kreal))
          !nbrterm = csurr*frac*
          nbrterm = nbrterm+csurr*frac*
     &        (abs(tempg(i,j,kreal)-tempg(i-1,j,kreal))
     &        +abs(tempg(i,j,kreal)-tempg(i+1,j,kreal))
     &        +abs(tempg(i,j,kreal)-tempg(i,j,kreal+1))
     &        +yterm )
          !localterm accounts for local influence on t4barsolid
          localterm = camb*(tempg(i,j,kreal)-tambientarray(i,j,kreal))
          wdth=localterm+nbrterm
          t4bar(i,j,k)=sqrt((tempg(i,j,kreal)**2+facp*wdth**2)
     &                     *(tempg(i,j,kreal)**2+facm*wdth**2))
          sourcegas(i,j,k)=4.0*papvgas(i,j,k)*
     &                     (t4bar(i,j,k)-tambientarray(i,j,kreal)**2)
     &                    *(t4bar(i,j,k)+tambientarray(i,j,kreal)**2)*
     &                     5.67051e-8
         endif  !if(tempg.gt.tambientarray)

        enddo !do i
       enddo !do j
      enddo !do k 
      
      ! BELLOW THAT LINE DEFINITION OF EAZ,EAZ2, Etot
      if (irad.eq.2) then
!        COMPUTATION OF SOURCE TERMS (EAZ, E2AZ, nphotemisAZ, nphotemisAZsplit)
      EAZ=0.0
      E2AZ = 0.0  !initialize full E array to zero
      Etot = 0.                 !Etot= Total energy
      if (iradeastflux.eq.1)
     +   eastFlux=0.0

!      IF(ibctopbot .EQ. 0) THEN
!       WRITE(6,*) 'Value ibctopbot - STOP'
!       STOP
!      ENDIF
      !calculate the energy from sources in each cell on a per process subdomain
      DO k=2,l  ! FPAZ all real cell except the first one
       DO j=1,mp
        ja=(mpos-1)*mp+j
        DO i=1,np
         ia=(npos-1)*np+i 
         EAZ(i,j,k+lphase) = (sourcegas(i,j,k) + sourcesol(i,j,k))*volumeAZ(i,j,k+lphase)
         Etot = Etot + EAZ(i,j,k+lphase)
        ENDDO !do i
       ENDDO !do j
      ENDDO !do k
      !FPAZ source computation in first real layer
      DO j=1,mp
        ja=(mpos-1)*mp+j
        DO i=1,np
          ia=(npos-1)*np+i 
          !gas + solid in the first cell
          if (nzonebelowcell1(i,j).eq.0) then
            EAZ(i,j,1) = (sourcegas(i,j,1)+sourcesol(i,j,1))*volumeAZ(i,j,1)
            !gas+solid in the second cell
            EAZ(i,j,2) = (sourcegas(i,j,1)+sourcesol(i,j,1))*volumeAZ(i,j,2)
            Etot = Etot + EAZ(i,j,1)+ EAZ(i,j,2)
          else
            !gas + solid in the first cell
            EAZ(i,j,1) = (sourcegas(i,j,1)+sourcesol(i,j,1))*volumeAZ(i,j,1)+
     +              sourcesol(i,j,1)*volumeAZ(i,j,2)
            !gas only in the second cell
            EAZ(i,j,2) = sourcegas(i,j,1)*volumeAZ(i,j,2)
            Etot = Etot + EAZ(i,j,1)+ EAZ(i,j,2)
          end if
        ENDDO !do i
      ENDDO !do j
     
      !gather the all subdomains to create full E array on all processes
       call allgather3d(EAZ,E2AZ,l+lphase)
      call mpi_allreduce(Etot,tmpval,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      Etot = tmpval
      IF(Etot.LT.Emin) THEN
!        if (mpi_rank.eq.0) WRITE(6,*) 'Etot=',Etot,'is neglected'
        if (mpi_rank.eq.0) WRITE(6,'(a,es17.10,a)') 
     +                          'Etot=',Etot,' is neglected'

      ELSE
      ! TODO check how to improve lmc computation
       lmc=6
       do k=1,l+lphase
        do j=1,m
         do i=1,n
!          if((papvtot2AZ(i,j,k).GT.papvthresh).or.
!     +       ((E2AZ(i,j,k)/real(nphotbatch)).ge.1)) then
          if((papvtot2AZ(i,j,k).GT.papvthresh2).or.
     +      ((E2AZ(i,j,k)/Etot*nphotbatch).ge.1)) then
           lmc = MAX(k,lmc)
          endif
         enddo !do i
        enddo !do j
       enddo !do k
      if(mpi_rank.eq.0)then
       WRITE(6,60) lmc  !,papvtot2AZ(72,47,1),papvtot2AZ(72,47,2),papvtot2AZ(72,47,3)
       !WRITE(6,*) 'kpapvgmax, kpapvmax: ',kpapvgmax,kpapvmax
      endif
 60   FORMAT('Maximum vertical cell index for Monte-Carlo',1X,I3)
      ENDIF
      endif !irad.eq.2
      end subroutine source_defs
!!!----------------------- end subroutine source defs -----------------------------------!!!


!!!*_*_*_*_*_*_*_*_*_*_*_*_*_ Monte Carlo method routines *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*!!!
!!!----------------------- subroutine firerad_MC  -----------------------------------!!!
      subroutine firerad_MC()

      use fireteca
      use metryic
      use turba
!      use constants,only: ss,rhomicrovalue
      Implicit None
      integer :: j!i, j, k !, kk, kf, it
      real, external :: zcart2

      double precision :: wt1, wt2, wtick !wt3 

      real:: Pi

      if(mpi_rank.eq.0) then
         wtick = mpi_Wtick()
         wt1=mpi_Wtime()
      endif

      !if(mpi_rank.EQ.0) print*,'Beginning of firerad_MC'

      !if((it.gt.1).or.(irst.eq.1))then


      ! TODO: check this section, that seems weird (FP 6/08/2012):
      ! I guess it should not permit to have a computation that does not
      ! depend of the processor number. Why using multiple seed
      if (ittot.ge.1) then
      !if(ittot.gt.itrestart+1) then
      ! iseed = iseed0   +ittot
      !endif
       !write(6,*) 'iseed = ',iseed
       call random_seed(size=seedsize)
       do j=1,numprocs
        seed(j) = (mpi_rank+1)*j*iseed
       enddo
       call random_seed(put=seed)
      endif


      Pi = ACOS(-1.)
      rstep=1.0

       nphotbatch = 1E6 !nproc*1e5 !CEILING(Etot/maxQphot)
      if (iradeastflux.eq.1) nphotbatch = 3E5
      !call rmaxmin1(papv,'papv',1,np,1,mp,l)
      !define papv, papvgas, papvtot 
      call papv_defs()

      !call rmaxmin1(actualfueldepth,'actualfueldepth',1,np,1,mp,l)
      !call rmaxmin1(papv,'papv',1,np,1,mp,l)
      !call rmaxmin1(papvgas,'papvgas',1,np,1,mp,l)
      !call rmaxmin1(temps,'temps',1,np,1,mp,l)
      
      !define sourcegas, sourcesol, EAZ,E2AZ,Etot, and lmc
      call source_defs()
      
      !call rmaxmin1(sourcegas,'sourcegas',1-ih,np+ih,1-ih,mp+ih,l)
      !call rmaxmin1(sourcesol,'sourcesol',1-ih,np+ih,1-ih,mp+ih,l)
      

       if (Etot.ge.Emin) call montecarlo_def()
      !if(mpi_rank.EQ.0) print*,'fin montecarlo'

      call rmaxmin1(-frhosiesrad, 'frhosiesrad', 1-ih, np+ih, 1-ih, mp+ih, l)
      call rmaxmin1(-firad, 'firad', 1-ih, np+ih, 1-ih, mp+ih, l)
      !JAS adding mpi wall timers to get some idea of performance gains.
      if(mpi_rank.eq.0)then
         wt2 = mpi_Wtime()
         !write(6,40) wt2 - wt1
! 40      FORMAT('Total time  in firerad Monte-Carlo method',1X,1pe15.5)
      endif

      end subroutine firerad_MC
!!!----------------------- End subroutine firerad_MC  -----------------------------------!!!

!!!-----------------------  subroutine allocate_photons  -----------------------------------!!!
      subroutine allocate_photons()
      USE gridsetup
      USE msga
      USE turba !FPAZ actualfueldepth
      USE metryic
!      USE fireteca,only :eastFlux
      IMPLICIT NONE
      INTEGER i, j, k, icpt, icptold, ia, ja !, j1,itmp,jtmp
      INTEGER  nphottotbis   !, ntot
      INTEGER  n1,n2, iphot
      real*8 Etot2,Qphotbatch
!      integer::tmpint
!      real*8:: tmpval
      INTEGER  mypstrt, mypstp
      Qphotbatch = Etot/dble(nphotbatch)
!      IF(mpi_rank .EQ. 0) THEN
!       WRITE(6,44) Qphot, Etot
! 44    FORMAT('Energy of a photon=',1X,1pd17.10,1X,'Total energy',1X,1pd17.10)
!      ENDIF
      !partition the total number of photons amongst all processes
      iphot = 0
      n1 = nphotbatch/nproc
      n2 = MOD(nphotbatch, nproc)
      DO i=0,nproc-1
       nperproc(i) = n1
       iphot = iphot + nperproc(i)
      ENDDO
      ! Remaining procs
      DO i=0,n2-1
       nperproc(i) = nperproc(i) + 1
       iphot = iphot + 1
      ENDDO
      IF(iphot.NE.nphotbatch) THEN
       IF(mpi_rank .EQ. 0) WRITE(6,*)' ERROR in distribution of photons'
      ENDIF
      !Fill an integer array with the number of photons to be emmitted in each cell
      !Do this on each processes' subdomian
      Etot2 = 0.
       DO k=1,lmc !FPAZ
        DO j=1,mp
         ja=(mpos-1)*mp+j
         DO i=1,np
          ia=(npos-1)*np+i 
           nphotemisAZtemp(ia,ja,k) = NINT(E2AZ(ia,ja,k)/Qphotbatch)
          Etot2 = Etot2 + nphotemisAZtemp(ia,ja,k)*Qphotbatch
         ENDDO !do i
        ENDDO !do j
       ENDDO !do k
      !gather the all subdomains to create full nphotemis array on all processes
      call mpi_allreduce(nphotemisAZtemp,nphotemisAZ,
     &      n*m*lmc,mpi_integer,mpi_sum,mpi_comm_world,ierror)
      !sum Etot2 over all processes
!      call mpi_allreduce(Etot2,tmpval,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
!      Etot2=tmpval
!      IF(mpi_rank.EQ.0) THEN
!         WRITE(6,46)  Etot, Etot2, 100.*ABS(Etot - Etot2)/Etot
!      ENDIF
! 46   FORMAT('Initial energy Etot= ',1X,1pd17.10,1X,
!     +     ' New energy Etot2=',1X,1pd17.10,1X,
!     +     ' % error on energy=',1X,1pd12.5)

      nphottotbis = SUM(nphotemisAZ)
      IF(mpi_rank .EQ. 0) then
        WRITE(6,45) nphottotbis
 45   FORMAT('There are',1X,I15,1X,' photons in the batch')
        if (real(nphotbatch-nphottotbis)/nphotbatch.ge.0.001)
     +     write(6,*) 'warning: nphotonbatch might be too low, error:',
     +          real(nphotbatch-nphottotbis)/nphotbatch*100, '%'
      ENDIF

145   FORMAT(I10,1X,'photon trajectories corrected')
      !define a start and stop point for the range of photons each process emits
      !*****-----------------JAS temporary test here!-------------------------******!
      !mypstrt = SUM(nperproc(0:mpi_rank-1)) !+1
      !mypstp = mypstrt + nperproc(mpi_rank)
      !icpt = 0
      !DO k=1,lmc
       !DO j=1,m
        !DO i=1,n 
        !icptold = icpt
        !icpt = icpt + nphotemis(i,j,k)
        !if((icptold.lt.mypstrt).and.(icpt.gt.mypstrt))then
         !This process must account for the second partition of the photons in cell i,j,k 
        ! nphot(i,j,k) = icpt-mypstrt
        !elseif((icptold.ge.mypstrt).and.(icpt.le.mypstp))then
         !This process must account for all of the photons in cell i,j,k 
        ! nphot(i,j,k) = nphotemis(i,j,k)
        !elseif((icptold.lt.mypstp).and.(icpt.gt.mypstp))then
         !This process must account for the first partition of the photons in cell i,j,k 
        ! nphot(i,j,k) = mypstp-icptold
        !endif
        !ENDDO !do i
       !ENDDO !do j
      !ENDDO !do k
      !*****-----------------JAS temporary test here!-------------------------******!

        mypstrt = SUM(nperproc(0:mpi_rank-1))+1
        mypstp = mypstrt + nperproc(mpi_rank)-1
        icpt = 0
        DO k=1,lmc
         DO j=1,m
          DO i=1,n
          icptold = icpt
          icpt = icpt + nphotemisAZ(i,j,k)
        !if((icptold.lt.mypstrt).and.(icpt.gt.mypstrt))then
          if(icptold.le.mypstrt)then
           if(icpt.ge.mypstp)then
          !all of my photons go here
            nphotemisAZsplit(i,j,k) = mypstp-mypstrt+1
         elseif(icpt.gt.mypstrt)then   !icpt is less than mypstp and greater than mypstrt
          nphotemisAZsplit(i,j,k) = icpt-mypstrt+1   !so only some of my photons go here
         endif
        elseif(icptold.le.mypstp)then !icptold>mypstrt and icptold<=mypstp
         if(icpt.ge.mypstp)then
          !the remainder of my photons go here 
          nphotemisAZsplit(i,j,k) = mypstp-icptold
         else   !icpt is less than mypstp so only some of my remiaining photons go here
          nphotemisAZsplit(i,j,k) = icpt-icptold
         endif
        endif
        ENDDO !do i
       ENDDO !do j
      ENDDO !do k
      end subroutine allocate_photons
!!!----------------------- end subroutine allocate_photons  -----------------------------------!!!


!!!----------------------- subroutine shoot_photons  -----------------------------------!!!
! this subroutine will shoot nphotemisAZsplit(ia,ja,k) and see where they arrive (number added in nphotonAZ)
!
      subroutine shoot_photons()
      
      USE gridsetup
      USE msga
      !USE turba !FPAZ actualfueldepth
      USE metryic
!      USE fireteca,only :eastFlux
      IMPLICIT NONE
      REAL p                    !optical thickness
      integer::i,j,k,iphot 
      REAL*8 p1, p2, p3
      REAL*8  zbottomcell, ztopcell
      REAL*8  zbottomreal, ztopreal
      REAL*8 z1top, z2top, z3top,z4top,z1bot,z2bot,z3bot,z4bot
      REAL*8 xf, yf, zf, dir(3) ! random direction
      INTEGER vox(3)            !voxel where photon stops
      EXTERNAL gdeform, zcart,zcart2
!      integer::tmpint
!       INTEGER ierr
      DO k=1,lmc
       DO j=1,m
        DO i=1,n ! Loop for photons emitted by cell (i,j,k)
         DO iphot=1,nphotemisAZsplit(i,j,k)
          p = RAND(1.)
          CALL RANDOMDIRNORM2(dir)
          z1top=zcartedge2AZ(1,i,j,k)
          z2top=zcartedge2AZ(2,i,j,k)
          z3top=zcartedge2AZ(3,i,j,k)
          z4top=zcartedge2AZ(4,i,j,k)
          z1bot=zcartedge2AZ(1,i,j,k-1)
          z2bot=zcartedge2AZ(2,i,j,k-1)
          z3bot=zcartedge2AZ(3,i,j,k-1)
          z4bot=zcartedge2AZ(4,i,j,k-1)
          !Determine initial photon starting positions (xf,yf,zf)
          zbottomcell=min(z1bot,z2bot,z3bot,z4bot)
          ztopcell=max(z1top,z2top,z3top,z4top)
          do
          !do while(zf.lt.zbottomreal.and.zf.gt.ztopreal)
          p1 = RAND(1.0)
          p2 = RAND(1.0)
          p3 = RAND(1.0)
          xf = (i - 1 + p1)*dx
          yf = (j - 1 + p2)*dy
          zf = zbottomcell + p3*(ztopcell-zbottomcell)
          zbottomreal=zsurface(p1,p2,z1bot,z2bot,z3bot,z4bot)
          ztopreal=zsurface(p1,p2,z1top,z2top,z3top,z4top)
          if (zf.ge.zbottomreal.and.zf.le.ztopreal) exit
          enddo
          !calcule la cellule ou le photon arrive
          CALL CELLPHOTARRIVE_DEFAZ(p, xf, yf, zf,i,j,k, dir, vox) 
 70               FORMAT('Maille',1X,7(I3,1X))
          IF( ((vox(1) .LT. 1) .OR. (vox(1) .GT. n)) .OR.
     &        ((vox(2) .LT. 1) .OR. (vox(2) .GT. m)) .OR.
     &        ((vox(3) .LT. 1) .OR. (vox(3) .GT. lmc))) THEN
         ELSE
           nphotonAZ(vox(1),vox(2),vox(3)) = nphotonAZ(vox(1),vox(2),vox(3)) + 1
          ENDIF
         ENDDO !do iphot
        ENDDO !do i
       ENDDO !do j
      ENDDO !do k
      !TODO: nphotwrong is supposed to be produced only when topo, can probably be removed
!            call mpi_allreduce(nphotwrong,tmpint,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)
!      nphotwrong=tmpint
!      IF(mpi_rank .EQ. 0) WRITE(6,145) nphotwrong
145   FORMAT(I10,1X,'photon trajectories corrected')
!      nphotwrong=0
       !TODO l+lphase could be reduced to lmc
!       CALL MPI_ALLREDUCE(nphotonAZ, nphotonAZgather, (l+lphase)*n*m, MPI_INTEGER, MPI_SUM,
!     &     MPI_COMM_WORLD, ierr)
      
      END SUBROUTINE shoot_photons

!!!----------------------- updateRadArrays  -----------------------------------!!!
! this subroutine update values of firad and frhosiesrad and nphotonEastFlux according to nphotonAZtot and Qphot
      SUBROUTINE updateRadArrays(Qphot)
       USE gridsetup
       USE msga
       USE turba !FPAZ actualfueldepth
       USE metryic
       USE fireteca,only :eastFlux,firad,frhosiesrad
       IMPLICIT NONE
       INTEGER i, j, k,   ia, ja !, j1,itmp,jtmp
       real :: papvsolreal
       real*8 :: Qphot
       real,EXTERNAL :: zcart
                
        DO k=1+lphase,lmc-1  !FPAZ
         DO j=1,mp
          ja=(mpos-1)*mp+j
          DO i=1,np
           ia=(npos-1)*np+i 
           IF(papvtot2AZ(ia,ja,k+lphase) .GT. papvthresh) THEN
             firad(i,j,k) = -sourcegas(i,j,k) +1.0/volumeAZ(i,j,k+lphase)*  !corrected by FP
     &         Qphot*nphotonAZgather(ia,ja,k+lphase)*papvgas(i,j,k)/papvtot2AZ(ia,ja,k+lphase)
             frhosiesrad(i,j,k) = -sourcesol(i,j,k) +1.0/volumeAZ(i,j,k+lphase)* !corrected by FP
     &         Qphot*nphotonAZgather(ia,ja,k+lphase)*papv(i,j,k)/papvtot2AZ(ia,ja,k+lphase)
           ENDIF !papvtot2.gt.papvthresh
           if (iradeastflux.eq.1)
     +       eastFlux(i,j,k)=Qphot*nphotonEastFlux(ia,ja,k+lphase)/(dy*
     +             (zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)))
          ENDDO !do i
         ENDDO !do j
        ENDDO !do k
 	!FPAZ FIRST CELL 
        DO j=1,mp
          ja=(mpos-1)*mp+j
          DO i=1,np
           ia=(npos-1)*np+i 
           if (iradeastflux.eq.1)
     +       eastFlux(i,j,1)=Qphot*(nphotonEastFlux(ia,ja,1)+
     +               nphotonEastFlux(ia,ja,1+lphase))/(dy*
     +             (zcart(zedge(2),i,j)-zcart(zedge(1),i,j)))
           firad(i,j,1) = -sourcegas(i,j,1) 
           frhosiesrad(i,j,1) = -sourcesol(i,j,1) 
       if (nzonebelowcell1(i,j).eq.0) then
       IF(papvtot2AZ(ia,ja,1) .GT. papvthresh) THEN
         firad(i,j,1) = firad(i,j,1)+1.0/(volumeAZ(i,j,1)+volumeAZ(i,j,2))*
     &         Qphot*nphotonAZgather(ia,ja,1)*papvgas(i,j,1)/papvtot2AZ(ia,ja,1)
         frhosiesrad(i,j,1) = frhosiesrad(i,j,1)+1.0/(volumeAZ(i,j,1)+volumeAZ(i,j,2))*
     &         Qphot*nphotonAZgather(ia,ja,1)*papv(i,j,1)/papvtot2AZ(ia,ja,1)
        ENDIF !papvtot2.gt.papvthresh
        IF(papvtot2AZ(ia,ja,2) .GT. papvthresh) THEN
         firad(i,j,1) = firad(i,j,1)+1.0/(volumeAZ(i,j,1)+volumeAZ(i,j,2))*
     &         Qphot*nphotonAZgather(ia,ja,1+lphase)*papvgas(i,j,1)/papvtot2AZ(ia,ja,1+lphase)
         frhosiesrad(i,j,1) = frhosiesrad(i,j,1)+1.0/(volumeAZ(i,j,1)+volumeAZ(i,j,2))*
     &         Qphot*nphotonAZgather(ia,ja,1+lphase)*papv(i,j,1)/papvtot2AZ(ia,ja,1+lphase)
         ENDIF !papvtot2.gt.papvthresh
       else !nzonebelowcell...
         IF(papvtot2AZ(ia,ja,1) .GT. papvthresh) THEN
           firad(i,j,1) = firad(i,j,1)+1.0/(volumeAZ(i,j,1)+volumeAZ(i,j,2))*
     &     Qphot*nphotonAZgather(ia,ja,1)*papvgas(i,j,1)/papvtot2AZ(ia,ja,1)
           papvsolreal=papv(i,j,1)*(volumeAZ(i,j,1)+volumeAZ(i,j,2))/
     +      volumeAZ(i,j,1)
           frhosiesrad(i,j,1) = frhosiesrad(i,j,1)+1.0/(volumeAZ(i,j,1)+volumeAZ(i,j,2))*
     &       Qphot*nphotonAZgather(ia,ja,1)*papvsolreal/papvtot2AZ(ia,ja,1)
         ENDIF !papvtot2.gt.1e-6
         IF(papvtot2AZ(ia,ja,2) .GT. papvthresh) THEN
            firad(i,j,1) = firad(i,j,1)+1.0/(volumeAZ(i,j,1)+volumeAZ(i,j,2))*
     &        Qphot*nphotonAZgather(ia,ja,1+lphase)*papvgas(i,j,1)/papvtot2AZ(ia,ja,1+lphase)
         ENDIF !papvtot2.gt.1e-6
       end if ! nzonebelow
         ENDDO !do i
         ENDDO !do j
       END SUBROUTINE updateRadArrays


!!!----------------------- subroutine montecarlo_def  -----------------------------------!!!

      SUBROUTINE MONTECARLO_DEF()

      USE gridsetup
      USE msga
      USE turba !FPAZ actualfueldepth
      USE metryic
!      USE fireteca,only :eastFlux
      IMPLICIT NONE
      !INCLUDE 'mpif.h'


      INTEGER ierr
      INTEGER i, j, k! ,ia, ja, j1,itmp,jtmp
      INTEGER ::nphottot=0
      INTEGER :: iternum  !,lnblnk
      !CHARACTER*30 filein
      !CHARACTER*5 s1
      !real :: papvsolreal
      real*8 :: maxpercabsdiff,maxpercabsdifftmp
      real*8 :: Qphot
      real,EXTERNAL :: zcart
       real::accuracythreshold
      !IF (mpi_rank.EQ.0) WRITE (*,*) 'Beginning of Montecarlo''s method'
      nphotonAZ = 0
      nphotonAZgather = 0
      nphotemisAZ = 0
      nphotemisAZsplit = 0
      nphotonEastFlux=0
      nphotontempEastFlux=0
      ! DEFINITION OF ARRAY nphotemisAZsplit that contains on the whole domain the number of photons
      ! each processor should shoot from ia,ja,k
      call allocate_photons()
      percabsAZold = 0.0
      percabsAZnew = 0.0
      maxpercabsdiff = 1.0
      maxpercabsdifftmp = 0.0
      iternum = 0
      nphottot=0
      accuracythreshold=5e-5
      if (iradeastflux.eq.1) accuracythreshold=1e-5
      do while(maxpercabsdiff.gt.accuracythreshold)  !loop for dynamic allocation of photons
      !do while(maxpercabsdiff.gt.1e-6)  !loop for dynamic allocation of photons
        call shoot_photons()
        nphottot=nphottot+nphotbatch 
        iternum = iternum+1
        !WRITE (6,48) MAXVAL(nphotonAZ(1:n,1:m,1:lmc)), mpi_rank
 48     FORMAT('Maximal number of absorbed photons',1X,I12,1X,
     &     'processor',1X,I4)
        ! TODO: the following sequence can be optimized now that nphotonAZtot is done,
        ! with arrays of size np,mp for percabsAZold and new
        do k=1,lmc !l+lphase
          do j=1,m
           do i=1,n
            percabsAZold(i,j,k) = percabsAZnew(i,j,k)
            percabsAZnew(i,j,k) = DBLE(nphotonAZ(i,j,k))/DBLE(nphottot)
            !percabsAZnew(i,j,k) = DBLE(nphotonAZgather(i,j,k))/DBLE(nphottot)
            enddo
          enddo
         enddo
       maxpercabsdiff = 0.0
       do k=1,lmc !l+lphase
        do j=1,m
         do i=1,n
          maxpercabsdifftmp = abs(percabsAZnew(i,j,k)-percabsAZold(i,j,k))
          if(maxpercabsdifftmp.gt.maxpercabsdiff)
     +     maxpercabsdiff = maxpercabsdifftmp
         enddo
        enddo
       enddo
       CALL MPI_ALLREDUCE(maxpercabsdiff,maxpercabsdifftmp,1, MPI_DOUBLE_PRECISION, MPI_MAX,
     &     MPI_COMM_WORLD, ierr)
        maxpercabsdiff=maxpercabsdifftmp        
       !if (mpi_rank.eq.0) WRITE(6,*) maxpercabsdiff
      enddo  !do while(maxpercabsdiff < 0.01) end of dynamical alllocation
    
      Qphot = Etot/dble(nphottot)
      IF(mpi_rank .EQ. 0) THEN
       WRITE(6,44) Qphot, Etot,iternum
 44    FORMAT('Energy of a photon=',1X,1pd17.10,1X,'Total energy',1X,1pd17.10,1X, 'Iteration',1X,I3)
      ENDIF

      CALL MPI_ALLREDUCE(nphotonAZ, nphotonAZgather, lmc*n*m, MPI_INTEGER, MPI_SUM,
     &     MPI_COMM_WORLD, ierr)
!        IF(mpi_rank .EQ. 0)
!     +     WRITE(6,*) sum(nphotonAZgather(:,:,:))
!     +     WRITE(6,*) 'last 47 86 1:',xxdble(nphotonAZgather(47,86,1))/dble(nphottot)
      if (iradeastflux.eq.1) then
         CALL MPI_ALLREDUCE(nphotonEastFlux, nphotontempEastFlux, lmc*n*m, MPI_INTEGER, MPI_SUM,
     &     MPI_COMM_WORLD, ierr)
          nphotonEastFlux = nphotontempEastFlux   
      endif 


      call updateRadArrays(Qphot)

 52   FORMAT('Calculation of total radiation',3(1X,1pe17.10))
      END SUBROUTINE MONTECARLO_DEF
!!!----------------------- End subroutine MONTECARLO_DEF-----------------------------------!!!

!!!----------------------- subroutine CELLPHOTOARRIVE_DEF -----------------------------------!!!
      SUBROUTINE CELLPHOTARRIVE_DEFAZ( p, xf, yf, zf,i0,j0,k0, dir, vox)
      USE gridsetup
      USE metryic
      USE msga
      USE restarta  !to have topo file
      IMPLICIT NONE
      REAL p
      REAL*8 x1, y1, z1, xf, yf, zf, lopt,x1old,y1old,z1old  ! position dans la flamme
      INTEGER suivvox(3),  i1,  j1, k1,i1old,j1old,k1old
      INTEGER suivvoxold(3), reversetop, reversebot
      INTEGER suivvoxold2(3)
      INTEGER :: i0,j0,k0
      INTEGER vox(3)            !voxel where photon stops
      !INTEGER iprocx, iprocy
      INTEGER itcount2
      REAL*8 t, tx, ty,tztop,tzbot  
      REAL*8 z1bot,z2bot,z3bot,z4bot,z1top,z2top,z3top,z4top   !bottom and top of the corners
      REAL*8 dir(3)             ! random direction
      REAL*8 acoeftop, bcoeftop, ccoeftop             ! monome of the nd degre equation on t for intersection
                                                      ! between photon trajectory and surface
      REAL*8 acoefbot, bcoefbot, ccoefbot             ! monome of the nd degre equation on t for intersection
                                                      ! between photon trajectory and surface
      REAL*8 tz1,tz2     ! root of the equation
      REAL*8 delta     ! b2-4ac
      REAL*8 p1,p2,ztopreal,zbottomreal !for tests 
      x1 = xf
      y1 = yf
      z1 = zf
      i1=i0
      j1=j0
      k1=k0
      itcount2=0
      suivvox=0
      if((dir(1).eq.0.0).or.(dir(2).eq.0.0).or.(dir(3).eq.0.0))then
       write(6,*) 'A trajectory vector component is zero!!!!'
       write(6,*) 'dir = ',dir
      endif
c      write (*,*)' suivvox(1), (2), (3)=', suivvox(1), suivvox(2), suivvox(3)

      lopt = -LOG(p)
      !DO WHILE((i1>=1).AND.(i1<=n).AND.(j1>=1).AND.(j1<=m).AND.(k1>=1).AND.(k1<=lmc).AND.(lopt>0.).AND.(itcount2.LE.10))
      DO WHILE((i1>=1).AND.(i1<=n).AND.(j1>=1).AND.(j1<=m).AND.(k1>=1).AND.(k1<=lmc).AND.(lopt>0.))
      if((dir(1).eq.0.0).or.(dir(2).eq.0.0).or.(dir(3).eq.0.0))then
       write(6,*) 'A trajectory vector component is zero!!!!'
       write(6,*) 'dir = ',dir
      endif
      itcount2=itcount2+1
       IF(dir(1)>0.) THEN
        tx=(i1*dx-x1)/dir(1)
       ELSE
        tx=(i1*dx-dx-x1)/dir(1)
       ENDIF
       IF(dir(2)>0.) THEN
        ty=(j1*dy-y1)/dir(2)
       ELSE
        ty=(j1*dy-dy-y1)/dir(2)
       ENDIF
        tztop=1e6
        tzbot=1e6
        delta=0.
        tz1=0.
        tz2=0.
          z1top=zcartedge2AZ(1,i1,j1,k1)
          z2top=zcartedge2AZ(2,i1,j1,k1)
          z3top=zcartedge2AZ(3,i1,j1,k1)
          z4top=zcartedge2AZ(4,i1,j1,k1)
          z1bot=zcartedge2AZ(1,i1,j1,k1-1)
          z2bot=zcartedge2AZ(2,i1,j1,k1-1)
          z3bot=zcartedge2AZ(3,i1,j1,k1-1)
          z4bot=zcartedge2AZ(4,i1,j1,k1-1)
    
          acoeftop=(z1top+z4top-z2top-z3top)*dir(1)*dir(2)/(dx*dy)
          bcoeftop=(-dir(3)+dir(2)/dy*(z3top-z1top)+dir(1)/dx*(z2top-z1top)
     +    +(z1top+z4top-z2top-z3top)*(dir(2)/dy*(x1/dx-(i1-1))+dir(1)/dx*(y1/dy-(j1-1))))
          ccoeftop=(z1top-z1+(z3top-z1top)*(y1/dy-(j1-1))+(z2top-z1top)*(x1/dx-(i1-1))
     +    +(z1top+z4top-z2top-z3top)*(x1/dx-(i1-1))*(y1/dy-(j1-1)))
           reversetop=1
           reversebot=1
          if (suivvox(3).ge.0.) then 
          if  (acoeftop.eq.0.) then    !a.eq.0
               if (ccoeftop*bcoeftop.lt.0.) tztop=-ccoeftop/bcoeftop  ! b.eq.0 correct intersection with tz>0
          else   ! the code never goes here in case of flat topo
             delta=bcoeftop**2-4*acoeftop*ccoeftop
             if (delta.gt.0.) then     !two solutions
             tz1=-0.5*(bcoeftop-dsqrt(delta))/acoeftop
             tz2=-0.5*(bcoeftop+dsqrt(delta))/acoeftop
                !if (tz1.ge.0.) then   ! choose the smallest positive root
                if (tz1.gt.0.) then   ! choose the smallest positive root
                     tztop=tz1 
                   !if (tz2.ge.0.)   tztop=min(tz1,tz2)
                   if (tz2.gt.0.)   tztop=min(tz1,tz2)
                !else if (tz2.ge.0.) then
                else if (tz2.gt.0.) then
                   tztop=tz2
                end if
             end if
          end if


c  this event seems to be very scarse, maybe not worse tracking it, at least with small cells
          else  if (acoeftop.ne.0.) then         ! look for potential second intersection
          ! here we use the fact that ccoef should be equal to 0 in this case (suivvox(3)=-1)
             if (bcoeftop*acoeftop.lt.0) then
                       tztop=-bcoeftop/acoeftop  ! choose the second root when positive (first should be 0)
                     reversetop=-1
             end if 
          end if   !end suivvox(3)
          
         acoefbot=(z1bot+z4bot-z2bot-z3bot)*dir(1)*dir(2)/(dx*dy)
          bcoefbot=(-dir(3)+dir(2)/dy*(z3bot-z1bot)+dir(1)/dx*(z2bot-z1bot)
     +    +(z1bot+z4bot-z2bot-z3bot)*(dir(2)/dy*(x1/dx-(i1-1))+dir(1)/dx*(y1/dy-(j1-1))))
          ccoefbot=(z1bot-z1+(z3bot-z1bot)*(y1/dy-(j1-1))+(z2bot-z1bot)*(x1/dx-(i1-1))
     +    +(z1bot+z4bot-z2bot-z3bot)*(x1/dx-(i1-1))*(y1/dy-(j1-1)))
         
          if (suivvox(3).le.0) then
          if  (acoefbot.eq.0.) then    !a.eq.0
               if (ccoefbot*bcoefbot.lt.0.) tzbot=-ccoefbot/bcoefbot ! correct intersection with tz>0
          else
             delta=bcoefbot**2-4*acoefbot*ccoefbot
             if (delta.gt.0.) then     !two solutions
             tz1=-0.5*(bcoefbot-dsqrt(delta))/acoefbot
             tz2=-0.5*(bcoefbot+dsqrt(delta))/acoefbot
                !if (tz1.ge.0.) then   ! choose the smallest positive root
                if (tz1.gt.0.) then   ! choose the smallest positive root
                     tzbot=tz1
                   !if (tz2.gt.0.) tzbot=min(tz1,tz2)
                   if (tz2.gt.0.) tzbot=min(tz1,tz2)
                !else if (tz2.ge.0.) then
                else if (tz2.gt.0.) then
                   tzbot=tz2
                end if
             end if
          end if
c  this event seems to be very scarse, maybe not worse tracking it, at least with small cells
          else  if (acoefbot.ne.0.) then         ! look for potential second intersection
          ! here we use the fact that ccoef should be equal to 0 in this case (suivvox(3)=-1)
             if (bcoefbot*acoefbot.lt.0) then
                       tzbot=-bcoefbot/acoefbot  ! choose the second root when positive (first should be 0)
                     reversebot=-1
             end if
          end if   !end suivvox(3)
c          if (tzbot.le.1000.and.mpi_rank.eq.0) then 
c            write (6,*) 'bot', tzbot,z1+dir(3)*tzbot,zsurface((x1+dir(1)*tzbot)/dx-(i1-1),
c     +       (y1+dir(2)*tzbot)/dy-(j1-1),z1bot,z2bot,z3bot,z4bot)
c          end if


       t = MIN(tx,ty,tzbot,tztop)
       
        suivvoxold2(1)=suivvoxold(1)
        suivvoxold2(2)=suivvoxold(2)
        suivvoxold2(3)=suivvoxold(3)
        suivvoxold(1)=suivvox(1)
        suivvoxold(2)=suivvox(2)
        suivvoxold(3)=suivvox(3)
       lopt = lopt - t*papvtot2AZ(i1,j1,k1)
       IF(lopt>=0)THEN
        suivvox=(/0,0,0/)
        IF(tx<=t)THEN
         IF(dir(1)>0) then
              suivvox = suivvox+(/1,0,0/)
              !compute flux leaving through the east face
              if (iradeastflux.eq.1)
     +       nphotonEastFlux(i1,j1,k1)=nphotonEastFlux(i1,j1,k1)+1
         ENDIF           
         IF(dir(1)<0) suivvox = suivvox+(/-1,0,0/)
        ENDIF
        IF(ty<=t) THEN
         IF (dir(2)>0) suivvox = suivvox+(/0,1,0/)
         IF (dir(2)<0) suivvox = suivvox+(/0,-1,0/)
        ENDIF
        IF (tztop<=t) then
            suivvox=suivvox+(/0,0,1/)
c        ELSE IF (tztop<=t*1.5) then
c              if (z1+dir(3)*t.ge.zsurface((x1+dir(1)*t)/dx-(i1-1),
c     +       (y1+dir(2)*t)/dy-(j1-1),z1top,z2top,z3top,z4top)) then
c              suivvox=suivvox+(/0,0,1/)
c           end if

        ENDIF
        IF (tzbot<=t) then
            suivvox=suivvox+(/0,0,-1/)
c       ELSE IF (tzbot<=t*1.5) then
c             if (z1+dir(3)*t.ge.zsurface((x1+dir(1)*t)/dx-(i1-1),
c    +       (y1+dir(2)*t)/dy-(j1-1),z1bot,z2bot,z3bot,z4bot)) then
c             suivvox=suivvox+(/0,0,-1/)
c          end if

        ENDIF
c    test before leaving cell   
        p1=((x1+t*dir(1))/dx-(i1-1))        
         p2=((y1+t*dir(2))/dy-(j1-1))
        zbottomreal=zsurface(p1,p2,z1bot,z2bot,z3bot,z4bot)
        ztopreal=zsurface(p1,p2,z1top,z2top,z3top,z4top)
          if (z1+t*dir(3).lt.zbottomreal-1e-5.or.z1+t*dir(3).gt.ztopreal+1e-5) then

c             if (suivvox(3).eq.0) nphotwrong=nphotwrong+1
           if (mpi_rank.eq.0) then
           write(6,*),mpi_rank,'wrong cell1 z1,k1,itcount',z1+t*dir(3),k1,itcount2
c           write(6,*),mpi_rank,'tx,ty,tztop,tzbot,t',tx,ty,tztop,tzbot,t
c           write(6,*),mpi_rank,'txold,tyold,tztopold,tzbotold,told',txold,tyold,tztopold,tzbotold,told
c           write(6,*),mpi_rank,'zbottomcell, ztop',zbottomreal,ztopreal
c           write(6,*),mpi_rank,'i0,j0,k0',i0,j0,k0
c           write(6,*),mpi_rank,'i1,j1,k1',i1,j1,k1
c           write(6,*),mpi_rank,'i1old,j1old,k1old',i1old,j1old,k1old
c           write(6,*),mpi_rank,'suivvox',suivvox
c           write(6,*),mpi_rank,'suivvoxold',suivvoxold
c           write(6,*),mpi_rank,'suivvoxold2',suivvoxold2
c           write(6,*),mpi_rank,'xf,yf,zf',xf,yf,zf
c           write(6,*),mpi_rank,'x1,y1,z1',x1,y1,z1
c           write(6,*),mpi_rank,'x1old,y1old,z1old',x1old,y1old,z1old
c           write(6,*),mpi_rank,'dir',dir
c           write(6,*),mpi_rank,'lopt',lopt
c           write(6,*),mpi_rank,'atop,btop,ctop',acoeftop,bcoeftop,ccoeftop
c           write(6,*),mpi_rank,'abot,bbot,cbot',acoefbot,bcoefbot,ccoefbot
c           write(6,*),mpi_rank,'atopold,btopold,ctopold',acoeftopold,bcoeftopold,ccoeftopold
c           write(6,*),mpi_rank,'abotold,bbotold,cbotold',acoefbotold,bcoefbotold,ccoefbotold
c           write(6,*),mpi_rank,'top',z1top,z2top,z3top,z4top
c           write(6,*),mpi_rank,'bot',z1bot,z2bot,z3bot,z4bot
c           write(6,*),mpi_rank,'p1,p2',p1,p2
c           write(6,*),mpi_rank,'reverse',reversetop,reversebot
c           write(6,*),mpi_rank,'tz1,tz2,delta',tz1,tz2,delta
c           stop
          end if
        end if
    
        i1old=i1
        j1old=j1
        k1old=k1
        i1 = i1 + suivvox(1)
        j1 = j1 + suivvox(2)
        k1 = k1 + suivvox(3)
       x1old=x1
       y1old=y1
       z1old=z1
       x1 = x1 + dir(1)*t
       y1 = y1 + dir(2)*t
       if(suivvox(3).eq.1)then
        z1=ztopreal
       elseif(suivvox(3).eq.-1)then
        z1=zbottomreal
       else
        z1 = z1 + dir(3)*t
       endif
c       cyclic boundary conditions implementation
        IF(ibcx .EQ. 1) THEN
         IF(i1 .LT. 1) THEN
          i1 = i1 + n
          x1 = x1 + n*dx
         ENDIF
         IF(i1 .GT. n) THEN
          i1 = i1 - n
          x1 = x1 - n*dx
         ENDIF
        ENDIF
        IF(ibcy .EQ. 1) THEN
         IF(j1 .LT. 1) THEN
          j1 = j1 + m
          y1 = y1 + m*dy
         ENDIF
         IF(j1 .GT. m) THEN
          j1 = j1 - m
          y1 = y1 - m*dy
         ENDIF
        ENDIF
c       end of cyclic bc
         if ((i1>=1).AND.(i1<=n).AND.(j1>=1).AND.(j1<=m).AND.(k1>=1).AND.(k1.le.2).AND.(lopt>0)) then
        ! check for fuel discontiuties in height of cell 1
    
        p1=(x1/dx-(i1-1))
        p2=(y1/dy-(j1-1))
           z1top=zcartedge2AZ(1,i1,j1,k1)
           z2top=zcartedge2AZ(2,i1,j1,k1)
           z3top=zcartedge2AZ(3,i1,j1,k1)
           z4top=zcartedge2AZ(4,i1,j1,k1)
           z1bot=zcartedge2AZ(1,i1,j1,k1-1)
           z2bot=zcartedge2AZ(2,i1,j1,k1-1)
           z3bot=zcartedge2AZ(3,i1,j1,k1-1)
           z4bot=zcartedge2AZ(4,i1,j1,k1-1)
           zbottomreal=zsurface(p1,p2,z1bot,z2bot,z3bot,z4bot)
            ztopreal=zsurface(p1,p2,z1top,z2top,z3top,z4top)
            if ((z1-zbottomreal).lt.(-1.0e-12)) then
                 k1=k1-1
            else if ((z1-ztopreal).gt.(1.0e-12)) then
                   k1=k1+1 
           end if
         end if


c this section is only needed in case of topo
        if (topofile(1:1).ne.'*'.and.topofile.ne.' ') then
         if ((i1>=1).AND.(i1<=n).AND.(j1>=1).AND.(j1<=m).AND.(k1>=1).AND.(k1.le.lmc).AND.(lopt>0)) then
c          if (x1.lt.(i1-1)*dx-1e-8.or.x1.gt.i1*dx+1e-8) then
c            write(6,*),mpi_rank,'wrong cell x1,i1,itcount',x1,i1,itcount2
c             nphotwrong=nphotwrong+1
c           end if
c          x1=min (x1,i1*dx)     ! not necessary usefull
c          x1=max(x1,(i1-1)*dx)
c          if (y1.lt.(j1-1)*dy-1e-8.or.y1.gt.j1*dy+1e-8) then
c            write(6,*),mpi_rank,'wrong cell y1,j1,itcount',y1,j1,itcount2
c             nphotwrong=nphotwrong+1
c           end if
c          y1=min (y1,j1*dy)     ! not necessary useful
c          y1=max(y1,(j1-1)*dy)
        p1=(x1/dx-(i1-1))
        p2=(y1/dy-(j1-1))
c        zold=zsurface((x1old+dir(1)*t)/dx-(i1old-1),
c     +       (y1old+dir(2)*t)/dy-(j1old-1),z1top,z2top,z3top,z4top)

        z1top=zcartedge2AZ(1,i1,j1,k1)
        z2top=zcartedge2AZ(2,i1,j1,k1)
        z3top=zcartedge2AZ(3,i1,j1,k1) 
        z4top=zcartedge2AZ(4,i1,j1,k1) 
        z1bot=zcartedge2AZ(1,i1,j1,k1-1)  
        z2bot=zcartedge2AZ(2,i1,j1,k1-1)
        z3bot=zcartedge2AZ(3,i1,j1,k1-1) 
        z4bot=zcartedge2AZ(4,i1,j1,k1-1)    
        zbottomreal=zsurface(p1,p2,z1bot,z2bot,z3bot,z4bot)
        ztopreal=zsurface(p1,p2,z1top,z2top,z3top,z4top)
          if (z1.lt.zbottomreal-1e-5.or.z1.gt.ztopreal+1e-5) then
             nphotwrong=nphotwrong+1
             ! this can happen 2 photons every millions when suivvox(3)=0
             ! this can happen 5 photons every millions when suivvox(3)<>0
   
c             if (suivvox(3).ne.0) nphotwrong=nphotwrong+1
c           if (mpi_rank.eq.0.) then
c           write(6,*),mpi_rank,'wrong cell2 z1,k1,itcount',z1,k1,itcount2
c           write(6,*),mpi_rank,'tx,ty,tztop,tzbot,t',tx,ty,tztop,tzbot,t
c           write(6,*),mpi_rank,'txold,tyold,tztopold,tzbotold,told',txold,tyold,tztopold,tzbotold,told
c           write(6,*),mpi_rank,'2:txold,tyold,tztopold,tzbotold,told',txold2,tyold2,tztopold2,tzbotold2,told2
c           write(6,*),mpi_rank,'tztopold,tzbotold',tztopold,tzbotold,itcount2
c           write(6,*),mpi_rank,'zbottomcell, ztop',zbottomreal,ztopreal
cc           write (6,*) 'top',z1old+dir(3)*t,zold
cc
c          write(6,*),mpi_rank,'i0,j0,k0',i0,j0,k0
c          write(6,*),mpi_rank,'i1,j1,k1',i1,j1,k1
c           write(6,*),mpi_rank,'i1old,j1old,k1old',i1old,j1old,k1old
c           write(6,*),mpi_rank,'suivvox',suivvox
c           write(6,*),mpi_rank,'suivvoxold',suivvoxold
c           write(6,*),mpi_rank,'suivvoxold2',suivvoxold2
c           write(6,*),mpi_rank,'xf,yf,zf',xf,yf,zf
c           write(6,*),mpi_rank,'x1,y1,z1',x1,y1,z1
c           write(6,*),mpi_rank,'x1old,y1old,z1old',x1old,y1old,z1old
c           write(6,*),mpi_rank,'dir',dir
c           write(6,*),mpi_rank,'lopt',lopt
c           write(6,*),mpi_rank,'atop,btop,ctop',acoeftop,bcoeftop,ccoeftop
c          write(6,*),mpi_rank,'abot,bbot,cbot',acoefbot,bcoefbot,ccoefbot
c           write(6,*),mpi_rank,'atopold,btopold,ctopold',acoeftopold,bcoeftopold,ccoeftopold
c           write(6,*),mpi_rank,'abotold,bbotold,cbotold',acoefbotold,bcoefbotold,ccoefbotold
c           write(6,*),mpi_rank,'top',z1top,z2top,z3top,z4top
c           write(6,*),mpi_rank,'bot',z1bot,z2bot,z3bot,z4bot
c           write(6,*),mpi_rank,'p1,p2',p1,p2

c            z1top=zcartedge2AZ(1,i1old,j1old,k1old)        
c            z2top=zcartedge2AZ(2,i1old,j1old,k1old)
c            z3top=zcartedge2AZ(3,i1old,j1old,k1old)
c            z4top=zcartedge2AZ(4,i1old,j1old,k1old)
c            z1bot=zcartedge2AZ(1,i1old,j1old,k1old-1)
c           z2bot=zcartedge2AZ(2,i1old,j1old,k1old-1)
c           z3bot=zcartedge2AZ(3,i1old,j1old,k1old-1)
c            z4bot=zcartedge2AZ(4,i1old,j1old,k1old-1)
c           write(6,*),mpi_rank,'top',z1top,z2top,z3top,z4top
c           write(6,*),mpi_rank,'bot',z1bot,z2bot,z3bot,z4bot
c           write(6,*),mpi_rank,'reverse',reversetop,reversebot
c           write(6,*),mpi_rank,'tz1,tz2,delta',tz1,tz2,delta
c           stop
c          end if
        end if
c        z1=min(z1,ztopreal)      !correction of photon's position on the frontier
c        z1=max(z1,zbottomreal)
         end if
        end if
c end of checking error in case of topo

         ! here we should have a test if suivvox(1) or suivvox(2) are not 0 and k1 le lphase+1, because the
         ! new k1 can be wrong when fuelheight change from one cell to the other
        END IF ! end of case lopt>=0





      ENDDO

      vox(1) = i1
      vox(2) = j1
      vox(3) = k1
      if((i1.lt.1).or.(i1.gt.n).or.(j1.lt.1).or.(j1.gt.m).or.(k1.lt.0).or.(k1.gt.lmc))then
       if(lopt.le.0)then
        write(6,*) 'i1,j1,k1 =',i1,j1,k1
       endif
      endif
      IF(lopt>0.) THEN         ! out of the domain or nounct2>10
       vox(1) = 0
       vox(2) = 0
       vox(3) = 0
      ENDIF

      END SUBROUTINE CELLPHOTARRIVE_DEFAZ

!!!----------------------- End subroutine CELLPHOTOARRIVE_DEF  -------------------------------!!!

!!!----------------------- subroutine randomdirnorm -------------------------------!!!
      SUBROUTINE randomdirnorm(dir)
      ! cette fonction calculeune direction d'Ã©mission pour une paroi caleii

      REAL Rtheta, sintheta, costheta, phi 
      REAL x, y, z
      REAL dir(3)               ! direction aleatoire
      
      Rtheta = RAND(1.)
      !theta = ASIN(SQRT(theta))
      sintheta = SQRT(Rtheta)
      costheta = SQRT(1.0 - Rtheta)
      phi = RAND(2.*3.1415927)
      x = sintheta * COS(phi)
      y = sintheta * SIN(phi)
      z = costheta

      ! horizontal
      dir(1) = z
      dir(2) = x
      dir(3) = y
      END SUBROUTINE randomdirnorm

!!!----------------------- End subroutine randomdirnorm -------------------------------!!!

!!!----------------------- subroutine randomdirnorm2 -------------------------------!!!
      SUBROUTINE RANDOMDIRNORM2(dir)
      ! cette fonction calcule une direction d'Ã©mission pour une source umique

      REAL Pi, Rtheta, Rphi, sintheta, costheta, phi
      REAL x, y, z
      REAL*8 dir(3)               ! direction aleatoire

      x = 0.0
      y = 0.0
      z = 0.0

      do while((abs(x).lt.1.0e-14).
     &       or.(abs(y).lt.1.0e-14).or.(abs(z).lt.1.0e-14))
      Rtheta = RAND(1.)
      Rphi = RAND(1.)
      !Need to check tp make sure we never get a direction component 
      !which is exactly zero!
      !if(Rtheta.lt.1.0e-14)then
      !   Rtheta = Rtheta+2.5e-10
      !elseif((1.0-Rtheta).lt.1.0e-14)then
      !   Rtheta = Rtheta-2.5e-10
      !endif 
      !if(Rphi.lt.1.0e-14)then
      !   Rphi = Rphi+2.5e-10
      !elseif((1.0-Rphi).lt.1.0e-14)then
      !   Rphi = Rphi-2.5e-10
      !endif 
      !Done with check

      Pi = ACOS(-1.)
      costheta = 1. - 2. * Rtheta
      sintheta = SQRT(1. - costheta**2)
      phi = 2. * pi * Rphi
      x = sintheta * COS(phi)
      y = sintheta * SIN(phi)
      z = costheta
      enddo  
    
      ! horizontal
      dir(1) = x                !z
      dir(2) = y                !x
      dir(3) = z                !y

      END SUBROUTINE randomdirnorm2

!!!----------------------- End subroutine randomdirnorm2 -------------------------------!!!

!!!----------------------- Utility Functions -------------------------------!!!
      FUNCTION rand(p)
      REAL rand
      REAL w, p
      CALL random_number(w)
      rand = w*p
      END FUNCTION rand
                                !
      FUNCTION zsurface(p1,p2,z1,z2,z3,z4)
      REAL*8  p1,p2  !(x-(i-1)dx)/dx,(y-(j-1)dy)/dy
      REAL*8 z1,z2,z3,z4
      REAL*8 zsurface
        zsurface=((z1+z4-z2-z3)*p2+z2-z1)*p1+(z3-z1)*p2+z1
      END FUNCTION zsurface


!!!----------------------- End Utility Functions -------------------------------!!!

!!!*_*_*_*_*_*_*_*_*_*_*_*_*_ Monte Carlo method routines *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*!!!
!!!----------------------- subroutine stretch_keff  -----------------------------------!!!
      subroutine stretch_keff()
      
      use metryic
      implicit none 
 
      integer :: i,j,k,kk
    
      keffx=0.0
      keffy=0.0
      keffz=0.0
      do k=1,LL
       do j=1-ih,mp+ih
        do i=1-ih,np+ih
         kk=k
         if(kk.lt.nz)kk=1
         if(k.ge.nz)kk=k-(nz-1)
         keff(i,j,k)=ckeff/gi0(i,j,k)
         if(i.gt.1-ih)then
          keffx(i,j,k)=0.5*(keff(i,j,k)+keff(i-1,j,k))
          if(i.eq.2-ih)then
           keffx(i-1,j,k)=keffx(i,j,k)
          endif
         endif
         if(j.gt.1-ih)then
          keffy(i,j,k)=0.5*(keff(i,j,k)+keff(i,j-1,k))
          if(j.eq.2-ih)then
           keffy(i,j-1,k)=keffy(i,j,k)
          endif
         endif
         if(k.gt.1)then
         keffz(i,j,k)=0.5*(keff(i,j,k)+keff(i,j,k-1))
         else
         keffz(i,j,k)=keff(i,j,k)
         endif
        enddo !do i
       enddo !do j
      enddo !do k
      !set halo cells for integer exponent vectors 
      do i=1,ih
       if(leftdedge.eq.1)then
        lix(1-i)=i
       endif
       if(rightdedge.eq.1)then
        lix(np+i)=i
       endif
      enddo
      do j=1,ih
       if(botdedge.eq.1)then
        liy(1-j)=j
       endif
       if(topdedge.eq.1)then
        liy(mp+j)=j
       endif
      enddo
      do k=1,ih
       liz(k)=ih-k+1
      enddo

      !stretch keff arrays 
      do k=1,LL
       do j=1-ih,mp+ih
        do i=1-ih,np+ih
         keff(i,j,k)=keff(i,j,k)
         if(i.gt.1-ih)then
          keffx(i,j,k)=keffx(i,j,k)*2/(cst**lix(i-1)+cst**lix(i))
          if(i.eq.2-ih)then
           keffx(i-1,j,k)=keffx(i,j,k)
          endif
         endif
         if(j.gt.1-ih)then
          keffy(i,j,k)=keffy(i,j,k)*2/(cst**liy(j-1)+cst**liy(j))
          if(j.eq.2-ih)then
           keffy(i,j-1,k)=keffy(i,j,k)
          endif
         endif
         if(k.gt.1)then
          keffz(i,j,k)=keffz(i,j,k)*2/(cst**liz(k-1)+cst**liz(k))
         endif
        enddo  !do i 
       enddo  !do j
      enddo  !do k 
      
      end subroutine stretch_keff
!!!-----------------------end subroutine stretch_keff -----------------------------------!!!

!!!*_*_*_*_*_*_*_*_*_*_*_*_*_ Diffusion method routines *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*!!!
!!!---------------------------- subroutine fire_radiation -----------------------------------!!!
c    
c diffusion method. from fireradiationpara1.f
c
      subroutine fire_radiation()

      use metryic
      use fireteca
      use turba
      use xvo

      Implicit None

      
      !real z(nz:Lnz),zs(1-ih:np+ih, 1-ih:mp+ih)
      
      !JAS 3/2/06 added explicits declaration to comply with implicit none
      integer :: i,j,k,kk,nsubdt,kf,kreal !nsub
      !real :: rstep
      real :: sinkterm
      !real :: keffbar,rl2,rl2old,rl22,rl2old2
      real,external :: zcart
      
      !double precision :: wt1,wt2,wtick   !JAS adding timing routine variables


      real surface(1-ih:n+ih,1-ih:m+ih)
      real surfacep(1-ih:np+ih,1-ih:mp+ih)
      real :: srcsolmax,srcgasmax,srctmp
      integer :: ia,ja  !JAS adding global indices

      logical implicit0 !, diril,edge
      !JAS adding mpi wall timers to get some idea of performance gains.
c     if(mpi_rank.eq.0)then
c      wtick = MPI_Wtick()
c      write(6,*) "Wtick = ",wtick
c      wt1=MPI_Wtime()
c      write(6,*) "Wall time 1 in fire_radiation = ",wt1
c     endif

      spval=99.e09
      implicit0=.true.  !Es will be solved by gcrke
      !implicit0=.false.  !Es will be solved using explcit method

      rnetsol = 0.0
      rnetgas = 0.0

      if(implicit0)nsubdt=1

      !define papv, papvgas, papvtot 
      call papv_defs()

      !define sourcegas, sourcesol
      call source_defs()

      call updated2(sourcegas,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(sourcesol,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
           !****************************** Ef energy field ******************************!
      ! setup src and sink terms for Ef virtual energy field      
      sinksol=0.0
      sinkgas=0.0
      srcsol=0.0
      srcgas=0.0
      do k=1,LL
       do j=1-ih,mp+ih
        do i=1-ih,np+ih
         kk=k
         if(k.lt.nz)kk=nz
         if(k.gt.Lnz)kk=Lnz
         kf=kfuel(i,j)
         if(k.gt.kf)srcgas(i,j,k)=sourcegas(i,j,k)
        enddo
       enddo
      enddo
      call updated2(srcgas,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(srcsol,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(sinkgas,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(sinksol,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)

       srcsolmax = MAXVAL(ABS(srcsol))
       srcgasmax = MAXVAL(ABS(srcgas))
       srctmp = 0.0
       call mpi_allreduce(srcsolmax,srctmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
       srcsolmax = srctmp
       srctmp = 0.0
       call mpi_allreduce(srcgasmax,srctmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
       srcgasmax = srctmp


         !PRINT*, 'MAXVAL(srcsol), MAXVAL(srcgas), mpi_rank',
      !+        MAXVAL(ABS(srcsol)), MAXVAL(ABS(srcgas)),mpi_rank

      if((srcsolmax.gt.0.0).or.(srcgasmax.gt.0.0))then
      if(implicit0)then
       call gcrke(E,5,LL)   !(E,5,20,LL)
       call lape(E,diffterm,LL)
      else
      !Solve for the Ef field explicitly
      call explicitE(E,LL)
      endif
      endif 

      !****************************** Ef energy field ******************************!
      ! setup src and sink terms for Ef virtual energy field      
      sinksol=0.0
      sinkgas=0.0
      srcsol=0.0
      srcgas=0.0
      do k=1,LL
       do j=1-ih,mp+ih
        do i=1-ih,np+ih
         kk=k
         if(k.lt.nz)kk=nz
         if(k.gt.Lnz)kk=Lnz
         kf=kfuel(i,j)
         if(k.gt.kf)srcgas(i,j,k)=sourcegas(i,j,k)
         if(k.le.kf)then
          sinksol(i,j,k)=papv(i,j,kk)*
     &                   (papvgas(i,j,kk)+papv(i,j,kk))*keff(i,j,k)
          sinkgas(i,j,k)=papvgas(i,j,kk)*
     &                   (papvgas(i,j,kk)+papv(i,j,kk))*keff(i,j,k)
         endif
        enddo
       enddo
      enddo
      call updated2(srcgas,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(srcsol,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(sinkgas,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(sinksol,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)

       srcsolmax = MAXVAL(ABS(srcsol))
       srcgasmax = MAXVAL(ABS(srcgas))
       srctmp = 0.0
       call mpi_allreduce(srcsolmax,srctmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
       srcsolmax = srctmp
       srctmp = 0.0
       call mpi_allreduce(srcgasmax,srctmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
       srcgasmax = srctmp


         !PRINT*, 'MAXVAL(srcsol), MAXVAL(srcgas), mpi_rank',
      !+        MAXVAL(ABS(srcsol)), MAXVAL(ABS(srcgas)),mpi_rank

      if((srcsolmax.gt.0.0).or.(srcgasmax.gt.0.0))then
      ! set cartesian z at fuel top
      sinkterm=0.
      if(nuemannbcs.eq.1)then
      do j=1,mp
       do i=1,np
        kf=kfuel(i,j)
         ! JMC surfacep is the top of the cell with non-zero density
         surfacep(i,j)=zcart(zedge(kf-(nz-1)+1),i,j)
       enddo
      enddo
      call mpi_allgather(surfacep(1:np,1:mp),np*mp,mpi_real,
     +                   surface(1:n,1:m),np*mp,mpi_real,
     +                   mpi_comm_world,ierror)
      call updated(surface,surface,n,m,1,1-ih,n+ih,1-ih,m+ih,1,0)

      ! mpi broadcast and compute cartesian fluxes
      call mpi_fluxes(surface,srcgas)
      !call mpi_fluxes(zs,surface,srcgas,sinkterm)


      hx=0.0
      hy=0.0
      hz=0.0
      do j=1-ih,mp+ih
       do i=1-ih,np+ih        ! Neumann BC at fuel top
        ia=(npos-1)*np+i !JAS introduced ia and ja to deal with xyzflux covering whole domain
        ja=(mpos-1)*mp+j

      ! JMC gi appears to be calculated at the cell center. It may need to be
      !     changed to reflect that hz is a cell face quantity.
        kf=kfuel(i,j)
        hz(i,j,kf+1)=xflux(ia,ja)*c13(i,j)*gmul0(kf+1)+
     +               yflux(ia,ja)*c23(i,j)*gmul0(kf+1)+
     +               zflux(ia,ja)*gi0(i,j,kf+1)
        hz(i,j,kf+1)=hz(i,j,kf+1)/gi0(i,j,kf+1)
        if(i.gt.1-ih)then
         hx(i,j,kf+1)=0.5*(xflux(ia,ja)+xflux(ia-1,ja))/(0.5*(gi0(i-1,j,kf+1)+gi0(i,j,kf+1)))
        else
         hx(i,j,kf+1)=xflux(ia,ja)/gi0(i,j,kf+1)
        endif
        if(j.gt.1-ih)then
         hy(i,j,kf+1)=0.5*(yflux(ia,ja)+yflux(ia,ja-j3))/(0.5*(gi0(i,j-j3,kf+1)+gi0(i,j,kf+1)))
        else
         hy(i,j,kf+1)=yflux(ia,ja)/gi0(i,j,kf+1)
        endif
       enddo
      enddo
      else    !using dirichlet bcs from E
      do k=1,LL
      do j=1-ih,mp+ih
       do i=1-ih,np+ih        ! Neumann BC at fuel top
        ia=(npos-1)*np+i !JAS introduced ia and ja to deal with xyzflux covering whole domain
        ja=(mpos-1)*mp+j
      
      ! JMC gi appears to be calculated at the cell center. It may need to be
      !     changed to reflect that hz is a cell face quantity.
        kf=kfuel(i,j)
        if(k.eq.kf+1)then
         Ef(i,j,k)=E(i,j,k)
        endif
       enddo
      enddo
      enddo
       srcsol = 0.0
       srcgas = 0.0
      endif !if neumannbcs .ne.1


      if(implicit0)then
       call gcrke(Ef,5,kf) !  (Ef,5,20,kf)
       call lape(Ef,diffterm,kf)
      else
      !Solve for the Ef field explicitly
       call explicitE(Ef,kf)
      endif

      call rmaxmin1(Ef,'Ef',1-ih,np+ih,1-ih,mp+ih,LL)
      !calculate net radiative change to gas or solids due to radiation
      do k=1,kf
       do j=1-ih,mp+ih
        do i=1-ih,np+ih
         if(nz.le.k.and.k.le.Lnz)then
          rnetgas(i,j,k)=rnetgas(i,j,k)+
     +               (srcgas(i,j,k)-sinkgas(i,j,k)*gi0(i,j,k)*Ef(i,j,k))
          rnetsol(i,j,k)=rnetsol(i,j,k)+
     +               (srcsol(i,j,k)-sinksol(i,j,k)*gi0(i,j,k)*Ef(i,j,k))
         endif !
        enddo !do i
       enddo !do j
      enddo !do k
      !******************************End  Ef energy field ******************************!
      endif  !srcsolmax or srcgasmax.gt.0.0
     
      !******************************  Es energy field ******************************!
      !setup src and sink terms for Es virtual energy field
      sinksol=0.0
      sinkgas=0.0
      srcsol=0.0
      srcgas=0.0
      do k=1,LL
       do j=1,mp
        do i=1,np
         kk=k
         if(k.lt.nz)kk=nz
         if(k.gt.Lnz)kk=Lnz
         kf=kfuel(i,j)
         sinksol(i,j,k)=papv(i,j,kk)*
     &                  (papvgas(i,j,kk)+papv(i,j,kk))*keff(i,j,k)
         sinkgas(i,j,k)=papvgas(i,j,kk)*
     &                  (papvgas(i,j,kk)+papv(i,j,kk))*keff(i,j,k)
         if(k.le.kf)srcgas(i,j,k)=sourcegas(i,j,k)
         if(k.le.kf)srcsol(i,j,k)=sourcesol(i,j,k)
        enddo
       enddo
      enddo

      call updated2(srcgas,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(srcsol,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(sinkgas,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(sinksol,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
       srcsolmax = MAXVAL(ABS(srcsol))
       srcgasmax = MAXVAL(ABS(srcgas))
       srctmp = 0.0
       call mpi_allreduce(srcsolmax,srctmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
       srcsolmax = srctmp
       srctmp = 0.0
       call mpi_allreduce(srcgasmax,srctmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
       srcgasmax = srctmp


         !PRINT*, 'MAXVAL(srcsol), MAXVAL(srcgas), mpi_rank',
      !+        MAXVAL(ABS(srcsol)), MAXVAL(ABS(srcgas)),mpi_rank

      if((srcsolmax.gt.0.0).or.(srcgasmax.gt.0.0))then

      if(implicit0)then
       call gcrke(Ef,5,kf) !  (Ef,5,20,kf)
       call lape(Es,diffterm,LL)
      else
       !Solve for the Ef field explicitly
       call explicitE(Es,LL)
      endif

      call rmaxmin1(Es,'Es',1-ih,np+ih,1-ih,mp+ih,l)
      !calculate net radiative change to gas or solids due to radiation
      do k=1,LL
       do j=1-ih,mp+ih
        do i=1-ih,np+ih
         if(nz.le.k.and.k.le.Lnz)then
          rnetgas(i,j,k)=rnetgas(i,j,k)+
     +               (srcgas(i,j,k)-sinkgas(i,j,k)*gi0(i,j,k)*Es(i,j,k))
          rnetsol(i,j,k)=rnetsol(i,j,k)+
     +               (srcsol(i,j,k)-sinksol(i,j,k)*gi0(i,j,k)*Es(i,j,k))
         endif
        enddo
       enddo
      enddo
      !******************************End  Es energy field ******************************!
      endif  !srcsolmax or srcgasmax.gt.0.0
      !update  solid and gas fields for rate of change due to radiative energy  exchange
      !call rmaxmin1(rnetsol, 'rnetsol', 1-ih, np+ih, 1-ih, mp+ih, l)
      !call rmaxmin1(rnetgas, 'rnetgas', 1-ih, np+ih, 1-ih, mp+ih, l)
      do k=nz,Lnz
       !adjust k for  real domain 
       kreal = k-nz+1
       do j=1,mp
        do i=1,np
         firad(i,j,kreal)=firad(i,j,kreal)-rnetgas(i,j,k)
         frhosiesrad(i,j,kreal)=frhosiesrad(i,j,kreal)-rnetsol(i,j,k)
        enddo !do i
       enddo !do j
      enddo !do k

      !JAS adding mpi wall timers to get some idea of performance gains.
c     if(mpi_rank.eq.0)then
c      wt2=MPI_Wtime()
c      write(6,*) "Wall time 2 in fire_radiation = ",wt2
c      write(6,*) "Total time  in fire_radiation = ",(wt2-wt1)
c     endif

      end subroutine fire_radiation
!!!------------------------end subroutine fire_radiation ----------------------------------------!!!


!!!---------------------------- subroutine explicitE ----------------------------------------!!!
      subroutine explicitE(Ee,Ltop)
      
      implicit none


      !JMC variable declaration
      integer :: i,j,k,it,Ltop
      real :: rl2,rl2old,eps,subdt,tmp
      real :: Ee(1-ih:np+ih,1-ih:mp+ih,LL)
      real :: rdt,rdzg,rdy0 

      eps=1.0e-6 !JMC convergence criteria
      rl2=3.0e38 !JMC initialize norm to be > eps
      rl2old=rl2+1.0 !JMC initialize old norm to be > norm
      it=0 !JMC counter for iterations
      
      !determine subintervals and sub-time step
      rdt=1.0/dt
      !rdzg=dzi/gimin
       rdzg=dzi*gimax
      rdy0=dyi
      if(j3.eq.0)rdy0=0.0

      subdt=0.49/((ckeff/gimin)* (dxi**2+rdy0**2+rdzg**2
     &        - (papvmax*(papvmax+papvgasmax)
     &        + papvgasmax*(papvmax+papvgasmax)) ) )
      write(6,*)'subdt=',subdt
      if(subdt.lt.0)then
       write(6,*)'In radiation, stability criteria failed for
     & explicit case, program terminated!'
       call mpi_finalize()
      endif

      do while(rl2.gt.eps)   !JMC Iterate until Ef converges
       it=it+1 !JMC count iterations
      
       !JMC get the residual from diffusion
       call lape(E,diffterm,Ltop)
      
       rl2old = rl2 !JMC store the norm from the last iteration
       rl2 = 0.0    !JMC initialize norm
       do k=1,Ltop
        do j=1,mp
         do i=1,np
          rl2 = rl2+(diffterm(i,j,k)-srcgas(i,j,k)-srcsol(i,j,k))**2  !L2 norm
         enddo !i
        enddo !j
       enddo !k
       tmp = 0.0
       call mpi_allreduce(rl2,tmp,1,mpi_real,mpi_sum,
     +                   mpi_comm_world,ierror)
       rl2 = tmp
       if(rl2.gt.rl2old)then
        write(6,*)'In explicitE: E is diverging,rl2,rl2old',rl2,rl2old
        exit
       endif

       call updated2(Ee,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
       do k=2,Ltop
        do j=1-ih+1,mp+ih-1
         do i=1-ih+1,np+ih-1
          !JMC update E
          Ee(i,j,k)=Ee(i,j,k)+subdt*(srcgas(i,j,k)+srcsol(i,j,k)+diffterm(i,j,k))
         enddo !i
        enddo !j
       enddo !k

      if(mod(it,10).eq.0)
     & write(6,*) 'In explicitE: rl2,iterations = ',rl2,it
      enddo !do while
       write(6,*) 'In explicitE: rl2,iterations = ',rl2,it
      return
      end subroutine explicitE
!!!---------------------------end subroutine explicitE----------------------------------------!!!

!!!---------------------------- subroutine lape ----------------------------------------!!!
      subroutine lape(Elp,r,Ltop)
      
      use metryic
      Implicit None
      
      real Elp(1-ih:np+ih,1-ih:mp+ih,LL),
     .     r(1-ih:np+ih,1-ih:mp+ih,LL) !,
      !logical Neumann,dirilicht,defined
      !JAS 3/2/06 added explicit declarations to comply with implicit none
      integer :: i,j,k, Ltop ! ,itest,im,jm,ihb,jhb,klpmax^M
      !real :: emx,emy,emz,epx,epy,epz^M
      real :: exbar,eybar,ezbar,coef,coef2,coef3 ! ee, det^M
      real :: g13,g23,g33,gii! ,gii2,ginv11,ginv12,ginv13,ginv22,ginv23,ginv33^M
      ! real :: hxflux,hyflux,hzflux^M
      !real :: hxz,hyzinteger :: i,j,k,itest,im,jm,ihb,jhb,Ltop,klpmax
      real :: xterm,yterm,zterm 
      
      !JMC initialize Laplacian residual to zero for Neumann BC
      r = 0.0
      Ex = 0.0
      Ey = 0.0
      Ez = 0.0
      !JMC Neumann BC on left edge for non-cyclic case
      if(leftdedge.eq.1)then
        Elp(1-ih,:,:)=0.0
      endif
      !JMC Neumann BC on right edge for non-cyclic case
      if(rightdedge.eq.1)then
        Elp(np+ih,:,:)=0.0
      endif
      !JMC Neumann BC on bottom edge for non-cyclic case
      if(botdedge.eq.1)then
          Elp(:,1-ih,:)=0.0
      endif
      !JMC Neumann BC on top edge for non-cyclic case
      if(topdedge.eq.1)then
        Elp(:,mp+ih,:)=0.0
      endif
      !JMC Neumann BC on the floor and ceiling of the domain
          !E(i,j,1)=2*E(i,j,2)-E(i,j,3)
           Elp(:,:,1)=0.0
      if(Ltop.eq.LL)then  !Es field, set dirichlet bc at domain ceiling
       do j=1-ih,mp+ih
        do i=1-ih,np+ih
         Elp(i,j,Ltop)=2*Elp(i,j,Ltop-1)-Elp(i,j,Ltop-2)
        enddo !do i
       enddo !do j
      elseif(nuemannbcs.eq.0)then
       !do j=1-ih,mp+ih
        !do i=1-ih,np+ih
         !Elp(i,j,Ltop)=0.5*(Elp(i,j,Ltop+1)+Elp(i,j,Ltop-1))
        !enddo !do i
       !enddo !do j
      endif

      !compute covariant derivatives in grid domain (capped by top of fuel bed)
      do k=1,Ltop
       do j=1-ih,mp+ih
        do i=1-ih,np+ih
         if(i.gt.1-ih)
     &   Ex(i,j,k)=dxi*(Elp(i,j,k)-Elp(i-1,j,k))
         if(j.gt.1-ih)
     &   Ey(i,j,k)=dyi*(Elp(i,j,k)-Elp(i,j-j3,k))
         if(k.ne.Ltop)then
          if(k.gt.1)then
           Ez(i,j,k)=dzi*(Elp(i,j,k)-Elp(i,j,k-1))
          endif
         else  !k.eq.Ltop
          if(Ltop.eq.LL)then
           !working on the Es field, so set constant z gradient at ceiling
           Ez(i,j,k)=dzi*(Elp(i,j,k)-Elp(i,j,k-1))
           !Ex(i,j,k+1)=Ex(i,j,k)
           !Ey(i,j,k+1)=Ey(i,j,k)
           Ez(i,j,k+1)=Ez(i,j,k)
          else
         !working on the Ef field, so back calculate Ex,Ey,Ez(i,j,k+1) from prescribed h(k+1) arrays
           if(nuemannbcs.eq.1)then
           if(i.gt.1-ih)then
             g13=0.5*(c13(i,j)+c13(i-1,j))*gmul0(k+1)
           else
             g13=c13(i,j)*gmul0(k+1)
           endif
           if(j.gt.1-ih)then
             g23=0.5*(c23(i,j)+c23(i,j-j3))*gmul0(k+1)
           else
             g23=c23(i,j)*gmul0(k+1)
           endif
           coef2=keffx(i,j,k+1)
           coef3=keffy(i,j,k+1)
           gii=0.5*(gi0(i,j,k)+gi0(i,j,k+1))
           g33=g13**2+g23**2+gii**2
           Ez(i,j,k)=dzi*(Elp(i,j,k)-Elp(i,j,k-1))
           coef=keffz(i,j,k+1)
           Ez(i,j,k+1)=(hz(i,j,k+1)/(-coef)
     &                 -g13*hx(i,j,k+1)/(-coef2)
     &                 -g23*hy(i,j,k+1)/(-coef3))
     &                 /(g33-g13**2-g23**2)
            Ex(i,j,k+1)=hx(i,j,k+1)/(-coef2)-g13*Ez(i,j,k+1)
            Ey(i,j,k+1)=hy(i,j,k+1)/(-coef3)-g23*Ez(i,j,k+1)
           else
            Ez(i,j,k+1)=dzi*(Elp(i,j,k+1)-Elp(i,j,k))
            Ez(i,j,k)=Ez(i,j,k+1)
           endif !nuemannbcs.eq.1
          endif !Ltop.eq.LL
         endif !k.ne.Ltop
        enddo !do i
       enddo !do j
      enddo !do k

      !compute contravariant x-flux
      do k=2,Ltop
       do j=1-ih,mp+ih
        do i=1-ih+1,np+ih
         g13=0.5*(c13(i,j)+c13(i-1,j))*gmul0(k)
         ezbar=0.25*(Ez(i,j,k+1)+Ez(i,j,k)+Ez(i-1,j,k+1)+Ez(i-1,j,k))
         coef=keffx(i,j,k)
         hx(i,j,k)=-coef*(Ex(i,j,k)+g13*ezbar)
        enddo !do i
       enddo !do j
      enddo !do k

      !compute contravariant y-flux
      do k=2,Ltop
       do j=1-ih+1,mp+ih
        do i=1-ih,np+ih
         g23=0.5*(c23(i,j)+c23(i,j-j3))*gmul0(k)
         ezbar=0.25*(Ez(i,j,k+1)+Ez(i,j,k)+Ez(i,j-j3,k+1)+Ez(i,j-j3,k))
         coef=keffy(i,j,k)
         hy(i,j,k)=-coef*(Ey(i,j,k)+g23*ezbar)
        enddo !do i
       enddo !do j
      enddo !do k

      !compute contravariant z-flux
      do k=2,Ltop
       do j=1-ih,mp+ih-1
        do i=1-ih,np+ih-1
         g13=0.5*c13(i,j)*(gmul0(k)+gmul0(k-1))
         g23=0.5*c23(i,j)*(gmul0(k)+gmul0(k-1))
         gii=0.5*(gi0(i,j,k)+gi0(i,j,k-1))
         g33=g13**2+g23**2+gii**2
         exbar=0.25*(Ex(i,j,k)+Ex(i+1,j,k)+Ex(i,j,k-1)+Ex(i+1,j,k-1))
         eybar=0.25*(Ey(i,j,k)+Ey(i,j+j3,k)+Ey(i,j,k-1)+Ey(i,j+j3,k-1))
         coef=keffz(i,j,k)
         hz(i,j,k)=-coef*(g13*exbar+g23*eybar+g33*Ez(i,j,k))
        enddo !do i
       enddo !do j 
      enddo !do k

      !compute laplacian and RHS for defined areas in the grid
      do k=2,Ltop
       do j=1-ih+1,mp+ih-1
        do i=1-ih+1,np+ih-1
         yterm=0.
         if(j3.eq.1)then
          yterm=dyi*(hy(i,j+j3,k)-hy(i,j,k))/(cst**liy(j))
         endif
         xterm = dxi*(hx(i+1,j,k)-hx(i,j,k))/(cst**lix(i))
         zterm = dzi*(hz(i,j,k+1)-hz(i,j,k))/(cst**liz(k)) 
         r(i,j,k)=-1.*(xterm+yterm+zterm)   
     &          - (sinkgas(i,j,k)+sinksol(i,j,k))*Elp(i,j,k)
         r(i,j,k)=r(i,j,k)*gi0(i,j,k)
        enddo !do i
       enddo !do j
      enddo !do k

      return
      end subroutine lape
!!!---------------------------end  subroutine lape ----------------------------------------!!!


!!!---------------------------- subroutine gcrke ----------------------------------------!!!
      subroutine gcrke(p,itmn,Ltop)  !(x,x, itmx, x)

      Implicit None

      real p(1-ih:np+ih,1-ih:mp+ih,LL)

      integer,parameter :: lord=3

      real r(1-ih:np+ih,1-ih:mp+ih,LL),
     &     ar(1-ih:np+ih,1-ih:mp+ih,LL),
     &     x(1-ih:np+ih,1-ih:mp+ih,LL,lord),
     &     ax(1-ih:np+ih,1-ih:mp+ih,LL,lord),
     &     ax2(lord),axar(lord),del(lord)


      integer :: i,j,k,ii,i0,it,itmn,itr,nlc,Ltop  !itmx
      real :: eps,eps0,epa,eer0,eem0,rl20,rax,beta,dvmx,rl2,tmp !eer,eem

      itr=6000/lord
      eps0=1.e-4
      eps=eps0*dti
      epa=1.e-30
      nlc=0
      r=0.0
      ar=0.0
      x=0.0
      ax=0.0

      !get initial r from initial p
      call lape(p,r,Ltop)
      call updated2(p,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(r,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      eer0=0.
      eem0=-1.e15
      rl20=0.
      do k=1,Ltop
       do j=1-ih+1,mp+ih-1
        do i0=1-ih+1,np+ih-1
         r(i0,j,k)=r(i0,j,k)+(srcgas(i0,j,k)+srcsol(i0,j,k))
         if((i0.ge.1).and.(j.ge.1).and.(i0.le.np).and.(j.le.mp))then
          eer0=eer0+r(i0,j,k)**2
          eem0=amax1(eem0,abs(r(i0,j,k)))
          rl20=rl20+r(i0,j,k)**2
         endif
        enddo !do i0
       enddo !do j
      enddo !do k
      call mpi_allreduce(rl20,tmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
      rl20=tmp
      call mpi_allreduce(eer0,tmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
      eer0=tmp
      call mpi_allreduce(eem0,tmp,1,mpi_real,mpi_max,mpi_comm_world,ierror)
      eem0=tmp
      eer0=amax1(eer0,epa)
      eem0=amax1(eem0,epa)
      rl20=amax1(rl20,epa)

      do k=1,Ltop
       do j=1-ih,mp+ih
        do i0=1-ih,np+ih
         x(i0,j,k,1)=r(i0,j,k)      ! (initial p0=r0 in notes)
        enddo !do i0
       enddo ! do j
      enddo !do k
      call lape(x,ax,Ltop)
      call updated2(x(1-ih,1-ih,1,1),np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      call updated2(ax(1-ih,1-ih,1,1),np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
      do 100 it=1,itr
       do i=1,lord
        ax2(i)=0.
        rax=0.
        do k=1,Ltop
         do j=1,mp
          do i0=1,np
           rax=rax+r(i0,j,k)*ax(i0,j,k,i)
           ax2(i)=ax2(i)+ax(i0,j,k,i)*ax(i0,j,k,i)
          enddo !do i0
         enddo !do j
        enddo !do k
        call mpi_allreduce(rax,tmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
        rax=tmp
        call mpi_allreduce(ax2(i),tmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
        ax2(i)=tmp
        ax2(i)=amax1(epa,ax2(i))
        beta=-rax/ax2(i)
        dvmx=-1.e15
        rl2=0.
        do k=1,Ltop
         do j=1-ih+1,mp+ih-1
          do i0=1-ih+1,np+ih-1
           p(i0,j,k)=p(i0,j,k)+beta* x(i0,j,k,i)
           r(i0,j,k)=r(i0,j,k)+beta*ax(i0,j,k,i)
           if((i0.ge.1).and.(j.ge.1).and.(i0.le.np).and.(j.le.mp))then
            dvmx=amax1(dvmx,abs(r(i0,j,k)))
            rl2=rl2+r(i0,j,k)*r(i0,j,k)
           endif
          enddo !do i0
         enddo !do j
        enddo !do k
        call updated2(p,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
        call updated2(r,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
        tmp = 0.0
        call mpi_allreduce(rl2,tmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
        rl2=tmp
        tmp = 0.0
        call mpi_allreduce(dvmx,tmp,1,mpi_real,mpi_max,mpi_comm_world,ierror)
        dvmx=tmp
        if(dvmx.le.eps.and.it.ge.itmn)then
         if(mpi_rank.eq.0)then
         write(6,*)'gcrke has converged,dvmx,eps,rl2,rl20,it' 
         write(6,*)dvmx,eps,rl2,rl20,it
         endif
         go to 200
        endif
        if(rl2.ge.rl20.and.it.ge.itmn)then
         if(mpi_rank.eq.0)then
         write(6,*)'gcrke is diverging,dvmx,eps,rl2,rl20, it' 
         write(6,*)dvmx,eps,rl2,rl20,it
         endif
         go to 200
        endif
        if(it.eq.itr)then
         if(mpi_rank.eq.0)then
         write(6,*)"gcrke didn't converge,dvmx,eps,rl2,rl20,it" 
         write(6,*)dvmx,eps,rl2,rl20,it
         endif
         go to 200  ! failed to converge
        endif
        rl20=amax1(rl2,epa)
        call lape(r,ar,Ltop)
        call updated2(r,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)
        call updated2(ar,np,mp,LL,1-ih,np+ih,1-ih,mp+ih,1)

        nlc=nlc+1
        do ii=1,i
         axar(ii)=0.
         do k=1,Ltop
          do j=1,mp
           do i0=1,np
            axar(ii)=axar(ii)+ax(i0,j,k,ii)*ar(i0,j,k)
           enddo !do i0
          enddo !do j
         enddo !do k
         call mpi_allreduce(axar(ii),tmp,1,mpi_real,mpi_sum,mpi_comm_world,ierror)
         axar(ii)=tmp
         del(ii)=-axar(ii)/ax2(ii)
        enddo !do ii
        if(i.lt.lord) then
         do k=1,Ltop
          do j=1-ih,mp+ih
           do i0=1-ih,np+ih
            x(i0,j,k,i+1)=r(i0,j,k)
            ax(i0,j,k,i+1)=ar(i0,j,k)
           enddo !do i0
          enddo !do j
         enddo !do k
         do ii=1,i
          do k=1,Ltop
           do j=1-ih,mp+ih
            do i0=1-ih,np+ih
             x(i0,j,k,i+1)= x(i0,j,k,i+1)+del(ii)* x(i0,j,k,ii)
             ax(i0,j,k,i+1)=ax(i0,j,k,i+1)+del(ii)*ax(i0,j,k,ii)
            enddo !do i0
           enddo !do j
          enddo !do k
         enddo !do ii
        else
         do k=1,Ltop
          do j=1-ih,mp+ih
           do i0=1-ih,np+ih
            x(i0,j,k,1)=r(i0,j,k)+del(1)* x(i0,j,k,1)
            ax(i0,j,k,1)=ar(i0,j,k)+del(1)*ax(i0,j,k,1)
           enddo !do i0
          enddo !do j
         enddo !do k
         do ii=2,i
          do k=1,Ltop
           do j=1-ih,mp+ih
            do i0=1-ih,np+ih
             x(i0,j,k,1 )= x(i0,j,k,1)+del(ii)* x(i0,j,k,ii)
             ax(i0,j,k,1 )=ax(i0,j,k,1)+del(ii)*ax(i0,j,k,ii)
             x(i0,j,k,ii)=0.
             ax(i0,j,k,ii)=0.
            enddo !do i0
           enddo !do j
          enddo !do k
         enddo !do ii
        endif !if(i.lt.lord)
       enddo !do i=1,lord
  100 continue   !enddo for it
  200 continue

      return
      end subroutine gcrke
!!!---------------------------end subroutine gcrke --------------------------------------!!!

!!!----------------------- subroutine mpi_fluxes -----------------------------------!!!
      subroutine mpi_fluxes(surface,src)  !(zs, x ,x, sinkterm)

!      use constants, only: rsourceht
      use metryic, only: x,y !zedge 
      Implicit None


      !JAS 3/6/06 added explicit declarations to comply wityh implicit none
      !real :: sinkterm
      integer :: i,j,k,ia,ja,ii,jj !iproc,npos0,mpos0
      !real :: zz,zcutoff,delx,dely,delz,rmag,rsquared,rcubed,expterm
      real :: srcterm,term1

      !real zs(1-ih:np+ih,1-ih:mp+ih)
      real src(1-ih:np+ih,1-ih:mp+ih,LL)
      real surface(1-ih:n+ih,1-ih:m+ih)
      !real tmpflux(1-ih:n+ih,1-ih:m+ih)
      real tmpxflux(1-ih:n+ih,1-ih:m+ih)
      real tmpyflux(1-ih:n+ih,1-ih:m+ih)
      real tmpzflux(1-ih:n+ih,1-ih:m+ih)
      !real zcar(np,mp,L)
      !real surfterm

      !JMC 6/11/7 variables used for cyclic bc
      integer :: ixdom,iydom,nxdom,nydom,numdomx,numdomy
      real :: rfalloff
      real :: xdomlen
      real :: ydomlen

      !JAS variables used for new source integral calculation
      real :: xt,yt,zt,xsl,xsu,ysl,ysu,zsl,zsu,srcintx,srcinty,srcintz
      real :: dxs,dys,a,b,c,d,f,g,tmp !signcd,signab
      real :: a0,b0,c0,d0,f0,g0
      integer :: nsubdx,nsubdy,iis,jjs,jsrc,isrc,fflag,abflag,cdflag

      tmpxflux = 0
      tmpyflux = 0
      tmpzflux = 0
      xflux = 0
      yflux = 0
      zflux = 0
      fflag = 1
      !JMC 6/11/7 Determine the number of extra domains used for cyclic bc
      rfalloff=200.0 !JMC 6/11/7 radius for negligible energy contribution
      numdomx=int(max(rfalloff/real(n*dx)+0.5,0.0))!number of x-domains
      numdomy=int(max(rfalloff/real(m*dy)+0.5,0.0))!number of y-domains
      nxdom=numdomx*ibcx!use extra domain space when cyclic x is on
      nydom=numdomy*ibcy!use extra domain space when cyclic y is on

      !adding high resolution loop variables for calculation of fluxes
      nsubdx = 21
      nsubdy = 21
      dxs = dx/nsubdx
      dys = dy/nsubdy

      do k=1,L    ! determine flux contributions on each process
       do j=1,mp
        ja=(mpos-1)*mp+j
        ysl=y(ja)-dy
        ysu=y(ja)
        do i=1,np
         ia=(npos-1)*np+i
         xsl=x(ia)-dx
         xsu=x(ia)
         zsl = zcoordsf(i,j,k)
         zsu = zcoordsf(i,j,k+1)

         srcterm=src(i,j,k+(nz-1))

         if(srcterm.gt.1e1)then
          jsrc = ja
          isrc = ia
          !******JMC 6/11/7 loop through extra domain space for cyclic*****!
          xdomlen=float(n)*dx
          ydomlen=float(m)*dy
          do iydom=(-1)*nydom,nydom
           do ixdom=(-1)*nxdom,nxdom
            do jj=1-ih,m+ih
             do jjs=1,nsubdy
              yt = y(jj)-dy+dys*(jjs-1)+0.5*dys
              c0 = (yt-ysu)
              d0 = (yt-ysl)
              c = abs(c0)
              d = abs(d0)
              if(d.gt.c)then
               c = abs(d0)
               d = abs(c0)
              endif
              if((d0*c0).lt.0.0)then
               cdflag = 1
              else
               cdflag = 0
              endif
              do ii=1-ih,n+ih     ! scan over surface
               zt = surface(ii,jj)-1e-6
               f0 = (zt-zsu)
               g0 = (zt-zsl)
               f = f0 !abs(f0)
               g =  g0 !abs(g0)
               if(g.gt.f)then
                !f = abs(g0)
                !g = abs(f0)
               endif
               do iis=1,nsubdx
                xt = x(ii)-dx+dxs*(iis-1)+0.5*dxs
                a0 = (xt-xsu)
                b0 = (xt-xsl)
                a = abs(a0)
                b = abs(b0)
                if(b.gt.a)then
                 a = abs(b0)
                 b = abs(a0)
                endif
                if((b0*a0).lt.0.0)then
                 abflag = 1
                else
                 abflag = 0
                endif

                call srcint(srcintx,f,g,c,d,a,b)
                call srcint(srcinty,f,g,a,b,c,d)
                if((cdflag.eq.0).and.(abflag.eq.0))then
                  call srcint(srcintz,a,b,c,d,f,g)
                elseif((cdflag.eq.0).and.(abflag.eq.1))then
                 tmp = 0.0
                 call srcint(srcintz,a,0.0,c,d,f,g)
                 tmp = srcintz
                 call srcint(srcintz,b,0.0,c,d,f,g)
                 srcintz = srcintz+tmp
                elseif((cdflag.eq.1).and.(abflag.eq.0))then
                 tmp = 0.0
                 call srcint(srcintz,a,b,c,0.0,f,g)
                 tmp = srcintz
                 call srcint(srcintz,a,b,d,0.0,f,g)
                 srcintz = srcintz+tmp
                else      !(cdflag.eq.1).and.(abflag.eq.1)
                  call srcint(srcintz,a0,b0,c0,d0,f0,g0)
                endif
                term1=srcterm/(4.0*3.14159)
                tmpxflux(ii,jj)=tmpxflux(ii,jj)+term1*srcintx/((nsubdx)*(nsubdy))
                tmpyflux(ii,jj)=tmpyflux(ii,jj)+term1*srcinty/((nsubdx)*(nsubdy))
                tmpzflux(ii,jj)=tmpzflux(ii,jj)+term1*srcintz/((nsubdx)*(nsubdy))
               enddo !do iis
              enddo !do ii
             enddo !do jjs
            enddo !do jj
           enddo !do ixdom
          enddo !do iydom
          !******JMC 6/11/7 loop through extra domain space for cyclic*****!
         endif !srcterm.gt.1e1
        enddo !do ia
       enddo !do ja
      enddo !do k

            call mpi_allreduce(tmpxflux,xflux,(n+2*ih)*(m+2*ih)
     &                      ,mpi_real,mpi_sum,mpi_comm_world,ierror)
            call mpi_allreduce(tmpyflux,yflux,(n+2*ih)*(m+2*ih)
     &                      ,mpi_real,mpi_sum,mpi_comm_world,ierror)
            call mpi_allreduce(tmpzflux,zflux,(n+2*ih)*(m+2*ih)
     &                      ,mpi_real,mpi_sum,mpi_comm_world,ierror)

      return
      end subroutine mpi_fluxes
!!!-----------------------end subroutine mpi_fluxes -----------------------------------!!!

!!!----------------------- subroutine srcint -----------------------------------!!!
      subroutine srcint(retval,a,b,c,d,f,g)

      real :: retval,a,b,c,d,f,g !,a56,b56,c34,d34
      !real :: sa,sb,sc,sd,sf,sg
      real :: term31,term32,term33,term34
      real :: term41,term42,term43,term44
      real :: term51,term52,term53,term54
      real :: term61,term62,term63,term64
      real :: term1,term2,term3,term4,term5,term6 !,tmps

      term1 = g*( atan(c*a/(g*sqrt(a**2+c**2+g**2)))
     &           -atan(c*b/(g*sqrt(b**2+c**2+g**2)))
     &           -atan(d*a/(g*sqrt(a**2+d**2+g**2)))
     &           +atan(d*b/(g*sqrt(b**2+d**2+g**2))) )
      term2 = f*( atan(d*a/(f*sqrt(a**2+d**2+f**2)))
     &           -atan(d*b/(f*sqrt(b**2+d**2+f**2)))
     &           -atan(c*a/(f*sqrt(a**2+c**2+f**2)))
     &           +atan(c*b/(f*sqrt(b**2+c**2+f**2))) )

      term31 =  c+sqrt(a**2+c**2+f**2)
      term32 =  c+sqrt(a**2+c**2+g**2)
      term33 =  d+sqrt(a**2+d**2+g**2)
      term34 =  d+sqrt(a**2+d**2+f**2)
      term3 = ((term31*term33)/(term32*term34))**a
      term41 =  c+sqrt(b**2+c**2+g**2)
      term42 =  c+sqrt(b**2+c**2+f**2)
      term43 =  d+sqrt(b**2+d**2+f**2)
      term44 =  d+sqrt(b**2+d**2+g**2)
      term4 = ((term41*term43)/(term42*term44))**b
      term51 =  a+sqrt(a**2+c**2+f**2)
      term52 =  a+sqrt(a**2+c**2+g**2)
      term53 =  b+sqrt(b**2+c**2+g**2)
      term54 =  b+sqrt(b**2+c**2+f**2)
      term5 = ((term51*term53)/(term52*term54))**c
      term61 =  a+sqrt(a**2+d**2+g**2)
      term62 =  a+sqrt(a**2+d**2+f**2)
      term63 =  b+sqrt(b**2+d**2+f**2)
      term64 =  b+sqrt(b**2+d**2+g**2)
      term6 = ((term61*term63)/(term62*term64))**d

      retval = term1+term2+log(term3)+log(term4)+log(term5)+log(term6)

      return
      end subroutine srcint
!!!-----------------------end subroutine srcint -----------------------------------!!!
!!!*_*_*_*_*_*_*_*_*_*_*_*_*_ End Diffusion method routines *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*!!!
      subroutine writeiorad(data,iunit,il,iu,jl,ju,nzdim)
      use gridsetup
      use msga

      Implicit None

      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,nzdim,iunit
c
      real data(il:iu,jl:ju,nzdim)
      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
      real :: chtemp

c
      real,allocatable::tmparray(:,:,:,:),outdata(:,:,:)
c
      if (mpi_rank.eq.0)
     +  allocate(tmparray(il:iu,jl:ju,nzdim,nproc),
     +            outdata(n,m,nzdim))
        nsize=(iu-il+1)*(ju-jl+1)*nzdim
        call mpi_gather(data,nsize,mpi_real,tmparray,nsize,mpi_real,
     +     0,mpi_comm_world,ierror)
      if (mpi_rank.eq.0) then
        do iprocx=1,nprocx
          do jprocy=1,nprocy
           iproc=1+(iprocx-1)+(jprocy-1)*nprocx
           do k=1,nzdim
             do j=1,mp
               do i=1,np
                 ia=(iprocx-1)*np + i
                 ja=(jprocy-1)*mp + j
                 outdata(ia,ja,k)=tmparray(i,j,k,iproc)
               enddo
             enddo
           enddo
         enddo
        enddo

        chtemp=sum(outdata)
        write (iunit)outdata
c
        deallocate(tmparray,outdata)
      endif
c
      return
      end subroutine writeiorad

      end module radiation
c2345678***************************************************
