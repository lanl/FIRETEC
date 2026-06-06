!-------------------------------------------------------------------------

!© 2026. Triad National Security, LLC. All rights reserved. This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

!-------------------------------------------------------------------------
          program compress

      use gridsetup
      use xvo
      use xve
      use xvi
      use workavg
      use turba
      use fireteca
      use weights
      use constants
      use metryic
      use restarta
      use msga
      use higrad   
      use advectionScheme
      use ignite
      use windfield !JMC
      use pres      !JMC
      use xvbin     !JMC
      use filesub   !JMC
      use radiation !, only: irad,icallrad ! KOO  
      use io
      use firebranda ! KOO
      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: nplwrites,i,j,k,it!kv 
      !integer :: idot,iend,it,ia,ja ! inv
      integer :: itname
      real :: tcv
      !real :: totalburnenergy,totalburnmass,xalocal,yalocal
      !real :: xalocal,yalocal
      real,external :: zcart
      integer::clength
       ! real:: sp2
      !JAS not used --common/rlookup/ ifa,ifb,ifc,ifd
      character(len=257)::pfilename ,fname
      integer :: ierr 
      integer :: iv ! FIXME
 
c     character end*10,fname*60
      data nplwrites /0/ !JAS not used -- nreswrites /0/
c
       namelist/compresslist/
     .irst,nt,nts,
     .isplit,nprocx,nprocy,n,m,l,aa1,nv,dx,dy,dz,dts,ntp,
     .ih,nr,ibcx,ibcy,irlx,irly,ibclatopen,ibctopopen,iab,
     .zab,zabt,tow,itheta,
     .irod,iturb,idrag,isa,ilapdo,rturbprandtl,ffparam,
     .irad,icallrad,crad,irandseed,iseed,
     .isootmodel,iradeastflux,
     .inonlocal,isubgridgas,  !FP
     .irhovapor,icorio,ilspgf,izlspgf,frqlspgf,islip,
     . tambient, pressground, zgroundref,iperturb,
     .u0,uramp,uramptime,uswitch,zu,
     .v0,vramp,vramptime,vswitch,
     .ius,iue,jus,jue,
     .slopeangle,slopeazimuth,
      !.st,
     .iord,isor,nonos,idiv,nfct,nonosold,
     .ifuelinra,fuelinranumber,ifuel,ivegread,
     .rhomicrovalue, cpwood,
     .kmax,kf77max,frqoutput,outname,frqfilstr,
     .restartfile,ioextra,topofile,ipotflow,
     .isoturb,ibctopbot,ignfile,icfmeflag,
     & windfieldstartfile,xvbdataname,
     +iwindfieldout,iwindfieldin,windspeedupfactor,itwindfield, !JMC 10/7/6
     &itinterp,ibcells,jbcells,  !JMC 10/7/6
     .is,ie,js,je,
     .iunstable,
     .iwallclock,                              !KOO 07/16/07
     .ifbrand,                                 !KOO
     .ifp



c     ! list of personal parameter for fp
       namelist/listfp/
     .iheatsource,hsros,hsint,ihsmass,lHeatSource,
     .outnamesub,frqoutputsub,issub,jssub,nisub,njsub,nksub  
       ! parameter for subdomain ouput 

! MMC 8/25/06  initialize MPI here rather than in rinitmsg
      call MPI_Init(ierror)
      call MPI_Comm_rank(mpi_comm_world,mpi_rank,ierror)
      call MPI_Comm_size(mpi_comm_world,numprocs,ierror)
c
! MMC 8/25/06; Read gridlist only by the rank=0 process, then MPI_BCAST gridlist variables
!       to the other processes. This allows the code to run on clusters where 
!       processes other than the rank=0 process may not have direct access to
!       gridlist, for example, if they don't all share the same PWD environment
!       variable.
!
      if (mpi_rank==0) then
        aa1=0.1  ! default value
        open(unit=15,file='gridlist',form='formatted',status='old')
        read (15,nml=compresslist)
        close(15)
        ! here fp open his personal set of parameters in file gridfp
        if (ifp.ge.1) then
          open(unit=115,file='gridlistfp',form='formatted',status='old')
          read (115,nml=listfp)
          close(115)
        endif
      endif
!mmc 6/28/06 added to prevent mixups with inonlocal and irhovapor, and nv in gridlist
      if (inonlocal==1) then    
        nv=nv+1
        irhovapor=0
      endif
!
! FIXME 
      iheatsource=0  

! MMC 8/25/06; Here we MPI_BCAST all the gridlist variables from the rank=0
!      process to all the other processes. If you modify gridlist, 
!      make sure that you also modify this series of MPI_BCAST calls
!      to reflect your changes to gridlist and ensure that all
!      gridlist variables are passed to the other processes:
!
      call MPI_Bcast(irst,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nt,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nts,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(isplit,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nprocx,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nprocy,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(n,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(m,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(l,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(aa1,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(nv,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(dx,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(dy,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(dz,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(dts,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(ntp,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ih,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nr,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ibcx,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ibcy,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(irlx,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(irly,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ibclatopen,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ibctopopen,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iab,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(zab,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(zabt,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(tow,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(irod,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iturb,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(idrag,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(isa,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ffparam,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(rturbprandtl,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(ilapdo,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(irad,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(icallrad,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(irandseed,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iseed,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(isootmodel,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(crad,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(iradeastflux,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(inonlocal,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(isubgridgas,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(irhovapor,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(tambient,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(pressground,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(zgroundref,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(itheta,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iperturb,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(icorio,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ilspgf,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(izlspgf,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(frqlspgf,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(islip,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(u0,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(uramp,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(uramptime,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(zu,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(uswitch,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(v0,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(vramp,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(vramptime,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(vswitch,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ius,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iue,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(jus,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(jue,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(slopeangle,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(slopeazimuth,1,mpi_real,0,mpi_comm_world,ierr)
      !call MPI_Bcast(st,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(iord,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(isor,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nonos,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(idiv,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nfct,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nonosold,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ifuel,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ivegread,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(cpwood,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(rhomicrovalue,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(ifuelinra,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(fuelinranumber,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(kmax,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(kf77max,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(frqoutput,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(frqfilstr,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iunstable,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ifbrand,1,mpi_integer,0,mpi_comm_world,ierr)  ! KOO
      if (mpi_rank==0) clength=len(outname)  !get length of outname so MPI_Bcast knows how long this character string is
      call MPI_Bcast(clength,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(outname,clength,mpi_character,0,mpi_comm_world,ierr)
      
      call MPI_Bcast(ioextra,1,mpi_integer,0,mpi_comm_world,ierr)
      
      if (mpi_rank==0) clength=len(restartfile)   !get length of restartfile so MPI_Bcast knows how long this character string is
      call MPI_Bcast(clength,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(restartfile,clength,mpi_character,0,mpi_comm_world,ierr)
      
      if (mpi_rank==0) clength=len(topofile)   !get length of topofile so MPI_Bcast knows how long this character string is
      call MPI_Bcast(clength,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(topofile,clength,mpi_character,0,mpi_comm_world,ierr)

      call MPI_Bcast(ipotflow,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(isoturb,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ibctopbot,1,mpi_integer,0,mpi_comm_world,ierr)
      
      if (mpi_rank==0) clength=len(ignfile)   !get length of ignfile so MPI_Bcast knows how long this character string is
      call MPI_Bcast(clength,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ignfile,clength,mpi_character,0,mpi_comm_world,ierr)
      
      call MPI_Bcast(icfmeflag,1,mpi_integer,0,mpi_comm_world,ierr)

!JMC added flags to gridlist for reading windfields
      if (mpi_rank==0) clength=len(windfieldstartfile)   
       !get length of windfieldstart so MPI_Bcast knows how long this character string is
      call MPI_Bcast(clength,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(windfieldstartfile,clength,mpi_character,0,mpi_comm_world,ierr)
      if (mpi_rank==0) clength=len(xvbdataname)   !get length of windfieldstart so MPI_Bcast knows how long this character string is
      call MPI_Bcast(clength,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(xvbdataname,clength,mpi_character,0,mpi_comm_world,ierr)
      call MPI_Bcast(iwindfieldout,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iwindfieldin,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(windspeedupfactor,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(itwindfield,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(itinterp,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ibcells,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(jbcells,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(is,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ie,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(js,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(je,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iwallclock,1,mpi_integer,0,mpi_comm_world,ierr)
! personal flag and personal parameters
      call MPI_Bcast(ifp,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iheatsource,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(hsros,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(hsint,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(lHeatSource,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(ihsmass,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(frqoutputsub,1,mpi_integer,0,mpi_comm_world,ierr)
      if (mpi_rank==0) clength=len(outnamesub)  !get length of outname
      !so MPI_Bcast knows how long this character string is
      call MPI_Bcast(clength,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(outnamesub,clength,mpi_character,0,mpi_comm_world,ierr)
      call MPI_Bcast(issub,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(jssub,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nisub,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(njsub,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nksub,1,mpi_integer,0,mpi_comm_world,ierr)

!JMC added flags to gridlist for reading windfields
      call MPI_BARRIER(mpi_comm_world,ierror)
! ------ Finished with initial MPI_BCAST of gridlist variables _____
c walltimer !KOO
      if(iwallclock.EQ.1) wtime0=MPI_Wtime()   
      call setup 
      call con
      call rinitmsg
      call definearray
      call rinit(tcv)
      if(ifbrand.eq.1) call rinitfbrand(irst,itrestart,nbmax,nb_imm) !KOO
      !totalburnenergy = 0
      !totalburnmass = 0
    
      itname=itrestart
      restarttime=itrestart*dts*nts
      
      !prepare to dump simulation parameters
      if(mpi_rank.eq.0)then
       pfilename = 'glist'
       call namefile(itrestart,pfilename,fname)
       open(unit=15,file=fname,form='formatted',status='unknown')
       write (15,nml=compresslist)
       close(15)
      endif
      
      if(iwallclock.EQ.1) call computeWallTime(wtime0,10,'init')
      
      do 10 it=1,nt
        if(iwallclock.EQ.1) wtime1=MPI_Wtime() 
        ittot=it+itrestart !FP
        ! if xe should evolve with time:  (windfieldin...)
        if (mpi_rank.eq.0) write(6,*)'it total= ',ittot,' it=', it
!        if (mpi_rank.eq.0) write (6,'(a12,f12.3)')'time(secs)=',time 

        if (ixevariation.eq.1) call xevariation()
       call rmaxmin(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)
       call rmaxmin1(temps,'temps',1-ih,np+ih,1-ih,mp+ih,l) 
c**************************************************************

      do k=1,l
         do j=1,mp
            do i=1,np

            xv(i,j,k,1)=xvb(i,j,k,1)
            xv(i,j,k,2)=xvb(i,j,k,2)
            xv(i,j,k,3)=xvb(i,j,k,3)
            xv(i,j,k,4)=xvb(i,j,k,4)
            xv(i,j,k,5)=xvb(i,j,k,nv)

!          if(it.eq.1) then 
!             open(7209,file='xvb.history',form='formatted',
!     +                 status='unknown') 
!             write(7209,182) it,xvb(i,j,k,1),xvb(i,j,k,2),xvb(i,j,k,3)
!     +            ,xvb(i,j,k,4),xvb(i,j,k,5),xvb(i,j,k,6)
!     +            ,xvb(i,j,k,7),xvb(i,j,k,8),pr(i,j,k)  
!          else ! it.eq.1
!             if(dts.eq.0.001) then 
!               if(MOD(it+1,2).eq.0) then  
!             open(7209,file='xvb.history',form='formatted',
!     +                 status='unknown',position='append') 
!             write(7209,182) it,xvb(i,j,k,1),xvb(i,j,k,2),xvb(i,j,k,3)
!     +            ,xvb(i,j,k,4),xvb(i,j,k,5),xvb(i,j,k,6)
!     +            ,xvb(i,j,k,7),xvb(i,j,k,8),pr(i,j,k)
!               endif   
!
!             else if(dts.eq.0.0005) then
!               if(MOD(it+3,4).eq.0) then
!             open(7209,file='xvb.history',form='formatted',
!     +                 status='unknown',position='append')
!             write(7209,182) it,xvb(i,j,k,1),xvb(i,j,k,2),xvb(i,j,k,3)
!     +            ,xvb(i,j,k,4),xvb(i,j,k,5),xvb(i,j,k,6)
!     +            ,xvb(i,j,k,7),xvb(i,j,k,8),pr(i,j,k)
!               endif
!       
!             else 
!             open(7209,file='xvb.history',form='formatted',
!     +                 status='unknown',position='append')         
!             write(7209,182) it,xvb(i,j,k,1),xvb(i,j,k,2),xvb(i,j,k,3)
!     +            ,xvb(i,j,k,4),xvb(i,j,k,5),xvb(i,j,k,6)
!     +            ,xvb(i,j,k,7),xvb(i,j,k,8),pr(i,j,k)
!          endif           
!
!182     format(I3,F9.5,F9.5,F9.5,F9.4,F9.5,F9.6,F9.6,F9.5,F13.2)  
!          endif 

          if(xvb(i,j,k,7)/xvb(i,j,k,8).GT.0.22
     +     .or. xvb(i,j,k,4)/xvb(i,j,k,8).LT.270) then

          print*,'ijkmpi',i,j,k,mpi_rank,it
          print*,'xvb123',xvb(i,j,k,1),xvb(i,j,k,2),xvb(i,j,k,3)
          print*,'xvb456',xvb(i,j,k,4),xvb(i,j,k,5),xvb(i,j,k,6)
          print*,'xvb78p',xvb(i,j,k,7),xvb(i,j,k,8),xvb(i,j,k,7)/xvb(i,j,k,8)
          print*,'temp',temps(i,j,k),tempg(i,j,k),xvb(i,j,k,4)/xvb(i,j,k,8),pr(i,j,k)

          STOP 

          endif 

          do iv=1, nv 
             if(xvb(i,j,k,iv).NE.xvb(i,j,k,iv) .or.
     +          xvb(i,j,k,iv)+1.EQ.xvb(i,j,k,iv)) then
           print*,'NaN',i,j,k,iv,xvb(i,j,k,iv)
           STOP 
             endif 
          enddo          


          enddo
        enddo
      enddo


c walltimer - JMC1 !KOO
      if(iwallclock.EQ.1) call computeWallTime(wtime1,2,'JMC1')
      if(iwallclock.EQ.1) wtime4=MPI_Wtime() 
      call rmoainner()

      if(iwallclock.EQ.1) call computeWallTime(wtime4,3,'rmoainner')
      if (mpi_rank.eq.0) write (6,'(a12,f12.3)')'time(secs)=',time

      call advec(xvb,1-ih,np+ih,1-ih,mp+ih,l)

      if(iwallclock.EQ.1) call computeWallTime(wtime2,5,'advec')
      if(iturb.ge.1) then
      do k=1,l
      do j=1,mp
      do i=1,np
      if (xvb(i,j,k,5).lt.rkmin) then
c        if(it.gt.1)write (*,*) 'rkmin violation 5, rkmin=',xvb(i,j,k,5)
         xvb(i,j,k,5)=rkmin
      endif
      if (xvb(i,j,k,6).lt.rkmin) then
c        if(it.gt.1)write (*,*) 'rkmin violation 6, rkmin=',xvb(i,j,k,6)
         xvb(i,j,k,6)=rkmin
      endif
      if (xvb(i,j,k,5).gt.rkmax) then
c        if(it.gt.1)write (*,*) 'rkmax violation 5, rkmax=',xvb(i,j,k,5)
         xvb(i,j,k,5)=rkmax
      endif
      if (xvb(i,j,k,6).gt.rkmax) then 
c        if(it.gt.1)write (*,*) 'rkmax violation 6, rkmax=',xvb(i,j,k,6)
         xvb(i,j,k,6)=rkmax
      endif
      enddo
      enddo
      enddo

      call updated(xvb(1-ih,1-ih,1,5),xe(1-ih,1-ih,1,5),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      call updated(xvb(1-ih,1-ih,1,6),xe(1-ih,1-ih,1,6),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      endif


      if(iwallclock.EQ.1) wtime2=MPI_Wtime()   
      call force

      if(iwallclock.EQ.1) call computeWallTime(wtime2,6,'force')
      call boundary(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)
      
      if(iwallclock.EQ.1) wtime2=MPI_Wtime()
      if(iturb.ge.1) then
         ! compute firetec and turb arrays with bc: u,v,w...
         call fieldUpdate()

         call turb()

         if(iwallclock.EQ.1) call computeWallTime(wtime2,7,'turb')
      endif

c complex burncode
      if(iwallclock.EQ.1) wtime4=MPI_Wtime()
      if(irod.eq.1) call firetec(it)
      if(iwallclock.EQ.1) call computeWallTime(wtime4,8,'firetec')
      if(ifbrand.EQ.1) call firebrand(it,itrestart)    !KOO

c call boundary conditions
      call boundary(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)

c  JMC: call a subroutine to write a wind profile.
c      if((it+itrestart).eq.10000)then
c
c        if(mpi_rank.eq.0)write(*,*)'going into icfmewind'
c        call icfmewind(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)
c        if(mpi_rank.eq.0)write(*,*)'exiting icfmewind'
c      endif
      if (iwallclock.EQ.1)wtime2=MPI_Wtime()

       if (ilspgf.ge.1.and.it/frqlspgf*frqlspgf.eq.it) call largeScalePGF(it)

          if (iwindfieldout.eq.1) then
            if(ittot.ge.itwindfield.and.mod(ittot,itinterp).eq.0) then
              itname=(ittot-itwindfield)
              call namefile(itname,xvbdataname,fxvbdataname)
              call windfld_frqwriteio(fxvbdataname)
            endif
            if(ittot.eq.itwindfield)then
              call frqwriteio(windfieldstartfile)
            endif
          endif
          if (frqfilstr>0) ! do filter every 10 time step
           ! to limit computational cost of filtering
!                +     if(it/10*10.eq.it) call filstr()
     +     call filstr()
          if(it/frqoutput*frqoutput.eq.it) then
            if (nplwrites.eq.0) then
              itname=it+itrestart
              call namefile(itname,outname,fname)
            endif
            if(mpi_rank.eq.0)write(6,*)trim(fname)
            call frqwriteio(fname)
          endif
          ! subdomain outputs
          if (frqoutputsub.gt.0)  then  ! if 0 no subdomain
          if(it/frqoutputsub*frqoutputsub.eq.it) then
             if (nplwrites.eq.0) then
               itname=it+itrestart
               call namefile(itname,outnamesub,fname)
             endif
             if(mpi_rank.eq.0)write(6,*)trim(fname)
             call frqwritesubio(fname)
          endif
          endif

c walltimer - IO and iteration !KOO 
      if(iwallclock.EQ.1) then
           call computeWallTime(wtime2,9,'io  ')
           call computeWallTime(wtime1,1,'it  ')
c           write(6,*) "cputimes at it=",it,",total cputimes"
c           write(6,*) "iteration/total:",cputime(1,1),cputime(1,2)
c           write(6,*) "JMC1:",cputime(2,1),cputime(2,2)
c           write(6,*) "rmoainner:",cputime(3,1),cputime(3,2)
c           write(6,*) "advvel:",cputime(4,1),cputime(4,2)
c           write(6,*) "advec:",cputime(5,1),cputime(5,2)
c           write(6,*) "force:",cputime(6,1),cputime(6,2)
c           write(6,*) "turb:",cputime(7,1),cputime(7,2)
c           write(6,*) "firetec:",cputime(8,1),cputime(8,2)
c           write(6,*) "I/O",cputime(9,1),cputime(9,2)
c      write(6,*) "mpdatanew3d in advec:",cputime(11,1),cputime(11,2)
c     write(6,*) "init. for firetec:",cputime(12,1),cputime(12,2)
c     write(6,*) "lapdo in firetec:",cputime(13,1),cputime(13,2)
c      write(6,*) "radiation in firetec:",cputime(14,1),cputime(14,2)
c      write(6,*) "fuel in firetec:",cputime(15,1),cputime(15,2)
c      write(6,*) "convec in firetec:",cputime(16,1),cputime(16,2)
c      write(6,*) "updates in firetec:",cputime(17,1),cputime(17,2)
       if (mpi_rank.eq.0) then  
      open (109,file='cputimes',form='formatted',status='unknown')
      write (109,11) it,cputime(1,1),cputime(2,1),cputime(3,1),
     +       cputime(4,1),cputime(5,1),cputime(6,1),cputime(7,1),
     +       cputime(8,1),cputime(9,1),cputime(11,1),cputime(12,1),
     +       cputime(13,1),cputime(14,1),cputime(15,2),cputime(16,2),
     +       cputime(17,2)
      write (6,*) " it#,  iterat, JMC1  , rmoain, advvel, advec
     + , force , turb  , fitec , I/O"
      write (6,12) it,cputime(1,1),cputime(2,1),cputime(3,1),
     +       cputime(4,1),cputime(5,1),cputime(6,1),cputime(7,1),
     +       cputime(8,1),cputime(9,1)
      write (6,*) "       mpdata, ftinit, lapdo , ftrad , ftfuel, convec
     +, ftupdt"
      write (6,13) cputime(11,1),cputime(12,1),
     +       cputime(13,1),cputime(14,1),cputime(15,2),cputime(16,2),
     +       cputime(17,2)

11    format(I5,F9.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4,
     +       F8.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4)

12    format(I5,F9.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4)

13    format(F14.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4)
        endif
      endif

10    continue


      if(iwallclock.EQ.1) then
        wtime2=MPI_Wtime()
        call mpi_reduce(wtime2-wtime0,wtime3,1,mpi_double_precision,
     +                  mpi_sum,0,mpi_comm_world,ierr)
        if(mpi_rank.EQ.0) then
          write(6,*) "total cputime",wtime3,"of iterations,",it
          close(109)
        endif
      endif

      call mpi_finalize(ierror)
      stop
      end

