!-------------------------------------------------------------------------

!© 2026. Triad National Security, LLC. All rights reserved. This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

!-------------------------------------------------------------------------

!----------------------------------------------------------------
! Primary driver program to run HIGRAD/FIRETEC
!----------------------------------------------------------------
program compress

  ! Global variables
  use gridlist_variables
  use gridsetup, only : itrestart,np,mp,ittot,time
  use msga_variables, only : mpi_rank,mpi_comm_world,ierror
  use xvall, only : ixevariation,xv,nv,xv_list,xe,relaxxv
  use higrad, only : rmoainner
  Implicit None

  ! Local variables
  integer :: it
  character(len=257):: fname
  namelist/compresslist/ &
    irst,nt,nts,dts,ntp, &
    n,m,l,dx,dy,dz,aa1, &
    nprocx,nprocy,ih, &
    nr,ibcx,ibcy,ibclatopen,ibctopopen,iab,islip, &
    zab,zabt,tow,itheta, &
    ifire,iturb,isa,rturbprandtl, &
    irad,icallrad,crad,irandseed,iseed, &
    isootmodel,iradeastflux, &
    inonlocal,irhovapor,iemissions,idiffsies,ifbrand, &
    icorio,ilspgf,izlspgf,frqlspgf, & 
    tambient,relativeHumidity,pressground,zgroundref,iperturb, &
    u0,uramp,uramptime,uswitch,zu, &
    v0,vramp,vramptime,vswitch, &
    ius,iue,jus,jue, &
    iord,nonos,idiv,nfct,nonosold, &
    ifuel,nfuel,ivegread, &
    rhomicro,cpwood,Water2WoodRatio, &
    frqoutput,outname,frqfilstr, &
    restartfile,topofile,ipotflow, &
    ignfile,ignVertExtent,icfmeflag, &
    windfieldstartfile,xvdataname, &
    iwindfieldout,iwindfieldin,windspeedupfactor,itwindfield, &
    itinterp,ibcells,jbcells, &
    is,ie,js,je

  ! Executable code
  
  ! Read gridlist
  open(unit=15,file='gridlist',form='formatted',status='old')
  read(15,nml=compresslist)
  close(15)

  ! Initialize Simulation
  call rinitmsg
  call xvpointers
  call definearray
  call rinit
   
  ! Dump simulation parameters
  if(mpi_rank.eq.0)then
    call namefile(itrestart,'glist',fname)
    open(unit=15,file=fname,form='formatted',status='unknown')
    write(15,nml=compresslist)
    close(15)
  endif

  it=1
  do while(ittot.lt.nt) 
    ittot=it+itrestart
    if(mpi_rank.eq.0)then
      print*,'it total=',ittot,' it=',it
      write(*,'(a12,f12.3)')'time(sec)=',time
    endif
    call rmaxmin(xv,xv_list,1-ih,np+ih,1-ih,mp+ih,l,nv)
    if(isnan(xv(1,1,1,1)))then
      print*,'NaNs in XV array'
      call mpi_finalize(ierror)
      STOP
    endif
    if(ixevariation.eq.1) call xevariation
    call rmoainner
    call advec(xv,1-ih,np+ih,1-ih,mp+ih,l)
    call forcing
    if(iturb.ge.1)then 
      call fieldUpdate()
      call turb
      if(ifire.eq.1)then 
        call firetec
        if(ifbrand.eq.1) call firebrand
      endif
    endif
    call boundary(xv,xe,relaxxv,1-ih,np+ih,1-ih,mp+ih,l,nv)
    if(ilspgf.ge.1.and.mod(ittot,frqlspgf).eq.0) call largeScalePGF
    if(iwindfieldout.eq.1)then
      if(ittot.ge.itwindfield.and.mod(ittot,itinterp).eq.0)then
        call namefile(ittot-itwindfield,xvdataname,fname)
        call windfld_write(fname)
      endif
      if(ittot.eq.itwindfield) call frqwriteio(windfieldstartfile)
    endif
    if(mod(it,frqoutput).eq.0)then
      call namefile(ittot,outname,fname)
      if(mpi_rank.eq.0) print*,trim(fname)
      call frqwriteio(fname)
    endif
    it=it+1
  enddo

  call mpi_finalize(ierror)

end program compress
