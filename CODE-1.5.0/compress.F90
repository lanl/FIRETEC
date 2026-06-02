!-----------------------------------------------------------------------
! Primary driver program to run HIGRAD/FIRETEC
!-----------------------------------------------------------------------
program compress

  ! Global variables
  use gridlist_variables, only : nt,ih,l,nfuel,windfieldout,frqoutput, &
    outname,isensor,ifbrand
  use gridsetup, only : itrestart,np,mp,ittot,time
  use msga_variables, only : mpi_rank,ierror
  use fuel_variables, only : lfuel
  use xvall, only : xvrho,nv,xv_list,xvfuel,xvfuel_list,nvfuel
  use sensor_gridlist_variables, only : se_frq_write
  Implicit None

  ! Local variables
  integer :: it,ift
  character(len=257) :: fname

  ! Executable code
  ! Read user inputs
  call input 

  ! Initialize Simulation
  call setup_MPI
  call definearray
  call rinit

  it=1
  do while(ittot.lt.nt) 
    ittot=it+itrestart
    if(mpi_rank.eq.0)then
      print*,'it total=',ittot,' it=',it
      write(*,'(a12,f12.3)')'time(sec)=',time
    endif
    
    call explicit_domain_updates ! Updates environmental variables (pressure, etc.)
    call rmaxmin_3D(xvrho(:,:,1:l,:),xv_list,1-ih,np+ih,1-ih,mp+ih,l,nv)
    do ift=1,nfuel
      if(mpi_rank.eq.0) print*,'Fuel Type',ift
      call rmaxmin_3D(xvfuel(ift,:,:,:,:),xvfuel_list,1,np,1,mp, &
        lfuel,nvfuel)
    enddo
    if(isnan(xvrho(1,1,1,1))) then
      print*,'NaNs in XV array'
      call mpi_finalize(ierror)
      STOP
    endif


    call small_explicit_forcings ! Calculates explicit large timestep forcings

!    if(implicit_methods.eq.1) then
      call MOA ! Mehod of Averaging for semi-implicit small timestep forcings
!    elseif(implicit_methods.eq.2) then
!      call Runge_Kutta ! Runge-Kutta for semi-implicit small timesteop forcings
!    endif

    call large_explicit_forcings

    call forcing ! Applies forcings, post-forcing effects, and advance timestep
    
    !!! IO
    if(windfieldout.eq.1) call windField(ittot)
    if((isensor.eq.1).and.(mod(it,se_frq_write).eq.0)) then
      call sensors_write
    endif
    if(mod(it,frqoutput).eq.0)then
      call namefile(ittot,outname,fname)
      if(mpi_rank.eq.0) print*,trim(fname)
      call frqwriteio(fname)
      if(ifbrand.ge.1) call writeio_fb(ittot)
    endif
    it=it+1
  enddo

  call mpi_finalize(ierror)

end program compress
