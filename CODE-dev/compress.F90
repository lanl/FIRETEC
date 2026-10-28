!-----------------------------------------------------------------------
! Primary driver program to run HIGRAD/FIRETEC
!-----------------------------------------------------------------------
program compress

  ! Global variables
  use gridlist_variables, only : nt,l,windfieldout,frqoutput, &
    isensor,ifbrand,dt
  use vegetation_gridlist_variables, only : nfuel
  use gridsetup, only : itrestart,ittot,time
  use msga_variables, only : mpi_rank,ierror
  use fuel_variables, only : temps
  use thermo_variables, only : pr,tempg
  use xvall, only : xv,xvrho,nv,xv_list,xvfuel,xvfuel_list,nvfuel, &
    isies
  use sensor_gridlist_variables, only : se_frq_write
  use rmaxmin_mod, only : rmaxmin
  Implicit None

  ! Local variables
  integer :: kv,it,ift

  ! Executable code
  ! Read user inputs
  call input 

  ! Initialize MPI
  call setup_MPI
  call definearray
  call rinit
  if(itrestart.eq.0) call binary_writeio(0)

  it=1
  do while(ittot.lt.nt) 
    ittot=it+itrestart
    if(mpi_rank.eq.0)then
      print*,'it total=',ittot,' it=',it
      write(*,'(a12,f12.3)')'time(sec)=',time
    endif
    
    call explicit_domain_updates ! Updates environmental variables (pressure, etc.)
    
    ! Write to the standard out state variables
    do kv=1,nv-1
      call rmaxmin(xvrho(:,:,1:l,kv),xv_list(kv))
    enddo
    call rmaxmin(xv(:,:,1:l,nv),xv_list(nv))
    call rmaxmin(tempg,'Tempg')
    call rmaxmin(pr,'pr')
    do ift=1,nfuel
      if(mpi_rank.eq.0) print*,'Fuel Type',ift
      do kv=1,nvfuel
        call rmaxmin(xvfuel(ift,:,:,:,kv),xvfuel_list(kv))
      enddo
      if(isies.ne.1) then
        call rmaxmin(temps(ift,:,:,:),'Temps')
      endif
    enddo

    if(isnan(xvrho(1,1,1,1)).or.abs(xvrho(1,1,1,1)).gt.1e15) then
      print*,'NaNs or Infinity in XV array'
      exit 
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
    time=time+dt
    if(windfieldout.eq.1) call windField(ittot)
    if((isensor.eq.1).and.(mod(ittot,se_frq_write).eq.0))then
      call sensors_write
    endif
    if(mod(ittot,frqoutput).eq.0)then
      call binary_writeio(ittot)
      if(ifbrand.ge.1) call writeio_fb(ittot)
    endif
    it=it+1
  enddo

  call mpi_finalize(ierror)

end program compress
