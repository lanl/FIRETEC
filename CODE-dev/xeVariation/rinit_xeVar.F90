!-----------------------------------------------------------------------
! rinitXeVariation initializes all arrays and variables related to 
! temporal evolutions in the environmental variables
!-----------------------------------------------------------------------
subroutine rinitXeVariation(xevariation)
  use gridlist_variables, only : l,irst,u0,v0
  use gridsetup, only : np,mp,itrestart
  use uvRamp_variables, only : uramp,utarget,vramp,vtarget, &
    uramptime,vramptime
  use sensor_variables, only : nSensors,senw,senDataOld,senDataNew, &
    senDataName
  use windfieldio, only : windfieldstartfile,nvwind,xvwind_list, &
    w2f_index,xvdataold,xvdatanew,itinterp,windspeedupfactor, &
    itabsnew,itabsold,xvdataname
  use msga_variables, only : mpi_rank
  use xvall, only : xv,nv,xv_list
  use constants, only : tolerance
  Implicit None

  ! Local Variables
  integer,intent(in) :: xevariation

  integer :: kv,kv2
  character(len=257) :: fxvdataname
  character(len=257) :: sensorDataName

  ! Executable Code
  select case(xevariation)
    case (1)
      if(utarget.lt.tolerance) &
        utarget=uramptime*uramp+u0
      if(vtarget.lt.tolerance) &
        vtarget=vramptime*vramp+v0
    case (2)
      ! access='stream' skips fortran-specific headers and trailers
      ! TODO: parallelize?
      open(48,file=windfieldstartfile, &
        form='unformatted',access='stream',status='old')  
      read(48) nvwind
      allocate(xvwind_list(nvwind),w2f_index(nvwind))
      allocate(xvdataold(np,mp,l,nvwind))
      allocate(xvdatanew(np,mp,l,nvwind))
      read(48) xvwind_list
      do kv=1,nvwind
        do kv2=1,nv
          if(xvwind_list(kv).eq.xv_list(kv2))then
            w2f_index(kv)=kv2
            cycle
          endif
        enddo
      enddo
      if(irst.eq.0)then ! restarting from windrun
        do kv=1,nvwind
          call readio(xv(1:np,1:mp,:,w2f_index(kv)),48,l)
        enddo
      endif
      close(48)
      
      ! Need to initiate xvdataold  
      itabsold=itrestart
      itabsnew=itrestart+itinterp
      if(irst.eq.1)then
        call namefile(int(itabsold/windspeedupfactor),xvdataname,fxvdataname)
        if(mpi_rank.eq.0) print*,'reading old windfield file ',fxvdataname
        call windfld_read(fxvdataname,xvdataold)
      else
        do kv=1,nvwind
          xvdataold(:,:,:,kv)=xv(1:np,1:mp,:,w2f_index(kv))
        enddo
      endif 
      call namefile(int(itabsnew/windspeedupfactor),xvdataname,fxvdataname)
      if(mpi_rank.eq.0) print*,'reading file ',fxvdataname
      call windfld_read(fxvdataname,xvdatanew)
 
    case (3) 
      allocate(senw(nSensors))
    case (4)
      allocate(senw(nSensors*2))
      allocate(senDataOld(nSensors,4),senDataNew(nSensors,4))
      do kv=1,nSensors
        if(mpi_rank.eq.0) print*,'Reading Sensor',kv,'datum'
        call namefile(kv,senDataName,sensorDataName)
        open(4884+kv,file=sensorDataName,form='formatted',status='old')
        read(4884+kv,*) ! File Headers 
        read(4884+kv,*) ! File Headers
        read(4884+kv,*) senDataOld(kv,1),senDataOld(kv,2), &
          senDataOld(kv,3),senDataOld(kv,4)
        read(4884+kv,*) senDataNew(kv,1),senDataNew(kv,2), &
          senDataNew(kv,3),senDataNew(kv,4)
      enddo
  end select

end subroutine rinitXeVariation
