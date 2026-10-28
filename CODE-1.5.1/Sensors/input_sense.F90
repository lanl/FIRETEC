!-----------------------------------------------------------------------
! Reads basic user input for emissions in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_sensors
  use gridlist_variables, only : n,m,l,dx,dy,dz
  use sensor_gridlist_variables
  use QueryGridlist_mod, only : QueryGridlist_variable, &
    QueryGridlist_array
  Implicit None
  
  ! Local Variables
  integer :: s, ios
  
  ! Executable Code
  call QueryGridlist_variable('se_frq_write',se_frq_write,48)
  if(se_frq_write.le.0)then
    print*,'ERROR: se_frq_write must be greater than zero'
    STOP
  endif
  call QueryGridlist_variable('numsensors',numsensors,48)
  call QueryGridlist_variable('ncycle',ncycle,48)
  call QueryGridlist_variable('locfile',locfile,48)

  allocate(se_xloc(numsensors))
  allocate(se_yloc(numsensors))
  allocate(se_zloc(numsensors))

  call QueryGridlist_variable('sensorfile',sensorfile,48)
  call ensure_directory_exists(sensorfile)

  if(locfile.eq.1) then
    call QueryGridlist_variable('sensorlocs',sensorlocs,48)
    open(unit=7825,file=sensorlocs,form='formatted',status='old', &
      iostat=ios)
    if(ios.ne.0) then
      write(*,*) 'ERROR: unable to open sensor location file: ', &
        trim(sensorlocs)
      stop 
    endif
    do s=1,numsensors
      read(7825,*,iostat=ios) se_xloc(s),se_yloc(s),se_zloc(s)
      if(ios.ne.0) then
        write(*,*) 'ERROR: failed reading sensor location ', s, &
          ' from ',trim(sensorlocs)
        stop 
      endif
    enddo
    close(7825)
  else
    call QueryGridlist_array('se_xloc',se_xloc,numsensors,48)
    call QueryGridlist_array('se_yloc',se_yloc,numsensors,48)
    call QueryGridlist_array('se_zloc',se_zloc,numsensors,48)
  endif
  !---------------- Below are warnings and stop functions ----------------!
  do s=1,numsensors
    if((se_xloc(s).gt.n*dx).or.(se_yloc(s).gt.m*dy).or.(se_zloc(s).gt.l*dz)) then
      print*,'SENSOR',s,'OUTSIDE OF DOMAIN'
      print*,'LOCATION x,y,z:',se_xloc(s),se_yloc(s),se_zloc(s)
      STOP
    endif
  enddo



end subroutine input_sensors
