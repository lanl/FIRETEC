!-----------------------------------------------------------------------
! Reads basic user input for emissions in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_sensors
  use sensor_gridlist_variables
  Implicit None
  
  ! Local Variables
  integer :: s 
  
  ! Executable Code
  call QueryGridlist_integer('se_frq_write',se_frq_write,48)
  call QueryGridlist_integer('numsensors',numsensors,48)
  call QueryGridlist_integer('ncycle',ncycle,48)
  call QueryGridlist_integer('locfile',locfile,48)

  allocate(se_xloc(numsensors))
  allocate(se_yloc(numsensors))
  allocate(se_zloc(numsensors))

  call QueryGridlist_string('sensorfile',sensorfile,48)

  if(locfile.eq.1) then
    call QueryGridlist_string('sensorlocs',sensorlocs,48)
    open(unit=7825,file=sensorlocs, form='formatted',status='old')
    do s=1,numsensors
      read(7825,*) se_xloc(s),se_yloc(s),se_zloc(s)
    enddo
    close(7825)
  else
    call QueryGridlist_real_array('se_xloc',se_xloc,numsensors,48)
    call QueryGridlist_real_array('se_yloc',se_yloc,numsensors,48)
    call QueryGridlist_real_array('se_zloc',se_zloc,numsensors,48)
  endif

end subroutine input_sensors
