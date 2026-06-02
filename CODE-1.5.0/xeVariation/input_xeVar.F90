!-----------------------------------------------------------------------
! Reads basic user input for xevariations in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_xevariation(ixevariation,windfieldout)
  use gridlist_variables, only : n,m
  use uvRamp_variables, only : uramp,vramp,uramptime,vramptime
  use windfieldio, only : windspeedupfactor,itwindfield,itinterp, &
    ibcells,jbcells,is,ie,js,je,windfieldstartfile,xvdataname
  use sensor_variables, only : nSensors,SensorTimes,SensorX,SensorY, &
    SensorZ,SensorUvel,SensorVvel,freqSensor
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: ixevariation,windfieldout

  ! Executable Code
  if(ixevariation.eq.1)then ! u/v ramp
    call QueryGridlist_real('uramp',uramp,48)
    call QueryGridlist_real('uramptime',uramptime,48)
    call QueryGridlist_real('vramp',vramp,48)
    call QueryGridlist_real('vramptime',vramptime,48)
  elseif(ixevariation.eq.2)then ! read-in windfield
    call QueryGridlist_integer('windspeedupfactor',windspeedupfactor,48)
    call QueryGridlist_integer('itwindfield',itwindfield,48)
    call QueryGridlist_integer('itinterp',itinterp,48)
    call QueryGridlist_integer('ibcells',ibcells,48)
    call QueryGridlist_integer('jbcells',jbcells,48)
    is=1; call QueryGridlist_integer('is',is,48)
    ie=n; call QueryGridlist_integer('ie',ie,48)
    js=1; call QueryGridlist_integer('js',js,48)
    je=m; call QueryGridlist_integer('je',je,48)
    call QueryGridlist_string ('windfieldstartfile', &
      windfieldstartfile,48)
    call QueryGridlist_string ('xvdataname',xvdataname,48)
  elseif(ixevariation.eq.3)then ! sparse sensors
    call QueryGridlist_integer('nSensors',nSensors,48)
    allocate(SensorX(nSensors),SensorY(nSensors),SensorZ(nSensors), &
      SensorTimes(nSensors),SensorUvel(nSensors),SensorVvel(nSensors)) 
    call QueryGridlist_real_array('SensorTimes',SensorTimes,nSensors,48)
    call QueryGridlist_real_array('SensorX',SensorX,nSensors,48)
    call QueryGridlist_real_array('SensorY',SensorY,nSensors,48)
    call QueryGridlist_real_array('SensorZ',SensorZ,nSensors,48)
    call QueryGridlist_real_array('SensorUvel',SensorUvel,nSensors,48)
    call QueryGridlist_real_array('SensorVvel',SensorVvel,nSensors,48)
  elseif(ixevariation.eq.4)then ! dense sensors
    call QueryGridlist_integer('nSensors',nSensors,48)
    allocate(SensorX(nSensors),SensorY(nSensors),SensorZ(nSensors))
    call QueryGridlist_real_array('SensorX',SensorX,nSensors,48)
    call QueryGridlist_real_array('SensorY',SensorY,nSensors,48)
    call QueryGridlist_real_array('SensorZ',SensorZ,nSensors,48)
    call QueryGridlist_integer('freqSensor',freqSensor,48)
  endif

  if(windfieldout.eq.1) then
    call QueryGridlist_integer('itwindfield',itwindfield,48)
    call QueryGridlist_integer('itinterp',itinterp,48)
    call QueryGridlist_integer('ibcells',ibcells,48)
    call QueryGridlist_integer('jbcells',jbcells,48)
    is=1; call QueryGridlist_integer('is',is,48)
    ie=n; call QueryGridlist_integer('ie',ie,48)
    js=1; call QueryGridlist_integer('js',js,48)
    je=m; call QueryGridlist_integer('je',je,48)
    call QueryGridlist_string ('windfieldstartfile', &
      windfieldstartfile,48)
    call QueryGridlist_string ('xvdataname',xvdataname,48)
  endif

end subroutine input_xevariation
