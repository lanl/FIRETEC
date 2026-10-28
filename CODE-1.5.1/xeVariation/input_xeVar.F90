!-----------------------------------------------------------------------
! Reads basic user input for xevariations in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_xevariation(ixevariation,windfieldout)
  use gridlist_variables, only : n,m,nprocx,nprocy,wp
  use uvRamp_variables, only : uramp,vramp,uramptime,vramptime, &
    utarget,vtarget
  use windfieldio, only : windspeedupfactor,itwindfield,itinterp, &
    ibcells,jbcells,is,ie,js,je,windfieldstartfile,xvdataname
  use sensor_variables, only : nSensors,SensorTimes,SensorX,SensorY, &
    SensorZ,SensorUvel,SensorVvel,senDataName
  use QueryGridlist_mod, only : QueryGridlist_variable, &
    QueryGridlist_array
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: ixevariation,windfieldout

  ! Executable Code
  if(ixevariation.eq.1)then ! u/v ramp
    call QueryGridlist_variable('uramp',uramp,48)
    call QueryGridlist_variable('utarget',utarget,48)
    call QueryGridlist_variable('uramptime',uramptime,48)
    call QueryGridlist_variable('vramp',vramp,48)
    call QueryGridlist_variable('vtarget',vtarget,48)
    call QueryGridlist_variable('vramptime',vramptime,48)
  elseif(ixevariation.eq.2)then ! read-in windfield
    call QueryGridlist_variable('windspeedupfactor',windspeedupfactor,48)
    call QueryGridlist_variable('itwindfield',itwindfield,48)
    call QueryGridlist_variable('itinterp',itinterp,48)
    call QueryGridlist_variable('ibcells',ibcells,48)
    call QueryGridlist_variable('jbcells',jbcells,48)
    is=1; call QueryGridlist_variable('is',is,48)
    ie=n; call QueryGridlist_variable('ie',ie,48)
    js=1; call QueryGridlist_variable('js',js,48)
    je=m; call QueryGridlist_variable('je',je,48)
    call QueryGridlist_variable('windfieldstartfile', &
      windfieldstartfile,48)
    call ensure_directory_exists(windfieldstartfile)
    call QueryGridlist_variable('xvdataname',xvdataname,48)
    call ensure_directory_exists(xvdataname)
  elseif(ixevariation.eq.3)then ! sparse sensors
    call QueryGridlist_variable('nSensors',nSensors,48)
    allocate(SensorX(nSensors),SensorY(nSensors),SensorZ(nSensors), &
      SensorTimes(nSensors),SensorUvel(nSensors),SensorVvel(nSensors)) 
    call QueryGridlist_array('SensorTimes',SensorTimes,nSensors,48)
    call QueryGridlist_array('SensorX',SensorX,nSensors,48)
    call QueryGridlist_array('SensorY',SensorY,nSensors,48)
    call QueryGridlist_array('SensorZ',SensorZ,nSensors,48)
    call QueryGridlist_array('SensorUvel',SensorUvel,nSensors,48)
    call QueryGridlist_array('SensorVvel',SensorVvel,nSensors,48)
  elseif(ixevariation.eq.4)then ! dense sensors
    call QueryGridlist_variable('nSensors',nSensors,48)
    allocate(SensorX(nSensors),SensorY(nSensors),SensorZ(nSensors))
    call QueryGridlist_array('SensorX',SensorX,nSensors,48)
    call QueryGridlist_array('SensorY',SensorY,nSensors,48)
    call QueryGridlist_array('SensorZ',SensorZ,nSensors,48)
    call QueryGridlist_variable('senDataName',senDataName,48)
  endif

  if(windfieldout.eq.1) then
    call QueryGridlist_variable('itwindfield',itwindfield,48)
    call QueryGridlist_variable('itinterp',itinterp,48)
    call QueryGridlist_variable('ibcells',ibcells,48)
    call QueryGridlist_variable('jbcells',jbcells,48)
    is=1; call QueryGridlist_variable('is',is,48)
    ie=n; call QueryGridlist_variable('ie',ie,48)
    js=1; call QueryGridlist_variable('js',js,48)
    je=m; call QueryGridlist_variable('je',je,48)
    call QueryGridlist_variable ('windfieldstartfile', &
      windfieldstartfile,48)
    call QueryGridlist_variable ('xvdataname',xvdataname,48)
    ! Create windfield output directory if it doesn't exist
    call ensure_directory_exists(xvdataname)
  endif

!-----------------------------------------------------------------------
! Everything below here is where others put all the warnings of garbage 
! that a user may try unintentionally to run with. It's put down here
! because in clunks up the code and makes it hard to read.
!-----------------------------------------------------------------------
  if(ixevariation.eq.1)then ! u/v ramp
  elseif(ixevariation.eq.2)then ! read-in windfield
    if(windspeedupfactor.le.0)then
      print*,'ERROR: windspeedupfactor must be greater than zero'
      STOP
    endif
    if(itinterp.le.0)then
      print*,'ERROR: itinterp must be greater than zero for windfield &
        &input'
      STOP
    endif
    if(nprocx.le.0.or.nprocy.le.0)then
      print*,'ERROR: nprocx and nprocy must be greater than zero'
      STOP
    endif
    if(ibcells.le.0.or.ibcells.gt.n/nprocx.or. &
       jbcells.le.0.or.jbcells.gt.m/nprocy)then
      print*,'ERROR: ibcells/jbcells must be positive and fit within &
        &each MPI subdomain'
      STOP
    endif
    if(2*ibcells.gt.n.or.2*jbcells.gt.m)then
      print*,'ERROR: windfield boundary strips overlap; reduce &
        &ibcells/jbcells'
      STOP
    endif
    if(is.ne.1.or.ie.ne.n.or.js.ne.1.or.je.ne.m)then
      print*,'ERROR: xeVariation windfield subdomain bounds are not &
        &implemented'
      STOP
    endif
  elseif(ixevariation.eq.3)then ! sparse sensors
    if(nSensors.le.0)then
      print*,'ERROR: nSensors must be greater than zero'
      STOP
    endif
    if(any(SensorZ.le.0.0_wp))then
      print*,'ERROR: all SensorZ values must be greater than zero'
      STOP
    endif
  elseif(ixevariation.eq.4)then ! dense sensors
    if(nSensors.le.0)then
      print*,'ERROR: nSensors must be greater than zero'
      STOP
    endif
    if(any(SensorZ.le.0.0_wp))then
      print*,'ERROR: all SensorZ values must be greater than zero'
      STOP
    endif
  endif

  if(windfieldout.eq.1) then
    if(itinterp.le.0)then
      print*,'ERROR: itinterp must be greater than zero for windfield &
        &output'
      STOP
    endif
    if(nprocx.le.0.or.nprocy.le.0)then
      print*,'ERROR: nprocx and nprocy must be greater than zero'
      STOP
    endif
    if(ibcells.le.0.or.ibcells.gt.n/nprocx.or. &
       jbcells.le.0.or.jbcells.gt.m/nprocy)then
      print*,'ERROR: ibcells/jbcells must be positive and fit within &
        &each MPI subdomain'
      STOP
    endif
    if(2*ibcells.gt.n.or.2*jbcells.gt.m)then
      print*,'ERROR: windfield boundary strips overlap; reduce &
        &ibcells/jbcells'
      STOP
    endif
    if(is.ne.1.or.ie.ne.n.or.js.ne.1.or.je.ne.m)then
      print*,'ERROR: xeVariation windfield subdomain bounds are not &
        &implemented'
      STOP
    endif
  endif

end subroutine input_xevariation
