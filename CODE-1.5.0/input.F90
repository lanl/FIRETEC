!----------------------------------------------------------------
! Reads basic user input for running HIGRAD/FIRETEC
!----------------------------------------------------------------
subroutine input
  use gridlist_variables
  Implicit None

  ! Local Variables

  ! Executable Code
  open(unit=48,file="gridlist",action="read")
  
  ! Temporal and spatial domain
  call QueryGridlist_integer('irst',irst,48)
  call QueryGridlist_integer('nt',nt,48)
  call QueryGridlist_integer('nts',nts,48)
  call QueryGridlist_real   ('dts',dts,48)
  call QueryGridlist_integer('ntp',ntp,48)
  call QueryGridlist_integer('n',n,48)
  call QueryGridlist_integer('m',m,48)
  call QueryGridlist_integer('l',l,48)
  call QueryGridlist_real   ('dx',dx,48)
  call QueryGridlist_real   ('dy',dy,48)
  call QueryGridlist_real   ('dz',dz,48)
  call QueryGridlist_real   ('aa1',aa1,48)
  call QueryGridlist_integer('nprocx',nprocx,48)
  call QueryGridlist_integer('nprocy',nprocy,48)
  call QueryGridlist_string ('topofile',topofile,48)
  
  ! General I/O
  call QueryGridlist_integer('frqoutput',frqoutput,48)
  call QueryGridlist_integer('frqfilstr',frqfilstr,48)
  call QueryGridlist_string ('outname',outname,48)
  call QueryGridlist_string ('restartfile',restartfile,48)
  call QueryGridlist_integer('icfmeflag',icfmeflag,48)
  call QueryGridlist_integer('ihdf',ihdf,48)

  ! Relaxation
  call QueryGridlist_integer('ih',ih,48)
  call QueryGridlist_integer('nr',nr,48)
  call QueryGridlist_integer('ibcx',ibcx,48)
  call QueryGridlist_integer('ibcy',ibcy,48)
  call QueryGridlist_integer('ibclatopen',ibclatopen,48)
  call QueryGridlist_integer('ibctopopen',ibctopopen,48)
  call QueryGridlist_integer('iab',iab,48)
  call QueryGridlist_real   ('zab',zab,48)
  call QueryGridlist_real   ('zabt',zabt,48)
  call QueryGridlist_real   ('tow',tow,48)
  call QueryGridlist_integer('itheta',itheta,48)
  
  ! Ambient Conditions
  call QueryGridlist_real   ('tambient',tambient,48)
  call QueryGridlist_real   ('pressground',pressground,48)
  call QueryGridlist_real   ('zgroundref',zgroundref,48)
  call QueryGridlist_integer('iperturb',iperturb,48)
  call QueryGridlist_real   ('u0',u0,48)
  call QueryGridlist_integer('uswitch',uswitch,48)
  call QueryGridlist_real   ('v0',v0,48)
  call QueryGridlist_integer('vswitch',vswitch,48)
  if(uswitch.eq.1.or.vswitch.eq.1) ixevariation=1
  call QueryGridlist_real   ('zu',zu,48)
  call QueryGridlist_integer('ius',ius,48)
  call QueryGridlist_integer('iue',iue,48)
  call QueryGridlist_integer('jus',jus,48)
  call QueryGridlist_integer('jue',jue,48)
  call QueryGridlist_real   ('slopeangle',slopeangle,48)
  call QueryGridlist_real   ('slopeazimuth',slopeazimuth,48)

  ! Atmospheric Transport
  call QueryGridlist_integer('icorio',icorio,48)
  call QueryGridlist_integer('ilspgf',ilspgf,48)
  call QueryGridlist_integer('izlspgf',izlspgf,48)
  call QueryGridlist_integer('frqlspgf',frqlspgf,48)
  call QueryGridlist_integer('islip',islip,48)
  call QueryGridlist_integer('iord',iord,48)
  call QueryGridlist_integer('nonos',nonos,48)
  call QueryGridlist_integer('idiv',idiv,48)
  call QueryGridlist_integer('nfct',nfct,48)
  call QueryGridlist_integer('nonosold',nonosold,48)

  ! Vegetation
  call QueryGridlist_integer('nfuel',nfuel,48)
  call QueryGridlist_integer('ivegread',ivegread,48)
  call QueryGridlist_real   ('rhoMicro',rhoMicro,48)
  call QueryGridlist_real   ('cpwood',cpwood,48)

  ! Burning
  call QueryGridlist_integer('ifire',ifire,48)
  call QueryGridlist_integer('inonlocal',inonlocal,48)
  call QueryGridlist_integer('idiffsies',idiffsies,48)
  call QueryGridlist_integer('ifuel',ifuel,48)
  call QueryGridlist_string ('ignfile',ignfile,48)

  ! Subgrid Turbulence
  call QueryGridlist_integer('iturb',iturb,48)
  if(iturb.gt.0) call input_subGridTurb(iturb)

  ! Radiation
  call QueryGridlist_integer('irad',irad,48)
  call QueryGridlist_integer('icallrad',icallrad,48)
  call QueryGridlist_real   ('crad',crad,48)
  call QueryGridlist_integer('irandseed',irandseed,48)
  call QueryGridlist_integer('iseed',iseed,48)
  call QueryGridlist_integer('isootmodel',isootmodel,48)
  call QueryGridlist_integer('iradeastflux',iradeastflux,48)

  ! XeVariation
  call QueryGridlist_integer('ixevariation',ixevariation,48)
  call QueryGridlist_integer('windfieldin',windfieldin,48)
  if(windfieldin.eq.1) ixevariation=2
  call QueryGridlist_integer('windfieldout',windfieldout,48)
  if(ixevariation.gt.0.or.windfieldout.eq.1) &
    call input_xevariation(ixevariation,windfieldout)

  ! Sensors
  call QueryGridlist_integer('isensor',isensor,48)
  if(isensor.gt.0) call input_sensors

  !Heat source
  call QueryGridlist_integer('iheatsource',iheatsource,48)
  if(iheatsource.ge.1) call input_heatsource(iheatsource)

  ! Firebrands
  call QueryGridlist_integer('ifbrand',ifbrand,48)
  if(ifbrand.eq.1) call input_firebrands(ifbrand)

  ! Emissions
  call QueryGridlist_integer('iemissions',iemissions,48)
  if(ifire.eq.1.and.iemissions.eq.0) iemissions=2
  if(iemissions.gt.0) call input_emissions(iemissions,ifire)
  call QueryGridlist_real('relativehumidity',relativeHumidity,48)

  close(48)

end subroutine input
