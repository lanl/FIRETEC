!----------------------------------------------------------------
! Reads basic user input for running HIGRAD/FIRETEC
!----------------------------------------------------------------
subroutine input
  use gridlist_variables
  use gridsetup, only : dts
  use constants, only : tolerance
  use vegetation_gridlist_variables, only : ifire
  use radiation_gridlist_variables, only : irad,isolar
  use QueryGridlist_mod, only: QueryGridlist_variable
  Implicit None

  ! Local Variables
  character(len=256) :: gridlist_file

  ! Executable Code
  gridlist_file="gridlist"
  if (command_argument_count().ge.1) &
    call get_command_argument(1,gridlist_file)
  open(unit=48,file=trim(gridlist_file),action="read")
  
  ! Temporal and spatial domain
  call QueryGridlist_variable('irst',irst,48)
  call QueryGridlist_variable('nt',nt,48)
  call QueryGridlist_variable('nts',nts,48)
  call QueryGridlist_variable('dt',dt,48)
  call QueryGridlist_variable('n',n,48)
  call QueryGridlist_variable('m',m,48)
  call QueryGridlist_variable('l',l,48)
  call QueryGridlist_variable('dx',dx,48)
  call QueryGridlist_variable('dy',dy,48)
  call QueryGridlist_variable('dz',dz,48)
  call QueryGridlist_variable('aa1',aa1,48)
  call QueryGridlist_variable('nprocx',nprocx,48)
  call QueryGridlist_variable('nprocy',nprocy,48)
  call QueryGridlist_variable ('topofile',topofile,48)
  
  ! General I/O
  call QueryGridlist_variable('frqoutput',frqoutput,48)
  call QueryGridlist_variable('frqfilstr',frqfilstr,48)
  call QueryGridlist_variable('outname',outname,48)
  call ensure_directory_exists(outname)
  call QueryGridlist_variable('restartfile',restartfile,48)
  call QueryGridlist_variable('icfmeflag',icfmeflag,48)
  call QueryGridlist_variable('ihdf',ihdf,48)

  ! Relaxation
  call QueryGridlist_variable('ih',ih,48)
  call QueryGridlist_variable('nr',nr,48)
  call QueryGridlist_variable('ibcx',ibcx,48)
  call QueryGridlist_variable('ibcy',ibcy,48)
  call QueryGridlist_variable('ibc_top',ibc_top,48)
  call QueryGridlist_variable('ibc_bot',ibc_bot,48)
  call QueryGridlist_variable('ibc_left',ibc_left,48)
  call QueryGridlist_variable('ibc_right',ibc_right,48)
  if(ibc_top.eq.7825)   ibc_top   = ibcy
  if(ibc_bot.eq.7825)   ibc_bot   = ibcy
  if(ibc_left.eq.7825)  ibc_left  = ibcx
  if(ibc_right.eq.7825) ibc_right = ibcx
  call QueryGridlist_variable('ibctopopen',ibctopopen,48)
  call QueryGridlist_variable('iab',iab,48)
  call QueryGridlist_variable('zab',zab,48)
  call QueryGridlist_variable('tow',tow,48)
  
  ! Ambient Conditions
  call QueryGridlist_variable('tambient',tambient,48)
  call QueryGridlist_variable('pressground',pressground,48)
  call QueryGridlist_variable('zgroundref',zgroundref,48)
  call QueryGridlist_variable('zStabBot',zStabBot,48)
  call QueryGridlist_variable('zStabMid',zStabMid,48)
  call QueryGridlist_variable('thetaStabMid',thetaStabMid,48)
  call QueryGridlist_variable('thetaStabTop',thetaStabTop,48)
  call QueryGridlist_variable('iperturb',iperturb,48)
  call QueryGridlist_variable('u0',u0,48)
  call QueryGridlist_variable('uswitch',uswitch,48)
  call QueryGridlist_variable('v0',v0,48)
  call QueryGridlist_variable('vswitch',vswitch,48)
  if(uswitch.eq.1.or.vswitch.eq.1) ixevariation=1
  call QueryGridlist_variable('zu',zu,48)
  call QueryGridlist_variable('ius',ius,48)
  call QueryGridlist_variable('iue',iue,48)
  call QueryGridlist_variable('jus',jus,48)
  call QueryGridlist_variable('jue',jue,48)
  call QueryGridlist_variable('slopeangle',slopeangle,48)
  call QueryGridlist_variable('slopeazimuth',slopeazimuth,48)

  ! Atmospheric Transport
  call QueryGridlist_variable('icorio',icorio,48)
  call QueryGridlist_variable('ilspgf',ilspgf,48)
  call QueryGridlist_variable('izlspgf',izlspgf,48)
  call QueryGridlist_variable('frqlspgf',frqlspgf,48)
  call QueryGridlist_variable('islip',islip,48)
  call QueryGridlist_variable('iord',iord,48)
  call QueryGridlist_variable('nonos',nonos,48)
  call QueryGridlist_variable('idiv',idiv,48)
  call QueryGridlist_variable('nfct',nfct,48)
  call QueryGridlist_variable('nonosold',nonosold,48)

  ! Vegetation
  call input_vegetation

  ! Subgrid Turbulence
  call QueryGridlist_variable('iturb',iturb,48)
  if(iturb.gt.0) call input_subGridTurb(iturb)

  ! Radiation
  call QueryGridlist_variable('irad',irad,48)
  if(irad.gt.0) call input_radiation(irad)

  ! XeVariation
  call QueryGridlist_variable('ixevariation',ixevariation,48)
  call QueryGridlist_variable('windfieldin',windfieldin,48)
  if(windfieldin.eq.1) ixevariation=2
  call QueryGridlist_variable('windfieldout',windfieldout,48)
  if(ixevariation.gt.0.or.windfieldout.eq.1) &
    call input_xevariation(ixevariation,windfieldout)

  ! Sensors
  call QueryGridlist_variable('isensor',isensor,48)
  if(isensor.gt.0) call input_sensors

  !Heat source
  call QueryGridlist_variable('iheatsource',iheatsource,48)
  if(iheatsource.ge.1) call input_heatsource(iheatsource)

  ! Firebrands
  call QueryGridlist_variable('ifbrand',ifbrand,48)
  if(ifbrand.eq.1) call input_firebrands

  ! Emissions
  call QueryGridlist_variable('iemissions',iemissions,48)
  call input_emissions(iemissions,ifire)

#ifdef LANL_MOD
  call input_LANL_MOD()
#endif
  
  if(abs(dt-0.01).lt.tolerance) then
    dts = dt/real(nts,kind=wp)
    call QueryGridlist_variable('dts',dts,48)
    dt = dts*real(nts,kind=wp)
  endif

  close(48)

  !-----------------------------------------------------------------------
  ! Everything below here is where others put all the warnings of garbage 
  ! that a user may try unintentionally to run with. It's put down here
  ! because in clunks up the code and makes it hard to read.
  !-----------------------------------------------------------------------
  if(frqoutput.le.0)then
    print*,'ERROR: frqoutput must be greater than zero'
    STOP
  endif
  if(icfmeflag.ne.0)then
    print*,'ERROR: icfmeflag is configured, but slice-grid input &
      &is not implemented in this tree'
    STOP
  endif

  ! Relaxation
  call validate_lateral_ibc('ibc_top',ibc_top)
  call validate_lateral_ibc('ibc_bot',ibc_bot)
  call validate_lateral_ibc('ibc_left',ibc_left)
  call validate_lateral_ibc('ibc_right',ibc_right)

  if(nr.lt.1)then
    print*,'ERROR: nr must be at least 1'
    STOP
  endif
  if(((ibc_top.ge.1.and.ibc_top.lt.5) .or. &
      (ibc_bot.ge.1.and.ibc_bot.lt.5) .or. &
      (ibc_left.ge.1.and.ibc_left.lt.5) .or. &
      (ibc_right.ge.1.and.ibc_right.lt.5)) .and. nr.lt.3)then
    print*,'ERROR: nr must be at least 3 for relaxation boundary  &
      &options 1-4'
    STOP
  endif
  if(((ibc_left.ge.1.and.ibc_left.lt.5).or. &
    (ibc_right.ge.1.and.ibc_right.lt.5)).and.nr.ge.n)then
    print*,'INVALID RELAXATION WIDTH: nr must be smaller than n'
    STOP
  endif
  if(((ibc_bot.ge.1.and.ibc_bot.lt.5).or. &
    (ibc_top.ge.1.and.ibc_top.lt.5)).and.nr.ge.m)then
    print*,'INVALID RELAXATION WIDTH: nr must be smaller than m'
    STOP
  endif
  if((ibc_left.ge.1.and.ibc_left.lt.5).and. &
    (ibc_right.ge.1.and.ibc_right.lt.5).and.2*nr.ge.n)then
    print*,'INVALID RELAXATION WIDTH: left/right sponge regions cover &
      &the x domain'
    STOP
  endif
  if((ibc_bot.ge.1.and.ibc_bot.lt.5).and. &
    (ibc_top.ge.1.and.ibc_top.lt.5).and.2*nr.ge.m)then
    print*,'INVALID RELAXATION WIDTH: bottom/top sponge regions cover &
      &the y domain'
    STOP
  endif

  if(ibc_top.eq.0.and.ibc_bot.ne.0) then
    print*,'CANNOT HAVE CYCLIC CONDITIONS ON TOP BUT NOT ON BOTTOM'
    print*,'Either both ibc_top and ibc_bot are 0 or both are > 0.'
    print*,'If both are > 0, they need not match.'
    STOP
  elseif(ibc_bot.eq.0.and.ibc_top.ne.0) then
    print*,'CANNOT HAVE CYCLIC CONDITIONS ON BOTTOM BUT NOT ON TOP'
    print*,'Either both ibc_top and ibc_bot are 0 or both are > 0.'
    print*,'If both are > 0, they need not match.'
    STOP
  endif

  if(ibc_left.eq.0.and.ibc_right.ne.0) then
    print*,'CANNOT HAVE CYCLIC CONDITIONS ON LEFT BUT NOT ON RIGHT'
    print*,'Either both ibc_left and ibc_right are 0 or both are > 0.'
    print*,'If both are > 0, they need not match.'
    STOP
  elseif(ibc_right.eq.0.and.ibc_left.ne.0) then
    print*,'CANNOT HAVE CYCLIC CONDITIONS ON RIGHT BUT NOT ON LEFT'
    print*,'Either both ibc_left and ibc_right are 0 or both are > 0.'
    print*,'If both are > 0, they need not match.'
    STOP
  endif
  if(tow.le.0.0_wp)then
    print*,'ERROR: tow must be positive'
    STOP
  endif
  if(ibctopopen.eq.0)then
    if(iab.lt.1.or.iab.gt.3)then
      print*,'ERROR: iab must be in the range 1:3 for the closed-top &
        &absorber'
      STOP
    endif
    if(iab.eq.2.and.tow.le.1.0_wp)then
      print*,'ERROR: iab=2 requires tow > 1'
      STOP
    endif
  endif
  
  ! Ambient Conditions

  ! Atmospheric Transport
  if(ilspgf.ge.1.and.frqlspgf.le.0)then
    print*,'ERROR: frqlspgf must be greater than zero when ilspgf &
      &is enabled'
    STOP
  endif

  ! Vegetation

  ! Subgrid Turbulence
  if(iturb.eq.1)then
    print*,'iturb=1, Smagorinski not smagging yet'
    STOP
  endif

  ! Radiation
  if(irad.lt.0.or.irad.gt.3)then
    print*,'ERROR: irad must be 0, 1, 2, or 3'
    STOP
  endif

  ! XeVariation
  if(windfieldin.gt.1)then
    print*,'ERROR: windfieldin=2 one-file mode is documented but not &
      &implemented'
    STOP
  endif
  if(windfieldout.gt.1)then
    print*,'ERROR: windfieldout=2 one-file mode is documented but not &
      &implemented'
    STOP
  endif
  if(windfieldin.eq.1)then
    if(uswitch.eq.1.or.vswitch.eq.1)then
      print*,'ERROR: windfieldin=1 cannot be combined with &
        &uswitch/vswitch=1 ramps'
      STOP
    endif
    if(ixevariation.gt.0.and.ixevariation.ne.2)then
      print*,'ERROR: windfieldin=1 requires ixevariation=0 or 2'
      STOP
    endif
  elseif(uswitch.eq.1.or.vswitch.eq.1)then
    if(ixevariation.gt.0.and.ixevariation.ne.1)then
      print*,'ERROR: uswitch/vswitch=1 requires ixevariation=0 or 1'
      STOP
    endif
    ixevariation=1
  endif

  ! Sensors
  if(isensor.ne.0.and.isensor.ne.1)then
    print*,'ERROR: isensor must be 0=off or 1=enabled'
    STOP
  endif

  !Heat source
  if(iheatsource.ne.0.and.iheatsource.ne.1.and.iheatsource.ne.2)then
    print*,'ERROR: iheatsource must be 0, 1 circular, or 2 rectangular'
    STOP
  endif

  ! Firebrands
  if(ifbrand.lt.0.or.ifbrand.gt.1) &
    STOP 'ERROR: ifbrand must be 0 or 1.'
  if(ifbrand.eq.1)then
    if(ifire.ne.1)then
      print*,'ERROR: ifbrand=1 requires ifire=1'
      STOP
    endif
  endif

  ! Emissions
  if(ifire.eq.1.and.iemissions.eq.0) iemissions=2
  if(iemissions.eq.3.and.ifire.ne.1)then
    print*,'ERROR: iemissions=3 (ZBEST) requires ifire=1'
    stop
  endif
  if(iemissions.lt.0.or.iemissions.gt.3)then
    print*,'ERROR: invalid iemissions value ',iemissions, &
    '. Valid values are 0, 1, 2, or 3.'
    stop
  endif

  if((ifire.eq.1.or.isolar.eq.1.or.irad.gt.0).and.iturb.lt.2)then
    print*,'Ya did wrong!!! ifire/isolar/irad fuel-energy paths &
      &require iturb=2 because they use ka and kb'
    STOP
  endif
end subroutine input

subroutine validate_lateral_ibc(name,ibc)
  implicit none

  character(len=*),intent(in) :: name
  integer,intent(in) :: ibc

  if(.not.(ibc.eq.0.or.(ibc.ge.1.and.ibc.le.5)))then
    print*,'INVALID BOUNDARY CONDITION: ',trim(name),' = ',ibc
    print*,'Allowed values are 0 cyclic, 1-4 relaxation, 5 symmetric.'
    STOP
  endif
end subroutine validate_lateral_ibc

!-----------------------------------------------------------------------
! ensure_directory_exists creates the directory path for a file if it
! doesn't already exist. Uses mkdir -p which is silent if dir exists.
!-----------------------------------------------------------------------
subroutine ensure_directory_exists(filepath)
  use mpi, only : MPI_Initialized,MPI_Barrier,MPI_COMM_WORLD
  use msga_variables, only : mpi_rank
  implicit none
  
  character(len=*), intent(in) :: filepath
  character(len=500) :: dirpath, mkdir_cmd
  integer :: last_slash, ierr
  logical :: dir_exists, mpi_is_initialized
  
  ! Find last slash to extract directory path
  last_slash = index(filepath, '/', back=.true.)

  if (last_slash > 0) then
    dirpath = filepath(1:last_slash-1)

    call MPI_Initialized(mpi_is_initialized, ierr)

    if (mpi_is_initialized) then
      ! Only rank 0 creates/checks/prints
      if (mpi_rank == 0) then
        inquire(file=trim(dirpath)//'/.', exist=dir_exists)

        if (.not. dir_exists) then
          mkdir_cmd = 'mkdir -p "' // trim(dirpath) // '"'
          call execute_command_line(trim(mkdir_cmd), wait=.true., &
            exitstat=ierr)

          if (ierr == 0) then
            print*, 'Created directory: ', trim(dirpath)
          else
            print*, 'ERROR creating directory: ', trim(dirpath)
          endif
        endif
      endif

      ! Make sure rank 0 creates the directory before others continue
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
    else
      ! Non-MPI case, or called before MPI_Init
      inquire(file=trim(dirpath)//'/.', exist=dir_exists)

      if (.not. dir_exists) then
        mkdir_cmd = 'mkdir -p "' // trim(dirpath) // '"'
        call execute_command_line(trim(mkdir_cmd), wait=.true., &
          exitstat=ierr)

        if (ierr == 0) then
          print*, 'Created directory: ', trim(dirpath)
        else
          print*, 'ERROR creating directory: ', trim(dirpath)
        endif
      endif

    endif
  endif

end subroutine ensure_directory_exists
