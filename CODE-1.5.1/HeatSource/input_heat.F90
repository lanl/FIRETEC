!-----------------------------------------------------------------------
! Reads basic user input for emissions in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_heatsource(iheatsource)
  use heatsource_gridlist_variables
  use gridlist_variables, only : wp,dt,nt
  use QueryGridlist_mod, only : QueryGridlist_variable, &
    QueryGridlist_array
  Implicit None
  
  ! Local Variables
  integer :: itmp
  integer,intent(in) :: iheatsource
  real(wp),allocatable :: tmp(:)

  ! Executable Code
  call QueryGridlist_variable('nsource_hs',nsource_hs,48)

  allocate(tmp(nsource_hs))
  tmp = 0.0_wp                 ! defaults: function=0, no mass, no MLR/HRR

  allocate(isource_function_hs(nsource_hs))
  allocate(starttime_hs(nsource_hs))
  allocate(endtime_hs(nsource_hs))
  allocate(imass_source_hs(nsource_hs))
  allocate(MLRHRRinputs_hs(nsource_hs))

  ! Handle case with no endtime_hs or starttime_hs
  endtime_hs = 1.0e20_wp  ! or HUGE(1.0_wp)
  starttime_hs = 0.0_wp

  call QueryGridlist_array('isource_function_hs',tmp,nsource_hs,48)
  isource_function_hs = nint(tmp)

  ! Legacy gridlists only have isource_function instead of isource_function_hs
  itmp = -1
  call QueryGridlist_variable('isource_function',itmp,48)
  if(itmp.ne.-1) isource_function_hs = itmp

  call QueryGridlist_array('starttime_hs',starttime_hs,nsource_hs,48)
  call QueryGridlist_array('endtime_hs',endtime_hs,nsource_hs,48)

  allocate(freq_hs(nsource_hs));  freq_hs=0._wp
  allocate(av_hs(nsource_hs));    av_hs=-1._wp
  allocate(flux_hs(nsource_hs));  flux_hs=-1._wp
  allocate(depth_hs(nsource_hs)); depth_hs=-1._wp
  call QueryGridlist_array('av_hs',av_hs,nsource_hs,48)
  call QueryGridlist_array('flux_hs',flux_hs,nsource_hs,48)
  call QueryGridlist_array('freq_hs',freq_hs,nsource_hs,48)
  call QueryGridlist_array('depth_hs',depth_hs,nsource_hs,48)

  if(iheatsource.eq.1)then
    allocate(xcen_hs(nsource_hs),ycen_hs(nsource_hs), &
      radius_hs(nsource_hs))
    call QueryGridlist_array('radius_hs',radius_hs,nsource_hs,48)
    call QueryGridlist_array('xcen_hs',xcen_hs,nsource_hs,48)
    call QueryGridlist_array('ycen_hs',ycen_hs,nsource_hs,48)
  elseif(iheatsource.eq.2)then
    allocate(xl_hs(nsource_hs),xu_hs(nsource_hs), &
      yl_hs(nsource_hs),yu_hs(nsource_hs))
    xl_hs = 0.0_wp; xu_hs = 0.0_wp
    yl_hs = 0.0_wp; yu_hs = 0.0_wp
    call QueryGridlist_array('xl_hs',xl_hs,nsource_hs,48)
    call QueryGridlist_array('xu_hs',xu_hs,nsource_hs,48)
    call QueryGridlist_array('yl_hs',yl_hs,nsource_hs,48)
    call QueryGridlist_array('yu_hs',yu_hs,nsource_hs,48)
  endif

  tmp = 0.0_wp
  call QueryGridlist_array('imass_source_hs',tmp,nsource_hs,48)
  imass_source_hs = nint(tmp)

  ! Legacy gridlists only include imass_source instead of imass_source_hs
  itmp = -1
  call QueryGridlist_variable('imass_source',itmp,48)
  if(itmp.ne.-1) imass_source_hs = itmp

  tmp = 0.0_wp
  call QueryGridlist_array('MLRHRRinputs_hs',tmp,nsource_hs,48)
  MLRHRRinputs_hs = nint(tmp)

  ! Legacy gridlists only include MLRHRRinputs instead of MLRHRRinputs_hs
  itmp = -1
  call QueryGridlist_variable('MLRHRRinputs',itmp,48)
  if(itmp.ne.-1) MLRHRRinputs_hs = itmp

  allocate(fm_hs(nsource_hs))
  fm_hs = 0.0_wp
  if(any(imass_source_hs.eq.2))then
    call QueryGridlist_array('fm_hs',fm_hs,nsource_hs,48)
  endif

  deallocate(tmp)

! -----------------------------------------------------------------------
!Instructions: -----------------------------------------------------------
!Multiple heat sources are specified using arrays of length nsource_hs.
!The ith entry of each array belongs to heat source i.
!
!All heat sources in a single simulation must use the same geometry type.
!The geometry type is selected once by iheatsource:
!
!  iheatsource = 1  circular heat sources
!  iheatsource = 2  rectangular heat sources
!
!This means circular and rectangular heat sources cannot be mixed in the
!same heat-source setup. If iheatsource = 1, every heat source uses the
!circle inputs radius_hs, xcen_hs, and ycen_hs. If iheatsource = 2, every
!heat source uses the rectangle inputs xl_hs, xu_hs, yl_hs, and yu_hs.
!
!Example with two circular heat sources:
!
!  nsource_hs = 2
!  iheatsource = 1
!
!  starttime_hs        = 0.0,  30.0
!  endtime_hs          = 60.0, 90.0
!  isource_function_hs = 0,    2
!  radius_hs           = r1,   r2
!  xcen_hs             = x1,   x2
!  ycen_hs             = y1,   y2
!  flux_hs             = ...
!  depth_hs            = ...
!
!Heat source 1 is active from 0 to 60 s with a constant function.
!Heat source 2 is active from 30 to 90 s with a sinusoidal function.
!starttime_hs defaults to 0 sec start time
!endtime_hs defaults to inf (never turns off)
!
!Example with two rectangular heat sources:
!
!  nsource_hs = 2
!  iheatsource = 2
!
!  xl_hs = xl1, xl2
!  xu_hs = xu1, xu2
!  yl_hs = yl1, yl2
!  yu_hs = yu1, yu2
!
!Both are rectangular because iheatsource = 2.
!
!If MLRHRRinputs_hs(i) = 1, heat source i reads its own MLR/HRR input file.
!During the simulation, HRR and MLR are linearly interpolated in time and
!distributed over that source's geometry.

!-----------------------------------------------------------------------
! Everything below here is where others put all the warnings of garbage 
! that a user may try unintentionally to run with. It's put down here
! because in clunks up the code and makes it hard to read.
!-----------------------------------------------------------------------
  if(nsource_hs.le.0)then
    print*,'ERROR: nsource_hs must be greater than zero when heat &
      &sources are enabled'
    STOP
  endif

  ! Handle case with no endtime_hs or starttime_hs
  if(any(isource_function_hs.lt.0).or.any(isource_function_hs.gt.3))then
    print*,'ERROR: isource_function_hs entries must be 0, 1, 2, or 3'
    STOP
  endif
  if(any(starttime_hs.le.real(nt,wp)*dt))then
    print*,'WARNING: heat source activates during windfield production.'
    print*,'This may affect how the windfield develops.'
    
  endif
  if(any(endtime_hs.le.starttime_hs))then
    print*,'ERROR: each endtime_hs must be greater than starttime_hs'
    STOP
  endif

  if(any((MLRHRRinputs_hs.ne.1).and. &
    ((av_hs.lt.0._wp).or.(flux_hs.lt.0._wp))))then
    print*,'ERROR: non-MLR/HRR heat sources require nonnegative av_hs and flux_hs'
    STOP
  endif  

  if(any(depth_hs.le.0._wp))then
    print*,'ERROR: depth_hs entries must be greater than zero'
    STOP
  endif
  if(any((isource_function_hs.eq.3).and.(freq_hs.le.0.0_wp)))then
    print*,'ERROR: pulse heat sources require freq_hs > 0'
    STOP
  endif

  if(iheatsource.eq.1)then
    if(any(radius_hs.le.0.0_wp))then
      print*,'ERROR: radius_hs entries must be greater than zero'
      STOP
    endif
  elseif(iheatsource.eq.2)then
    if(any(xu_hs.le.xl_hs).or.any(yu_hs.le.yl_hs))then
      print*,'ERROR: rectangular heat sources require xl_hs < xu_hs &
        &and yl_hs < yu_hs'
      STOP
    endif
  endif

  if(any(imass_source_hs.lt.0).or.any(imass_source_hs.gt.2))then
    print*,'ERROR: imass_source_hs entries must be 0, 1, or 2'
    STOP
  endif

  if(any(MLRHRRinputs_hs.lt.0).or.any(MLRHRRinputs_hs.gt.1))then
    print*,'ERROR: MLRHRRinputs_hs entries must be 0 or 1'
    STOP
  endif

end subroutine input_heatsource
