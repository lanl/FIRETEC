!-----------------------------------------------------------------------
! Reads basic user input for emissions in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_heatsource(iheatsource)
  use heatsource_gridlist_variables
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: iheatsource

  ! Executable Code
  call QueryGridlist_integer('isource_function',isource_function,48)
  call QueryGridlist_integer('nsource_hs',nsource_hs,48)
  allocate(freq_hs(nsource_hs),av_hs(nsource_hs),flux_hs(nsource_hs), &
    depth_hs(nsource_hs))
  call QueryGridlist_real_array('av_hs',av_hs,nsource_hs,48)
  call QueryGridlist_real_array('flux_hs',flux_hs,nsource_hs,48)
  call QueryGridlist_real_array('freq_hs',freq_hs,nsource_hs,48)
  call QueryGridlist_real_array('depth_hs',depth_hs,nsource_hs,48)
  if(iheatsource.eq.1)then
    allocate(xcen_hs(nsource_hs),ycen_hs(nsource_hs), &
      radius_hs(nsource_hs))
    call QueryGridlist_real_array('radius_hs',radius_hs,nsource_hs,48)
    call QueryGridlist_real_array('xcen_hs',xcen_hs,nsource_hs,48)
    call QueryGridlist_real_array('ycen_hs',ycen_hs,nsource_hs,48)
  elseif(iheatsource.eq.2)then
    allocate(xl_hs(nsource_hs),xu_hs(nsource_hs), &
      yl_hs(nsource_hs),yu_hs(nsource_hs))
    call QueryGridlist_real_array('xl_hs',xl_hs,nsource_hs,48)
    call QueryGridlist_real_array('xu_hs',xu_hs,nsource_hs,48)
    call QueryGridlist_real_array('yl_hs',yl_hs,nsource_hs,48)
    call QueryGridlist_real_array('yu_hs',yu_hs,nsource_hs,48)
  endif
  call QueryGridlist_integer('imass_source',imass_source,48)
  allocate(fm_hs(nsource_hs))
  if(imass_source.eq.2)then
    call QueryGridlist_real_array('fm_hs',fm_hs,nsource_hs,48)
  endif
  call QueryGridlist_integer('MLRHRRinputs',MLRHRRinputs,48)

end subroutine input_heatsource
