!-----------------------------------------------------------------------
! Reads basic user input for emissions in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_firebrands
  use firebrand_gridlist_variables
  use QueryGridlist_mod, only : QueryGridlist_variable
  Implicit None
  
  ! Local Variables

  ! Executable Code
  call QueryGridlist_variable('shapeFB',shapeFB,48)
  call QueryGridlist_variable('dsizeFB',dsizeFB,48)
  call QueryGridlist_variable('radiusFB',radius,48)
  call QueryGridlist_variable('heightFB',height,48)
  call QueryGridlist_variable('rhoFB',rhoFB,48)
  call QueryGridlist_variable('TempThreshFB',TempHot,48)
  call QueryGridlist_variable ('outname_fb',outname_fb,48)
  call ensure_directory_exists(outname_fb)
  call QueryGridlist_variable ('restartfile_fb',restartfile_fb,48)

  if(radius.le.0.0_wp) STOP 'ERROR: radiusFB must be positive.'
  if(height.le.0.0_wp) STOP 'ERROR: heightFB must be positive.'
  if(rhoFB.le.0.0_wp) STOP 'ERROR: rhoFB must be positive.'

end subroutine input_firebrands
