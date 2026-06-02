!-----------------------------------------------------------------------
! Reads basic user input for emissions in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_firebrands
  use firebrand_gridlist_variables
  Implicit None
  
  ! Local Variables

  ! Executable Code
  call QueryGridlist_integer('shapeFB',shapeFB,48)
  call QueryGridlist_integer('dsizeFB',dsizeFB,48)
  call QueryGridlist_real('radiusFB',radius,48)
  call QueryGridlist_real('heightFB',height,48)
  call QueryGridlist_real('rhoFB',rhoFB,48)
  call QueryGridlist_real('TempThreshFB',TempHot,48)
  call QueryGridlist_string ('outname_fb',outname_fb,48)
  call QueryGridlist_string ('restartfile_fb',restartfile_fb,48)

end subroutine input_firebrands
