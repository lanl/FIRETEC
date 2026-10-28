!-----------------------------------------------------------------------
! Reads basic user input for vegetation in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_radiation(irad)
  use radiation_gridlist_variables, only : icallrad,crad,isootmodel, &
    isolar,unixTimeStart,plotLong,plotLat
  use gridlist_variables, only : wp
  use QueryGridlist_mod, only : QueryGridlist_variable
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: irad

  ! Executable Code
  call QueryGridlist_variable('icallrad',icallrad,48)
  if(irad.eq.2)then
    call QueryGridlist_variable('crad',crad,48)
    call QueryGridlist_variable('isootmodel',isootmodel,48)
  endif
  
  ! Solar Heating
  call QueryGridlist_variable('isolar',isolar,48)
  call QueryGridlist_variable('unixT',unixTimeStart,48)
  call QueryGridlist_variable('plotLong',plotLong,48)
  call QueryGridlist_variable('plotLat',plotLat,48)

!-----------------------------------------------------------------------
! Everything below here is where others put all the warnings of garbage 
! that a user may try unintentionally to run with. It's put down here
! because in clunks up the code and makes it hard to read.
!-----------------------------------------------------------------------
  if(irad.eq.2.and.icallrad.le.0)then
    print*,'ERROR: icallrad must be greater than zero for Monte Carlo radiation'
    STOP
  endif
  if(irad.eq.2)then
    if(isootmodel.lt.0.or.isootmodel.gt.1)then
      print*,'ERROR: isootmodel must be 0 or 1'
      STOP
    endif
    if(isootmodel.eq.0.and.crad.lt.0.0_wp)then
      print*,'ERROR: crad must be nonnegative'
      STOP
    endif
  endif

end subroutine input_radiation
