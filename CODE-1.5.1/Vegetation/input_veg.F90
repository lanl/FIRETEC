!-----------------------------------------------------------------------
! Reads basic user input for vegetation in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_vegetation()
  use gridlist_variables, only : wp
  use vegetation_gridlist_variables, only : nfuel,ivegread,rhoMicro, &
    cpwood,ifire,inonlocal,idiffsies,ifuel,ignfile
  use QueryGridlist_mod, only : QueryGridlist_variable
  Implicit None
  
  ! Local Variables

  ! Executable Code
  call QueryGridlist_variable('nfuel',nfuel,48)
  call QueryGridlist_variable('ivegread',ivegread,48)
  call QueryGridlist_variable('rhoMicro',rhoMicro,48)
  call QueryGridlist_variable('cpwood',cpwood,48)

  ! Burning
  call QueryGridlist_variable('ifire',ifire,48)
  call QueryGridlist_variable('inonlocal',inonlocal,48)
  call QueryGridlist_variable('idiffsies',idiffsies,48)
  call QueryGridlist_variable('ifuel',ifuel,48)
  call QueryGridlist_variable('ignfile',ignfile,48)

!-----------------------------------------------------------------------
! Everything below here is where others put all the warnings of garbage 
! that a user may try unintentionally to run with. It's put down here
! because in clunks up the code and makes it hard to read.
!-----------------------------------------------------------------------
  if(nfuel.lt.1)then
    print*,'ERROR: nfuel must be at least 1; received ',nfuel
    STOP
  endif
  if(rhoMicro.le.0.0_wp)then
    print*,'ERROR: rhoMicro must be positive; received ',rhoMicro
    STOP
  endif
  if(cpwood.le.0.0_wp)then
    print*,'ERROR: cpwood must be positive; received ',cpwood
    STOP
  endif

  ! Burning
  if(ifire.ne.0.and.ifire.ne.1)then
    print*,'ERROR: ifire must be 0 or 1; received ',ifire
    STOP
  endif
  if(inonlocal.ne.0.and.inonlocal.ne.1)then
    print*,'ERROR: inonlocal must be 0 or 1; received ',inonlocal
    STOP
  endif
  if(idiffsies.ne.0.and.idiffsies.ne.1)then
    print*,'ERROR: idiffsies must be 0 or 1; received ',idiffsies
    STOP
  endif
  if(idiffsies.eq.1.and.ifire.ne.1)then
    print*,'ERROR: idiffsies=1 requires ifire=1'
    STOP
  endif
  if(ifuel.lt.1.or.ifuel.gt.3)then
    print*,'ERROR: ifuel must be 1, 2, or 3; received ',ifuel
    STOP
  endif

end subroutine input_vegetation
