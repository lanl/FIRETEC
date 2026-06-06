!----------------------------------------------------------------
! startFromFile initializes the higrad/firetec domain based on a
! specified file, either a windfield start file or a restart file
!----------------------------------------------------------------
subroutine startFromFile
  use gridlist_variables, only : irst,restartfile,windfieldstartfile, &
    iwindfieldin,itinterp,windspeedupfactor,xvdataname
  use gridsetup, only : itrestart
  use msga_variables, only : mpi_rank
  use xvall, only : xvdataold,xv
  Implicit None

  ! Local Variables
  integer :: idot,iend
  character(len=257) :: fxvdataname

  ! Executable Code
  if(irst.eq.1)then
    idot=index(restartfile,'out.')+4
    iend=len_trim(restartfile)
    read(restartfile(idot:iend),*) itrestart
    if(mpi_rank.eq.0) print*,'restarting with ',itrestart
    call irstreadio(restartfile)
  else ! restart from a windfieldstart (irst.eq.0 and iwindfieldin.eq.1)
    if(mpi_rank.eq.0) print*,'starting with windfieldstart ', &
      windfieldstartfile
    call irstreadio(windfieldstartfile)
  endif

  ! Restart in the case of iwindfieldin (old comes from restart and new from xvdata)
  if(iwindfieldin.eq.1)then
    xvdataold=xv
    call namefile((itrestart+itinterp)/windspeedupfactor,xvdataname,fxvdataname)
    if(mpi_rank.eq.0) print*,'reading file ',fxvdataname
    call windfld_read(fxvdataname)
  endif

end subroutine startFromFile
