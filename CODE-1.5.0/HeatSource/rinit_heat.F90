!-----------------------------------------------------------------------
! rinitEmission initializes all arrays and variables related to 
! smoke, emissions, and plume evolution throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitHeatSource()
  use heatsource_gridlist_variables, only : MLRHRRinputs,TME,MLR,HRR, &
    inputName,nsource_hs
  use gridlist_variables, only : irst
  use gridsetup, only : time
  use msga_variables, only : mpi_rank
  Implicit None

  ! Local Variables
  real :: tmetemp,mlrtemp,hrrtemp
  integer :: iostat,itsource
  character(len=257) :: filename

  ! Executable Code
  if(MLRHRRinputs.eq.1) then
    allocate(TME(nsource_hs,2))
    allocate(MLR(nsource_hs,2))
    allocate(HRR(nsource_hs,2))
    do itsource = 1,nsource_hs
      call namefile(itsource,inputName,filename)
      open(7825+itsource,file=filename,form='formatted',status='old')
      read(7825+itsource,*) ! headers
      if(irst.eq.1) then
        read(7825+itsource,*) tmetemp, mlrtemp, hrrtemp
        do while (tmetemp.lt.time)
          read(7825+itsource,*,iostat=iostat) tmetemp, mlrtemp, hrrtemp ! Get to appropriate line
          if(iostat.lt.0) exit 
        enddo
        backspace(7825+itsource) ! Go back one line since tmetemp.lt.time
        backspace(7825+itsource)
        backspace(7825+itsource)
      endif
      read(7825+itsource,*) TME(itsource,1),MLR(itsource,1),HRR(itsource,1)
      read(7825+itsource,*,iostat=iostat) TME(itsource,2),MLR(itsource,2),HRR(itsource,2)
      if(mpi_rank.eq.0) print*,'Reading MLR/HRR data'
      if(iostat.lt.0) then ! iostat<0 means end of file
        TME(itsource,2) = TME(itsource,1)
        MLR(itsource,2) = MLR(itsource,1)
        HRR(itsource,2) = HRR(itsource,1)
      endif
    enddo
  endif
end subroutine rinitHeatSource
