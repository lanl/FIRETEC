!-----------------------------------------------------------------------
! This package provides a passive sensor capability for writing
! high frequency sensor data to a text file. For now, sensors
! are fixed in space. Future capability will provide sensors
! that advect with the flow field. 
!-----------------------------------------------------------------------
subroutine sensors_write()
  use sensor_gridlist_variables, only : numsensors,ncycle,sensorfile
  use sensor_general_variables, only : se_proc,se_i,se_j,se_k, &
    xv_values
  use msga_variables, only : mpi_rank
  use gridsetup, only : ittot,time
  use xvall, only : xv,xv_list,irho
  Implicit none

  ! Local Variables
  integer :: unitnum
  character(len=200) :: sensorFileName
  integer :: ni,i

  ! Executable Code
  ! Locate the xv values for each sensor and store them
  do ni=1,numsensors
    unitnum = 16400+ni
    call namefile(ni,sensorfile,sensorFileName)
    if(mpi_rank.eq.se_proc(ni)) then
      do i=1,size(xv_list)
        if(i.eq.irho) then
          xv_values(ni,i) = xv(se_i(ni),se_j(ni),se_k(ni),i)
        else
          xv_values(ni,i) = xv(se_i(ni),se_j(ni),se_k(ni),i) &
                       / xv(se_i(ni),se_j(ni),se_k(ni),irho)
        endif
      enddo
      ! Write sensor variables to file
      open(unit=unitnum,file=trim(sensorFileName)//'.csv', &
          form='formatted',status='old',position='append')
      if(ncycle.lt.ittot) then
        write(unitnum,'(A1)') ''
        write(unitnum,'(I8,A1)',advance='no')ittot,','
        write(unitnum,'(F14.6,A1)',advance='no')time,','
        do i=1,size(xv_list)
          write(unitnum,'(F23.16)',advance='no')xv_values(ni,i)
          if(i.ne.size(xv_list)) then
            write(unitnum,'(A1)',advance='no')','
          endif
        enddo
      endif
      close(unitnum)
    endif
  enddo

end subroutine sensors_write
