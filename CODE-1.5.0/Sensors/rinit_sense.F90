!-----------------------------------------------------------------------
! rinitSense initializes all arrays and variables related to 
! sensor evolution throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitSense()
  use gridlist_variables, only : l,dx,dy,prec
  use sensor_gridlist_variables, only : numsensors,se_xloc,se_yloc, &
    se_zloc,ncycle,sensorfile
  use sensor_general_variables, only : se_proc,se_i,se_j,se_k,se_x, &
    se_y,se_z,xv_values
  use gridsetup, only : itrestart,np,mp
  use msga_variables, only : mpi_rank,npos,mpos
  use metric_variables, only : zedge,zs
  use xvall, only : xv_list
  use zcart_function, only : zcart
  implicit none

  ! Local Variables
  integer :: i,j,k,ia,ja,ni
  integer :: unitnum
  character(len=200) :: sensorFileName
  logical :: exists
  real(prec) :: zmh,zph
  real(prec) :: fuzz=0.0001

  ! Executable Code
  allocate(se_proc(numsensors))
  allocate(se_i(numsensors),se_j(numsensors),se_k(numsensors))
  allocate(se_x(numsensors),se_y(numsensors),se_z(numsensors))
  allocate(xv_values(numsensors,size(xv_list)))

  !Open sensor files
  if(mpi_rank.eq.0)then
    do i=1,numsensors
      unitnum = 16400+i

      !Check if file exists:
      call namefile(i,sensorfile,sensorFileName)
      inquire(file=trim(sensorFileName)//'.csv', exist = exists)
      if(exists)then ! if file does exists
        open(unit=unitnum,file=trim(sensorFileName)//'.csv', &
          form='formatted',status='old',position='append')
        if(ncycle.lt.itrestart) then
          write(*,*)'WARNING: itrestart > ncycle in sensors.'
          write(*,*)'WARNING: sensor times may be out of order.'
        endif
      else ! if file does not exist
        open(unit=unitnum,file=trim(sensorFileName)//'.csv', &
          form='formatted',status='unknown')
        write(unitnum,'(A16)',advance='no')'cycle, time, '
        do ni=1,size(xv_list)
          write(unitnum,'(A8)',advance='no') xv_list(ni)
          if(ni.ne.size(xv_list)) then
            write(unitnum,'(A1)',advance='no')','
          endif
        enddo
      endif
      close(unitnum)
    enddo
  endif

  ! Search and find the sensor processor, cell indices and physical locations
  se_proc(:)=-1
  do ni=1,numsensors
    do k=1,l
      do j=1,mp
        ja = (mpos-1)*mp + j - 1
        do i=1,np
          ia = (npos-1)*np + i - 1
          zmh = zcart(zedge(k),i,j)-zs(i,j)
          zph = zcart(zedge(k+1),i,j)-zs(i,j)

          if((se_xloc(ni).ge.real(ia)*dx + fuzz).and. &
            (se_xloc(ni).lt.(real(ia)+1.)*dx + fuzz).and. &
            (se_yloc(ni).ge.real(ja)*dy + fuzz).and. &
            (se_yloc(ni).lt.(real(ja)+1.)*dy + fuzz).and. &
            (se_zloc(ni).ge.zmh) .and. &
            (se_zloc(ni).lt.zph) ) then
            print*,'sensor',ni,'is in processor',mpi_rank

            se_proc(ni) = mpi_rank
            se_x(ni) = (real(ia)+0.5)*dx
            se_y(ni) = (real(ja)+0.5)*dy
            se_z(ni) = zmh
            se_i(ni) = i
            se_j(ni) = j
            se_k(ni) = k
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine rinitSense
