module rmaxmin_mod

  contains
  !-------------------------------------------------------------------
  ! rmaxmin will display the min/max of each array variable
  !-------------------------------------------------------------------
  subroutine rmaxmin(array,string)
    use gridlist_variables, only : wp
    use msga_variables, only : mpi_rank,mpi_comm_world,ierror, &
      mpi_max,mpi_min
    use msga_variables, only : mpi_wp  
    Implicit None
  
    ! Local Variables
    real(wp),intent(in) :: array(..)
    character(*),intent(in) :: string
    
    real(wp) :: rmax,rmin,xmax,xmin
  
    ! Executable Code
    select rank(array)
      rank (0)
        rmax=array
        rmin=array
      rank (1)
        rmax=maxval(array)
        rmin=minval(array)
      rank (2)
        rmax=maxval(array)
        rmin=minval(array)
      rank (3)
        rmax=maxval(array)
        rmin=minval(array)
      rank (4)
        rmax=maxval(array)
        rmin=minval(array)
      rank (5)
        rmax=maxval(array)
        rmin=minval(array)
      rank default
        print*,'Unsupported rank',rank(array)
    end select
    call mpi_reduce(rmax,xmax,1,mpi_wp,mpi_max,0,mpi_comm_world,ierror)
    call mpi_reduce(rmin,xmin,1,mpi_wp,mpi_min,0,mpi_comm_world,ierror)
    if(mpi_rank.eq.0)then
      print*,' RMAXMIN for ',string,' : max = ',xmax,' ,min = ',xmin
    endif
  end subroutine rmaxmin

end module rmaxmin_mod
