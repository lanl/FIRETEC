!-----------------------------------------------------------------------
! rmaxmin_3D will display the min/max of each higrad array variable
!-----------------------------------------------------------------------
subroutine rmaxmin_1D(array,string,il,iu,nvp)
  use gridlist_variables, only : prec
  use msga_variables, only : mpi_rank,mpi_comm_world,ierror,mpi_max, &
    mpi_min
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif  
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,nvp
  real(prec),intent(in) :: array(il:iu,nvp)
  character(*) :: string(nvp)
  
  integer :: kv
  real(prec) :: rmax,rmin,xmax,xmin
  real :: xmaxprint,xminprint

  ! Executable Code
  do kv=1,nvp
    if(iu.gt.0)then
      rmax=maxval(array(:,kv))
      rmin=minval(array(:,kv))
    else
      rmax=1e-20
      rmin=1e20
    endif
#ifdef DBL_PREC
    call mpi_reduce(rmax,xmax,1,mpi_double_precision,mpi_max,0, &
      mpi_comm_world,ierror)
    call mpi_reduce(rmin,xmin,1,mpi_double_precision,mpi_min,0, &
      mpi_comm_world,ierror)
    xmaxprint = real(xmax)
    xminprint = real(xmin)
#else
    call mpi_reduce(rmax,xmax,1,mpi_real,mpi_max,0, &
      mpi_comm_world,ierror)
    call mpi_reduce(rmin,xmin,1,mpi_real,mpi_min,0, &
      mpi_comm_world,ierror)
    xmaxprint = xmax
    xminprint = xmin
#endif
    if(mpi_rank.eq.0) print*,' RMAXMIN for ',string(kv), &
      ' : max = ',xmaxprint,' ,min = ',xminprint
  enddo

end subroutine rmaxmin_1D

!-----------------------------------------------------------------------
! rmaxmin_3D will display the min/max of each higrad array variable
!-----------------------------------------------------------------------
subroutine rmaxmin_3D(array,string,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : prec
  use msga_variables, only : mpi_rank,mpi_comm_world,ierror,mpi_max, &
    mpi_min
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif  
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(prec),intent(in) :: array(il:iu,jl:ju,lls,nvp)
  character(*) :: string(nvp)
  
  integer :: kv
  real(prec) :: rmax,rmin,xmax,xmin
  real :: xmaxprint,xminprint

  ! Executable Code
  do kv=1,nvp
    rmax=maxval(array(:,:,:,kv))
    rmin=minval(array(:,:,:,kv))
#ifdef DBL_PREC
    call mpi_reduce(rmax,xmax,1,mpi_double_precision,mpi_max,0, &
      mpi_comm_world,ierror)
    call mpi_reduce(rmin,xmin,1,mpi_double_precision,mpi_min,0, &
      mpi_comm_world,ierror)
    xmaxprint = real(xmax)
    xminprint = real(xmin)
#else
    call mpi_reduce(rmax,xmax,1,mpi_real,mpi_max,0, &
      mpi_comm_world,ierror)
    call mpi_reduce(rmin,xmin,1,mpi_real,mpi_min,0, &
      mpi_comm_world,ierror)
    xmaxprint = xmax
    xminprint = xmin
#endif
    if(mpi_rank.eq.0) print*,' RMAXMIN for ',string(kv), &
      ' : max = ',xmaxprint,' ,min = ',xminprint
  enddo

end subroutine rmaxmin_3D
