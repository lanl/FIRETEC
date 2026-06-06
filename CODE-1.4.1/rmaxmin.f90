!----------------------------------------------------------------
! rmaxmin will display the min/max of each higrad array variable
!----------------------------------------------------------------
subroutine rmaxmin(array,string,il,iu,jl,ju,lls,nvp)
  use msga_variables, only : mpi_rank,mpi_real,mpi_comm_world, &
    ierror,mpi_max,mpi_min
  use gridsetup, only : np,mp
  use xvall, only : irho,nv
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real,intent(in) :: array(il:iu,jl:ju,lls,nvp)
  character(*) :: string(nvp)
  
  integer :: kv
  real :: rmax,rmin,xmax,xmin
  real,allocatable :: tmp(:,:,:)

  ! Executable Code
  allocate(tmp(np,mp,lls))
  do kv=1,nvp
    tmp=array(1:np,1:mp,1:lls,kv)
    if(nvp.eq.nv.and.kv.ne.irho) &
      tmp=tmp/array(1:np,1:mp,1:lls,irho)
    rmax=maxval(tmp)
    rmin=minval(tmp)
    call mpi_reduce(rmax,xmax,1,mpi_real,mpi_max,0, &
      mpi_comm_world,ierror)
    call mpi_reduce(rmin,xmin,1,mpi_real,mpi_min,0, &
      mpi_comm_world,ierror)
    if(mpi_rank.eq.0) print*,' RMAXMIN for ',string(kv), &
      ' : max = ',xmax,' ,min = ',xmin
  enddo
  deallocate(tmp)

end subroutine rmaxmin
