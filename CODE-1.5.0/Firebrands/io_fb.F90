!-----------------------------------------------------------------------
! writeio_fb writes a single lagrangian array to a specified output file
!-----------------------------------------------------------------------
subroutine writeio_fb(ittot)
  use gridlist_variables, only : prec
  use firebrand_gridlist_variables, only : outname_fb
  use firebrand_general_variables, only : xvbrand,numBrands,nvbrand
  use msga_variables, only : numprocs,mpi_integer,mpi_comm_world, &
    mpi_rank,ierror
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  Implicit None

  ! Local Variables
  integer,intent(in) :: ittot

  integer :: it,iproc
  integer :: totNumBrands(numprocs),displs(numprocs)
  real(prec),allocatable :: xvtotBrands(:,:)
  character(len=257) :: fname

  ! Executable Code
  ! Gather all brand information to processor 0
  call mpi_gather(numBrands,1,mpi_integer,totNumBrands,1,mpi_integer, &
    0,mpi_comm_world,ierror)
  allocate(xvtotBrands(sum(totNumBrands),nvbrand))
  if(mpi_rank.eq.0)then
    !allocate(xvtotBrands(sum(totNumBrands),nvbrand))
    displs(1) = 0 ! compute displacements
    if(numprocs.gt.1)then
      do iproc = 2,numprocs
        displs(iproc)=displs(iproc-1)+totNumBrands(iproc-1)
      enddo
    endif
  endif
  do it=1,nvbrand
#ifdef DBL_PREC
    call mpi_gatherv(xvbrand(:,it),numBrands,mpi_double_precision, &
      xvtotBrands(:,it),totNumBrands,displs,mpi_double_precision,0, &
      mpi_comm_world,ierror)
#else
    call mpi_gatherv(xvbrand(:,it),numBrands,mpi_real, &
      xvtotBrands(:,it),totNumBrands,displs,mpi_real,0, &
      mpi_comm_world,ierror)
#endif
  enddo

  ! Write brand information to new file
  if(mpi_rank.eq.0) then
    call namefile(ittot,outname_fb,fname)
    print*,trim(fname)

    open(48,file=fname,form='unformatted',status='unknown')
    write(48) sum(totNumBrands)
    write(48) xvtotBrands
    close(48)
    deallocate(xvtotBrands)
  endif

end subroutine writeio_fb

!-----------------------------------------------------------------------
! irstreadio_fb reads the firebrand output file and initiates and 
! partitions firebrands throughout the domain
!-----------------------------------------------------------------------
subroutine irstreadio_fb
  use gridlist_variables, only : prec,dx,dy,nprocx
  use firebrand_gridlist_variables, only : restartfile_fb
  use firebrand_general_variables, only : ix,iy, &
    xvbrand,numBrands,nvbrand
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank
  Implicit None
  
  ! Local Variables
  integer :: it,count
  integer :: npos,mpos,mpi_rank_fb
  integer :: numBrand
  real(prec),allocatable :: xvtmpbrand(:,:)
  logical,allocatable :: retain(:)

  ! Executable Code
  ! Read-in restart file
  open(48,file=restartfile_fb,form='unformatted',status='old')
  read(48) numBrands
  allocate(xvbrand(numBrands,nvbrand))
  read(48) xvbrand

  ! Identify the brands to be removed from processor
  allocate(retain(numBrands)); retain = .true.
  numBrand=numBrands
  do it=1,numBrands
    npos=floor(xvbrand(it,ix)/dx/np)
    mpos=floor(xvbrand(it,iy)/dy/mp)
    mpi_rank_fb=mpos*nprocx+npos
    if(mpi_rank_fb.ne.mpi_rank)then
      retain(it)=.false.
      numBrand=numBrand-1
    endif
  enddo
  
  ! Remove irrelevant brands
  allocate(xvtmpbrand(numBrand,nvbrand))
  count=1
  do it=1,numBrands
    if(retain(it))then
      xvtmpbrand(count,:)=xvbrand(it,:)
      count=count+1
    endif
  enddo
  deallocate(xvbrand)
  allocate(xvbrand(numBrand,nvbrand))
  numBrands=numBrand
  xvbrand=xvtmpbrand
  deallocate(xvtmpbrand)

end subroutine irstreadio_fb
