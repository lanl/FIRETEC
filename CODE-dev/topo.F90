!-----------------------------------------------------------------------
! topo reads in topography data and populates the zs domain with said 
! data for all the processors
!-----------------------------------------------------------------------
subroutine topo
  use gridlist_variables, only : topofile,n,m,ih,wp
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank,npos,mpos,leftedge, &
    rightedge,botedge,topedge
  use metric_variables, only : zs 
  Implicit None

  ! Local Variables
  integer :: i,j
  real(wp),allocatable :: zsio(:,:)

  ! Executable Code
  allocate(zsio(n,m))
  open(25,file=trim(topofile),form='unformatted',status='old')
    read(25) zsio
  close(25)
  zs(1:np,1:mp)=zsio((npos-1)*np+1:(npos-1)*np+np, &
    (mpos-1)*mp+1:(mpos-1)*mp+mp)
  if(mpi_rank.eq.0) print*,'Read topography file ',trim(topofile)
  deallocate(zsio)

  ! Initialize topo values for ghost points
  call update(zs,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0,0)
  
  ! Set edges of the domain
  if(leftedge.eq.1) then
    do j=1-ih,mp+ih
      do i=1,ih
        zs(1-i,j)=zs(1,j)
      enddo
    enddo
  endif
  if(rightedge.eq.1) then
    do j=1-ih,mp+ih
      do i=1,ih
        zs(np+i,j)=zs(np,j)
      enddo
    enddo
  endif
  if(botedge.eq.1) then
    do i=1-ih,np+ih
      do j=1,ih
        zs(i,1-j)=zs(i,1)
      enddo
    enddo
  endif
  if(topedge.eq.1) then
    do i=1-ih,np+ih
      do j=1,ih
        zs(i,mp+j)=zs(i,mp)
      enddo
    enddo
  endif

end subroutine topo
