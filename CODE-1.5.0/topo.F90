!----------------------------------------------------------------
! topo reads in topography data and populates the zs domain with
! said data for all the processors
!----------------------------------------------------------------
subroutine topo
  use gridlist_variables, only : topofile,n,m,ih,prec
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank,npos,mpos,leftdedge, &
    rightdedge,botdedge,topdedge
  use metric_variables, only : zs 
  Implicit None

  ! Local Variables
  integer :: i,j
  real(prec),allocatable :: zsio(:,:)

  ! Executable Code
  ! Read topo data into master if present
  allocate(zsio(n,m)); zsio = 0.0
  open(25,file=topofile,form='unformatted',status='old')
    read(25) zsio
  close(25)
  zs(1:np,1:mp)=zsio((npos-1)*np+1:(npos-1)*np+np,(mpos-1)*mp+1:(mpos-1)*mp+mp)
  if(mpi_rank.eq.0) print*,'Read topography file ',topofile
  deallocate(zsio)

  ! Initialize topo values for ghost points
  call update(zs,zs,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0)
  
  ! Set edges of the domain (corners?)
  if(leftdedge.eq.1) then
    do j=1-ih,mp+ih
      zs(1,j)=zs(2,j) !zero gradient
      do i=1,ih
        zs(1-i,j)=zs(1,j)
      enddo
    enddo
  endif
  if(rightdedge.eq.1) then
    do j=1-ih,mp+ih
      zs(np,j)=zs(np-1,j) !zero gradient
      do i=1,ih
        zs(np+i,j)=zs(np,j)
      enddo
    enddo
  endif
  if(botdedge.eq.1) then
    do i=1-ih,np+ih
      zs(i,1)=zs(i,2)
      do j=1,ih
        zs(i,1-j)=zs(i,1)
      enddo
    enddo
  endif
  if(topdedge.eq.1) then
    do i=1-ih,np+ih
      zs(i,mp)=zs(i,mp-1)
      do j=1,ih
        zs(i,mp+j)=zs(i,mp)
      enddo
    enddo
  endif

end subroutine topo
