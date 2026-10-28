!-----------------------------------------------------------------------
! binary_writeio is the common function called for the standard writing
! of HIGRAD output and restart files to a binary fortran file
!-----------------------------------------------------------------------
subroutine binary_writeio(ittot)
  use gridlist_variables, only : l,outname
  use vegetation_gridlist_variables, only : nfuel
  use msga_variables, only : mpi_rank
  use fuel_variables, only : lfuel
  use gridsetup, only : np,mp
  use xvall, only : nv,xv,nvfuel,xvfuel
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: ittot
  integer :: it,ift
  character(len=257) :: fname 

  ! Executable Code
  call namefile(ittot,outname,fname)
  if(mpi_rank.eq.0) print*,trim(fname)
  
  ! TODO: parallelize?
  if(mpi_rank.eq.0) open(48,file=fname,form='unformatted', &
    status='unknown')
  do it=1,nv
    call writeio(xv(1:np,1:mp,:,it),48,l)
  enddo
  do it=1,nvfuel
    do ift=1,nfuel
      call writeio(xvfuel(ift,:,:,:,it),48,lfuel)
    enddo
  enddo
  if(mpi_rank.eq.0) close(48)
end subroutine binary_writeio

!-----------------------------------------------------------------------
! irstreadio is the common function called for the standard reading of
! HIGRAD output and restart files
!-----------------------------------------------------------------------
subroutine irstreadio
  use gridlist_variables, only : l,restartfile
  use vegetation_gridlist_variables, only : nfuel
  use msga_variables, only : mpi_rank
  use fuel_variables, only : lfuel
  use gridsetup, only : np,mp,itrestart
  use xvall, only : nv,xv,nvfuel,xvfuel
  Implicit None
  
  ! Local Variables
  integer :: it,ift
  integer :: idot,iend

  ! Executable Code
  iend=len_trim(restartfile)
  idot=index(restartfile(:iend),'.',back=.true.)
  read(restartfile(idot+1:iend),*) itrestart

  if(mpi_rank.eq.0) then
    print*,'restarting with it = ',itrestart
    open(48,file=restartfile,form='unformatted',status='old')
  endif
  do it=1,nv
    call readio(xv(1:np,1:mp,:,it),48,l)
  enddo
  do it=1,nvfuel
    do ift=1,nfuel
      call readio(xvfuel(ift,:,:,:,it),48,lfuel)
    enddo
  enddo
  if(mpi_rank.eq.0) close(48)

end subroutine irstreadio

!-----------------------------------------------------------------------
! writeio writes a single array to a specified output file
!-----------------------------------------------------------------------
subroutine writeio(indata,iunit,nzdim)
  use gridlist_variables, only : n,m,nprocx,nprocy,wp
  use msga_variables, only : numprocs,mpi_comm_world,mpi_rank,ierror
  use msga_variables, only : mpi_wp
  use gridsetup, only : np,mp
  Implicit None
  
  ! Local Variables    
  integer,intent(in) :: nzdim,iunit
  real(wp),intent(in) :: indata(np,mp,nzdim)

  integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
  real(wp),allocatable :: tmparray(:,:,:,:),outdata(:,:,:)
     
  ! Executable Code
  allocate(tmparray(np,mp,nzdim,numprocs),outdata(n,m,nzdim))
  nsize=np*mp*nzdim
  call mpi_gather(indata,nsize,mpi_wp,tmparray,nsize,mpi_wp, &
    0,mpi_comm_world,ierror)
  if (mpi_rank.eq.0) then
    do iprocx=1,nprocx
      do jprocy=1,nprocy
        iproc=iprocx+(jprocy-1)*nprocx
        do k=1,nzdim
          do j=1,mp
            do i=1,np
              ia=(iprocx-1)*np+i
              ja=(jprocy-1)*mp+j
              outdata(ia,ja,k)=tmparray(i,j,k,iproc)
            enddo
          enddo
        enddo
      enddo
    enddo
    write(iunit) outdata
  endif
  deallocate(tmparray,outdata)
end subroutine writeio

!-----------------------------------------------------------------------
! readio reads a single array from a specified output file
!-----------------------------------------------------------------------
subroutine readio(array,iunit,nzdim)
  use gridlist_variables, only : n,m,nprocx,nprocy,wp
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_comm_world,mpi_rank,numprocs,ierror
  use msga_variables, only : mpi_wp
  Implicit None

  ! Local Variables
  integer,intent(in) :: nzdim,iunit
  real(wp),intent(out) :: array(np,mp,nzdim)
      
  integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
  real(wp),allocatable:: tmparray(:,:,:,:),indata(:,:,:)

  ! Executable Code
  allocate(indata(n,m,nzdim),tmparray(np,mp,nzdim,numprocs))
  if(mpi_rank.eq.0)then
    indata = 0.0_wp
    tmparray = 0.0_wp
    read(iunit) indata
    do iprocx=1,nprocx
      do jprocy=1,nprocy
        iproc=iprocx+(jprocy-1)*nprocx
        do k=1,nzdim
          do j=1,mp
            do i=1,np
              ia=(iprocx-1)*np+i
              ja=(jprocy-1)*mp+j
              tmparray(i,j,k,iproc)=indata(ia,ja,k)
            enddo
          enddo
        enddo
      enddo
    enddo
  endif

  nsize=np*mp*nzdim
  call mpi_scatter(tmparray,nsize,mpi_wp,array,nsize,mpi_wp,0, &
    mpi_comm_world,ierror)
  deallocate(tmparray,indata)

end subroutine readio
