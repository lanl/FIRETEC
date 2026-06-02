!-----------------------------------------------------------------------
! frqwriteio is the common function called for the standard writing of
! HIGRAD output and restart files
!-----------------------------------------------------------------------
subroutine frqwriteio(fname)
  use gridlist_variables, only : nfuel,l
  use msga_variables, only : mpi_rank
  use fuel_variables, only : lfuel
  use gridsetup, only : np,mp
  use xvall, only : nv,xv,nvfuel,xvfuel
  Implicit None
  
  ! Local Variables
  character(*),intent(in) :: fname 
  integer :: it,ift

  ! Executable Code
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
end subroutine frqwriteio

!-----------------------------------------------------------------------
! irstreadio is the common function called for the standard reading of
! HIGRAD output and restart files
!-----------------------------------------------------------------------
subroutine irstreadio
  use gridlist_variables, only : nfuel,l,restartfile
  use msga_variables, only : mpi_rank
  use fuel_variables, only : lfuel
  use gridsetup, only : np,mp,itrestart
  use xvall, only : nv,xv,nvfuel,xvfuel
  Implicit None
  
  ! Local Variables
  integer :: it,ift
  integer :: idot,iend

  ! Executable Code
  idot=index(restartfile,'out.')+4
  iend=len_trim(restartfile)
  read(restartfile(idot:iend),*) itrestart

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
  use gridlist_variables, only : n,m,nprocx,nprocy,prec
  use msga_variables, only : numprocs,mpi_comm_world,mpi_rank,ierror
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use gridsetup, only : np,mp
  Implicit None
  
  ! Local Variables    
  integer,intent(in) :: nzdim,iunit
  real(prec),intent(in) :: indata(np,mp,nzdim)

  integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
  real(prec),allocatable :: tmparray(:,:,:,:),outdata(:,:,:)
     
  ! Executable Code
  if (mpi_rank.eq.0) &
    allocate(tmparray(np,mp,nzdim,numprocs),outdata(n,m,nzdim))
  nsize=np*mp*nzdim
#ifdef DBL_PREC
  call mpi_gather(indata,nsize,mpi_double_precision,tmparray, &
    nsize,mpi_double_precision,0,mpi_comm_world,ierror)
#else
  call mpi_gather(indata,nsize,mpi_real,tmparray,nsize,mpi_real, &
    0,mpi_comm_world,ierror)
#endif
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
    deallocate(tmparray,outdata)
  endif
end subroutine writeio

!-----------------------------------------------------------------------
! readio reads a single array from a specified output file
!-----------------------------------------------------------------------
subroutine readio(array,iunit,nzdim)
  use gridlist_variables, only : n,m,nprocx,nprocy,prec
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_comm_world,mpi_rank,numprocs,ierror
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  Implicit None

  ! Local Variables
  integer,intent(in) :: nzdim,iunit
  real(prec),intent(in) :: array(np,mp,nzdim)
      
  integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
  real(prec),allocatable:: tmparray(:,:,:,:),indata(:,:,:)

  ! Executable Code
  if(mpi_rank.eq.0)then
    allocate(indata(n,m,nzdim),tmparray(np,mp,nzdim,numprocs))
    indata = 0.0
    tmparray = 0.0
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
#ifdef DBL_PREC
  call mpi_scatter(tmparray,nsize,mpi_double_precision,array, &
    nsize,mpi_double_precision,0,mpi_comm_world,ierror)
#else
  call mpi_scatter(tmparray,nsize,mpi_real,array,nsize,mpi_real,0, &
    mpi_comm_world,ierror)
#endif
  if(mpi_rank.eq.0) deallocate(tmparray,indata)

end subroutine readio



