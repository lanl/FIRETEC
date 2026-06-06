!----------------------------------------------------------------
! frqwriteio is the common function called for the standard
! writing of HIGRAD output and restart files
!----------------------------------------------------------------
subroutine frqwriteio(fname)
  use gridlist_variables, only : nfuel,windfieldstartfile,l, &
    nfuel
  use msga_variables, only : mpi_rank
  use fuel_variables, only : rhoFuel,rhoWater,sizeScale,sies,psiwmax,lfuel
  use gridsetup, only : np,mp
  use xvall, only : nv,xv
  Implicit None
  
  ! Local Variables
  character(*),intent(in) :: fname 
  integer :: it

  ! Executable Code
  ! TODO: parallelize?
  if(mpi_rank.eq.0) open(48,file=fname,form='unformatted',status='unknown')
  do it=1,nv
    call writeio(xv(1:np,1:mp,:,it),48,l)
  enddo
  if(fname.eq.windfieldstartfile) then
    if(mpi_rank.eq.0) close(48)
    return
  endif
  do it=1,nfuel
    call writeio(rhoFuel(it,1:np,1:mp,:),48,lfuel)
    call writeio(rhoWater(it,1:np,1:mp,:),48,lfuel)
    call writeio(sizeScale(it,1:np,1:mp,:),48,lfuel)
    call writeio(sies(it,1:np,1:mp,:),48,lfuel)
    call writeio(psiwmax(it,1:np,1:mp,:),48,lfuel) ! TODO: derive this from rhoWater?
  enddo
  if(mpi_rank.eq.0) close(48)
end subroutine frqwriteio

!----------------------------------------------------------------
! irstreadio is the common function called for the standard
! reading of HIGRAD output and restart files
!----------------------------------------------------------------
subroutine irstreadio(fname)
  use gridlist_variables, only : nfuel,windfieldstartfile,l, &
    nfuel
  use msga_variables, only : mpi_rank
  use fuel_variables, only : rhoFuel,rhoWater,sizeScale,sies,psiwmax,lfuel
  use gridsetup, only : np,mp
  use xvall, only : nv,xv
  Implicit None
  
  ! Local Variables
  character(*),intent(in) :: fname 
  integer :: it

  ! Executable Code
  if(mpi_rank.eq.0) open(48,file=fname,form='unformatted',status='old')
  do it=1,nv
    call readio(xv(1:np,1:mp,:,it),48,l)
  enddo
  if(fname.eq.windfieldstartfile) return
  do it=1,nfuel
    call readio(rhoFuel(it,1:np,1:mp,:),48,lfuel)
    call readio(rhoWater(it,1:np,1:mp,:),48,lfuel)
    call readio(sizeScale(it,1:np,1:mp,:),48,lfuel)
    call readio(sies(it,1:np,1:mp,:),48,lfuel)
    call readio(psiwmax(it,1:np,1:mp,:),48,lfuel)
  enddo
  if(mpi_rank.eq.0) close(48)

end subroutine irstreadio

!----------------------------------------------------------------
! writeio writes a single array to a specified output file
!----------------------------------------------------------------
subroutine writeio(indata,iunit,nzdim)
  use gridlist_variables, only : n,m,nprocx,nprocy
  use msga_variables, only : numprocs,mpi_real,mpi_comm_world, &
    mpi_rank,ierror
  use gridsetup, only : np,mp
  Implicit None
  
  ! Local Variables    
  integer,intent(in) :: nzdim,iunit
  real,intent(in) :: indata(np,mp,nzdim)

  integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
  real,allocatable :: tmparray(:,:,:,:),outdata(:,:,:)
     
  ! Executable Code
  if (mpi_rank.eq.0) &
    allocate(tmparray(np,mp,nzdim,numprocs),outdata(n,m,nzdim))
  nsize=np*mp*nzdim
  call mpi_gather(indata,nsize,mpi_real,tmparray,nsize,mpi_real, &
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
    deallocate(tmparray,outdata)
  endif
end subroutine writeio

!----------------------------------------------------------------
! readio reads a single array from a specified output file
!----------------------------------------------------------------
subroutine readio(array,iunit,nzdim)
  use gridlist_variables, only : n,m,nprocx,nprocy
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_real,mpi_comm_world,mpi_rank, &
    numprocs,ierror
  Implicit None

  ! Local Variables
  integer,intent(in) :: nzdim,iunit
  real,intent(in) :: array(np,mp,nzdim)
      
  integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
  real,allocatable:: tmparray(:,:,:,:),indata(:,:,:)

  ! Executable Code
  if(mpi_rank.eq.0)then
    allocate(indata(n,m,nzdim),tmparray(np,mp,nzdim,numprocs))
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
  call mpi_scatter(tmparray,nsize,mpi_real,array,nsize,mpi_real,0, &
    mpi_comm_world,ierror)
  if(mpi_rank.eq.0) deallocate(tmparray,indata)

end subroutine readio
