!----------------------------------------------------------------
! windfld_read reads xv data along domain boundaries from a wind
! field to buffer those boundaries
!----------------------------------------------------------------
subroutine windfld_read(fname)
  use gridlist_variables, only : n,m,l,ih,ibcells,jbcells, &
    nprocx,nprocy
  use msga_variables, only : mpi_rank,numprocs,mpi_real,ierror, &
    mpi_comm_world
  use gridsetup, only : np,mp
  use xvall, only : xvdatanew,nv
  Implicit None

  ! Local Variables
  character(*),intent(in) :: fname

  integer :: kv,i,j,k,ii,jj
  integer :: iproc,iprocx,jprocy,ia,ja
  integer :: nsize
  real,allocatable :: tmparray(:,:,:,:,:),indata(:,:,:,:)
  real,allocatable :: indatai(:,:,:,:),indataj(:,:,:,:)

  ! Executable Code
  if(mpi_rank.eq.0)then
    allocate(indatai(ibcells*2,m,l,nv),indataj(n,jbcells*2,l,nv), &
      indata(n,m,l,nv),tmparray(1-ih:np+ih,1-ih:mp+ih,l,nv,numprocs))
    open(41,file=fname,form='unformatted',status='old')
    read(41) indatai,indataj
    close(41)
    do k=1,l
      do j=1,m
        do i=1,n
          if(i.le.ibcells)then
            indata(i,j,k,:)=indatai(i,j,k,:)
          elseif(i.gt.n-ibcells)then
            ii=i-(n-ibcells*2)
            indata(i,j,k,:)=indatai(ii,j,k,:)
          else
            indata(i,j,k,:)=0.
          endif
          if(j.le.jbcells)then
            indata(i,j,k,:)=indataj(i,j,k,:)
          elseif(j.gt.m-jbcells)then
            jj=j-(m-jbcells*2)
            indata(i,j,k,:)=indataj(i,jj,k,:)
          endif
        enddo
      enddo
    enddo
    do iprocx=1,nprocx
      do jprocy=1,nprocy
        iproc=iprocx+(jprocy-1)*nprocx
        do kv=1,nv
          do k=1,l
            do j=1,mp
              do i=1,np
                ia=(iprocx-1)*np+i
                ja=(jprocy-1)*mp+j
                tmparray(i,j,k,kv,iproc)=indata(ia,ja,k,kv)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif

  nsize=(np+2.*ih)*(mp+2.*ih)*l
  call mpi_scatter(tmparray,nsize,mpi_real,xvdatanew, &
    mpi_real,0,mpi_comm_world,ierror)
  
  if(mpi_rank.eq.0) deallocate(indatai,indataj,indata,tmparray)

end subroutine windfld_read

!----------------------------------------------------------------
! windfld_write writes xv data along domain boundaries from a 
! wind field to buffer those boundaries in future runs
!----------------------------------------------------------------
subroutine windfld_write(fname)
  use gridlist_variables, only : n,m,l,ibcells,jbcells, &
    nprocx,nprocy
  use msga_variables, only : mpi_rank,numprocs,mpi_real,ierror, &
    mpi_comm_world
  use gridsetup, only : np,mp
  use xvall, only : xv,nv
  Implicit None

  ! Local Variables
  character(*),intent(in) :: fname

  integer :: kv,i,j,k,ii,jj
  integer :: iproc,iprocx,jprocy,ia,ja
  integer :: nsize
  real,allocatable :: tmparray(:,:,:,:,:)
  real,allocatable :: outdatai(:,:,:,:),outdataj(:,:,:,:)

  ! Executable Code
  if(mpi_rank.eq.0) allocate(tmparray(np,mp,l,nv,numprocs))
  nsize=np*mp*l*nv
  call mpi_gather(xv(1:np,1:mp,:,:),nsize,mpi_real,tmparray,nsize, &
    mpi_real,0,mpi_comm_world,ierror)

  if(mpi_rank.eq.0)then
    allocate(outdatai(ibcells*2,m,l,nv),outdataj(n,jbcells*2,l,nv))
    open(unit=21,file=fname,form='unformatted',status='unknown')
    do iprocx=1,nprocx
      do jprocy=1,nprocy
        iproc=iprocx+(jprocy-1)*nprocx
        do kv=1,nv
          do k=1,l
            do j=1,mp
              do i=1,np
                ia=(iprocx-1)*np+i
                ja=(jprocy-1)*mp+j
                if(ia.le.ibcells)then
                  outdatai(ia,ja,l,:)=tmparray(i,j,k,:,iproc)
                elseif(ia.gt.n-ibcells)then
                  ii=2*ibcells-n+ia
                  outdatai(ii,ja,l,:)=tmparray(i,j,k,:,iproc)
                endif
                if(ja.le.jbcells)then
                  outdataj(ia,ja,l,:)=tmparray(i,j,k,:,iproc)
                elseif(ja.gt.m-jbcells)then
                  jj=2*jbcells-m+ja
                  outdataj(ia,jj,l,:)=tmparray(i,j,k,:,iproc)
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    write(21) outdatai,outdataj
    deallocate(outdatai,outdataj,tmparray)
  endif
  
end subroutine windfld_write
