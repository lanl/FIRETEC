!-----------------------------------------------------------------------
! Windfield master file that governs the outputting of windfields only
!-----------------------------------------------------------------------
subroutine windField(ittot)
  use windfieldio, only : itwindfield,itinterp,xvdataname
  Implicit None

  ! Local Variables
  integer,intent(in) :: ittot

  character(len=257) :: fname

  ! Executable Code
  if(ittot.ge.itwindfield.and.mod(ittot-itwindfield,itinterp).eq.0)then
    call namefile(ittot-itwindfield,xvdataname,fname)
    call windfld_write(fname)
  endif
  if(ittot.eq.itwindfield) call windFldRstWrite

end subroutine windField

!-----------------------------------------------------------------------
! windfld_write writes xv data along domain boundaries from a wind field
! to buffer those boundaries in future runs
!-----------------------------------------------------------------------
subroutine windfld_write(fname)
  use gridlist_variables, only : n,m,l,nprocx,nprocy,wp
  use windfieldio, only : ibcells,jbcells
  use msga_variables, only : mpi_rank,numprocs,ierror,mpi_comm_world
  use msga_variables, only : mpi_wp
  use gridsetup, only : np,mp
  use xvall, only : xv,nv
  Implicit None

  ! Local Variables
  character(*),intent(in) :: fname

  integer :: i,j,ii,jj
  integer :: iproc,iprocx,jprocy,ia,ja
  integer :: nsize
  real(wp),allocatable :: tmparray(:,:,:,:,:)
  real(wp),allocatable :: outdatai(:,:,:,:),outdataj(:,:,:,:)

  ! Executable Code
  allocate(tmparray(np,mp,l,nv,numprocs))
  nsize=np*mp*l*nv
  call mpi_gather(xv(1:np,1:mp,:,:),nsize,mpi_wp,tmparray,nsize, &
    mpi_wp,0,mpi_comm_world,ierror)

  if(mpi_rank.eq.0)then
    allocate(outdatai(ibcells*2,m,l,nv),outdataj(n,jbcells*2,l,nv))
    open(unit=21,file=fname,form='unformatted',status='replace', &
      action='write')
    do iprocx=1,nprocx
      do jprocy=1,nprocy
        iproc=iprocx+(jprocy-1)*nprocx
        do j=1,mp
          do i=1,np
            ia=(iprocx-1)*np+i
            ja=(jprocy-1)*mp+j
            if(ia.le.ibcells)then
              outdatai(ia,ja,:,:)=tmparray(i,j,:,:,iproc)
            elseif(ia.gt.n-ibcells)then
              ii=2*ibcells-n+ia
              outdatai(ii,ja,:,:)=tmparray(i,j,:,:,iproc)
            endif
            if(ja.le.jbcells)then
              outdataj(ia,ja,:,:)=tmparray(i,j,:,:,iproc)
            elseif(ja.gt.m-jbcells)then
              jj=2*jbcells-m+ja
              outdataj(ia,jj,:,:)=tmparray(i,j,:,:,iproc)
            endif
          enddo
        enddo
      enddo
    enddo
    write(21) outdatai,outdataj
    close(21)
    deallocate(outdatai,outdataj)
  endif
  deallocate(tmparray)
  
end subroutine windfld_write

!-----------------------------------------------------------------------
! windfld_read reads xv data along domain boundaries from a wind field
! to buffer those boundaries
!-----------------------------------------------------------------------
subroutine windfld_read(fname,xvdata)
  use gridlist_variables, only : n,m,l,nprocx,nprocy,wp
  use windfieldio, only : ibcells,jbcells,nvwind
  use msga_variables, only : npos,mpos
  use gridsetup, only : np,mp
  Implicit None

  ! Local Variables
  character(*),intent(in) :: fname
  real(wp),intent(inout) :: xvdata(np,mp,l,nvwind)
  real(wp),allocatable :: indatai(:,:,:,:),indataj(:,:,:,:)

  ! Executable Code
  allocate(indatai(ibcells*2,m,l,nvwind),indataj(n,jbcells*2,l,nvwind))
  open(41,file=fname,form='unformatted',status='old')
  read(41) indatai,indataj
  close(41)
  if(npos.eq.1) xvdata(1:ibcells,:,:,:)= &
    indatai(1:ibcells,mp*(mpos-1)+1:mp*mpos,:,:)
  if(npos.eq.nprocx) xvdata(np-ibcells+1:np,:,:,:)= &
    indatai(ibcells+1:ibcells*2,mp*(mpos-1)+1:mp*mpos,:,:)
  if(mpos.eq.1) xvdata(:,1:jbcells,:,:)= &
    indataj(np*(npos-1)+1:np*npos,1:jbcells,:,:)
  if(mpos.eq.nprocy) xvdata(:,mp-jbcells+1:mp,:,:)= &
    indataj(np*(npos-1)+1:np*npos,jbcells+1:jbcells*2,:,:)
  deallocate(indatai,indataj)

end subroutine windfld_read

!-----------------------------------------------------------------------
! windFldRstWrite writes unique xv data for a fire run restarted from a 
! wind run
!-----------------------------------------------------------------------
subroutine windFldRstWrite
  use gridlist_variables, only : l
  use windfieldio, only : windfieldstartfile
  use msga_variables, only : mpi_rank
  use gridsetup, only : np,mp
  use xvall, only : xv,xv_list,nv
  Implicit None

  ! Local Variables
  integer :: kv

  ! Executable Code
  ! access='stream' skips fortran-specific headers and trailers
  ! TODO: parallelize?
  if(mpi_rank.eq.0)then
    open(48,file=windfieldstartfile, &
      form='unformatted',access='stream',status='unknown')
    write(48) nv
    write(48) xv_list
  endif
  do kv=1,nv  
    call writeio(xv(1:np,1:mp,:,kv),48,l)
  enddo
  if(mpi_rank.eq.0) close(48)

end subroutine windFldRstWrite
