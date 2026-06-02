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
  if(ittot.ge.itwindfield.and.mod(ittot,itinterp).eq.0)then
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
  use gridlist_variables, only : n,m,l,nprocx,nprocy,prec
  use windfieldio, only : ibcells,jbcells
  use msga_variables, only : mpi_rank,numprocs,ierror,mpi_comm_world
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use gridsetup, only : np,mp
  use xvall, only : xv,nv
  Implicit None

  ! Local Variables
  character(*),intent(in) :: fname

  integer :: i,j,ii,jj
  integer :: iproc,iprocx,jprocy,ia,ja
  integer :: nsize
  real(prec),allocatable :: tmparray(:,:,:,:,:)
  real(prec),allocatable :: outdatai(:,:,:,:),outdataj(:,:,:,:)

  ! Executable Code
  if(mpi_rank.eq.0) allocate(tmparray(np,mp,l,nv,numprocs))
  nsize=np*mp*l*nv
#ifdef DBL_PREC
  call mpi_gather(xv(1:np,1:mp,:,:),nsize,mpi_double_precision, &
    tmparray,nsize,mpi_double_precision,0,mpi_comm_world,ierror)
#else
  call mpi_gather(xv(1:np,1:mp,:,:),nsize,mpi_real,tmparray,nsize, &
    mpi_real,0,mpi_comm_world,ierror)
#endif

  if(mpi_rank.eq.0)then
    allocate(outdatai(ibcells*2,m,l,nv),outdataj(n,jbcells*2,l,nv))
    outdatai = 0.0
    outdataj = 0.0
    open(unit=21,file=fname,form='unformatted',status='unknown')
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
    deallocate(outdatai,outdataj,tmparray)
  endif
  
end subroutine windfld_write

!-----------------------------------------------------------------------
! windfld_read reads xv data along domain boundaries from a wind field
! to buffer those boundaries
!-----------------------------------------------------------------------
subroutine windfld_read(fname)
  use gridlist_variables, only : n,m,l,nprocx,nprocy,prec
  use windfieldio, only : ibcells,jbcells,xvdatanew,nvwind
  use msga_variables, only : mpi_rank,numprocs,ierror,mpi_comm_world
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use gridsetup, only : np,mp
  Implicit None

  ! Local Variables
  character(*),intent(in) :: fname

  integer :: i,j,ii,jj
  integer :: iproc,iprocx,jprocy,ia,ja
  integer :: nsize
  real(prec),allocatable :: tmparray(:,:,:,:,:),indata(:,:,:,:)
  real(prec),allocatable :: indatai(:,:,:,:),indataj(:,:,:,:)

  ! Executable Code
  if(mpi_rank.eq.0)then
    allocate(indata(n,m,l,nvwind)); indata=0.
    allocate(tmparray(np,mp,l,nvwind,numprocs), &
      indatai(ibcells*2,m,l,nvwind),indataj(n,jbcells*2,l,nvwind))
    tmparray = 0.0
    indatai = 0.0
    indataj = 0.0
    open(41,file=fname,form='unformatted',status='old')
    read(41) indatai,indataj
    close(41)
    do j=1,m
      do i=1,n
        if(i.le.ibcells)then
          indata(i,j,:,:)=indatai(i,j,:,:)
        elseif(i.gt.n-ibcells)then
          ii=i-(n-ibcells*2)
          indata(i,j,:,:)=indatai(ii,j,:,:)
        else
          indata(i,j,:,:)=0.
        endif
        if(j.le.jbcells)then
          indata(i,j,:,:)=indataj(i,j,:,:)
        elseif(j.gt.m-jbcells)then
          jj=j-(m-jbcells*2)
          indata(i,j,:,:)=indataj(i,jj,:,:)
        endif
      enddo
    enddo
    do iprocx=1,nprocx
      do jprocy=1,nprocy
        iproc=iprocx+(jprocy-1)*nprocx
        do j=1,mp
          do i=1,np
            ia=(iprocx-1)*np+i
            ja=(jprocy-1)*mp+j
            tmparray(i,j,:,:,iproc)=indata(ia,ja,:,:)
          enddo
        enddo
      enddo
    enddo
  endif

  nsize=np*mp*l*nvwind
#ifdef DBL_PREC
  call mpi_scatter(tmparray,nsize,mpi_double_precision, &
    xvdatanew,nsize,mpi_double_precision,0,mpi_comm_world,ierror)
#else
  call mpi_scatter(tmparray,nsize,mpi_real,xvdatanew,nsize, &
    mpi_real,0,mpi_comm_world,ierror)
#endif
  
  if(mpi_rank.eq.0) deallocate(indatai,indataj,indata,tmparray)

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
