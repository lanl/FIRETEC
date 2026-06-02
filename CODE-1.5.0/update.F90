!-----------------------------------------------------------------------
!  This routine takes an arraya of dim (il:iu,jl:ju,nzp)
!  send for i=[1,ih] and [nxp-ih+1,nxp] and
!  receive for i=[1-ih,0] and i=[nxp+1,nxp+ih]
!     if icorner.eq.0 : do i and j ghost cells
!     if icorner.eq.1 : do i and j ghost cells and corner pieces
!     if icorner.eq.2 : do i ghost cells only
!     if icorner.eq.3 : do j ghost cells only
!
!  It also included boundary conditions when relaxation to xe is used
!  otherwhise, bc are hard coded (pr, u1,u2,u3, turbulence...)
!
!  NB1: FP introduced a fix here for arrrays of period np defined on
!  np+1 arrays (for ex advective velocities). the flag iperiod is
!  defined:
!     if iperiod = 0 i0=j0=0
!     if iperiod = 1 i0=1 and j0=0
!     if iperiod = 2 i0=0 and j0=1
!  The use of i0 (or j0 with nyp) is the following:
!  For arrays defined on cell face dim [1-ih, np+1+ih], 
!  il=1-ih,iu=np+1+ih, nxp=np+1 and i0=1
!  the block of data [1+i0,ih+i0] and [np+1+ih-1-i0, np+1-i0] is send
!  receive for i=[1-ih,0] and i=[np+1+1,np+1+ih] which is correct
!  because a(np+1) of current proc is equal to a(1) of peright and
!  a(1) of current proc to a(np+1) of peleft. (period is np and not
!  np+1)
!
!  NB2: when we are at a boundary of the domain and NOT with cyclic bc
!  (topdedge, ...), an additionnal treatment is done using the env
!  field : the ghost cells are forced with xe values in
!  [1,nxp][1,nyp][1,nzp]
!-----------------------------------------------------------------------
subroutine update(a,ae,nxp,nyp,nzp,il,iu,jl,ju,icorner,iperiod)
  use gridlist_variables, only : ih,nprocy,ibclatopen,prec
  use msga_variables, only : ierror,mpi_status_size,pebelow,peabove, &
    peleft,peright,mpi_comm_world,botdedge,topdedge,leftdedge, &
    rightdedge,peleftbelow,perightbelow,peleftabove,perightabove
#if DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use gridsetup, only : timestep
  Implicit None

  ! Local Variables
  integer,intent(in) :: nxp,nyp,nzp,il,iu,jl,ju,icorner,iperiod
  real(prec),intent(in) :: ae(il:iu,jl:ju,nzp)
  real(prec),intent(inout) :: a(il:iu,jl:ju,nzp)

  integer :: i,j,k
  integer :: i0,j0
  integer :: nxbuf,nybuf,ncorbuf
  integer :: icnt,itag,iref,jref,ncount
  real(prec),allocatable :: tmpxsnd(:),tmpxrcv(:),tmpysnd(:), &
    tmpyrcv(:),tmpcorsnd(:),tmpcorrcv(:),sndrcvStatus(:)

  ! Executable Code
  i0=0
  j0=0
  if(iperiod.eq.1) i0=1
  if(iperiod.eq.2) j0=1
  nxbuf=nxp*ih*nzp
  nybuf=nyp*ih*nzp
  ncorbuf=ih*ih*nzp

  allocate(tmpxsnd(nxbuf)); tmpxsnd = 0.0
  allocate(tmpxrcv(nxbuf)); tmpxrcv = 0.0
  allocate(tmpysnd(nybuf)); tmpysnd = 0.0
  allocate(tmpyrcv(nybuf)); tmpyrcv = 0.0
  allocate(tmpcorsnd(ncorbuf)); tmpcorsnd = 0.0
  allocate(tmpcorrcv(ncorbuf)); tmpcorrcv = 0.0
  allocate(sndrcvStatus(mpi_status_size)); sndrcvStatus = 0.0

  if(icorner.ne.2) then
!  send bottom points to processor below, receive top ghost 
!  points from processor above
    if(nprocy.gt.1)then
      icnt=1
      do k=1,nzp
        do j=1,ih
          do i=1,nxp
            tmpxsnd(icnt)=a(i,j+j0,k)
            icnt=icnt+1
          enddo
        enddo
      enddo
      ncount=nxp*ih*nzp
      itag=(timestep*10000)+1000
#ifdef DBL_PREC
      call MPI_Sendrecv(tmpxsnd,ncount,mpi_double_precision,pebelow, &
        itag,tmpxrcv,ncount,mpi_double_precision,peabove,itag, &
        mpi_comm_world,sndrcvStatus,ierror)
#else
      call MPI_Sendrecv(tmpxsnd,ncount,mpi_real,pebelow,itag,tmpxrcv, &
        ncount,mpi_real,peabove,itag,mpi_comm_world,sndrcvStatus,ierror)
#endif
      icnt=1
      do k=1,nzp
        do j=1,ih
          do i=1,nxp
            a(i,nyp+j,k)=tmpxrcv(icnt)
            icnt=icnt+1
          enddo
        enddo
      enddo
!  send top points to processor above, receive bottom ghost 
!  points from processor below
      icnt=1
      do k=1,nzp
        do j=1,ih
          do i=1,nxp
            tmpxsnd(icnt)=a(i,nyp-j0-ih+j,k)
            icnt=icnt+1
          enddo
        enddo
      enddo
      ncount=nxp*ih*nzp
      itag=(timestep*10000)+1100

#ifdef DBL_PREC
      call MPI_Sendrecv(tmpxsnd,ncount,mpi_double_precision,peabove, &
        itag,tmpxrcv,ncount,mpi_double_precision,pebelow,itag, &
        mpi_comm_world,sndrcvStatus,ierror)
#else
      call MPI_Sendrecv(tmpxsnd,ncount,mpi_real,peabove,itag,tmpxrcv, &
        ncount,mpi_real,pebelow,itag,mpi_comm_world,sndrcvStatus,ierror)
#endif

      icnt=1
      do k=1,nzp
        do j=1,ih
          do i=1,nxp
            a(i,-ih+j,k)=tmpxrcv(icnt)
            icnt=icnt+1
          enddo
        enddo
      enddo
    elseif(nprocy.eq.1)then !JMC case where there is 1 y-processor
      do j=1,ih
        a(:,nyp+j,:)=a(:,j+j0,:)
        a(:,1-j,:)=a(:,nyp-j0+1-j,:)
        
      enddo
    endif
  endif

  if(icorner.ne.3) then
!  send left points to processor on left, receive right ghost 
!  points from processor on right
    icnt=1
    do k=1,nzp
      do j=1,nyp
        do i=1,ih
          tmpysnd(icnt)=a(i+i0,j,k)
          icnt=icnt+1
        enddo
      enddo
    enddo
    ncount=nyp*ih*nzp
    itag=(timestep*10000)+1200

#ifdef DBL_PREC
    call MPI_Sendrecv(tmpysnd,ncount,mpi_double_precision,peleft, &
      itag,tmpyrcv,ncount,mpi_double_precision,peright,itag, &
      mpi_comm_world,sndrcvStatus,ierror)
#else
    call MPI_Sendrecv(tmpysnd,ncount,mpi_real,peleft,itag,tmpyrcv, &
      ncount,mpi_real,peright,itag,mpi_comm_world,sndrcvStatus,ierror)
#endif

    icnt=1
    do k=1,nzp
      do j=1,nyp
        do i=1,ih
            a(nxp+i,j,k)=tmpyrcv(icnt)
            icnt=icnt+1
         enddo
       enddo
    enddo
!  send right points to processor on right, receive left ghost 
!  points from processor on left
    icnt=1
    do k=1,nzp
      do j=1,nyp
        do i=1,ih
          tmpysnd(icnt)=a(nxp-i0-ih+i,j,k)
          icnt=icnt+1
        enddo
      enddo
    enddo
    ncount=nyp*ih*nzp
    itag=(timestep*10000)+1300

#ifdef DBL_PREC
    call MPI_Sendrecv(tmpysnd,ncount,mpi_double_precision,peright, &
      itag,tmpyrcv,ncount,mpi_double_precision,peleft,itag, &
      mpi_comm_world,sndrcvStatus,ierror)
#else
    call MPI_Sendrecv(tmpysnd,ncount,mpi_real,peright,itag,tmpyrcv, &
      ncount,mpi_real,peleft,itag,mpi_comm_world,sndrcvStatus,ierror)
#endif

    icnt=1
    do k=1,nzp
      do j=1,nyp
        do i=1,ih
          a(-ih+i,j,k)=tmpyrcv(icnt)
          icnt=icnt+1
        enddo
      enddo
    enddo
  endif
!  now send and receive corner pieces if requested
  if(icorner.eq.1)then
!  send right above points to processor on above right, receive  
!  left below points from processor on below left
    icnt=1
    do k=1,nzp
      do j=1,ih
        do i=1,ih
          tmpcorsnd(icnt)=a(nxp-i0-ih+i,nyp-j0-ih+j,k)
          icnt=icnt+1
        enddo
      enddo
    enddo
    ncount=ih*ih*nzp
    itag=(timestep*10000)+1400

#ifdef DBL_PREC
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_double_precision,perightabove, &
      itag,tmpcorrcv,ncount,mpi_double_precision,peleftbelow,itag, &
      mpi_comm_world,sndrcvStatus,ierror)
#else
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,perightabove,itag,tmpcorrcv, &
      ncount,mpi_real,peleftbelow,itag,mpi_comm_world,sndrcvStatus,ierror)
#endif

    if(botdedge.ne.1.or.leftdedge.ne.1)then
      icnt=1
      do k=1,nzp
        do j=1,ih
          do i=1,ih
            a(i-ih,j-ih,k)=tmpcorrcv(icnt)
            icnt=icnt+1
          enddo
        enddo
      enddo
    endif
!  send right below points to processor on below right, receive  
!  left above points from processor on above left
    icnt=1
    do k=1,nzp
      do j=1,ih
        do i=1,ih
          tmpcorsnd(icnt)=a(nxp-i0-ih+i,j+j0,k)
          icnt=icnt+1
        enddo
      enddo
    enddo
    ncount=ih*ih*nzp
    itag=(timestep*10000)+1500

#ifdef DBL_PREC
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_double_precision,perightbelow, &
      itag,tmpcorrcv,ncount,mpi_double_precision,peleftabove,itag, &
      mpi_comm_world,sndrcvStatus,ierror)
#else
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,perightbelow,itag,tmpcorrcv, &
      ncount,mpi_real,peleftabove,itag,mpi_comm_world,sndrcvStatus,ierror)
#endif

    if (topdedge.ne.1.or.leftdedge.ne.1) then
      icnt=1
      do k=1,nzp
        do j=1,ih
          do i=1,ih
            a(i-ih,nyp+j,k)=tmpcorrcv(icnt)
            icnt=icnt+1
          enddo
        enddo
      enddo
    endif
!  send left below points to processor on below left, receive  
!  right above points from processor on above right
    icnt=1
    do k=1,nzp
      do j=1,ih
        do i=1,ih
          tmpcorsnd(icnt)=a(i+i0,j+j0,k)
          icnt=icnt+1
        enddo
      enddo
    enddo
    ncount=ih*ih*nzp
    itag=(timestep*10000)+1600
#ifdef DBL_PREC
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_double_precision,peleftbelow, &
      itag,tmpcorrcv,ncount,mpi_double_precision,perightabove,itag, &
      mpi_comm_world,sndrcvStatus,ierror)
#else
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,peleftbelow,itag,tmpcorrcv, &
      ncount,mpi_real,perightabove,itag,mpi_comm_world,sndrcvStatus,ierror)
#endif
    if (topdedge.ne.1.or.rightdedge.ne.1) then
      icnt=1
      do k=1,nzp
        do j=1,ih
          do i=1,ih
            a(nxp+i,nyp+j,k)=tmpcorrcv(icnt)
            icnt=icnt+1
          enddo
        enddo
      enddo
    endif
!  send left above points to processor on above left, receive  
!  right below points from processor on below right
    icnt=1
    do k=1,nzp
      do j=1,ih
        do i=1,ih
          tmpcorsnd(icnt)=a(i+i0,nyp-j0-ih+j,k)
          icnt=icnt+1
        enddo
      enddo
    enddo
    ncount=ih*ih*nzp
    itag=(timestep*10000)+1700
#ifdef DBL_PREC
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_double_precision,peleftabove, &
      itag,tmpcorrcv,ncount,mpi_double_precision,perightbelow,itag, &
      mpi_comm_world,sndrcvStatus,ierror)
#else
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,peleftabove,itag,tmpcorrcv, &
      ncount,mpi_real,perightbelow,itag,mpi_comm_world,sndrcvStatus,ierror)
#endif
    if (botdedge.ne.1.or.rightdedge.ne.1) then
      icnt=1
      do k=1,nzp
        do j=1,ih
          do i=1,ih
            a(nxp+i,j-ih,k)=tmpcorrcv(icnt)
            icnt=icnt+1
          enddo
        enddo
      enddo
    endif
  endif
  !Go back and fill the global domain edge halos as appropriate for non cyclic at edges
  if(topdedge.eq.1)then
    if(ibclatopen.eq.0)then
      do k=1,nzp
        do j=1,ih
          do i=1-ih,nxp+ih
            iref=max(1,min(i,nxp))
            a(i,nyp+j,k)=ae(iref,nyp,k)
          enddo
        enddo
      enddo
    else  ! ibclatopen>0, not implemented here
      stop
    endif
  endif
  if(botdedge.eq.1)then
    if(ibclatopen.eq.0)then
      do k=1,nzp
        do j=1,ih
          do i=1-ih,nxp+ih
            iref=max(1,min(i,nxp))
            a(i,1-j,k)=ae(iref,1,k)
          enddo
        enddo
      enddo
    else  ! ibclatopen>0, not implemented here
      stop
    endif
  endif
  if(rightdedge.eq.1)then
    if(ibclatopen.eq.0)then
      do k=1,nzp
        do j=1-ih,nyp+ih
          do i=1,ih
            jref=max(1,min(j,nyp))
            a(nxp+i,j,k)=ae(nxp,jref,k)
          enddo
        enddo
      enddo
    else  ! ibclatopen>0, not implemented here
      stop
    endif
  endif  
  if(leftdedge.eq.1)then
    do k=1,nzp
      do j=1-ih,nyp+ih
        do i=1,ih
          jref=max(1,min(j,nyp))
          a(1-i,j,k)=ae(1,jref,k)
        enddo
      enddo
    enddo
  endif


  deallocate (tmpxsnd) 
  deallocate (tmpxrcv) 
  deallocate (tmpysnd) 
  deallocate (tmpyrcv) 
  deallocate (tmpcorsnd) 
  deallocate (tmpcorrcv)   
  deallocate (sndrcvStatus)

end subroutine update
