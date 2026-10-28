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
!  NB1: FP introduced a fix here for arrays of period np defined on
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
!  (topedge, ...), an additionnal treatment is done using the env
!  field : the ghost cells are forced with xe values in
!  [1,nxp][1,nyp][1,nzp]
!  
!  "iuv" inputs in the subroutine dictate 
!  Dirichlet(1) or Neumann(0) for the top and sides if the boundary 
!  conditions are symmetric (ibc_* = 5).  For symmetric or mirror 
!  purposes, all variables will be Neumann except for the velocities 
!  normal to the particular boundary.
!
!-----------------------------------------------------------------------
subroutine update(a,nxp,nyp,nzp,il,iu,jl,ju,icorner,iperiod,iuv)
  use gridlist_variables, only : ih,wp, &
    ibc_top,ibc_bot,ibc_left,ibc_right
  use msga_variables, only : ierror,mpi_status_size,pebelow,peabove, &
    peleft,peright,mpi_comm_world,botedge,topedge,leftedge, &
    rightedge,peleftbelow,perightbelow,peleftabove,perightabove
  use msga_variables, only : mpi_wp
  Implicit None

  ! Local Variables
  integer,intent(in) :: nxp,nyp,nzp,il,iu,jl,ju,icorner,iperiod,iuv
  real(wp),intent(inout) :: a(il:iu,jl:ju,nzp)

  integer :: i,j,k
  integer :: i0,j0,imirror,jmirror
  integer :: itag,ncount
  real(wp) :: tmpxsnd(nxp*ih*nzp),tmpxrcv(nxp*ih*nzp), &
    tmpysnd(nyp*ih*nzp),tmpyrcv(nyp*ih*nzp), &
    tmpcorsnd(ih*ih*nzp),tmpcorrcv(ih*ih*nzp)
  integer,allocatable :: sndrcvStatus(:)

  ! Executable Code
  i0=0
  j0=0
  if(iperiod.eq.1) i0=1
  if(iperiod.eq.2) j0=1
  imirror=1
  jmirror=1
  if(iuv.eq.1) imirror=-1
  if(iuv.eq.2) jmirror=-1
  allocate(sndrcvStatus(mpi_status_size))

  if(icorner.ne.2)then
!  send bottom points to processor below, receive top ghost 
!  points from processor above
    ncount=0
    do k=1,nzp
      do j=1,ih
        do i=1,nxp
          ncount=ncount+1
          tmpxsnd(ncount)=a(i,j+j0,k)
        enddo
      enddo
    enddo
    itag=1000
    call MPI_Sendrecv(tmpxsnd,ncount,mpi_wp,pebelow,itag,tmpxrcv, &
      ncount,mpi_wp,peabove,itag,mpi_comm_world,sndrcvStatus,ierror)
    if(topedge.ne.1)then
      ncount=0
      do k=1,nzp
        do j=1,ih
          do i=1,nxp
            ncount=ncount+1
            a(i,nyp+j,k)=tmpxrcv(ncount)
          enddo
        enddo
      enddo
    else
      if(ibc_top.eq.5)then ! Mirror
        do k=1,nzp
          do j=1,ih
            do i=1,nxp
              a(i,nyp+j,k)=jmirror*a(i,nyp-j+1,k)
            enddo
          enddo
        enddo
      else  ! Neuman
        do k=1,nzp
          do j=1,ih
            do i=1,nxp
              a(i,nyp+j,k)=a(i,nyp,k)
            enddo
          enddo
        enddo
      endif
    endif

!  send top points to processor above, receive bottom ghost 
!  points from processor below
    ncount=0
    do k=1,nzp
      do j=1,ih
        do i=1,nxp
          ncount=ncount+1
          tmpxsnd(ncount)=a(i,nyp-j0-ih+j,k)
        enddo
      enddo
    enddo
    itag=1100
    call MPI_Sendrecv(tmpxsnd,ncount,mpi_wp,peabove,itag,tmpxrcv, &
      ncount,mpi_wp,pebelow,itag,mpi_comm_world,sndrcvStatus,ierror)
    if(botedge.ne.1)then
      ncount=0
      do k=1,nzp
        do j=1,ih
          do i=1,nxp
            ncount=ncount+1
            a(i,-ih+j,k)=tmpxrcv(ncount)
          enddo
        enddo
      enddo
    else
      if(ibc_bot.eq.5)then ! Mirror
        do k=1,nzp
          do j=1,ih
            do i=1,nxp
              a(i,1-j,k)=jmirror*a(i,j,k)
            enddo
          enddo
        enddo
      else ! Neuman
        do k=1,nzp
          do j=1,ih
            do i=1,nxp
              a(i,1-j,k)=a(i,1,k)
            enddo
          enddo
        enddo
      endif
    endif
  endif
  
  if(icorner.ne.3)then
!  send left points to processor on left, receive right ghost 
!  points from processor on right
    ncount=0
    do k=1,nzp
      do j=1,nyp
        do i=1,ih
          ncount=ncount+1
          tmpysnd(ncount)=a(i+i0,j,k)
        enddo
      enddo
    enddo
    itag=1200
    call MPI_Sendrecv(tmpysnd,ncount,mpi_wp,peleft,itag,tmpyrcv, &
      ncount,mpi_wp,peright,itag,mpi_comm_world,sndrcvStatus,ierror)
    if(rightedge.ne.1)then
      ncount=0
      do k=1,nzp
        do j=1,nyp
          do i=1,ih
            ncount=ncount+1
            a(nxp+i,j,k)=tmpyrcv(ncount)
          enddo
        enddo
      enddo
    else
      if(ibc_right.eq.5)then ! Mirror
        do k=1,nzp
          do j=1,nyp
            do i=1,ih
              a(nxp+i,j,k)=imirror*a(nxp-i+1,j,k)
            enddo
          enddo
        enddo
      else ! Neuman
        do k=1,nzp
          do j=1,nyp
            do i=1,ih
              a(nxp+i,j,k)=a(nxp,j,k)
            enddo
          enddo
        enddo
      endif
    endif  

!  send right points to processor on right, receive left ghost 
!  points from processor on left
    ncount=0
    do k=1,nzp
      do j=1,nyp
        do i=1,ih
          ncount=ncount+1
          tmpysnd(ncount)=a(nxp-i0-ih+i,j,k)
        enddo
      enddo
    enddo
    itag=1300
    call MPI_Sendrecv(tmpysnd,ncount,mpi_wp,peright,itag,tmpyrcv, &
      ncount,mpi_wp,peleft,itag,mpi_comm_world,sndrcvStatus,ierror)
    if(leftedge.ne.1)then
      ncount=0
      do k=1,nzp
        do j=1,nyp
          do i=1,ih
            ncount=ncount+1
            a(-ih+i,j,k)=tmpyrcv(ncount)
          enddo
        enddo
      enddo
    else
      if(ibc_left.eq.5)then ! Mirror
        do k=1,nzp
          do j=1,nyp
            do i=1,ih
              a(1-i,j,k)=imirror*a(i,j,k)
            enddo
          enddo
        enddo
      else ! Neuman
        do k=1,nzp
          do j=1,nyp
            do i=1,ih
              a(1-i,j,k)=a(1,j,k)
            enddo
          enddo
        enddo
      endif
    endif
  endif

!  now send and receive corner pieces if requested
  if(icorner.eq.1)then
!  send right above points to processor on above right, receive  
!  left below points from processor on below left
    ncount=0
    do k=1,nzp
      do j=1,ih
        do i=1,ih
          ncount=ncount+1
          tmpcorsnd(ncount)=a(nxp-i0-ih+i,nyp-j0-ih+j,k)
        enddo
      enddo
    enddo
    itag=1400
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_wp,perightabove,itag, &
      tmpcorrcv,ncount,mpi_wp,peleftbelow,itag,mpi_comm_world, &
      sndrcvStatus,ierror)
    if(botedge.eq.1.and.leftedge.eq.1)then ! Domain bot left corner
      if(ibc_left.eq.5)then ! Mirror left edge, both are edges
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,1-j,k)=imirror*a(i,1-j,k)
            enddo
          enddo
        enddo
      else if(ibc_bot.eq.5)then ! Mirror bot edge, both are edges
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,1-j,k)=jmirror*a(1-i,j,k)
            enddo
          enddo
        enddo
      else ! Both edges are Neuman
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,1-j,k)=a(1,1,k)
            enddo
          enddo
        enddo
      endif
    else if(leftedge.eq.1)then
      if(ibc_left.eq.5)then ! Left side mirror, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,1-j,k)=imirror*a(i,1-j,k)
            enddo
          enddo
        enddo
      else ! Left side Neuman, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,1-j,k)=a(1,1-j,k)
            enddo
          enddo
        enddo
      endif
    else if(botedge.eq.1)then
      if(ibc_bot.eq.5)then ! Bot side mirror, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,1-j,k)=jmirror*a(1-i,j,k)
            enddo
          enddo
        enddo
      else ! Bot side Neuman, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,1-j,k)=a(1-i,1,k)
            enddo
          enddo
        enddo
      endif
    else ! Not an edge of the domain (or cyclic in both directions)
      ncount=0
      do k=1,nzp
        do j=1,ih
          do i=1,ih
            ncount=ncount+1
            a(i-ih,j-ih,k)=tmpcorrcv(ncount)
          enddo
        enddo
      enddo
    endif

!  send right below points to processor on below right, receive  
!  left above points from processor on above left
    ncount=0
    do k=1,nzp
      do j=1,ih
        do i=1,ih
          ncount=ncount+1
          tmpcorsnd(ncount)=a(nxp-i0-ih+i,j+j0,k)
        enddo
      enddo
    enddo
    itag=1500
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_wp,perightbelow,itag, &
      tmpcorrcv,ncount,mpi_wp,peleftabove,itag,mpi_comm_world, &
      sndrcvStatus,ierror)
    if(topedge.eq.1.and.leftedge.eq.1)then ! Domain top left corner
      if(ibc_left.eq.5)then
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,nyp+j,k)=imirror*a(i,nyp+j,k) ! Mirror velocity along left edge
            enddo
          enddo
        enddo
      else if(ibc_top.eq.5)then
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,nyp+j,k)=jmirror*a(1-i,nyp-j+1,k) ! Mirror velocity along top edge
            enddo
          enddo
        enddo
      else
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,nyp+j,k)=a(1,nyp,k) ! Neuman Condition along both edges
            enddo
          enddo
        enddo
      endif
    else if(leftedge.eq.1)then
      if(ibc_left.eq.5)then
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,nyp+j,k)=imirror*a(i,nyp+j,k) ! Mirror velocity along left edge
            enddo
          enddo
        enddo
      else
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,nyp+j,k)=a(1,nyp+j,k) ! Neuman Condition along left edge
            enddo
          enddo
        enddo
      endif
    else if(topedge.eq.1)then
      if(ibc_top.eq.5)then
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,nyp+j,k)=jmirror*a(1-i,nyp-j+1,k) ! Mirror velocity along top edge
            enddo
          enddo
        enddo
      else
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(1-i,nyp+j,k)=a(1-i,nyp,k) ! Neuman Condition along top edge
            enddo
          enddo
        enddo
      endif
    else ! Not an edge of the domain (or cyclic in both directions)
      ncount=0
      do k=1,nzp
        do j=1,ih
          do i=1,ih
            ncount=ncount+1
            a(i-ih,nyp+j,k)=tmpcorrcv(ncount)
          enddo
        enddo
      enddo
    endif

!  send left below points to processor on below left, receive  
!  right above points from processor on above right
    ncount=0
    do k=1,nzp
      do j=1,ih
        do i=1,ih
          ncount=ncount+1
          tmpcorsnd(ncount)=a(i+i0,j+j0,k)
        enddo
      enddo
    enddo
    itag=1600
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_wp,peleftbelow,itag, &
      tmpcorrcv,ncount,mpi_wp,perightabove,itag,mpi_comm_world, &
      sndrcvStatus,ierror)
    if(topedge.eq.1.and.rightedge.eq.1)then ! Domain top right corner
      if(ibc_right.eq.5)then ! Mirror right edge, both are edges
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,nyp+j,k)=imirror*a(nxp-i+1,nyp+j,k)
            enddo
          enddo
        enddo
      else if(ibc_top.eq.5)then ! Mirror top edge, both are edges
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,nyp+j,k)=jmirror*a(nxp+i,nyp-j+1,k)
            enddo
          enddo
        enddo
      else ! Both edges are Neuman
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,nyp+j,k)=a(nxp,nyp,k)
            enddo
          enddo
        enddo
      endif
    else if(rightedge.eq.1)then
      if(ibc_right.eq.5)then ! Right side mirror, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,nyp+j,k)=imirror*a(nxp-i+1,nyp+j,k)
            enddo
          enddo
        enddo
      else ! Right side Neuman, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,nyp+j,k)=a(nxp,nyp+j,k)
            enddo
          enddo
        enddo
      endif
    else if(topedge.eq.1)then ! Top side mirror, one side edge
      if(ibc_top.eq.5)then
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,nyp+j,k)=jmirror*a(nxp+i,nyp-j+1,k)
            enddo
          enddo
        enddo
      else ! Top side Neuman, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,nyp+j,k)=a(nxp+i,nyp,k)
            enddo
          enddo
        enddo
      endif
    else ! Not an edge of the domain (or cyclic in both directions)
      ncount=0
      do k=1,nzp
        do j=1,ih
          do i=1,ih
            ncount=ncount+1
            a(nxp+i,nyp+j,k)=tmpcorrcv(ncount)
          enddo
        enddo
      enddo
    endif

!  send left above points to processor on above left, receive  
!  right below points from processor on below right
    ncount=0
    do k=1,nzp
      do j=1,ih
        do i=1,ih
          ncount=ncount+1
          tmpcorsnd(ncount)=a(i+i0,nyp-ih-j0+j,k)
        enddo
      enddo
    enddo
    itag=1700
    call MPI_Sendrecv(tmpcorsnd,ncount,mpi_wp,peleftabove,itag, &
      tmpcorrcv,ncount,mpi_wp,perightbelow,itag,mpi_comm_world, &
      sndrcvStatus,ierror)
    if(botedge.eq.1.and.rightedge.eq.1)then ! Domain bot right corner
      if(ibc_right.eq.5)then
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,1-j,k)=imirror*a(nxp-i+1,1-j,k)
            enddo
          enddo
        enddo
      else if(ibc_bot.eq.5)then ! Mirror bot edge, both are edges
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,1-j,k)=jmirror*a(nxp+i,j,k)
            enddo
          enddo
        enddo
      else ! Both edges are Neuman
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,1-j,k)=a(nxp,1,k)
            enddo
          enddo
        enddo
      endif
    else if(rightedge.eq.1)then
      if(ibc_right.eq.5)then ! Right side mirror, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,1-j,k)=imirror*a(nxp-i+1,1-j,k)
            enddo
          enddo
        enddo
      else ! Right side Neuman, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,1-j,k)=a(nxp,1-j,k)
            enddo
          enddo
        enddo
      endif
    else if(botedge.eq.1)then
      if(ibc_bot.eq.5)then ! Bot side mirror, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,1-j,k)=jmirror*a(nxp+i,j,k)
            enddo
          enddo
        enddo
      else ! Bot side Neuman, one side edge
        do k=1,nzp
          do j=1,ih
            do i=1,ih
              a(nxp+i,1-j,k)=a(nxp+i,1,k)
            enddo
          enddo
        enddo
      endif
    else ! Not an edge of the domain (or cyclic in both directions)
      ncount=0
      do k=1,nzp
        do j=1,ih
          do i=1,ih
            ncount=ncount+1
            a(nxp+i,j-ih,k)=tmpcorrcv(ncount)
          enddo
        enddo
      enddo
    endif
  endif

  deallocate (sndrcvStatus)

end subroutine update

!-----------------------------------------------------------------------
! This is a helper code for when one has to loop update through the xv
! array.
!-----------------------------------------------------------------------

subroutine loop_update_flags(kv,iuv)
  use xvall, only : iuvel, ivvel
  implicit none

  integer,intent(in)  :: kv
  integer,intent(inout) :: iuv

  if(kv.eq.iuvel)then
    iuv = 1
  elseif(kv.eq.ivvel)then
    iuv = 2
  else
    iuv = 0
  endif
end subroutine loop_update_flags
