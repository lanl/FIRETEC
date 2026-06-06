c***********************************************************
!      This routine takes an arraya of dim (il:iu,jl:ju,nzp)
!      send for i=[1,ih] and [nxp-ih+1,nxp] and
!      receive for i=[1-ih,0] and i=[nxp+1,nxp+ih]
!     if icorner.eq.0 : do i and j ghost cells
!     if icorner.eq.1 : do i and j ghost cells and corner pieces
!     if icorner.eq.2 : do i ghost cells only
!     if icorner.eq.3 : do j ghost cells only
!
!
!     It also included boundary conditions when relaxation to xe is used
!     otherwhise, bc are hard coded (pr, u1,u2,u3, turbulence...)
!
!     NB1: FP introduced a fix here for arrrays of period np defined on
!     np+1 arrays (for ex advective velocities). the flag iperiod is
!     defined:
!      if iperiod = 0 i0=j0=0
!      if iperiod = 1 i0=1 and j0=0
!      if iperiod = 0 i0=0 and j0=1
!      the use of i0 (or j0 with nyp) is the following:For arrays defined on cell face dim [1-ih, np+1+ih], 
!     il=1-ih,iu=np+1+ih, nxp=np+1 and i0=1
!     the block of data [1+i0,ih+i0] and [np+1+ih-1-i0, np+1-i0] is send
!      receive for i=[1-ih,0] and i=[np+1+1,np+1+ih] which is correct
!     because a(np+1) of current proc is equal to a(1) of peright and
!     a(1) of current proc to a(np+1) of peleft. (period is np and not
!     np+1)
!
!     NB2: when we are at a boundary of the domain and NOT with cyclic bc
!     (topdedge, ...), an additionnal treatment is done using the env
!     field : the ghost cells are forced with xe values in
!     [1,nxp][1,nyp][1,nzp]
!     FP fixed the corner pieces that were not done propertly with xe
!
!     NOTES FOR 2D RUNS:
!     works with nprocy.eq.1 and mp.eq.1 (2d), 
!     but not efficiently when with nprocx.eq.1 (sendrecev to himself)
!**********************************************************
      subroutine updated(a,ae,nxp,nyp,nzp,il,iu,jl,ju,icorner, iperiod)
      use gridsetup
      use msga
      implicit none
      integer:: nxp,nyp,nzp,il,iu,jl,ju,icorner, iperiod
      real a(il:iu,jl:ju,nzp) ,ae(il:iu,jl:ju,nzp)
      integer,allocatable:: status(:)
      real,allocatable:: tmpxsnd(:),tmpxrcv(:),
     +                   tmpysnd(:),tmpyrcv(:),
     +                   tmpcorsnd(:),tmpcorrcv(:)
      integer:: i,j,k,icnt,itag,iref,jref
      integer:: nxbuf,nybuf,ncorbuf,ncount,ierr
      integer:: i0, j0
      i0=0
      j0=0
      if (iperiod.eq.1) i0=1
      if (iperiod.eq.2) j0=1
!     parameter (nxbuf=(n+2)*ih*(l+2), nybuf=(m+2)*ih*(l+2),
!    +           ncorbuf=ih*ih*(l+2))

      ! FPFIX: lnz is radiation grid!!!!
      !nxbuf=(n+2)*ih*(lnz+2) 
      !nybuf=(m+2)*ih*(lnz+2) !wss
      !ncorbuf=ih*ih*(lnz+2)
      
      nxbuf = nxp*ih*nzp
      nybuf = nyp*ih*nzp
      ncorbuf = ih*ih*nzp

      allocate (tmpxsnd(nxbuf)) 
      allocate (tmpxrcv(nxbuf)) 
      allocate (tmpysnd(nybuf)) 
      allocate (tmpyrcv(nybuf)) 
      allocate (tmpcorsnd(ncorbuf)) 
      allocate (tmpcorrcv(ncorbuf)) 
      allocate (status(mpi_status_size))
      if(icorner.ne.2) then
!  send bottom points to processor below, receive top ghost 
!  points from processor above
        if(nprocy.gt.1)then !JMC Only call mpi routines when necessary
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
          call MPI_Sendrecv(tmpxsnd,ncount,mpi_real,pebelow,itag,
     +                  tmpxrcv,ncount,mpi_real,peabove,itag,
     +                  mpi_comm_world,status,ierr)
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
          call MPI_Sendrecv(tmpxsnd,ncount,mpi_real,peabove,itag,
     +                  tmpxrcv,ncount,mpi_real,pebelow,itag,
     +                  mpi_comm_world,status,ierr)
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
          if(mp.ne.1)then   !JMC 1 y-processor but m .gt.1
            do j=1,ih
              a(:,nyp+j,:)=a(:,j+j0,:)
              a(:,1-j,:)=a(:,nyp-j0+1-j,:)
            enddo
          else   !JMC 1 y-processor but m .eq.1 or 2-D case
            do j=1,ih
              a(:,1+j,:)=a(:,1,:)
              a(:,1-j,:)=a(:,1,:)
            enddo
          endif
        else
          write(6,*)'In updated, nprocy is less than 0'
          write(6,*)'Your job is antiparallel: Terminated'
          call mpi_finalize()
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
          call MPI_Sendrecv(tmpysnd,ncount,mpi_real,peleft,itag,
     +                  tmpyrcv,ncount,mpi_real,peright,itag,
     +                  mpi_comm_world,status,ierr)
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
          call MPI_Sendrecv(tmpysnd,ncount,mpi_real,peright,itag,
     +                  tmpyrcv,ncount,mpi_real,peleft,itag,
     +                  mpi_comm_world,status,ierr)

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
      if (icorner.eq.1) then
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
        call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,perightabove,itag,
     +                     tmpcorrcv,ncount,mpi_real,peleftbelow,itag,
     +                     mpi_comm_world,status,ierr)
        if (botdedge.ne.1.or.leftdedge.ne.1) then
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
        call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,perightbelow,itag,
     +                    tmpcorrcv,ncount,mpi_real,peleftabove,itag,
     +                    mpi_comm_world,status,ierr)
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
        call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,peleftbelow,itag,
     +                     tmpcorrcv,ncount,mpi_real,perightabove,itag,
     +                     mpi_comm_world,status,ierr)
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
        call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,peleftabove,itag,
     +                     tmpcorrcv,ncount,mpi_real,perightbelow,itag,
     +                     mpi_comm_world,status,ierr)
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
       !Go back and fill the global domain edge halos as appropriate for
       !non cyclic at edges
          if (topdedge.eq.1) then
           if (ibclatopen.eq.0) then
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
          if (botdedge.eq.1) then
           if (ibclatopen.eq.0) then
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
          if (rightdedge.eq.1) then
           if (ibclatopen.eq.0) then
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
          if (leftdedge.eq.1) then
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
      deallocate (status)
      return

      end
