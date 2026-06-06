c***********************************************************
      subroutine updated2(a,nxp,nyp,nzp,il,iu,jl,ju,icorner)  !ae removed by KOO added d adn ae rrl
      use gridsetup
      use msga
c
      implicit none
c     this routine updates the halo (or ghost) cells surrounding
c     each processors subgrid
c
      integer:: nxp,nyp,nzp,il,iu,jl,ju,icorner
      real a(il:iu,jl:ju,nzp),abc !,ae(il:iu,jl:ju,nzp),abc  !ae removed by KOO added ae and abc rrl
      integer,allocatable:: status(:)
      real,allocatable:: tmpxsnd(:),tmpxrcv(:),
     +                   tmpysnd(:),tmpyrcv(:),
     +                   tmpcorsnd(:),tmpcorrcv(:)
      integer:: i,j,k,icnt,itag ! ,i2d3d,iref,jref
      integer:: nxbuf,nybuf,ncorbuf,ncount,ierr
      abc=1.
      nxbuf=(nxp)*ih*(nzp) 
      nybuf=(nyp)*ih*(nzp) 
      ncorbuf=ih*ih*(nzp)
      allocate (tmpxsnd(nxbuf)) 
      allocate (tmpxrcv(nxbuf)) 
      allocate (tmpysnd(nybuf)) 
      allocate (tmpyrcv(nybuf)) 
      allocate (tmpcorsnd(ncorbuf)) 
      allocate (tmpcorrcv(ncorbuf)) 
      allocate (status(mpi_status_size))
      if(icorner.ne.2) then
c
c  send bottom points to processor below, receive top ghost 
c  points from processor above
        if(nprocy.gt.1)then !JMC Only call mpi routines when necessary
          icnt=1
          do k=1,nzp
            do j=1,ih
              do i=1,nxp
                 tmpxsnd(icnt)=a(i,j,k)
                 icnt=icnt+1
              enddo
            enddo
          enddo
c
          ncount=nxp*ih*nzp
          itag=(timestep*10000)+1000
          call MPI_Sendrecv(tmpxsnd,ncount,mpi_real,pebelow,itag,
     +                  tmpxrcv,ncount,mpi_real,peabove,itag,
     +                  mpi_comm_world,status,ierr)
c
           if(topdedge.ne.1)then
            icnt=1
            do k=1,nzp
              do j=1,ih
                do i=1,nxp
                  a(i,nyp+j,k)=tmpxrcv(icnt)
                  icnt=icnt+1
                enddo
              enddo
            enddo
           endif
c
c  send top points to processor above, receive bottom ghost 
c  points from processor below
          icnt=1
          do i=1,nxp
            do j=1,ih
              do k=1,nzp
                tmpxsnd(icnt)=a(i,nyp-ih+j,k)
                icnt=icnt+1
              enddo
            enddo
          enddo
c
          ncount=nxp*ih*nzp
          itag=(timestep*10000)+1100
          call MPI_Sendrecv(tmpxsnd,ncount,mpi_real,peabove,itag,
     +                  tmpxrcv,ncount,mpi_real,pebelow,itag,
     +                  mpi_comm_world,status,ierr)
c
           if(botdedge.ne.1)then 
            icnt=1
            do i=1,nxp
              do j=1,ih
                do k=1,nzp
                  a(i,-ih+j,k)=tmpxrcv(icnt)
                  icnt=icnt+1
                enddo
              enddo
            enddo
           endif
        
        elseif(nprocy.eq.1)then !JMC case where there is 1 y-processor
          if(mp.ne.1)then   !JMC 1 y-processor but m .gt.1
           if(topdedge.eq.0)then
            do j=1,ih
              a(:,nyp+j,:)=a(:,j,:)
              a(:,1-j,:)=a(:,nyp+1-j,:)
            enddo
           endif
          else   !JMC 1 y-processor but m .eq.1 or 2-D case
           if(topdedge.eq.0)then
            do j=1,ih
              a(:,1+j,:)=a(:,1,:)
              a(:,1-j,:)=a(:,1,:)
            enddo
           endif
          endif
        else
          write(6,*)'In updated, nprocy is less than 0'
          write(6,*)'Your job is antiparallel: Terminated'
          call mpi_finalize()
        endif
      endif

      if(icorner.ne.3) then
c
c  send left points to processor on left, receive right ghost 
c  points from processor on right
          icnt=1
          do i=1,ih
            do j=1,nyp
              do k=1,nzp
                tmpysnd(icnt)=a(i,j,k)
                icnt=icnt+1
              enddo
            enddo
          enddo
c
          ncount=nyp*ih*nzp
          itag=(timestep*10000)+1200
          call MPI_Sendrecv(tmpysnd,ncount,mpi_real,peleft,itag,
     +                  tmpyrcv,ncount,mpi_real,peright,itag,
     +                  mpi_comm_world,status,ierr)
ci
           if(rightdedge.ne.1)then
            icnt=1
            do i=1,ih
              do j=1,nyp
                do k=1,nzp
                  a(nxp+i,j,k)=tmpyrcv(icnt)
                  icnt=icnt+1
                enddo
              enddo
            enddo
           endif
c
c  send right points to processor on right, receive left ghost 
c  points from processor on left
          icnt=1
          do i=1,ih
            do j=1,nyp
              do k=1,nzp
                tmpysnd(icnt)=a(nxp-ih+i,j,k)
                icnt=icnt+1
              enddo
            enddo
          enddo
c
          ncount=nyp*ih*nzp
          itag=(timestep*10000)+1300
          call MPI_Sendrecv(tmpysnd,ncount,mpi_real,peright,itag,
     +                  tmpyrcv,ncount,mpi_real,peleft,itag,
     +                  mpi_comm_world,status,ierr)
c
           if(leftdedge.ne.1)then
            icnt=1
            do i=1,ih
              do j=1,nyp
                do k=1,nzp
                  a(-ih+i,j,k)=tmpyrcv(icnt)
                  icnt=icnt+1
                enddo
              enddo
            enddo
           endif
      endif
c
c  now send and receive corner pieces if requested
      if (icorner.eq.1) then
c
c  send right above points to processor on above right, receive  
c  left below points from processor on below left
        icnt=1
        do i=1,ih
          do j=1,ih
            do k=1,nzp
              tmpcorsnd(icnt)=a(nxp-ih+i,nyp-ih+j,k)
              icnt=icnt+1
            enddo
          enddo
        enddo
c
        ncount=ih*ih*nzp
        itag=(timestep*10000)+1400
        call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,perightabove,itag,
     +                     tmpcorrcv,ncount,mpi_real,peleftbelow,itag,
     +                     mpi_comm_world,status,ierr)
c
        if (botdedge.ne.1.or.leftdedge.ne.1) then
          icnt=1
          do i=1,ih
            do j=1,ih
              do k=1,nzp
                a(i-ih,j-ih,k)=tmpcorrcv(icnt)
                icnt=icnt+1
              enddo
            enddo
          enddo
        endif
c
c  send right below points to processor on below right, receive  
c  left above points from processor on above left
        icnt=1
        do i=1,ih
          do j=1,ih
            do k=1,nzp
              tmpcorsnd(icnt)=a(nxp-ih+i,j,k)
              icnt=icnt+1
            enddo
          enddo
        enddo
c
        ncount=ih*ih*nzp
        itag=(timestep*10000)+1500
        call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,perightbelow,itag,
     +                    tmpcorrcv,ncount,mpi_real,peleftabove,itag,
     +                    mpi_comm_world,status,ierr)
c
        if (topdedge.ne.1.or.leftdedge.ne.1) then
          icnt=1
          do i=1,ih
            do j=1,ih
              do k=1,nzp
                a(i-ih,nyp+j,k)=tmpcorrcv(icnt)
                icnt=icnt+1
              enddo
            enddo
          enddo
        endif
c
c  send left below points to processor on below left, receive  
c  right above points from processor on above right
        icnt=1
        do i=1,ih
          do j=1,ih
            do k=1,nzp
              tmpcorsnd(icnt)=a(i,j,k)
              icnt=icnt+1
            enddo
          enddo
        enddo
c
        ncount=ih*ih*nzp
        itag=(timestep*10000)+1600
        call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,peleftbelow,itag,
     +                     tmpcorrcv,ncount,mpi_real,perightabove,itag,
     +                     mpi_comm_world,status,ierr)
c
        if (topdedge.ne.1.or.rightdedge.ne.1) then
          icnt=1
          do i=1,ih
            do j=1,ih
              do k=1,nzp
                a(nxp+i,nyp+j,k)=tmpcorrcv(icnt)
                icnt=icnt+1
              enddo
            enddo
          enddo
        endif
c 
c  send left above points to processor on above left, receive  
c  right below points from processor on below right
        icnt=1
        do i=1,ih
          do j=1,ih
            do k=1,nzp
              tmpcorsnd(icnt)=a(i,nyp-ih+j,k)
              icnt=icnt+1
            enddo
          enddo
        enddo
c
        ncount=ih*ih*nzp
        itag=(timestep*10000)+1700
        call MPI_Sendrecv(tmpcorsnd,ncount,mpi_real,peleftabove,itag,
     +                     tmpcorrcv,ncount,mpi_real,perightbelow,itag,
     +                     mpi_comm_world,status,ierr)
c
        if (botdedge.ne.1.or.rightdedge.ne.1) then
          icnt=1
          do i=1,ih
            do j=1,ih
              do k=1,nzp
                a(nxp+i,j-ih,k)=tmpcorrcv(icnt)
                icnt=icnt+1
              enddo
            enddo
          enddo
        endif
c
      endif
c


      deallocate (tmpxsnd)
      deallocate (tmpxrcv)
      deallocate (tmpysnd)
      deallocate (tmpyrcv)
      deallocate (tmpcorsnd)
      deallocate (tmpcorrcv)
      deallocate (status)
      return

      end
