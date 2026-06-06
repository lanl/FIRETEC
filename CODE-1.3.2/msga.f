c2345678***************************************************
      module msga

      use gridsetup
      Implicit None
      include 'mpif.h'

      integer mpi_rank,timestep,npos,mpos,ierror,numprocs,blocktype,blocktype0
      integer middle,rightedge,leftedge,botedge,topedge
      integer rightdedge,leftdedge,botdedge,topdedge                         !added d rrl
      integer peleft,peright,peabove,pebelow
      integer perightabove,perightbelow,peleftbelow,peleftabove
      integer iwallclock         !Koo: flag for the walltimers 07/16/07
      double precision wtime0,wtime1,wtime2,wtime3,wtime4, wtdelta
      double precision, allocatable:: cputime(:,:)  !cputimes 07/16/07

      save


      contains

!!!--------------------------- subroutine rinitmsg --------------------------------!!!
      subroutine rinitmsg
      npos = mod((mpi_rank+nprocx),nprocx)+1
      mpos = mpi_rank/nprocx+1
c
c  set up edge information (0=>not an edge, 1=>edge)
      rightedge=0
      if (mod((mpi_rank+1+nprocx),nprocx).eq.0) rightedge=1
c
      leftedge=0
      if (mod((mpi_rank+1+nprocx),nprocx).eq.1.or.nprocx.eq.1)
     +                                                    leftedge=1
c
      botedge=0
      if ((mpi_rank+1).le.nprocx) botedge=1
c
      topedge=0
      if (((nproc)-(mpi_rank+1)).lt.nprocx) topedge=1
c
      middle=0
      if (rightedge.eq.0.and.leftedge.eq.0.and.botedge.eq.0
     +                                    .and.topedge.eq.0) middle=1
c
c  set up neighbor information
c
      peleft=mpi_rank-1
      if (peleft.lt.((mpos-1)*nprocx)) peleft=mpi_rank+(nprocx-1)
c
      peright=mpi_rank+1
      if (peright.gt.(mpos*nprocx-1)) peright=mpi_rank-(nprocx-1)
c
      peabove=mpi_rank+nprocx
      if (peabove.gt.(nproc-1)) peabove=npos-1
c
      pebelow=mpi_rank-nprocx
      if (pebelow.lt.0) pebelow=mpi_rank+((nprocy-1)*nprocx)
c
      if (npos.lt.nprocx) then
         perightabove=peabove+1
         perightbelow=pebelow+1
      else
         perightabove=peabove-(nprocx-1)
         perightbelow=pebelow-(nprocx-1)
      end if
      if (npos.ne.1) then
         peleftabove=peabove-1
         peleftbelow=pebelow-1
      else
         peleftabove=peabove+(nprocx-1)
         peleftbelow=pebelow+(nprocx-1)
      end if
c
      if (mpi_rank.eq.0)
     +print *,' RINITMSG ',mpi_rank,peleft,peright,pebelow,peabove,
     +          perightabove,perightbelow,peleftabove,peleftbelow
c
      timestep=0.

       rightdedge=rightedge*(1-ibcx)                   !added d rrl
       leftdedge=leftedge*(1-ibcx)                  !added d rrl
       topdedge=topedge*(1-ibcy)                  !added d rrl
       botdedge=botedge*(1-ibcy)                  !added d rrl



      return
      end subroutine rinitmsg

!!!--------------------------- end subroutine rinitmsg --------------------------------!!!

!!!--------------------------- subroutine allgather3d --------------------------------!!!
      subroutine allgather3d(procdata, domdata,llim)

      Implicit None

      !real :: procdata(1-ih:np+ih,1-ih:mp+ih,l)
      real :: procdata(1:np,1:mp,llim)
      real :: domdata(1:n,1:m,llim)
      integer :: llim
      integer :: i,j,k,iproc,source,tag,stat(MPI_STATUS_SIZE)
      integer ::tmpnpos,tmpmpos,ia,ja
       
      !Setup the MPI derived datatypes for use in scatter gather of full domain arrays
      call MPI_TYPE_VECTOR(mp*llim, np, np, MPI_REAL, blocktype, ierror)
      call MPI_TYPE_COMMIT(blocktype, ierror)
     
      if(mpi_rank.ne.0)then
       tag = mpi_rank
       !call MPI_SEND(procdata(1,1,1), np*mp*l, mpi_real, 0, tag,
       call MPI_SEND(procdata(1,1,1), np*mp*llim, mpi_real, 0, tag,
     &                       MPI_COMM_WORLD, ierror)  !FP
       else
       do k=1,llim
        do j=1,mp
         do i=1,np
        !domdata(1:np,1:mp,1:l)=procdata(1:np,1:mp,1:l)
        domdata(i,j,k)=procdata(i,j,k)
         enddo
        enddo
       enddo
       endif
      if(mpi_rank.eq.0)then
       do iproc=1,numprocs-1
        tmpnpos = mod((iproc+nprocx),nprocx)+1
        tmpmpos = iproc/nprocx+1
        ia=(tmpnpos-1)*np+1
        ja=(tmpmpos-1)*mp+1
        source = iproc
       !call MPI_RECV(domdata(ia:ia+np-1,ja:ja+mp-1,1:l), 1, blocktype, source, iproc,
       call MPI_RECV(domdata(ia:ia+np-1,ja:ja+mp-1,1:llim), 1, blocktype, source, iproc,
     &                   MPI_COMM_WORLD, stat, ierror)  !FP
       enddo
      endif
       !call MPI_Bcast(domdata,n*m*l,mpi_real,0,mpi_comm_world,ierror)
       call MPI_Bcast(domdata,n*m*llim,mpi_real,0,mpi_comm_world,ierror)   !FP
       !call MPI_Bcast(domdata(1:n,1:m,1:l),n*m*l,mpi_real,0,mpi_comm_world,ierror)

      call MPI_TYPE_FREE(blocktype,ierror)
      return
      end subroutine allgather3d
!!!--------------------------- end subroutine allgather3d --------------------------------!!!
!!!--------------------------- subroutine allgather3d --------------------------------!!!
      subroutine allgather3d0(procdata, domdata,llim)
       !start from 0 on k
      Implicit None

      !real :: procdata(1-ih:np+ih,1-ih:mp+ih,l)
      real :: procdata(1:np,1:mp,0:llim)
      real :: domdata(1:n,1:m,0:llim)
      integer :: llim
      integer :: i,j,k,iproc,source,tag,stat(MPI_STATUS_SIZE)
      integer ::tmpnpos,tmpmpos,ia,ja
      
      !JAS at the moment I will define blocktype0 to have mp*l+1 element but need to revisit
      ! for the case when lphase is bigger than 1  
      call MPI_TYPE_VECTOR(mp*(llim+1), np, np, MPI_REAL, blocktype0, ierror)
      call MPI_TYPE_COMMIT(blocktype0, ierror)
      
      if(mpi_rank.ne.0)then
       tag = mpi_rank
       !call MPI_SEND(procdata(1,1,1), np*mp*l, mpi_real, 0, tag,
       call MPI_SEND(procdata(1,1,0), np*mp*(llim+1), mpi_real, 0, tag,
     &                       MPI_COMM_WORLD, ierror)  !FP
       else
       do k=0,llim
        do j=1,mp
         do i=1,np
        !domdata(1:np,1:mp,1:l)=procdata(1:np,1:mp,1:l)
        domdata(i,j,k)=procdata(i,j,k)
         enddo
        enddo
       enddo
       endif
      if(mpi_rank.eq.0)then
       do iproc=1,numprocs-1
        tmpnpos = mod((iproc+nprocx),nprocx)+1
        tmpmpos = iproc/nprocx+1
        ia=(tmpnpos-1)*np+1
        ja=(tmpmpos-1)*mp+1
        source = iproc
       !call MPI_RECV(domdata(ia:ia+np-1,ja:ja+mp-1,1:l), 1, blocktype, source, iproc,
       call MPI_RECV(domdata(ia:ia+np-1,ja:ja+mp-1,0:llim), 1, blocktype0, source, iproc,
     &                   MPI_COMM_WORLD, stat, ierror)  !FP
       enddo
      endif
       !call MPI_Bcast(domdata,n*m*l,mpi_real,0,mpi_comm_world,ierror)
       call MPI_Bcast(domdata,n*m*(llim+1),mpi_real,0,mpi_comm_world,ierror)   !FP
       !call MPI_Bcast(domdata(1:n,1:m,1:l),n*m*l,mpi_real,0,mpi_comm_world,ierror)
     
       call MPI_Barrier(mpi_comm_world,ierror)
 
      call MPI_TYPE_FREE(blocktype0,ierror)

      return
      end subroutine allgather3d0
!!!--------------------------- end subroutine allgather3d --------------------------------!!!
!!!--------------------------- subroutine finalizemsg --------------------------------!!!
      subroutine finalize_msg()

      Implicit None 

      call mpi_finalize(ierror) 
      
      end subroutine finalize_msg
!!!--------------------------- end subroutine finalizemsg --------------------------------!!!
      subroutine computeWallTime(wtini,i,name)
      double precision wtini
      integer i
      character(len=4)::name
      integer :: ierr 
        wtime3=MPI_Wtime();
        wtdelta=wtime3-wtini
        call mpi_reduce(wtdelta,cputime(i,1),1,mpi_double_precision,
     +                  mpi_sum,0,mpi_comm_world,ierr)
        if(mpi_rank.EQ.0) then
           write(6,*) "cputime for ",name," ",cputime(i,1)
           cputime(i,2)=cputime(i,1)+cputime(i,2)  !cumulative 
        endif 
        wtime2=wtime3
       end subroutine computeWallTime 
      end module msga
c2345678***************************************************

