c***********************************************************
      subroutine topo
      use gridsetup
      use metryic
      use restarta
      use msga

      Implicit None

      
      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer :: i,j,ia,ja

c
c  read topo input data into the master or set to flat (ASCII for
c  benchmarks)
      if (mpi_rank.eq.0) then
        if(topofile(1:1).ne.'*'.and.topofile.ne.' ') then
c         open (25,file=topofile,form='formatted',status='old')
          open (25,file=topofile,form='unformatted',status='old')
c         do j=1,m
c           read (25,*) (zsio(i,j),i=1,n)
c         enddo
          read (25) zsio
          close (25)
          print *,' Read topography file ',topofile
        else
        zsio=0.0
        endif
      endif
c
c  broadcast the topo data to each node
      call MPI_Bcast(zsio,n*m,mpi_real,0,mpi_comm_world,ierror)
c
c  loop to extract topo values for each node
      do j=1,mp
        do i=1,np
          ia = (npos-1)*np + i
          ja = (mpos-1)*mp + j
          zs(i,j)=zsio(ia,ja)
        enddo
      enddo
c
c  initialize topo values for ghost points in case of cyclic
      call updated(zs,zs,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0)
c
c  initialize topo around edges depending on whether the data is
c  cyclic or not

c topo needs some careful consideration with respect to dedges!!!!!Judy please help here!!!!
c  FP and rrl (06/2012) decided that in case of non cyclic (ibcx=0) and edge: ibclefdege.eq.1:
c - 0 gradient should be implemented between cell 1 and 2
c halo cells should be updated with the 0 gradient

      !if(leftedge.eq.1) then
        !do j=1,mp
        !  zs(1,j)=zs(2,j)*(1-ibcx)+zs(0,j)*ibcx
        !enddo
      if(leftdedge.eq.1) then
        do j=1-ih,mp+ih
         zs(1,j)=zs(2,j) !zero gradient
         do i=1,ih
          zs(1-i,j)=zs(1,j)
         enddo
        enddo
      endif
      !if(rightedge.eq.1) then
        !do j=1,mp
        !  zs(np,j)=zs(np-1,j)*(1-ibcx)+zs(np+1,j)*ibcx
        !enddo
      if(rightdedge.eq.1) then
        do j=1-ih,mp+ih
         zs(np,j)=zs(np-1,j) !zero gradient
         do i=1,ih
          zs(np+i,j)=zs(np,j)
         enddo
        enddo
      endif
c
      if(j3.ne.0) then
        !if(botedge.eq.1) then
          !do i=1,np
          !  zs(i,1)=zs(i,2)*(1-ibcy)+zs(i,0)*ibcy
          !enddo
        if(botdedge.eq.1) then
          do i=1-ih,np+ih
          zs(i,1)=zs(i,2)
         do j=1,ih
            zs(i,1-j)=zs(i,1)
          enddo
         enddo
        endif
        !if(topedge.eq.1) then
          !do i=1,np
          !  zs(i,mp)=zs(i,mp-1)*(1-ibcy)+zs(i,mp+1)*ibcy
          !enddo
        if(topdedge.eq.1) then
          do i=1-ih,np+ih
          zs(i,mp)=zs(i,mp-1)
         do j=1,ih
            zs(i,mp+j)=zs(i,mp)
          enddo
         enddo
        endif
      endif
c
      !call updated(zs,zs,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1)
      !call updated(zs,zs,np,mp,1,1-ih,np+ih,1-ih,mp+ih,0)
      if (mpi_rank.eq.0) write(6,*) 'topo.f returning!'
      return
      end
