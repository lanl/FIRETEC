c***********************************************************
      subroutine rmaxmin1(array,string,il,iu,jl,ju,lls)
      use gridsetup
      use msga

      implicit none

      integer :: il,iu,jl,ju,lls
      real array(il:iu,jl:ju,lls)
      real   tmp(il:iu,jl:ju,lls)
      character*(*) string
      real,allocatable::tmparray(:,:,:,:)
      integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,imaxtemp
     &          ,jmaxtemp,kmaxtemp,imaxlocaltemp,jmaxlocaltemp
     &          ,iprocmax,imintemp,jmintemp,kmintemp,iminlocaltemp
     &          ,jminlocaltemp,iprocmin
      integer :: ierr,nsize
      real    :: rmin,rmax,xmin,xmax
      integer :: loc_prt=1

      if (mpi_rank.eq.0)
     +  allocate(tmparray(il:iu,jl:ju,lls,nproc))
c
c
c  compute max and mim of data that is local to the process
      tmp=array
      rmax=array(1,1,1)
      rmin=array(1,1,1)
      do k=1,l
        do j=1,mp
          do i=1,np
            rmax=amax1(array(i,j,k),rmax)
            rmin=amin1(array(i,j,k),rmin)
          enddo
        enddo
      enddo
c
c  compute global max and min and put it on the master
      call MPI_Reduce(rmax,xmax,1,mpi_real,mpi_max,0,
     +                   mpi_comm_world,ierr)
      call MPI_Reduce(rmin,xmin,1,mpi_real,mpi_min,0,
     +                   mpi_comm_world,ierr)
c
c  print max and min values on the master
c     if (mpi_rank.eq.0) print *,' RMAXMIN for ',string,
c    +                           ' - max = ',xmax,
c    +                           ' ,min = ',xmin
c
      if(loc_prt.EQ.1) then 
        nsize=(iu-il+1)*(ju-jl+1)*lls
        call mpi_gather(tmp,nsize,mpi_real,tmparray,nsize,mpi_real,
     +     0,mpi_comm_world,ierror)
         if(mpi_rank.eq.0) then
        do iprocx=1,nprocx
          do jprocy=1,nprocy
          iproc=1+(iprocx-1)+(jprocy-1)*nprocx
         do i=1,np
            do j=1,mp
               do k=1,lls
c                 zl=zcart(z(k),i,j,zb)
c                 zla=zcart(z(k),i,j,zb)-zs(i,j)
                  if (xmax.eq.tmparray(i,j,k,iproc)) then
                     ia=(iprocx-1)*np+i
                     ja=(jprocy-1)*mp+j
                     imaxtemp=ia
                     jmaxtemp=ja
                     kmaxtemp=k
                     imaxlocaltemp=i
                     jmaxlocaltemp=j
                     iprocmax=iproc
                  endif
                  if (xmin.eq.tmparray(i,j,k,iproc)) then
                     ia=(iprocx-1)*np+i
                     ja=(jprocy-1)*mp+j
                     imintemp=ia
                     jmintemp=ja
                     kmintemp=k
                     iminlocaltemp=i
                     jminlocaltemp=j
                     iprocmin=iproc
                  endif
               enddo
            enddo
         enddo
         enddo
         enddo
         endif
        endif 

         if(mpi_rank.eq.0)then
!                 write(6,*)'-RMAXMIN for ',string,
!     +                              ' - max = ',xmax,' ,min = ',xmin
          write(6,'(a,1x,a," - max = ",es13.4," , min = ",es13.4)')
     +      '  RMAXMIN for',trim(string),xmax,xmin
         if(loc_prt.eq.1) then         
                 write(6,*)' :cell of maximum above is '
     &                              ,imaxtemp,jmaxtemp,kmaxtemp
                 write(6,*)' :cell of minimum above is '
     &                              ,imintemp,jmintemp,kmintemp
         endif 

                 deallocate(tmparray)  ! bug fix Koo 08/25/2011
         endif

      return
      end
