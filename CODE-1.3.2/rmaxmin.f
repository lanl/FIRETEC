c***********************************************************
      !subroutine rmaxmin(array,string,il,iu,jl,ju,lls,nvp)
      !JAS 3/6/06 changed routine signature her and in calling routines
      !(compress,firetec,rinit) because "string" is not used!
      subroutine rmaxmin(array,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use constants
      use msga
      use nonlocal

      Implicit None
 

      !JAS 3/6/06 added e3xplicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp
      
c
      real::array(il:iu,jl:ju,lls,nvp)
      real::tmp(il:iu,jl:ju,lls),pressure(il:iu,jl:ju,lls)
      real,allocatable::tmparray(:,:,:,:)
      !  JAS not used     --character*(*) string
      character*(5) names(10) /'u','v','w','temp','ka','kb',
     +                         'rhoo','dens',' ',' '/
      character*(5) namesnonlocal(10) /'u','v','w','temp','ka','kb',
     +                         'rhoo','rhoh','dens',' '/
      !JAS 3/6/06 added e3xplicit declarations to comply with implicit none
! KOO
!      character(len=13) :: srcfile = "  RMAXMIN for"

      integer :: kv
      real :: rmax,rmin,xmax,xmin,rmaxtheta,rmintheta,xmaxth,xminth

! KOO
      integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,imaxtemp
     &          ,jmaxtemp,kmaxtemp,imaxlocaltemp,jmaxlocaltemp
     &          ,iprocmax,imintemp,jmintemp,kmintemp,iminlocaltemp
     &          ,jminlocaltemp,iprocmin
      integer :: nsize ! ierr
      integer :: loc_prt=0  

      if (inonlocal==1) names=namesnonlocal
      
      if (mpi_rank.eq.0)
     +  allocate(tmparray(il:iu,jl:ju,lls,nproc))
c
c  compute max and mim of data that is local to the process
      do kv=1,nvp
         tmp(1:np,1:mp,1:lls)=array(1:np,1:mp,1:lls,kv)
         if(kv.ne.4.and.kv.ne.nv)then
         !if(kv.ne.4.and.kv.ne.nv.and.kv.ne.7)then    !JASTEST
           tmp(1:np,1:mp,1:lls)=tmp(1:np,1:mp,1:lls)/
     +       array(1:np,1:mp,1:lls,nv)
         endif
         if(kv.eq.4)then
           pressure(1:np,1:mp,1:lls)=                     ! tmp is rho theta
     +       (tmp(1:np,1:mp,1:lls)*rg/prrcp)**(cp/cv)
           tmp(1:np,1:mp,1:lls)=tmp(1:np,1:mp,1:lls)/     ! tmp is theta
     +       array(1:np,1:mp,1:lls,nv)
           rmaxtheta   =maxval(tmp(1:np,1:mp,1:lls))
           rmintheta   =minval(tmp(1:np,1:mp,1:lls))

           call mpi_reduce(rmaxtheta,xmaxth,1,mpi_real,mpi_max,0,
     +                     mpi_comm_world,ierror)
           call mpi_reduce(rmintheta,xminth,1,mpi_real,mpi_min,0,
     +                     mpi_comm_world,ierror)
           if(mpi_rank.eq.0)write(6,*)' RMAXMIN for theta',
     +                      ' - max = ',xmaxth,' ,min = ',xminth
!           if(mpi_rank.eq.0)
!     +      write(6,'(a,1x,a," - max = ",es13.4," , min = ",es13.4)')
!     +      srcfile,'theta',xmaxth,xminth 

           tmp(1:np,1:mp,1:lls)=tmp(1:np,1:mp,1:lls)*     ! tmp is T (K)
     +                  (pressure(1:np,1:mp,1:lls)/1.e5)**cap
         endif

         rmax   =maxval(tmp(1:np,1:mp,1:lls))
         rmin   =minval(tmp(1:np,1:mp,1:lls))

         call mpi_reduce(rmax,xmax,1,mpi_real,mpi_max,0,
     +                   mpi_comm_world,ierror)
         call mpi_reduce(rmin,xmin,1,mpi_real,mpi_min,0,
     +                   mpi_comm_world,ierror)
         if(mpi_rank.eq.0) write(6,*)' RMAXMIN for ',names(kv),
     +                     ' - max = ',xmax,' ,min = ',xmin
!           if(mpi_rank.eq.0)
!     +      write(6,'(a,1x,a," - max = ",es13.4," , min = ",es13.4)')
!     +      srcfile,names(kv),xmax,xmin

! KOO for finding location 
        if(kv.LT.4 .and. loc_prt.EQ.1) then
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
                 write(6,*)' :cell of maximum above is '
     &                              ,imaxtemp,jmaxtemp,kmaxtemp
                 write(6,*)' :cell of minimum above is '
     &                              ,imintemp,jmintemp,kmintemp

        endif 
       endif
! KOO for finding location

      enddo

      if(mpi_rank.eq.0) deallocate(tmparray)
      end


