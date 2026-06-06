c234567****************************************************
      module io

      implicit none

      contains
!----------------------------------------------------------
      subroutine frqwriteio(fname)

      use gridsetup
      use xvo
      use xve
      use xvi
      use workavg
      use turba
      use fireteca
      use weights
      use constants
      use metryic
      use radiation
      use msga
      use pres
      use lspgf
      use nonlocal
      Implicit None
      
      character(*) fname

      if (mpi_rank.eq.0) open(unit=21,file=fname,form='unformatted',
     +                         status='unknown')
      ! nb: the last index between 1 and 4  is for specific useful
      ! outputs in the log file (see writeio for details)
      call writeio(xvb(1-ih,1-ih,1,1),21,1-ih,np+ih,1-ih,mp+ih,l,2) ! field 1
      call writeio(xvb(1-ih,1-ih,1,2),21,1-ih,np+ih,1-ih,mp+ih,l,4) ! field 2
      call writeio(xvb(1-ih,1-ih,1,3),21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 3
      call writeio(xvb(1-ih,1-ih,1,4),21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 4

      if(iturb.ge.1)then
        call writeio(xvb(1-ih,1-ih,1,5),21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 5
        call writeio(xvb(1-ih,1-ih,1,6),21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 6
      endif

      call writeio(xvb(1-ih,1-ih,1,nv),21,1-ih,np+ih,1-ih,mp+ih,l,3) ! field nv
      if (ilspgf.ge.1) call writeio(sinthetaf,21,1-ih,np+ih,1-ih,mp+ih,l,0) 
      if(irod.eq.1) then
        call writeio(rhof,21,1-ih,np+ih,1-ih,mp+ih,l,0)               ! field 7
        call writeio(rhofinitial,21,1-ih,np+ih,1-ih,mp+ih,l,0)        ! field 8
        call writeio(rhowater,21,1-ih,np+ih,1-ih,mp+ih,l,0)        ! field 15
        call writeio(actualfueldepth,21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 17
        call writeio(sizescale,21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 17
        call writeio(temps,21,1-ih,np+ih,1-ih,mp+ih,l,1)              ! field 9
        call writeio(xvb(1-ih,1-ih,1,7),21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 10
c        call writeio(tempg,21,1-ih,np+ih,1-ih,mp+ih,l)              ! inserted as field 11 rrl

        if(inonlocal.eq.1) then
          call writeio(xvb(1-ih,1-ih,1,8),21,1-ih,np+ih,1-ih,mp+ih,l,0)
        endif

        if(irhovapor.eq.1) then
          call writeio(xvb(1-ih,1-ih,1,8),21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 11
        endif
        call writeio(sies,21,1-ih,np+ih,1-ih,mp+ih,l,0)            ! field 12
      !  call writeio(rhos,21,1-ih,np+ih,1-ih,mp+ih,l,0)            ! field 13
!        call rmaxmin1(rhos,'rhos in io',1-ih,np+ih,1-ih,mp+ih,l,0)
        if (idirt.eq.1) call writeio(rhodirt,21,1-ih,np+ih,1-ih,mp+ih,l,0)         ! field 14
        call writeio(psiwmax,21,1-ih,np+ih,1-ih,mp+ih,l,0)         ! field 16
        if (irad.GE.1.and.irod.eq.1) then  !KOO eq->GE
          call writeio(firad,21,1-ih,np+ih,1-ih,mp+ih,l,0)       ! field 23
          call writeio(frhosiesrad,21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 24
         endif

! below this point things will not be read during restart
        if (ioextra.eq.1) then
          if(irad.GE.1) then ! KOO eq->GE
         if(irad.eq.1)then
          call writeio(Ef,21,1-ih,np+ih,1-ih,mp+ih,LL,0)     ! field 18
          call writeio(Es,21,1-ih,np+ih,1-ih,mp+ih,LL,0)     ! field 19
          call writeio(rnetsol,21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 20
          call writeio(rnetgas,21,1-ih,np+ih,1-ih,mp+ih,l,0) ! field 21
         endif
        endif
      endif
        call writeio(convht,21,1-ih,np+ih,1-ih,mp+ih,l,0)  ! field 22

        call writeio(tempg,21,1-ih,np+ih,1-ih,mp+ih,l,0)  ! field 22
      if (iradeastflux.eq.1) 
     +      call writeio(eastFlux,21,1-ih,np+ih,1-ih,mp+ih,l,0)
      if (inonlocal.eq.1) then
        call writeio(fg,21,1-ih,np+ih,1-ih,mp+ih,l,0) !jjc
      endif
       end if
      if (mpi_rank.eq.0) close(21) ! JAS 10/4/06  JLW 10/12/06
 
      return 
      end subroutine frqwriteio
!----------------------------------------------------------
!    writeio on a subdomain of the whole domain
!----------------------------------------------------------
      subroutine frqwritesubio(fname)

      use gridsetup
      use xvo
      use xve
      use xvi
      use workavg
      use turba
      use fireteca
      use weights
      use constants
      use metryic
      use radiation
      use msga
      use pres
      use nonlocal
      Implicit None

      character(*) fname
      if (mpi_rank.eq.0) open(unit=21,file=fname,form='unformatted',
     +                         status='unknown')
      call writesubio(xvb(1-ih,1-ih,1,1),21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)
      call writesubio(xvb(1-ih,1-ih,1,2),21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)
      call writesubio(xvb(1-ih,1-ih,1,3),21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)
      if(iturb.ge.1) then
        call writesubio(xvb(1-ih,1-ih,1,5),21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)
        call writesubio(xvb(1-ih,1-ih,1,6),21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)
      endif
      call writesubio(xvb(1-ih,1-ih,1,7),21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)
      call writesubio(xvb(1-ih,1-ih,1,nv),21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)
      if(irod.eq.1) then
        call writesubio(temps,21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)!
      endif
       call writesubio(tempg,21,1-ih,np+ih,1-ih,mp+ih,l,issub,jssub,nisub,njsub,nksub)!
      if (mpi_rank.eq.0) close(21)

      return
      end subroutine frqwritesubio

!----------------------------------------------------------------------------
      subroutine irstreadio(fname)

      use gridsetup
      use xvo
      use xve
      use xvi
      use workavg
      use turba
      use fireteca
      use weights
      use constants
      use metryic
      use restarta
      use radiation
      use msga
      use lspgf
      use pres
      use nonlocal
      Implicit None


      !JAS 3/6/06 added explicit declarationis to comply with implicit none 
      !integer :: ir
      character(len=60)::fname  !JMC added filename that is passed in
      real, allocatable :: unused(:,:,:)

      if(mpi_rank.eq.0)
     + open(41,file=fname,form='unformatted',status='old')
      !do ir=1,nrst
        call    readio(xvb(1-ih,1-ih,1,1),41,l,1-ih,np+ih,1-ih,mp+ih)
        call    readio(xvb(1-ih,1-ih,1,2),41,l,1-ih,np+ih,1-ih,mp+ih)
        call    readio(xvb(1-ih,1-ih,1,3),41,l,1-ih,np+ih,1-ih,mp+ih)
        call    readio(xvb(1-ih,1-ih,1,4),41,l,1-ih,np+ih,1-ih,mp+ih)

        if(iturb.ge.1)then
          call    readio(xvb(1-ih,1-ih,1,5),41,l,1-ih,np+ih,1-ih,mp+ih)
          call    readio(xvb(1-ih,1-ih,1,6),41,l,1-ih,np+ih,1-ih,mp+ih)
        endif

        call    readio(xvb(1-ih,1-ih,1,nv),41,l,1-ih,np+ih,1-ih,mp+ih)
        if (ilspgf.ge.1) call readio(sinthetaf,41,l,1-ih,np+ih,1-ih,mp+ih) 

        if(irod.eq.1) then
         if (irst.eq.1) then
          call    readio(rhof,41,l,1-ih,np+ih,1-ih,mp+ih)
          call    readio(rhofinitial,41,l,1-ih,np+ih,1-ih,mp+ih)
          call    readio(rhowater,41,l,1-ih,np+ih,1-ih,mp+ih) !mv fp
          call    readio(actualfueldepth,41,l,1-ih,np+ih,1-ih,mp+ih) !mv fp
          call    readio(sizescale,41,l,1-ih,np+ih,1-ih,mp+ih) !mv fp
         else if (irst.eq.0.or.irst.eq.2) then !irst.eq.2 or 0 (for iwindfiedin)
           allocate(unused(1-ih:np+ih,1-ih:mp+ih,l))
           call    readio(unused,41,l,1-ih,np+ih,1-ih,mp+ih)
           call    readio(unused,41,l,1-ih,np+ih,1-ih,mp+ih)
           call    readio(unused,41,l,1-ih,np+ih,1-ih,mp+ih)
           call    readio(unused,41,l,1-ih,np+ih,1-ih,mp+ih)
           call    readio(unused,41,l,1-ih,np+ih,1-ih,mp+ih)
         end if !irst.eq.2.or.0
          call    readio(temps,41,l,1-ih,np+ih,1-ih,mp+ih)
          call    readio(xvb(1-ih,1-ih,1,7),41,l,1-ih,np+ih,1-ih,mp+ih)

          if(inonlocal.eq.1) then
            call    readio(xvb(1-ih,1-ih,1,8),41,l,1-ih,np+ih,1-ih,mp+ih)
          endif

          if(irhovapor.eq.1) then
            call    readio(xvb(1-ih,1-ih,1,8),41,l,1-ih,np+ih,1-ih,mp+ih)
          endif

          call    readio(sies,41,l,1-ih,np+ih,1-ih,mp+ih)
          !call    readio(rhos,41,l,1-ih,np+ih,1-ih,mp+ih)    !rm fp
          if (idirt.eq.1) call    readio(rhodirt,41,l,1-ih,np+ih,1-ih,mp+ih)  !rm fp
          call    readio(psiwmax,41,l,1-ih,np+ih,1-ih,mp+ih)

c          if(irad.GE.1) then   !KOO eq->GE
c           if(irad.eq.1)then
c            call readio(Ef,41,LL,1-ih,np+ih,1-ih,mp+ih)
c            call readio(Es,41,LL,1-ih,np+ih,1-ih,mp+ih)
c           endif
c            call readio(rnetsol,41,l,1-ih,np+ih,1-ih,mp+ih)
c            call readio(rnetgas,41,l,1-ih,np+ih,1-ih,mp+ih)
c          endif
c            call readio(convht,41,l,1-ih,np+ih,1-ih,mp+ih)
        if (irad.GE.1) then  !KOO eq->GE
          call readio(firad,41,l,1-ih,np+ih,1-ih,mp+ih) !field 23
          call readio(frhosiesrad,41,l,1-ih,np+ih,1-ih,mp+ih) !field 24
        endif
       endif
            !call readio(tempg,41,l,1-ih,np+ih,1-ih,mp+ih)  ! moved by FP

      !enddo

      if (mpi_rank.eq.0) close(41) ! JAS 10/4/06  JLW 10/12/06
      

      return
      end subroutine irstreadio

!----------------------------------------------------------
      subroutine irstreadio2(fname)

      use gridsetup
      use xvo
      use xve
      use xvi
      use workavg
      use turba
      use fireteca
      use weights
      use constants
      use metryic
      use restarta
      use radiation
      use msga
      use pres
      use nonlocal
      use lspgf
      Implicit None


      !JAS 3/6/06 added explicit declarationis to comply with implicit none
      !integer :: ir
      character(len=60)::fname  !JMC added filename that is passed in

      open(41,file=fname,form='unformatted',status='old')
      !do ir=1,nrst
        call    readio(xvb(1-ih,1-ih,1,1),41,l,1-ih,np+ih,1-ih,mp+ih)
        call    readio(xvb(1-ih,1-ih,1,2),41,l,1-ih,np+ih,1-ih,mp+ih)
        call    readio(xvb(1-ih,1-ih,1,3),41,l,1-ih,np+ih,1-ih,mp+ih)
        call    readio(xvb(1-ih,1-ih,1,4),41,l,1-ih,np+ih,1-ih,mp+ih)

        if(iturb.ge.1)then
          call    readio(xvb(1-ih,1-ih,1,5),41,l,1-ih,np+ih,1-ih,mp+ih)
          call    readio(xvb(1-ih,1-ih,1,6),41,l,1-ih,np+ih,1-ih,mp+ih)
        endif

        call    readio(xvb(1-ih,1-ih,1,nv),41,l,1-ih,np+ih,1-ih,mp+ih)
        if (ilspgf.ge.1) call readio(sinthetaf,41,l,1-ih,np+ih,1-ih,mp+ih) 
      !enddo
      close(41)

      return

      end subroutine irstreadio2
!--------------------------------------------------------------------------------------------
      subroutine writeio(data,iunit,il,iu,jl,ju,nzdim,iflag)
      use gridsetup
      use msga
      use metryic
      Implicit None
      ! FP (10/2012) added an iflag that entails to do specific computations on fields that can be printed in the
      ! log such as fireposition, mean wind at a a given height...
      ! iflag = 1 deals with temps for fire position
      ! iflag = 2, 3, 4 deals with mean u, rhog, v at the reference
      ! height zu for check of wind field evolution          
      integer,intent(in) :: il,iu,jl,ju,nzdim,iunit,iflag
      !FP added iflag in sept 2012 to track the firefront, mean wind at a given height...
      
c
      real data(il:iu,jl:ju,nzdim)
      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
      real :: chtemp
      integer:: front=1
      real::z1,z2,var1,var2,varmean
      real,external::zcart
c
      real,allocatable::tmparray(:,:,:,:),outdata(:,:,:)
c
      if (mpi_rank.eq.0)
     +  allocate(tmparray(il:iu,jl:ju,nzdim,nproc),
     +            outdata(n,m,nzdim))
        nsize=(iu-il+1)*(ju-jl+1)*nzdim
        call mpi_gather(data,nsize,mpi_real,tmparray,nsize,mpi_real,
     +     0,mpi_comm_world,ierror)
      if (mpi_rank.eq.0) then
        do iprocx=1,nprocx
          do jprocy=1,nprocy
           iproc=1+(iprocx-1)+(jprocy-1)*nprocx
           do k=1,nzdim
             do j=1,mp
               do i=1,np
                 ia=(iprocx-1)*np + i
                 ja=(jprocy-1)*mp + j
                 outdata(ia,ja,k)=tmparray(i,j,k,iproc)
               enddo
             enddo
           enddo
         enddo
        enddo

        chtemp=sum(outdata)
        write (iunit)outdata
c
       ! if variable is temps we track the firefront
      ! iflag = 1 deals with temps for fire position
       if(iflag.eq.1) then
         do ia=2,n-1
         do ja=2,m-1
            if (outdata(ia,ja,1).gt.600) 
     .         front=real(i)
         enddo
         enddo
         if (front.gt.1) write(6,*) 'firefront has reached cell :',front
       else if (iflag.le.4.and.iflag.ne.0) then
      ! iflag = 2, 3, 4 deals with mean variable (u, rhog, v) at the reference
      ! height zu for check of wind field evolution. NB: this does not
      ! account for topo          
       var1=0.0
       var2=0.0
       varmean=0.0
       do k=1,l-1
        z1=zcart(z(k),1,1)-zs(1,1)
        z2=zcart(z(k+1),1,1)-zs(1,1)
        !write (6,*) 'k,z1,z2,zu',k,z1,z2,zu
        if (z1<=zu.and.zu<z2) then
        do ja=1,m
         do ia=1,n
           var1=var1+outdata(ia,ja,k)
           var2=var2+outdata(ia,ja,k+1)
         enddo
        enddo
         varmean=(var1*(z2-zu)+var2*(zu-z1))
     +    /(n*m*(z2-z1))
        endif
       enddo
       if (iflag.eq.2) write(6,*) 'rho*u at height ',zu, ' is',varmean
       if (iflag.eq.4) write(6,*) 'rho*v at height ',zu, ' is ',varmean
       if (iflag.eq.3) write(6,*) 'target rhou is',u0*varmean
       if (iflag.eq.3) write(6,*) 'target rhov is',v0*varmean
       endif ! end if iflag

        deallocate(tmparray,outdata)
      endif
      return
      end subroutine writeio
!----------------------------------------------------------
      subroutine writesubio(data,iunit,il,iu,jl,ju,nzdim,iis,jjs,ni,nj,nk)
      use gridsetup
      use msga
      use metryic
      Implicit None
      integer,intent(in) :: il,iu,jl,ju,nzdim,iunit,iis,jjs,ni,nj, nk

c
      real data(il:iu,jl:ju,nzdim)
      !JAS 3/6/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,ia2,ja2,nsize
      real :: chtemp
      real,external::zcart
      real,allocatable::tmparray(:,:,:,:),outdata(:,:,:)
      if (mpi_rank.eq.0)
     +  allocate(tmparray(il:iu,jl:ju,nzdim,nproc),
     +            outdata(ni,nj,nk))
        nsize=(iu-il+1)*(ju-jl+1)*nzdim
        call mpi_gather(data,nsize,mpi_real,tmparray,nsize,mpi_real,
     +     0,mpi_comm_world,ierror)
      if (mpi_rank.eq.0) then
        do iprocx=1,nprocx
          do jprocy=1,nprocy
           iproc=1+(iprocx-1)+(jprocy-1)*nprocx
           do k=1,nk
             do j=1,mp
               do i=1,np
                 ia=(iprocx-1)*np + i
                 ja=(jprocy-1)*mp + j
                 ia2 = ia - iis + 1
                 ja2 = ja - jjs + 1
                 if ((ia2 >= 1).and.(ja2 >= 1).and.(ia2 <= ni).and.(ja2<= nj))
     +              outdata(ia2,ja2,k)=tmparray(i,j,k,iproc)
               enddo
             enddo
           enddo
         enddo
        enddo

        chtemp=sum(outdata)
        write (iunit)outdata

        deallocate(tmparray,outdata)
      endif
      return

      end subroutine writesubio
!----------------------------------------------------------

      subroutine readio(array,iunit,nzdim,il,iu,jl,ju)
      use gridsetup
      use msga

      Implicit None


      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,nzdim,iunit

c
      real array(il:iu,jl:ju,nzdim)
      real,allocatable:: tmparray(:,:,:,:),indata(:,:,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,iproc,iprocx,jprocy,ia,ja,nsize
      real :: chtmp
c
      if (mpi_rank.eq.0) then
        allocate(indata(n,m,nzdim),
     +           tmparray(il:iu,jl:ju,nzdim,nproc))
        read (iunit)indata
        chtmp=sum(indata)
        do iprocx=1,nprocx
          do jprocy=1,nprocy
           iproc=1+(iprocx-1)+(jprocy-1)*nprocx
           do k=1,nzdim
             do j=1,mp
               do i=1,np
                 ia=(iprocx-1)*np + i
                 ja=(jprocy-1)*mp + j
                 tmparray(i,j,k,iproc)=indata(ia,ja,k)
               enddo
             enddo
           enddo
         enddo
        enddo
      endif

      nsize=(iu-il+1)*(ju-jl+1)*nzdim
      call mpi_scatter(tmparray,nsize,mpi_real,array,nsize,mpi_real,
     +   0,mpi_comm_world,ierror)
      if(mpi_rank.eq.0)deallocate(tmparray,indata)

      return
      end subroutine readio
!----------------------------------------------------------

      end module io
