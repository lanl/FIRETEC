! 03/19/15 make option of clear fb_imm after imfb file written, and no
!          read from restart.
!
! 07/24/14 burnout added in fbrandlist as burnout limit length scale 
!
! 07/14/14 nfb_imm is changed to be real to have weighted numbers of
!          effective firebreands. 
!
! 03/01/13 fb_igniite will be added based on effective landed firebrand
!          numbers stored in nfb_imm(np,mp,l)  

!  01/13/11 ishape=8 is added for massless particles
!
!
!234567890123456789012345678901234567890123456789012345678901234567890**
      subroutine rinitfbrand(irst1,it1,nbmax1,nb_imm1)
!
      use gridsetup
      use metryic
      use msga
      use fireteca
      use firebranda

      Implicit None

      !include 'mpif.h'
!
      character base*5 /'fbvo.'/, fname1*7
      character basei*5 /'imfb.'/,fnamei1*7
      character fname2*8,fnamei2*8
      character fname3*9,fnamei3*9
      character fname4*10,fnamei4*10
      character fname5*11,fnamei5*11
      character fname6*12,fnamei6*12

      real,allocatable::g_brand(:,:),g_fb_imm(:,:)
      integer,allocatable::g_ibrand(:,:),g_ifb_imm(:,:)
      logical,allocatable::iproc_fb(:),iproc_fb_imm(:)
      integer i,j,k,nb,ia,ib,ja,jb,kb,ik,ierr
!      integer iscan_start,iscan_end,jscan_start,jscan_end
      integer nbmax_tl,nb_imm_tl,fdim,ifdim
      integer it1,irst1,nbmax1,nb_imm1,ibmax1,ib1,ib_imm
      real bx,by,bz 
      real,external :: zcart   


      namelist/firebrandlist/
     .w_ter,Cd_dn,Cd_cn,Cd_sp,tck,radi,aka,rho_fb,temphot,fden_limit,
     .tck_limit_d,rad_limit_c,rad_limit_sp,ivar_fbsize,nb_per_cell,
     .fb_start,lau_frq,lau_low_limit,lau_up_limit,iunit1,iunit2,tx_out,
     .ishape,lrrat_d,lrrat_c,v_in_rat,nfb_io,iland_fuel,burnout

      call MPI_Comm_rank(mpi_comm_world,mpi_rank,ierr)
!
! Read fbrandlist file only by the rank=0 process, then MPI_BCAST gridlist variables
!       to the other processes. This allows the code to run on clusters where 
!       processes other than the rank=0 process may not have direct access to
!       gridlist, for example, if they don't all share the same PWD environment
!       variable.
!
      burnout=0.00005

      if (mpi_rank.EQ.0) then
        open(unit=701,file='fbrandlist',form='formatted',status='old')
        read (701,nml=firebrandlist)
        close(701)
      endif

      allocate (nbrand(n,m,l))
      allocate (nfb_imm(0:n+1,0:m+1,0:l+1))
      allocate (zposition(0:np+1,0:mp+1,0:l+1))
      allocate (cl_hgt(0:np+1,0:mp+1,0:l))
      allocate (dzk(0:np+1,0:mp+1,0:l))
      allocate (nbproc(0:nproc-1))

      if (mpi_rank.EQ.0) then 
          print*,"rinitfbrand"
             if (w_ter) print*, " w_ter=ture"
             if (.NOT.w_ter) print*, " w_ter=false","v_in_rat=",v_in_rat
          print*," Cd_dn=",Cd_dn
          print*," Cd_sp=", Cd_sp
          print*," tck=",tck
          print*," radi=",radi
          print*," aka=",aka
          print*," rho_fb=",rho_fb
          print*," temphot=",temphot
          print*," fden_limit=",fden_limit
          print*," tck_limit_d=",tck_limit_d
          print*," rad_limit_c=",rad_limit_c
          print*," rad_limit_sp=",rad_limit_sp
          print*," ivar_fbsize=", ivar_fbsize     
          print*," nb_per_cell=",nb_per_cell
          print*," fb_start=",fb_start
          print*," lau_frq=", lau_frq
          print*," lau_low_limit=",lau_low_limit
          print*," lau_up_limit=",lau_up_limit
          print*," iunit1=",iunit1
          print*," iunit2=",iunit2
            if(tx_out) print*," tx_out=true"
            if(.NOT.tx_out) print*," tx_put=false"
            if(nfb_io) print*," nfb_io=true"
            if(.NOT.nfb_io) print*," nfb_io=false"
          print*," ishape=",ishape
c  temporarily one shape allowed
c   ishape:  0=dsk_ub, 1=dsk_dh/dt, 2=dsk_dr/dt
c            3=cyl_ub, 4=cyl_dh/dt, 5=cyl_dr/dt
c            6=sph_ub, 7=sph_dr/dt
         print*," lrrat_d=",lrrat_d
         print*," lrrat_c=",lrrat_c
         print*," v_in_rat=",v_in_rat
      endif 
      call MPI_Bcast(w_ter,1,mpi_logical,0,mpi_comm_world,ierr)
      call MPI_Bcast(Cd_dn,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(Cd_cn,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(Cd_sp,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(tck,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(radi,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(aka,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(rho_fb,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(temphot,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(fden_limit,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(tck_limit_d,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(rad_limit_c,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(rad_limit_sp,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(ivar_fbsize,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(nb_per_cell,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(fb_start,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(lau_frq,1,mpi_integer,0,mpi_comm_world,ierr)     

      call MPI_Bcast(lau_low_limit,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(lau_up_limit,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iunit1,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(iunit2,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(tx_out,1,mpi_logical,0,mpi_comm_world,ierr)
      call MPI_Bcast(nfb_io,1,mpi_logical,0,mpi_comm_world,ierr)
      call MPI_Bcast(ishape,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(lrrat_d,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(lrrat_c,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(v_in_rat,1,mpi_real,0,mpi_comm_world,ierr)
      call MPI_Bcast(iland_fuel,1,mpi_integer,0,mpi_comm_world,ierr)
      call MPI_Bcast(burnout,1,mpi_real,0,mpi_comm_world,ierr)

cell hgt calculations
c
cc for topo stuff
c      if(topofile(1:1).NE.'*' .AND. topoile.NE.' ') then
c      print*,"topo init for brands"
      do j=0,mp+1
         do i=0,np+1
             zposition(i,j,0)=zs(i,j)
             cl_hgt(i,j,0)=zposition(i,j,0)
             zposition(i,j,1)=zcart(z(1),i,j) 
             dzk(i,j,0)=zposition(i,j,1)-zs(i,j)
             cl_hgt(i,j,1)=2*dzk(i,j,0)+zs(i,j) 
             do k=2,l
                zposition(i,j,k)=zcart(z(k),i,j)
                dzk(i,j,k-1)=zposition(i,j,k)-zposition(i,j,k-1)
                cl_hgt(i,j,k)=2*zposition(i,j,k)-cl_hgt(i,j,k-1)
             enddo
             dzk(i,j,l)=cl_hgt(i,j,l)-zposition(i,j,l)
             zposition(i,j,l+1)=cl_hgt(i,j,l)
             
         enddo
       enddo
c      print*,"topo init done for brands"

      nbrand=0
      nfb_imm=0
      nbproc(0:nproc-1)=0
       
       nbmax1=0
       nb_imm1=0
       
       if (irst1.GE.1 . AND. it1.GT.fb_start) then
c234567

! reading imfb file 
       if (mpi_rank.eq.0) then
            if (it1.LT.100) then
        write (fnamei1,'(a,i2.2)') basei,it1
        open (iunit2,file=fnamei1,form='unformatted',status='unknown')

           elseif (it1.LT.1000) then
        write (fnamei2,'(a,i3.3)') basei,it1
        open (iunit2,file=fnamei2,form='unformatted',status='unknown')

            elseif (it1.LT.10000) then
        write (fnamei3,'(a,i4.4)') basei,it1
        open (iunit2,file=fnamei3,form='unformatted',status='unknown')

            elseif (it1.LT.100000) then
        write (fnamei4,'(a,i5.5)') basei,it1
        open (iunit2,file=fnamei4,form='unformatted',status='unknown')

            elseif (it1.LT.1000000) then
        write (fnamei5,'(a,i6.6)') basei,it1
        open (iunit2,file=fnamei5,form='unformatted',status='unknown')

             elseif (it1.LT.10000000) then
        write (fnamei6,'(a,i7.7)') basei,it1
        open (iunit2,file=fnamei6,form='unformatted',status='unknown')

            endif   ! it1
        endif ! mpi_rank  

        if(mpi_rank.EQ.0) then
! reading fbvo file - flying firebrands 
            if (it1.LT.100) then
        write (fname1,'(a,i2.2)') base,it1
        open (iunit1,file=fname1,form='unformatted',status='unknown')

           elseif (it1.LT.1000) then
        write (fname2,'(a,i3.3)') base,it1
        open (iunit1,file=fname2,form='unformatted',status='unknown')

            elseif (it1.LT.10000) then
        write (fname3,'(a,i4.4)') base,it1
        open (iunit1,file=fname3,form='unformatted',status='unknown')

            elseif (it1.LT.100000) then
        write (fname4,'(a,i5.5)') base,it1
        open (iunit1,file=fname4,form='unformatted',status='unknown')

            elseif (it1.LT.1000000) then
        write (fname5,'(a,i6.6)') base,it1
        open (iunit1,file=fname5,form='unformatted',status='unknown')

             elseif (it1.LT.10000000) then
        write (fname6,'(a,i7.7)') base,it1
        open (iunit1,file=fname6,form='unformatted',status='unknown')

             endif   ! it1

          read (iunit1) nbmax_tl,fdim,ifdim
          allocate(g_brand(nbmax_tl+1,fdim))
          allocate(g_ibrand(nbmax_tl+1,ifdim))
          read (iunit1) g_brand
          read (iunit1) g_ibrand
          print*, "reading brand: it=",it1,"nbmax=",nbmax_tl

          if(.NOT. nfb_io) then
          read (iunit2) nb_imm_tl,fdim,ifdim
          allocate(g_fb_imm(nb_imm_tl+1,fdim))
          allocate(g_ifb_imm(nb_imm_tl+1,ifdim))
          read (iunit2) g_fb_imm
          read (iunit2) g_ifb_imm
          print*, "reading im_brand: it=",it1,"nb_imm=",nb_imm_tl
          endif 
cc print brands
c       do nb=1,nbmax_tl
c       print*,"read_fb",nb,g_brand(nb,1),g_brand(nb,2),g_brand(nb,3),
c     +                  g_ibrand(nb,1)
c       enddo
c       do nb=1,nb_imm_tl
c       print*,"read_im",nb,g_fb_imm(nb,1),g_fb_imm(nb,2),g_fb_imm(nb,3),
c     +                  g_ifb_imm(nb,1)
c       enddo
       endif 

      call MPI_Comm_rank(mpi_comm_world,mpi_rank,ierr)

      call mpi_bcast(nbmax_tl,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(fdim,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(ifdim,1,mpi_integer,0,mpi_comm_world,ierror)
      if(mpi_rank.NE.0) then
        allocate(g_brand(nbmax_tl+1,fdim))
        allocate(g_ibrand(nbmax_tl+1,ifdim))
      endif
c234567
       call mpi_bcast(g_brand,(nbmax_tl+1)*fdim,mpi_real,
     +               0,mpi_comm_world,ierror)
       call mpi_bcast(g_ibrand,(nbmax_tl+1)*ifdim,mpi_integer,
     +               0,mpi_comm_world,ierror)

       allocate(iproc_fb(nbmax_tl))

      if(.NOT. nfb_io) then
       call mpi_bcast(nb_imm_tl,1,mpi_integer,0,mpi_comm_world,ierror)
      if(mpi_rank.NE.0) then
         allocate(g_fb_imm(nb_imm_tl+1,fdim))
         allocate(g_ifb_imm(nb_imm_tl+1,ifdim))
       endif ! mpi_rank.NE.0

       call mpi_bcast(g_fb_imm,(nb_imm_tl+1)*fdim,mpi_real,
     +               0,mpi_comm_world,ierror)
       call mpi_bcast(g_ifb_imm,(nb_imm_tl+1)*ifdim,mpi_integer,
     +               0,mpi_comm_world,ierror)

       allocate(iproc_fb_imm(nb_imm_tl)) 
  
      endif 

       do nb=1,nbmax_tl
        bx=g_brand(nb,1)
        by=g_brand(nb,2)
        bz=g_brand(nb,3)

        ia=nint(bx/dx+0.5*n+0.5)
        ja=nint(by/dy+0.5*m+0.5)

        ib=ia-((npos-1)*np)
        jb=ja-((mpos-1)*mp)
        kb=l
       if ( ((ib.GE.1) .AND. (ib.LE.np))
     +    .AND. ((jb.GE.1) .AND. (jb.LE.mp)) ) then

          ik=1
          do while(bz.GE.cl_hgt(ib,jb,ik))
             ik=ik+1
          enddo 
          kb=ik

          nbmax1=nbmax1+1
          iproc_fb(nb)=.true.
!       print*,"fb ID#",g_ibrand(nb,1),"is in",mpi_rank
       else
          iproc_fb(nb)=.false.
       endif
       enddo
 
!      print*,"nbmax1=",nbmax1,"in",mpi_rank
       if(nbmax1.GT.0) then
          allocate(brand(nbmax1,FBDIM))
          allocate(ibrand(nbmax1,IFBDIM))

          ibmax1=0
          do ib1=1,nbmax_tl
            if(iproc_fb(ib1)) then
              ibmax1=ibmax1+1
              brand(ibmax1,1:FBDIM)=g_brand(ib1,1:FBDIM)
             ibrand(ibmax1,1:IFBDIM)=g_ibrand(ib1,1:IFBDIM)
            endif
          enddo  
!          print*,"# of fb=",ibmax1,"in",mpi_rank
       endif  
     
       if( .NOT. nfb_io) then

       do nb=1,nb_imm_tl
         bx=g_fb_imm(nb,1)
         by=g_fb_imm(nb,2)
         bz=g_fb_imm(nb,3)

         ia=nint(bx/dx+0.5*n+0.5)
         ja=nint(by/dy+0.5*m+0.5)

         ib=ia-((npos-1)*np)
         jb=ja-((mpos-1)*mp)
         kb=l
         if ( ((ib.GE.1) .AND. (ib.LE.np))
     +    .AND. ((jb.GE.1) .AND. (jb.LE.mp)) ) then
            ik=1
            do while(bz.GE.cl_hgt(ib,jb,ik))
               ik=ik+1
            enddo
            kb=ik

            nb_imm1=nb_imm1+1
            iproc_fb_imm(nb)=.true.
         else
            iproc_fb_imm(nb)=.false.
         endif
       enddo

!       print*,"nb_imm1=",nb_imm1,"in",mpi_rank
       if(nb_imm1.GT.0) then 
         allocate(fb_imm(nb_imm1,FBDIM))
         allocate(ifb_imm(nb_imm1,IFBDIM))
       
         ib_imm=0
         do ib1=1,nb_imm_tl
           if(iproc_fb_imm(ib1)) then
             ib_imm=ib_imm+1
             fb_imm(ib_imm,1:FBDIM)=g_fb_imm(ib1,1:FBDIM)
            ifb_imm(ib_imm,1:IFBDIM)=g_ifb_imm(ib1,1:IFBDIM)
           endif
         enddo
!         print*,"# of imfb=",ib_imm,"in",mpi_rank
       endif

       deallocate(g_fb_imm)
       deallocate(g_ifb_imm)
       deallocate(iproc_fb_imm)

! nfb_imm recount
       do nb=1,nb_imm1
        bx=fb_imm(nb,1)
        by=fb_imm(nb,2)
        bz=fb_imm(nb,3)
           ia=nint(bx/dx+0.5*n+0.5)
           ja=nint(by/dy+0.5*m+0.5)
        ib=ia-((npos-1)*np)
        jb=ja-((mpos-1)*mp)
!ccc it is possible that change mpi_rank on landing. 
      if(ib.LT.0) ib = 0
      if(ib.GT.np) ib= np
      if(jb.LT.0) jb=0
      if(jb.GT.mp) jb = mp

        ik=1
        do while(bz.GE.cl_hgt(ib,jb,ik))
           ik=ik+1
        enddo
        kb=ik
        if(ifb_imm(nb,2).EQ.5) nfb_imm(ib,jb,kb)=nfb_imm(ib,jb,kb)+1
       enddo !nb_1,nb_imm1

       endif ! nfb_io        

       deallocate(g_brand) 
       deallocate(g_ibrand)
       deallocate(iproc_fb)

       endif 
       
       return
       end subroutine rinitfbrand
!234567890123456789012345678901234567890123456789012345678901234567890**

       subroutine firebrand(itfb,itr)

       use metryic
       use gridsetup
       use xvo  
       use workavg
       use msga
       use firebranda
       !use ignite
       use fireteca
       !include 'mpif.h'

       integer itfb,ierr,itabs,itr
!       integer iproc 
       integer,allocatable::nbmaxs1(:),nb_imms1(:)       
       integer,allocatable::nbcrts1(:),nb_igns1(:)
  
       itabs = itfb + itr
!       if(nbmax.LT.0) nbmax=0
!       if(nb_imm.LT.0) nb_imm=0
       call MPI_Comm_rank(mpi_comm_world,mpi_rank,ierr)

       if ((itabs.GE.fb_start) .AND.
     + MOD(itabs-fb_start,lau_frq).EQ.0) then 
         call fb_launch(nbmax,itabs)
       endif 

       if (itabs.GE.fb_start) then

         if (nbmax.GE.1) then
            call fb_move(nbmax)
         endif
 
         call fb_update(nbmax,nb_imm) 
         if (mpi_rank.EQ.0) then
            print*,"fb_updated"
            allocate(nbmaxs1(0:nproc-1))
            allocate(nb_imms1(0:nproc-1))
            allocate(nbcrts1(0:nproc-1))
            allocate(nb_igns1(0:nproc-1))
         endif

!         call fb_ignite(nb_crt,nb_ign)

      call mpi_gather(nbmax,1,mpi_integer,nbmaxs1,1,mpi_integer,
     +                0,mpi_comm_world,ierr)
      call mpi_gather(nb_imm,1,mpi_integer,nb_imms1,1,mpi_integer,
     +                0,mpi_comm_world,ierr)

!         call mpi_gather(nb_crt,1,mpi_integer,nbcrts1,1,mpi_integer,
!     +                0,mpi_comm_world,ierr)
!         call mpi_gather(nb_ign,1,mpi_integer,nb_igns1,1,mpi_integer,
!     +                0,mpi_comm_world,ierr)

!        print*,"nbmax,nb_imm,mpi_rank",nbmax,nb_imm,mpi_rank
        if(mpi_rank.EQ.0) then
!           do iproc=0,nproc-1
!              print*,"In",iproc,"nbmax=",nbmaxs1(iproc),
!     +               "nb_imm=",nb_imms1(iproc)
!           enddo
           print*,"total nbmax=,",sum(nbmaxs1)
           print*,"total nb_imm=,",sum(nb_imms1)
!           print*,"total nbcrt=,",sum(nbcrts1)
!           print*,"total nb_ign=,",sum(nb_igns1)

           deallocate(nbmaxs1)
           deallocate(nb_imms1)
           deallocate(nbcrts1)
           deallocate(nb_igns1)

        endif 

        if (MOD(itabs,frqoutput).EQ.0) then  
          call fb_writio(iunit1+2,iunit2+2,itabs)
          if(nb_imm.GT.0) then
              deallocate(fb_imm)
              deallocate(ifb_imm)
              nb_imm=0
           endif
        endif

!      call rmaxmin_int(nbrand,'nbrand',1,np,1,mp,l)
!       call rmaxmin1(nfb_imm,'nfb_im',1,np,1,mp,l)

      endif
      return
      end subroutine firebrand
!234567890123456789012345678901234567890123456789012345678901234567890**
!
      subroutine srand1 (i)
!***********************************************************************
!  This routine sets a seed for a random number generator.
!  The seed must be a positive integer.
!***********************************************************************
      implicit none
!
      integer*8 i
!
      integer*4 is1, is2
      common /randomii/ is1, is2
!
      is1=i
      is2=is1*1.8
!
      return
      end
!234567890123456789012345678901234567890123456789012345678901234567890**
!
      real*8 function rand1 ()
!***********************************************************************
!   Uniform random number generator in range [0,1] based on an
!   algorithm by Lecuyer.
!
!   P. L'Ecuyer, "Efficient and Portable Combined Random Number
!   Generators", Communications of the ACM, 31 (1988), 742--749 and 774
!
!   Programmer: Henry J. Happ III, Applied Research Assoc.  Feb 2006
!   Updated Feb 2011 to make sure there are good seeds.
!***********************************************************************
!
      implicit none
!
      integer*4 is1, is2
      common /randomii/ is1, is2
!
      integer*4 k, iz
      real*8 factor, z
      logical first
!
      integer*4 n1, n2, n3, n4, n5, n6, n7, n8, n9
      integer*4 zero, one
!
      data n1 /53668/
      data n2 /40014/
      data n3 /40692/
      data n4 /52774/
      data n5 /12211/
      data n6 /2147483563/
      data n7 /2147483399/
      data n8 /2147483562/
      data n9 /3791/
!
      data zero /0/
      data one  /1/
      data first /.true./
!
      save factor
!
      if (first) then
         first=.false.
         factor=1.0/2147483563.0
!
!   If the random number seeds are not set, set them.
!
         if (is1 .eq. 0 .or. is2 .eq. 0) then
            is1=37
            is2=66
         endif
      endif
!
      k=is1/n1
      is1=n2*mod(is1, n1)-k*n5
      if (is1 .lt. zero) is1=is1+n6
      k=is2/n4
      is2=n3*mod(is2, n4)-k*n9
      if (is2 .lt. zero) is2=is2+n7
!
      iz=(is1-n6)+is2
      if (iz .lt. one) iz=iz+n8
!
      z=iz
      rand1=z*factor
!
      return
      end

!234567890123456789012345678901234567890123456789012345678901234567890**
!
      subroutine fb_launch(nbmax2,ita)
     
      use metryic
      use gridsetup
      use xvo
      use msga
      use turba
      use fireteca
      use firebranda

      Implicit None
      !include 'mpif.h'

      integer timeA 
      real*8, external::rand1
      integer*8 :: iseed 
!      real, external::rand
      real arad,brad,crad
      real bx_new,by_new,bz_new,ux_fb,uy_fb,uz_fb
      real force_drag,force_grav,bh_max,br_max
      integer i,j,k,ia,ja,ii,ita,jseed
      integer ia_new,ja_new,ib_new,jb_new
      integer iscan_start,jscan_start,iscan_end,jscan_end
      integer nbmax1,nbmax2,nbmax3,nb_new,ibmax1,nb_cell
      logical newbrand(np,mp,lau_low_limit:lau_up_limit) 
      real,allocatable:: new_fb(:,:)
      integer,allocatable:: inew_fb(:,:)    

!     determine shape
!
!      ishape=4
!     
! ishape in rinitfbrand.f, temporarily one shape allowed
!  b_shape:  0=dsk_ub, 1=dsk_dh/dt, 2=dsk_dr/dt
!            3=cyl_ub, 4=cyl_dh/dt, 5=cyl_dr/dt
!            6=sph_ub, 7=sph_dr/dt

! scan nbrand(n,m,l) for (npos-1)*np+1 to npos*np 
!                        (mpos-1)*mp+1 to mpos*mp
!
      if(mpi_rank.EQ.0) print*,'fb_launch'
!      print*,'fb_launch,mpi=',mpi_rank
      jscan_start = (mpos-1)*mp+1
      jscan_end = mpos*mp
      iscan_start = (npos-1)*np+1
      iscan_end = npos*np

      nbmax1=nbmax2
      nb_new=0

      if (MOD(ita-fb_start,lau_frq).EQ.0) then

      do k=lau_low_limit,lau_up_limit
      do j=1,mp                           ! scan to find launching cell
      do i=1,np

         ia=(npos-1)*np+i
         ja=(mpos-1)*mp+j

         if ( ( (temps(i,j,k).GE.temphot)
     +      .AND.(rhos(i,j,k).GE.fden_limit) )
     +      .AND.(xvb(i,j,k,3).GT.0) 
!     +     .OR. (ishape.NE.8)
     +    ) then
!
!calculate cell-centered values for launching criteria
!
          if (ishape.LE.2) then                                 ! DISK
            force_drag=0.5*Cd_dn*sin(aka)
!c     +         *(xvb(i,j,k,3)**2)
!ccccc allow vertical launching
     +        *((xvb(i,j,k,1)**2)+(xvb(i,j,k,2)**2)+(xvb(i,j,k,3)**2))
     +        /xvb(i,j,k,nv)
            force_grav=rho_fb*9.81*tck
            bh_max = Cd_dn *(xvb(i,j,k,3)**2)
     +                /(2*9.81*rho_fb)/xvb(i,j,k,nv)
!c            br_max=0.0 
          elseif ((ishape.GT.2) .AND. (ishape.LE.5)) then    ! CYLINDER
            force_drag=Cd_cn*(sin(aka)**3)
!c     +         *(xvb(i,j,k,3)**2)
!cccc allow vertical launching
     +        *((xvb(i,j,k,1)**2)+(xvb(i,j,k,2)**2)+(xvb(i,j,k,3)**2))
     +         /xvb(i,j,k,nv)
            force_grav=3.141592*rho_fb*9.81*radi
            br_max = Cd_cn*(xvb(i,j,k,3)**2)
     +         /(3.141592*9.81*rho_fb)/xvb(i,j,k,nv)
          elseif (ishape.GT.5 .AND. ishape.LT.8) then          ! SPHERE
            force_drag=0.5*Cd_sp
!c     +         *(xvb(i,j,k,3)**2)
!cccc allow vertical launching
     +        *((xvb(i,j,k,1)**2)+(xvb(i,j,k,2)**2)+(xvb(i,j,k,3)**2))
     +         /xvb(i,j,k,nv)
            force_grav=4*rho_fb*9.81*radi/3
            br_max = 3*Cd_sp*(xvb(i,j,k,3)**2)
     +         /(8*9.81*rho_fb)/xvb(i,j,k,nv)
          endif
c
cc   criteria of launching CELL
c
          if ( ((ivar_fbsize.EQ.0)
     +            .AND.(force_drag.GT.force_grav))
     +       .OR.
     +         ((ivar_fbsize.EQ.1)
     +           .AND.
     +          ((((bh_max.GT.(1.1*tck_limit_d)).AND.(ishape.LE.2))
     +       .OR. ((br_max.GT.(1.1*rad_limit_c)).AND.(ishape.LE.5)
     +                .AND.(ishape.GT.2)))
     +          .OR. ((br_max.GT.rad_limit_sp).AND.(ishape.GE.6))))
     +       ) then
             newbrand(i,j,k)=.true.
             nb_new=nb_new+nb_per_cell
c             if(mpi_rank.Eq.43) 
c              print*,"launching cell=",mpi_rank,i,j,k
             nbrand(i,j,k)=nbrand(i,j,k)+nb_per_cell
          else
              newbrand(i,j,k)=.false.
          endif
       else !       
          newbrand(i,j,k)=.false.
         endif ! firebrand launching cell     

          if(ishape.EQ.8 .AND. ia.EQ.25) then !where we put there massless particles 
             if(j.EQ.3 .AND. k.EQ.10) print*,"ishape8 launch",mpi_rank
             newbrand(i,j,k)=.true.
             nb_new=nb_new+nb_per_cell
          endif 

        enddo !x
       enddo !y
      enddo !z
!      print*,"nb_new",nb_new,mpi_rank
c12     continue
     
      if (nb_new.GT.0) then  

      nbmax2=nbmax1+nb_new   
      allocate(new_fb(nbmax2,FBDIM))       
      allocate(inew_fb(nbmax2,IFBDIM))
      if (nbmax1.GT.0) then
        do ibmax1=1,nbmax1
            new_fb(ibmax1,1:FBDIM)=brand(ibmax1,1:FBDIM)
          inew_fb(ibmax1,1:IFBDIM)=ibrand(ibmax1,1:IFBDIM)
        enddo
      endif 
c      if(mpi_rank.Eq.43)  print*,"nbmax1",nbmax1,"nbmax2",
c     + nbmax2,"mpi",mpi_rank  
      nbmax3=nbmax1+1

      !call itime(timeArray)
      call system_clock(count=timeA)
      iseed=timeA+ita+mpi_rank
      call srand1(iseed)  
    
      do k=lau_low_limit,lau_up_limit
      do j=1,mp
      do i=1,np

        if(newbrand(i,j,k)) then
        ia=(npos-1)*np+i
        ja=(mpos-1)*mp+j
        nb_cell=0
        do while(nb_cell.LT.nb_per_cell)
             arad=0.000001+rand1()*0.999998
             bx_new= (ia-0.5*n+arad-1)*dx
             brad=0.000001+rand1()*0.999998
             by_new= (ja-0.5*m+brad-1)*dy
             crad=0.000001+rand1()*0.999998
             bz_new=cl_hgt(i,j,k-1)
     +             +crad*(cl_hgt(i,j,k)-cl_hgt(i,j,k-1))
          call fb_getwind(bx_new,by_new,bz_new,ux_fb,uy_fb,uz_fb)

          if (ishape.LE.2) then                                 !  DISK
               force_drag=0.5*xvb(i,j,k,nv)*Cd_dn*sin(aka)
!     +            *(uz_fb**2)
!cccc allow vertical launching
     +            *(ux_fb**2+uy_fb**2+uz_fb**2)
                bh_max = xvb(i,j,k,nv)*Cd_dn*(uz_fb*uz_fb)
     +                 /(2*9.81*rho_fb)
!cc
          elseif ((ishape.GT.2) .AND. (ishape.LE.5)) then    ! CYLINDER
               force_drag=Cd_cn*xvb(i,j,k,nv)*((sin(aka))**3)
!c     +          *(uz_fb**2)
!cccc allow vertical launching
     +            *(ux_fb*ux_fb+uy_fb*uy_fb+uz_fb*uz_fb)
            br_max = xvb(i,j,k,nv)*Cd_cn*(uz_fb**2)
     +              /(3.141592*9.81*rho_fb)
            bh_max = 0.0
!c
          elseif ((ishape.GT.5 .AND. ishape.LT.8))then           ! SPHERE
               force_drag=0.5*xvb(i,j,k,nv)*Cd_sp
!c     +         *(uz_fb**2)
!ccccc allow vertical launching
     +        *(ux_fb**2+uy_fb**2+uz_fb**2)
               br_max = 3*xvb(i,j,k,nv)*Cd_sp*(uz_fb**2)/(8*9.81*rho_fb)
               bh_max=0.0
          endif
c
cc   criteria of launching FIREBRAND
c
        if (((ivar_fbsize.EQ.0)
     +            .AND.(force_drag.GT.force_grav))
     +     .OR.
     +       ((ivar_fbsize.EQ.1).AND.
     +       ((((bh_max.GE.tck_limit_d).AND.(ishape.LE.2))
     +         .OR. ((br_max.GE.rad_limit_c).AND.(ishape.LE.5)
     +         .AND.(ishape.GT.2)))
     +         .OR. ((br_max.GE.rad_limit_sp).AND.(ishape.GE.6))))
     +     .OR. (ishape.EQ.8)
     +       ) then

ccc  store information for new firebrand.
             new_fb(nbmax3,1)=bx_new
             new_fb(nbmax3,2)=by_new 
             new_fb(nbmax3,3)=bz_new
cc check new_fb is in the domain - for debugging..
           ia_new=nint(bx_new/dx+0.5*n+0.5)
           ja_new=nint(by_new/dy+0.5*m+0.5)
           ib_new=ia_new-((npos-1)*np)
           jb_new=ja_new-((mpos-1)*mp)
!      if(ib_new.LT.1 .OR. ib_new.GT.np) then 
!         print*,"newfb ib wrong",ib_new,jb_new,ia_new,ja_new
!         print*,"i,j,k,mpi=",i,j,k,mpi_rank
!         print*,"ia,ja,x,y=",ia,ja,dx*(ia-0.5*n-0.5),dy*(ja-0.5*m-0.5)
!         print*,"k,cl_hgts",k,cl_hgt(i,j,k-1),cl_hgt(i,j,k)
!         print*,"arad,brad,crad=",arad,brad,crad         
!         print*,"x,y,z_new",bx_new,by_new,bz_new
!         print*,"bhmax,tck_limit_d",bh_max,tck_limit_d
!         print*,"u_c=",xvb(i,j,k,1)/xvb(i,j,k,nv),
!     +   xvb(i,j,k,2)/xvb(i,j,k,nv),xvb(i,j,k,3)/xvb(i,j,k,nv)     
!         print*,"u=",ux_fb,uy_fb,uz_fb
!      endif
!      if(jb_new.LT.1 .OR. jb_new.GT.mp) then
!         print*,"newfb jb wrong",ib_new,jb_new,ia_new,ja_new
!         print*,"i,j,k,mpi=",i,j,k,mpi_rank
!         print*,"ia,ja,x,y=",ia,ja,dx*(ia-0.5*n-0.5),dy*(ja-0.5*m-0.5)
!         print*,"k,cl_hgts",k,cl_hgt(i,j,k-1),cl_hgt(i,j,k)
!         print*,"arad,brad,crad=",arad,brad,crad
!         print*,"x,y,z_new",bx_new,by_new,bz_new
!         print*,"bhmax,tck_limit_d",bh_max,tck_limit_d
!         print*,"u_c=",xvb(i,j,k,1)/xvb(i,j,k,nv),
!     +   xvb(i,j,k,2)/xvb(i,j,k,nv),xvb(i,j,k,3)/xvb(i,j,k,nv)
!         print*,"u=",ux_fb,uy_fb,uz_fb
!      endif
             if(ishape.EQ.8) w_ter=.true.
             if(w_ter) v_in_rat=1

             new_fb(nbmax3,4)=v_in_rat*ux_fb
             new_fb(nbmax3,5)=v_in_rat*uy_fb
             new_fb(nbmax3,6)=v_in_rat*uz_fb
!ccc
             new_fb(nbmax3,7)=ita*dts*nts
             if(ivar_fbsize.EQ.1) then
                if(ishape.LE.2) then
                     new_fb(nbmax3,8)=bh_max/lrrat_d
                     new_fb(nbmax3,9)=bh_max
                endif
                if(ishape.GE.3) then
                     new_fb(nbmax3,8)=br_max
                     new_fb(nbmax3,9)=br_max*lrrat_c
                endif
                if(ishape.EQ.8) then
                     new_fb(nbmax3,8)=0.05
                     new_fb(nbmax3,9)=0.05
                endif
             else
                new_fb(nbmax3,8)=radi
                new_fb(nbmax3,9)=tck
             endif

             new_fb(nbmax3,10)=temps(i,j,k) 
             if (ishape.EQ.8) then
                 if(sqrt((bx_new**2)+(bz_new-90.0)**2).LT.50.0) then 
                        new_fb(nbmax3,10)=1200
                 else
                        new_fb(nbmax3,10)=300
                 endif
             endif 
             new_fb(nbmax3,11)=bx_new
             new_fb(nbmax3,12)=by_new 
             new_fb(nbmax3,13)=bz_new 
             new_fb(nbmax3,14)=rho_fb 
             new_fb(nbmax3,15)=rho_fb 
             new_fb(nbmax3,16)=aka 
             new_fb(nbmax3,17)=new_fb(nbmax3,7)
             new_fb(nbmax3,18)=new_fb(nbmax3,8) 
             new_fb(nbmax3,19)=new_fb(nbmax3,9)
             new_fb(nbmax3,20)=new_fb(nbmax3,10)      
             inew_fb(nbmax3,1)=10000*mpi_rank+nbmax3 
             inew_fb(nbmax3,2)=1
             inew_fb(nbmax3,3)=ishape 
             inew_fb(nbmax3,4)=mpi_rank
             inew_fb(nbmax3,5)=0
             inew_fb(nbmax3,6)=0
             nbmax3=nbmax3+1
             nb_cell=nb_cell+1
          endif !lauching
          enddo !(do while)
        endif ! newbrand
      enddo !x
      enddo !y
      enddo !z
     
      if (nbmax1.GT.0) then
        deallocate(brand)
        deallocate(ibrand)
      endif

      allocate(brand(nbmax2,FBDIM))
      allocate(ibrand(nbmax2,IFBDIM))

      do ii=1,nbmax2
        brand(ii,1:FBDIM)=new_fb(ii,1:FBDIM)
        ibrand(ii,1:IFBDIM)=inew_fb(ii,1:IFBDIM)
      enddo

      deallocate(new_fb)
      deallocate(inew_fb)

      endif ! formally line 13 2 ! nb_new > 0
      endif ! formally line 13 1 ! MOD

13    continue
c      print*,"end-launch,mpi=",mpi_rank,nbmax2-nbmax1,"brands lched"
      return
 
      end subroutine fb_launch

 
c234567890123456789012345678901234567890123456789012345678901234567890**

      subroutine fb_move(nbmax_m)

      use metryic
      use gridsetup
      use xvo
      use msga
      use fireteca
      use firebranda

      Implicit None

      real ak,A_c,A_ek,fb_mass,fb_dm,fb_mo,rho_g
      real u_bx,u_by,u_bz,u_abs,v_bx,v_by,v_bz,u_bh
      real w_bx,w_by,w_bz,w_bh,w_abs,v_bxo,v_byo,v_bzo,r_o,h_o
      integer ib,jb,kb,iscan_start,jscan_start,iscan_end,jscan_end
      integer nbmax_m,nb,ia,ja,ik  
    
      jscan_start = (mpos-1)*mp+1
      jscan_end = mpos*mp
      iscan_start = (npos-1)*np+1
      iscan_end = npos*np

c      if(mpi_rank.EQ.0) print*,'fb_move',nbmax_m
     
      if (nbmax_m.GT.0) then
  
      do nb=1,nbmax_m
   
      r_o=brand(nb,8)
      h_o=brand(nb,9)

!      if(ibrand(nb,4).NE.mpi_rank) then
!        print*,"there is non-mpi_rank ibrand(nb,4) in proc#."
!        print*,ibrand(nb,4),mpi_rank
!      endif       
      if(ibrand(nb,2).NE.1) then
        print*,"immo. firebrand in fb_move routine."
      endif   
c
!      if(brand(nb,8).LE.0.00002 .OR. brand(nb,9).LE.0.00002) then
      if(brand(nb,8).GT.burnout .AND. brand(nb,9).GT.burnout) then

      call fb_getwind(brand(nb,1),brand(nb,2),brand(nb,3),
     +                 u_bx,u_by,u_bz)
c      if(u_bx.GT.10) print*,"too fast after getwind",u_bx,mpi_rank
      if(ibrand(nb,3).NE.8) then 

      v_bx=brand(nb,4)
      v_by=brand(nb,5)
      v_bz=brand(nb,6)

      v_bxo=v_bx
      v_byo=v_by
      v_bzo=v_bz

      w_bx=u_bx-v_bx
      w_by=u_by-v_by
      w_bz=u_bz-v_bz

      w_abs=sqrt(w_bx**2+w_by**2+w_bz**2)
      w_bh =sqrt(w_bx**2+w_by**2)
      u_bh =sqrt(u_bx**2+u_by**2)
      u_abs=sqrt(u_bx**2+u_by**2+u_bz**2)

           ia=nint(brand(nb,1)/dx+0.5*n+0.5)
           ja=nint(brand(nb,2)/dy+0.5*m+0.5)

         ib=ia-((npos-1)*np)
         jb=ja-((mpos-1)*mp)
         kb=l
      if(ib.LT.0 .or. ib.GT.np) then
         print*,"[mv cl_hgt] ib=",ib
         print*,"x,y,z",brand(nb,1),brand(nb,2),brand(nb,3)
         print*,"ijb,ija,mpi",ib,jb,ia,ja,mpi_rank
      endif
      if(jb.LT.0 .or. jb.GT.mp) then
         print*,"[mv cl_hgt] jb=",jb
         print*,"x,y,z",brand(nb,1),brand(nb,2),brand(nb,3)
         print*,"ijb,ija,mpi",ib,jb,ia,ja,mpi_rank
      endif

         ik=1
         do while(brand(nb,3).GE.cl_hgt(ib,jb,ik))
            ik=ik+1
         enddo
         kb=ik

c      print*,"before move posi",brand(nb,1),brand(nb,2),brand(nb,3)
c      print*,"before move cell",ia,ja,ib,jb,kb
      if(kb.GT.0) rho_g=xvb(ib,jb,kb,nv)
      if(kb.LE.0) rho_g=xvb(ib,jb,1,nv)
      ak=brand(nb,16)
c
ccc property change - combusion model for firebrand
c
       call fb_property(ibrand(nb,3),fb_mo,fb_dm,brand(nb,8)
     +                 ,brand(nb,9),brand(nb,14),rho_g,w_abs) 

       fb_mass=fb_mo+fb_dm           

       if(fb_mo.LT.0) then 
          fb_mo=0
       endif
       if(fb_mass.LT.0) then
          fb_mass=0  
       endif
      if(fb_mass.LE.0.000000001) then
         fb_mass=0
      endif
      if ((fb_mo.GT.0)
     +    .OR. (fb_mass.GT.0)
     +    .OR. (brand(nb,8).GT.burnout .AND. brand(nb,9).GT.burnout)
     +    ) then
ccc smaller than 1 micgogram    

       if (ibrand(nb,3).LE.2) then                             
c DISK 
          A_c =3.141692*(brand(nb,8)**2)
          A_ek=0.5*A_c*rho_g*Cd_dn
c          print*,"A_c,A_ek",A_c,A_ek
       elseif ((ibrand(nb,3).GE.3) .AND. (ibrand(nb,3).LE.5)) then 
c CYL
          A_c =2*brand(nb,8)*brand(nb,9)
          A_ek=0.5*A_c*rho_g*Cd_cn
       else                                                    
c SPH
          A_c =3.141692*(brand(nb,8)**2)
          A_ek=0.5*A_c*rho_g*Cd_sp
       endif

       if (A_c.LT.0) then
           A_c=0
           A_ek=10000
c           write(6,*) "A_c<0",ibrand(nb,1),brand(nb,8),brand(nb,9)
       endif
       if (w_ter) then
          w_abs=sqrt(9.81*fb_mo/A_ek)
          v_bz=u_bz-sin(ak)*w_abs 
          v_bx=u_bx-cos(ak)*w_abs*u_bx/u_bh
          v_by=u_by-cos(ak)*w_abs*u_by/u_bh
       else
c          w_bh =sqrt(w_bx**2+w_by**2)
          if(w_bh.LE.0) w_bh=0.000001
          v_bz= 
     +(fb_mo*(v_bz-9.81*dt)+A_ek*dt*w_abs*(sin(ak)*u_bz+cos(ak)*w_bh) )
     +/ ( fb_mass          +A_ek*dt*w_abs *sin(ak) )

          w_bz=u_bz-v_bz
          w_abs=sqrt(w_bx**2+w_by**2+w_bz**2)

          v_bx= 
     +  ( (fb_mo*v_bx)+u_bx*A_ek*dt*w_abs*(sin(ak)-cos(ak)*w_bz/w_bh) )
     + /(   fb_mass   +     A_ek*dt*w_abs*(sin(ak)-cos(ak)*w_bz/w_bh) )
          v_by= 
     +  ( (fb_mo*v_by)+u_by*A_ek*dt*w_abs*(sin(ak)-cos(ak)*w_bz/w_bh) )
     + /(   fb_mass   +     A_ek*dt*w_abs*(sin(ak)-cos(ak)*w_bz/w_bh) )

c       print*,"w_bh=",w_bh,"w_bz=",w_bz
c       print*,"vb_xyz=",v_bx,v_by,v_bz
c       print*,"ub_xyz=",u_bx,u_by,u_bz
c       print*,"wb_xyz=",w_bx,w_by,w_bz
       endif

       else  ! ibrand(nb,3).NE.8
         v_bx=u_bx
         v_by=u_by
         v_bz=u_bz
       endif ! ibrand(nb,3).eq.8 (formally 7571)  
 
c       write(6,*) "vel.",v_bx,v_by,v_bz,ibrand(nb,1)
       brand(nb,1)=brand(nb,1)+dt*v_bx
       brand(nb,2)=brand(nb,2)+dt*v_by
       brand(nb,3)=brand(nb,3)+dt*v_bz 
       brand(nb,4)=v_bx
       brand(nb,5)=v_by
       brand(nb,6)=v_bz
       brand(nb,7)=brand(nb,7)+dt
cccc test 
           ia=nint(brand(nb,1)/dx+0.5*n+0.5)
           ja=nint(brand(nb,2)/dy+0.5*m+0.5)

         ib=ia-((npos-1)*np)
         jb=ja-((mpos-1)*mp)

      if(ib.LT.0 .OR. ib.GT.np+1) then
          print*,"too fast, ib,jb =",ib,jb,ibrand(nb,2)
          print*,"vo: ",v_bxo,v_byo,v_bzo
          print*,"mo:",fb_mo,fb_mass,ak
          print*,"shape",ibrand(nb,3),"r:",brand(nb,8),"h:",brand(nb,9)
          print*,"r_o=",r_o,"h_o=",h_o
          print*,"v: ",v_bx,v_by,v_bz
          print*,"u: ",u_bx,u_by,u_bz
          print*,"w: ",w_bx,w_by,w_bz
          print*,"pos: ",brand(nb,1),brand(nb,2),brand(nb,3)
          print*,"mpi,nb,time0",mpi_rank,nb,brand(nb,17),time
      endif  
      if(jb.LT.0 .OR. jb.GT.mp+1) then
          print*,"too fast, ib,jb=",ib,jb,ibrand(nb,2)
          print*,"vo: ",v_bxo,v_byo,v_bzo
          print*,"mo:",fb_mo,fb_mass,ak
          print*,"shape",ibrand(nb,3),"r:",brand(nb,8),"h:",brand(nb,9)
          print*,"r_o=",r_o,"h_o=",h_o
          print*,"v: ",v_bx,v_by,v_bz
          print*,"u: ",u_bx,u_by,u_bz
          print*,"w: ",w_bx,w_by,w_bz
          print*,"pos: ",brand(nb,1),brand(nb,2),brand(nb,3)
          print*,"mpi,nb,time0",mpi_rank,nb,brand(nb,17),time
      endif  
      if(ib.LT.0) ib=0
      if(ib.GT.np+1) ib=np+1
      if(jb.LT.0) jb=0
      if(jb.GT.mp+1) jb=mp+1

         ik=1
         do while(brand(nb,3).GE.cl_hgt(ib,jb,ik))
            ik=ik+1
         enddo
         kb=ik
      endif ! if not burnned out, formally line 757 2 3  4
      endif ! if not burnned out, formally line 757 1 

c      print*,"after move posi",brand(nb,1),brand(nb,2),brand(nb,3)
c      print*,"after move cell",ia,ja,ib,jb,kb
      enddo

      endif ! nbmax > 0

      return
 
      end subroutine fb_move


!234567890123456789012345678901234567890123456789012345678901234567890**

      subroutine fb_update(nbmax_o,nbmax_imm)

      use metryic
      use gridsetup
      use xvo
      use msga
      use turba
      use fireteca
      use firebranda 
      Implicit None
      !include 'mpif.h'

      real,allocatable::fb_imm_new(:,:)
      real,allocatable::send_fb(:,:,:),recv_fb(:,:,:)
      real,allocatable::s_buff1(:),r_buff1(:),s_buff2(:),r_buff2(:)
      real,allocatable::s_buff3(:),r_buff3(:),s_buff4(:),r_buff4(:)
      real,allocatable::s_buff5(:),r_buff5(:),s_buff6(:),r_buff6(:)
      real,allocatable::s_buff7(:),r_buff7(:),s_buff8(:),r_buff8(:)

      integer,allocatable::ifb_imm_new(:,:)
      integer,allocatable::isend_fb(:,:,:),irecv_fb(:,:,:)
      integer,allocatable::isbuff1(:),ir_buff1(:)      
      integer,allocatable::isbuff2(:),ir_buff2(:)       
      integer,allocatable::isbuff3(:),ir_buff3(:)       
      integer,allocatable::isbuff4(:),ir_buff4(:)       
      integer,allocatable::isbuff5(:),ir_buff5(:)       
      integer,allocatable::isbuff6(:),ir_buff6(:)       
      integer,allocatable::isbuff7(:),ir_buff7(:)       
      integer,allocatable::isbuff8(:),ir_buff8(:)       

      real,allocatable::upd_fb(:,:),upd_fb_imm(:,:)
      integer,allocatable::iupd_fb(:,:),iupd_fb_imm(:,:)
      integer,allocatable::status(:)

      integer imm_new,nupd,nsend(8),nrecv(8),nrecv_oa
      integer iimm_new,inupd,isend(8),irecv(8),nsend_oa
      integer nbmax_o,nbmax_imm,nb,i,j,k
      integer imf,ifb,lnb,ifs,ir,ierr,ibsd,nbmax_new_imm
      integer nbmax_new,nbmax1,nbmax_imm1,nbmax2,maxnsend,maxnrecv

      integer ia,ja,ib,jb,kb,ik,itag,ns_buff,nr_buff
      integer ial,jal,ibl,jbl
      integer iscan_start,jscan_start,iscan_end,jscan_end
      integer ns_fb,nr_fb,irecv_tl,nrecv_total,scount,rcount
      real bx,by,bz,bmass,zs_real,ddx,ddy
      real travel_dist,fuel_lnd_hgt ! added for landing criteria
      integer idn
c
      jscan_start = (mpos-1)*mp+1
      jscan_end = mpos*mp
      iscan_start = (npos-1)*np+1
      iscan_end = npos*np

!      if(mpi_rank.EQ.0) print*,'fb_update'

      nbmax1=nbmax_o
      nbmax_imm1=nbmax_imm
      imm_new=0
      nsend(1:8)=0
      nrecv(1:8)=0
      allocate (status(mpi_status_size))
   
!      print*,"fb_update begin, mpi=",mpi_rank
     
      if (nbmax1.GT.0) then
!
      do nb=1,nbmax1
!
       idn=0

        bx=brand(nb,1)    
        by=brand(nb,2)
        bz=brand(nb,3)
!
        ia=nint(bx/dx+0.5*n+0.5)
        ja=nint(by/dy+0.5*m+0.5)

        ib=ia-((npos-1)*np)
        jb=ja-((mpos-1)*mp)

      if(ib.LT.0 .OR. ib.GT.np+1) then
          print*,"early in update"
          print*,"too fast, ib,jb =",ib,jb,ibrand(nb,2)
          print*,"shape",ibrand(nb,3),"r:",brand(nb,8),"h:",brand(nb,9)
          print*,"mpi,nb,time0",mpi_rank,nb,brand(nb,17),time
      endif
      if(jb.LT.0 .OR. jb.GT.mp+1) then
          print*,"early in update"
          print*,"too fast, ib,jb=",ib,jb,ibrand(nb,2)
          print*,"shape",ibrand(nb,3),"r:",brand(nb,8),"h:",brand(nb,9)
          print*,"pos: ",brand(nb,1),brand(nb,2),brand(nb,3)
          print*,"mpi,nb,time0",mpi_rank,nb,brand(nb,17),time
      endif

      if(ib.LT.0) ib=0
      if(ib.GT.np+1) ib=np+1
      if(jb.LT.0) jb=0
      if(jb.GT.mp+1) jb=mp+1

        ik=1
        do while(bz.GE.cl_hgt(ib,jb,ik))
           ik=ik+1
        enddo
        kb=ik
         
!      print*,"upd posi",brand(nb,1),brand(nb,2),brand(nb,3)
!      print*,"upd cell",ia,ja,ib,jb,kb,mpi_rank

! 
!   burnout?
!       if (brand(nb,8).LE.0.00001 .OR. brand(nb,9).LE.0.00001) then
       if (brand(nb,8).LE.burnout .OR. brand(nb,9).LE.burnout) then
           brand(nb,8)=max(brand(nb,8),0.0)        
           brand(nb,9)=max(brand(nb,9),0.0)
          imm_new=imm_new+1
          ibrand(nb,2)=0 
!      print*,"burnout",nb,mpi_rank,ibrand(nb,3),brand(nb,8),brand(nb,9)
          idn=1
       endif 

c    another burnout - mass threshold
       if (ibrand(nb,3).LT.6 .and. idn.eq.0) then
         bmass=brand(nb,14)*3.141592*brand(nb,9)*brand(nb,8)*brand(nb,8)
       else
         bmass=brand(nb,14)*4*3.141592*brand(nb,8)**3/3
       endif 
       if (bmass.LT.0.000000001 .and. idn.eq.0) then
!          print*,"burnout, bmass=",bmass,mpi_rank,nb
          imm_new=imm_new+1
          ibrand(nb,2)=0
          idn=1
       endif
!
!   
!   too high?
       if (kb.EQ.l .and. idn.eq.0) then
          imm_new=imm_new+1
c       print*,"too high",nb,mpi_rank,brand(nb,1),brand(nb,2),brand(nb,3aa)
c     +       ,brand(nb,4),brand(nb,5),brand(nb,6)
          ibrand(nb,2)=2 
          idn=1
       endif
c
c   landed?  
       if(idn.eq.0) then  

!       if (kb.EQ.1) then
          ial=nint(bx/dx+0.5*n)
          jal=nint(by/dy+0.5*m)
          ibl=ial-((npos-1)*np)
          jbl=jal-((mpos-1)*mp)
          ddx=bx-(ial-0.5*n-0.5)*dx
          ddy=by-(jal-0.5*m-0.5)*dy
          zs_real=( (dx-ddx)*(dy-ddy)*zs(ibl,jbl)
     +                 + ddx*(dy-ddy)*zs(ibl+1,jbl)
     +            + (dx-ddx)*ddy*zs(ibl,jbl+1)
     +                  + ddx*ddx*zs(ibl+1,jbl+1) )/(dx*dy) 

      endif

       if (iland_fuel.EQ.1.and. idn.eq.0) then

! effective firebrand
         if(kb.EQ.1) then
           fuel_lnd_hgt=zs_real+0.5*actualfueldepth(ib,jb,kb)
         else
           fuel_lnd_hgt=zs_real+0.5*actualfueldepth(ib,jb,kb)
     +                         +actualfueldepth(ib,jb,kb-1)
         endif !kb.eq.1    


         if(rhof(ib,jb,kb).GT.fden_limit
     +   .AND. bz.LE.fuel_lnd_hgt) then
            travel_dist=sqrt((bx-brand(nb,11))**2
     +                      +(by-brand(nb,11))**2)
!            print*,"travel dist:",travel_dist,temps(1,ib,jb,kb)

            if(travel_dist.GT.2*sqrt(dx**2+dy**2)
     +      .AND.temps(ib,jb,kb).LT.temphot) then
              imm_new=imm_new+1
              ibrand(nb,2)=5
              brand(nb,3)= fuel_lnd_hgt
              nfb_imm(ib,jb,kb)=nfb_imm(ib,jb,kb)+1
!              print*,"landed on fuel cell!!!",nb,mpi_rank 
              idn=1
            elseif (travel_dist.GT.2*sqrt(dx**2+dy**2)
     +      .AND.temps(ib,jb,kb).GE.temphot) then ! ibrand(nb,2)=5 

              imm_new=imm_new+1
              brand(nb,3)= fuel_lnd_hgt
              ibrand(nb,2)=6
!              print*,"landed on fire",nb,mpi_rank 
              idn=1
            endif ! ibrand(nb,2)=6
         endif ! rhof > fden_limit .AND....

         if ((bz.LE.zs_real .OR. bz.LT.0.0).and. idn.eq.0 ) then
             brand(nb,3)=zs(ib,jb)
             imm_new=imm_new+1
             if (
     +          (temps(ib,jb,kb).LT.temphot)
     +           .AND.
     +          (rhof(ib,jb,kb).GT.fden_limit)
     +          ) then
                ibrand(nb,2)=5
c                print*,"landed on fuel",nb,mpi_rank
                 nfb_imm(ib,jb,kb)=nfb_imm(ib,jb,kb)+1
                idn=1
             elseif (temps(ib,jb,kb).GE.temphot) then
               ibrand(nb,2)=6
c                print*,"landed on fire",nb,mpi_rank
                idn=1
             else
                ibrand(nb,2)=3
c               print*,"landed on no fuel",nb,mpi_rank
                idn=1
             endif ! catergorize
         endif ! hitting the ground.. 
       endif ! iland_fuel.EQ.1 .and. idn.eq.0

c  hitting ground as landing criteria     
!
       if (iland_fuel.NE.1 .and. idn.eq.0)then
        if (kb.EQ.1) then

          if (bz.LE.zs_real .OR. bz.LT.0.0) then
             brand(nb,3)=zs(ib,jb)
             imm_new=imm_new+1
             if (
     +          (temps(ib,jb,kb).LT.temphot) 
     +           .AND.
     +          (rhos(ib,jb,kb).GT.fden_limit) 
     +           .AND. iland_fuel.EQ.0) then
                ibrand(nb,2)=5
c                print*,"landed on fuel",nb,mpi_rank
                 nfb_imm(ib,jb,kb)=nfb_imm(ib,jb,kb)+1
                idn=1
             elseif (temps(ib,jb,kb).GE.temphot) then
               ibrand(nb,2)=6
c                print*,"landed on fire",nb,mpi_rank
                idn=1
             else
                ibrand(nb,2)=3
c               print*,"landed on no fuel",nb,mpi_rank
                idn=1
             endif
         endif
      endif
      endif ! iland_fuel.NE.1 .and. idn.eq.0 
c
c  hit boundary - send or immobilize?
c  left
      if (ib.LT.1 .and. idn.eq.0) then
         if ((leftedge.EQ.1 .OR. ia.LT.1).AND.ibcx.EQ.0) then
            imm_new=imm_new+1
            ibrand(nb,2)=4 
            brand(nb,1)=-dx*(0.5*n+1) 
            idn=1
         else
            if (jb.LT.1) then
c leftbelow 777777777777777
                if ((botedge.EQ.1 .OR. ja.LT.1).AND.ibcy.EQ.0) then
                   imm_new=imm_new+1
                   ibrand(nb,2)=4
                   brand(nb,2)=-dy*(0.5*m+1)
                   idn=1
                else 
                   nsend(7)=nsend(7)+1
                   ibrand(nb,2)=17
                   if(ia.LT.1) then 
                       brand(nb,1)=brand(nb,1)+dx*n
                      ibrand(nb,5)=ibrand(nb,5)-1
                   endif
                   if(ja.LT.1) then
                       brand(nb,2)=brand(nb,2)+dy*m
                      ibrand(nb,6)=ibrand(nb,6)-1
                   endif
                   idn=1
                endif 
c 
             elseif (jb.GT.mp .and. idn.eq.0) then
c leftabove 666666666666666
                if ((topedge.EQ.1 .OR. ja.GT.m).AND.ibcy.EQ.0) then
                   imm_new=imm_new+1
                   ibrand(nb,2)=4
                   brand(nb,2)=dy*m-dy*(0.5*m+1)
                   idn=1
                else
                   nsend(6)=nsend(6)+1
                   ibrand(nb,2)=16
                   if(ia.LT.1) then 
                       brand(nb,1)=brand(nb,1)+dx*n
                      ibrand(nb,5)=ibrand(nb,5)-1
                   endif
                   if(ja.GT.m) then
                       brand(nb,1)=brand(nb,1)-dy*m
                      ibrand(nb,6)=ibrand(nb,6)+1
                   endif
                   idn=1
                endif
            else
c just left 1111111111111111
                nsend(1)=nsend(1)+1
                ibrand(nb,2)=11
                if(ia.LT.1) then
                    brand(nb,1)=brand(nb,1)+dx*n
                   ibrand(nb,5)=ibrand(nb,5)+1
                endif
                idn=1
            endif
         endif
       endif
c  right
      if (ib.GT.np .and.idn.eq.0) then
         if ((rightedge.EQ.1 .OR. ia.GT.n).AND.ibcx.EQ.0)  then
            imm_new=imm_new+1
            ibrand(nb,2)=4
            brand(nb,1)=dx*n-dx*(0.5*n+1)
            idn=1
         else
            if (jb.LT.1) then
c rightbelow 888888888888888888    
                if ((botedge.EQ.1 .OR. ja.LT.1).AND.ibcy.EQ.0) then
                   imm_new=imm_new+1
                   ibrand(nb,2)=4
                   brand(nb,2)=-dy*(0.5*m+1)
                   idn=1
                else
                   nsend(8)=nsend(8)+1
                   ibrand(nb,2)=18
                   if(ia.GT.n) then
                       brand(nb,1)=brand(nb,1)-dx*n
                      ibrand(nb,5)=ibrand(nb,5)+1
                   endif
                   if(ja.LT.1) then
                       brand(nb,2)=brand(nb,2)+dy*m
                      ibrand(nb,6)=ibrand(nb,6)-1
                   endif
                   idn=1
                endif
            elseif (jb.GT.mp) then
c rightabove 5555555555555555555    
                if ((topedge.EQ.1 .OR. ja.GT.m).AND.ibcy.EQ.0)  then
                   imm_new=imm_new+1
                   ibrand(nb,2)=4
                   brand(nb,2)=dt*m-dy*(0.5*m+1)
                   idn=1
                else
                   nsend(5)=nsend(5)+1
                   ibrand(nb,2)=15
                   if(ia.GT.n) then
                       brand(nb,1)=brand(nb,1)-dx*n
                      ibrand(nb,5)=ibrand(nb,5)+1
                   endif
                   if(ja.GT.m) then 
                       brand(nb,2)=brand(nb,2)-dy*m
                      ibrand(nb,6)=ibrand(nb,6)+1
                   endif
                   idn=1
                endif
            else
c just right 22222222222222222222
                nsend(2)=nsend(2)+1
                ibrand(nb,2)=12
                if(ia.GT.n) then
                    brand(nb,1)=brand(nb,1)-dx*n
                   ibrand(nb,5)=ibrand(nb,5)+1
                endif
                idn=1
            endif
         endif
       endif
c below
      if (jb.LT.1 .and. idn.eq.0) then
         if ((botedge.EQ.1 .OR. ja.LT.1).AND.ibcy.EQ.0) then
            imm_new=imm_new+1
            ibrand(nb,2)=4
            brand(nb,2)=-dy*(0.5*m+1) 
            idn=1
         else
c just below 44444444444444444444
             nsend(4)=nsend(4)+1
             ibrand(nb,2)=14
             if(ja.LT.1) then
                 brand(nb,2)=brand(nb,2)+dy*m
                ibrand(nb,6)=ibrand(nb,6)-1
              endif
             idn=1
         endif
       endif
c above 
      if (jb.GT.mp) then
         if ((topedge.EQ.1 .OR. ja.GT.m).AND.ibcy.EQ.0) then
            imm_new=imm_new+1
            ibrand(nb,2)=4
            brand(nb,2)=dt*m-dy*(0.5*m+1) 
            idn=1
         else
c just above 333333333333333333333
             nsend(3)=nsend(3)+1
             ibrand(nb,2)=13
             if(ja.GT.m) then
                 brand(nb,2)=brand(nb,2)-dy*m
                ibrand(nb,6)=ibrand(nb,6)+1
             endif
             idn=1
         endif
       endif
!21    continue         
      enddo
      endif  ! if nbmax1 > 0  
!23    continue
c
      nupd=imm_new+sum(nsend)
      inupd=0
      iimm_new=0
      isend=0
      irecv=0
c
c      print*,"ready to send numbers,mpi=",mpi_rank
c
ccc prepare to receive
       call MPI_Comm_rank(mpi_comm_world,mpi_rank,ierr)
c
          itag=(timestep*10000)+7901
          ns_fb=nsend(1)
          call MPI_Sendrecv(ns_fb,1,mpi_integer,peleft,itag,
     +                      nr_fb,1,mpi_integer,peright,itag,
     +                      mpi_comm_world,status,ierr) 
          nrecv(1)=nr_fb

          itag=(timestep*10000)+7906
          ns_fb=nsend(6)
          call MPI_Sendrecv(ns_fb,1,mpi_integer,peleftabove,itag,
     +                      nr_fb,1,mpi_integer,perightbelow,itag,
     +                      mpi_comm_world,status,ierr)
          nrecv(6)=nr_fb

          itag=(timestep*10000)+7907
          ns_fb=nsend(7)
          call MPI_Sendrecv(ns_fb,1,mpi_integer,peleftbelow,itag,
     +                      nr_fb,1,mpi_integer,perightabove,itag,
     +                      mpi_comm_world,status,ierr)
          nrecv(7)=nr_fb

          itag=(timestep*10000)+7902
          ns_fb=nsend(2)
          call MPI_sendrecv(ns_fb,1,mpi_integer,peright,itag,
     +                      nr_fb,1,mpi_integer,peleft,itag,
     +                      mpi_comm_world,status,ierr)
          nrecv(2)=nr_fb

          itag=(timestep*10000)+7905
          ns_fb=nsend(5)
          call MPI_sendrecv(ns_fb,1,mpi_integer,perightabove,itag,
     +                      nr_fb,1,mpi_integer,peleftbelow,itag,
     +                      mpi_comm_world,status,ierr)
          nrecv(5)=nr_fb

          itag=(timestep*10000)+7908
          ns_fb=nsend(8)
          call MPI_sendrecv(ns_fb,1,mpi_integer,perightbelow,itag,
     +                      nr_fb,1,mpi_integer,peleftabove,itag,
     +                      mpi_comm_world,status,ierr)
          nrecv(8)=nr_fb

          itag=(timestep*10000)+7903
          ns_fb=nsend(3)
          call MPI_sendrecv(ns_fb,1,mpi_integer,peabove,itag,
     +                      nr_fb,1,mpi_integer,pebelow,itag,
     +                      mpi_comm_world,status,ierr)
          nrecv(3)=nr_fb

          itag=(timestep*10000)+7904
          ns_fb=nsend(4)
          call MPI_sendrecv(ns_fb,1,mpi_integer,pebelow,itag,
     +                      nr_fb,1,mpi_integer,peabove,itag,
     +                      mpi_comm_world,status,ierr)
          nrecv(4)=nr_fb
cc allocate arrays
      nrecv_total=sum(nrecv) 
      maxnsend=maxval(nsend)
      maxnrecv=maxval(nrecv)

         allocate (fb_imm_new(imm_new+1,FBDIM))
         allocate (ifb_imm_new(imm_new+1,IFBDIM))
         allocate (send_fb(maxnsend+1,FBDIM,8))
         allocate (isend_fb(maxnsend+1,IFBDIM,8))
         allocate (recv_fb(maxnrecv+1,FBDIM,8))
         allocate (irecv_fb(maxnrecv+1,IFBDIM,8))
c
c     print*,"numbers sent, prepare to send brands,mpi=",mpi_rank 
c   
c234567890123456789012345678901234567890123456789012345678901234567890**
c  updates brand(:,:) & prepare to send
c  
      nbmax2=nbmax1-nupd
      if(nupd.GT.0) then
      do nb=1,nbmax2
       if (ibrand(nb,2).NE.1) then
         inupd=inupd+1
         if (ibrand(nb,2).LT.10) then
            iimm_new=iimm_new+1
            fb_imm_new(iimm_new,1:FBDIM)=brand(nb,1:FBDIM)
            ifb_imm_new(iimm_new,1:IFBDIM)=ibrand(nb,1:IFBDIM)
         else 
            ibsd=ibrand(nb,2)-10
            isend(ibsd)=isend(ibsd)+1
            do ifb=1,FBDIM
               send_fb(isend(ibsd),ifb,ibsd)=brand(nb,ifb)
               if (ifb.LE.IFBDIM)
     +         isend_fb(isend(ibsd),ifb,ibsd)=ibrand(nb,ifb)
            enddo
         endif
24       continue
          do while (ibrand(nbmax1-inupd+1,2).NE.1)
             lnb=nbmax1-inupd+1
             if (ibrand(lnb,2).LT.10) then
                iimm_new=iimm_new+1
                fb_imm_new(iimm_new,1:FBDIM)=brand(lnb,1:FBDIM)
                ifb_imm_new(iimm_new,1:IFBDIM)=ibrand(lnb,1:IFBDIM)
             else
                ibsd=ibrand(lnb,2)-10
                isend(ibsd)=isend(ibsd)+1
               do ifb=1,FBDIM
                  send_fb(isend(ibsd),ifb,ibsd)=brand(lnb,ifb)
                  if (ifb.LE.IFBDIM)
     +            isend_fb(isend(ibsd),ifb,ibsd)=ibrand(lnb,ifb)
               enddo
             endif
             brand(lnb,1:FBDIM)=0.0
             ibrand(lnb,1:IFBDIM)=0
             inupd=inupd+1
          enddo
        brand(nb,1:FBDIM)=brand(nbmax1-inupd+1,1:FBDIM)    
        ibrand(nb,1:IFBDIM)=ibrand(nbmax1-inupd+1,1:IFBDIM)
        brand(nbmax1-inupd+1,1:FBDIM)=0.0
        ibrand(nbmax1-inupd+1,1:IFBDIM)=0
       endif   
      enddo

      do nb=nbmax2+1,nbmax1
         if (ibrand(nb,1).GT.0) then
          if (ibrand(nb,2).EQ.1) then 
           print*,"there is a flying brand not updated.proc:",mpi_rank
           print*,"nb",nb,brand(nb,1),brand(nb,2),brand(nb,3),mpi_rank
           print*,"ibrands",ibrand(nb,1),ibrand(nb,2),ibrand(nb,3),mpi_rank
           print*,"nbmax1,2,nupd,inupd",nbmax1,nbmax2,nupd,inupd,mpi_rank

         elseif (ibrand(nb,2).NE.1 .AND. ibrand(nb,2).LE.10) then 
            iimm_new=iimm_new+1
            fb_imm_new(iimm_new,1:FBDIM)=brand(nb,1:FBDIM)
            ifb_imm_new(iimm_new,1:IFBDIM)=ibrand(nb,1:IFBDIM)

         elseif (ibrand(nb,2).GT.10) then
            ibsd=ibrand(nb,2)-10
            isend(ibsd)=isend(ibsd)+1
            do ifb=1,FBDIM
               send_fb(isend(ibsd),ifb,ibsd)=brand(nb,ifb)
               if(ifb.LE.IFBDIM)
     +       isend_fb(isend(ibsd),ifb,ibsd)=ibrand(nb,ifb)
            enddo
         endif
      endif
      enddo

      endif
!26    continue
c
c      print*,"before send/recieve,mpi=",mpi_rank
c
c234567890123456789012345678901234567890123456789012345678901234567890**
c Send-Receive 
c
        itag=(timestep*10000)+7801
        scount=nsend(1)*FBDIM
        rcount=nrecv(1)*FBDIM
        allocate(s_buff1(scount))
        allocate(r_buff1(rcount))
        ns_buff=1
        do ifb=1,FBDIM
           do ifs=1,nsend(1)
              s_buff1(ns_buff)=send_fb(ifs,ifb,1)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(s_buff1,scount,mpi_real,peleft,itag,
     +                    r_buff1,rcount,mpi_real,peright,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,FBDIM
          do ir=1,nrecv(1)
             recv_fb(ir,ifb,1)=r_buff1(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7701
        scount=nsend(1)*IFBDIM
        rcount=nrecv(1)*IFBDIM
        allocate(isbuff1(scount))
        allocate(ir_buff1(rcount))
        ns_buff=1
        do ifb=1,IFBDIM
           do ifs=1,nsend(1)
              isbuff1(ns_buff)=isend_fb(ifs,ifb,1)
              ns_buff=ns_buff+1
           enddo
        enddo

        call MPI_sendrecv(isbuff1,scount,mpi_integer,peleft,itag,
     +                    ir_buff1,rcount,mpi_integer,peright,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,IFBDIM
          do ir=1,nrecv(1)
             irecv_fb(ir,ifb,1)=ir_buff1(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7806
        scount=nsend(6)*FBDIM
        rcount=nrecv(6)*FBDIM
        allocate(s_buff6(scount))
        allocate(r_buff6(rcount))
        ns_buff=1
        do ifb=1,FBDIM
           do ifs=1,nsend(6)
              s_buff6(ns_buff)=send_fb(ifs,ifb,6)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(s_buff6,scount,mpi_real,peleftabove,itag,
     +                     r_buff6,rcount,mpi_real,perightbelow,itag,
     +                     mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,FBDIM
          do ir=1,nrecv(6)
             recv_fb(ir,ifb,6)=r_buff6(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7706
        scount=nsend(6)*IFBDIM
        rcount=nrecv(6)*IFBDIM
        allocate(isbuff6(scount))
        allocate(ir_buff6(rcount))
        ns_buff=1
        do ifb=1,IFBDIM
           do ifs=1,nsend(6)
              isbuff6(ns_buff)=isend_fb(ifs,ifb,6)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(isbuff6,scount,mpi_integer,peleftabove,itag,
     +                    ir_buff6,rcount,mpi_integer,perightbelow,itag,
     +                     mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,IFBDIM
          do ir=1,nrecv(6)
             irecv_fb(ir,ifb,6)=ir_buff6(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7807
        scount=nsend(7)*FBDIM
        rcount=nrecv(7)*FBDIM
        allocate(s_buff7(scount))
        allocate(r_buff7(rcount))
        ns_buff=1
        do ifb=1,FBDIM
           do ifs=1,nsend(7)
              s_buff7(ns_buff)=send_fb(ifs,ifb,7)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(s_buff7,scount,mpi_real,peleftbelow,itag,
     +                    r_buff7,rcount,mpi_real,perightabove,itag,
     +                     mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,FBDIM
          do ir=1,nrecv(7)
             recv_fb(ir,ifb,7)=r_buff7(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7707
        scount=nsend(7)*IFBDIM
        rcount=nrecv(7)*IFBDIM
        allocate(isbuff7(scount))
        allocate(ir_buff7(rcount))
        ns_buff=1
        do ifb=1,IFBDIM
           do ifs=1,nsend(7)
              isbuff7(ns_buff)=isend_fb(ifs,ifb,7)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(isbuff7,scount,mpi_integer,peleftbelow,itag,
     +                    ir_buff7,rcount,mpi_integer,perightabove,itag,
     +                     mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,IFBDIM
          do ir=1,nrecv(7)
             irecv_fb(ir,ifb,7)=ir_buff7(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7802
        scount=nsend(2)*FBDIM
        rcount=nrecv(2)*FBDIM
        allocate(s_buff2(scount))
        allocate(r_buff2(rcount))
        ns_buff=1
        do ifb=1,FBDIM
           do ifs=1,nsend(2)
              s_buff2(ns_buff)=send_fb(ifs,ifb,2)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(s_buff2,scount,mpi_real,peright,itag,
     +                    r_buff2,rcount,mpi_real,peleft,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,FBDIM
          do ir=1,nrecv(2)
             recv_fb(ir,ifb,2)=r_buff2(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7702
        scount=nsend(2)*IFBDIM
        rcount=nrecv(2)*IFBDIM
        allocate(isbuff2(scount))
        allocate(ir_buff2(rcount))
        ns_buff=1
        do ifb=1,IFBDIM
           do ifs=1,nsend(2)
              isbuff2(ns_buff)=isend_fb(ifs,ifb,2)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(isbuff2,scount,mpi_integer,peright,itag,
     +                    ir_buff2,rcount,mpi_integer,peleft,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,IFBDIM
          do ir=1,nrecv(2)
             irecv_fb(ir,ifb,2)=ir_buff2(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7805
        scount=nsend(5)*FBDIM
        rcount=nrecv(5)*FBDIM
        allocate(s_buff5(scount))
        allocate(r_buff5(rcount))
        ns_buff=1
        do ifb=1,FBDIM
           do ifs=1,nsend(5)
              s_buff5(ns_buff)=send_fb(ifs,ifb,5)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(s_buff5,scount,mpi_real,perightabove,itag,
     +                    r_buff5,rcount,mpi_real,peleftbelow,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,FBDIM
          do ir=1,nrecv(5)
             recv_fb(ir,ifb,5)=r_buff5(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7705
        scount=nsend(5)*IFBDIM
        rcount=nrecv(5)*IFBDIM
        allocate(isbuff5(scount))
        allocate(ir_buff5(rcount))
        ns_buff=1
        do ifb=1,IFBDIM
           do ifs=1,nsend(5)
              isbuff5(ns_buff)=isend_fb(ifs,ifb,5)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(isbuff5,scount,mpi_integer,perightabove,itag,
     +                    ir_buff5,rcount,mpi_integer,peleftbelow,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,IFBDIM
          do ir=1,nrecv(5)
             irecv_fb(ir,ifb,5)=ir_buff5(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7808
        scount=nsend(8)*FBDIM
        rcount=nrecv(8)*FBDIM
        allocate(s_buff8(scount))
        allocate(r_buff8(rcount))
        ns_buff=1
        do ifb=1,FBDIM
           do ifs=1,nsend(8)
              s_buff8(ns_buff)=send_fb(ifs,ifb,8)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(s_buff8,scount,mpi_real,perightbelow,itag,
     +                    r_buff8,rcount,mpi_real,peleftabove,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,FBDIM
          do ir=1,nrecv(8)
             recv_fb(ir,ifb,8)=r_buff8(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7708
        scount=nsend(8)*IFBDIM
        rcount=nrecv(8)*IFBDIM
        allocate(isbuff8(scount))
        allocate(ir_buff8(rcount))
        ns_buff=1
        do ifb=1,IFBDIM
           do ifs=1,nsend(8)
              isbuff8(ns_buff)=isend_fb(ifs,ifb,8)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(isbuff8,scount,mpi_integer,perightbelow,itag,
     +                    ir_buff8,rcount,mpi_integer,peleftabove,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,IFBDIM
          do ir=1,nrecv(8)
             irecv_fb(ir,ifb,8)=ir_buff8(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7803
        scount=nsend(3)*FBDIM
        rcount=nrecv(3)*FBDIM
        allocate(s_buff3(scount))
        allocate(r_buff3(rcount))
        ns_buff=1
        do ifb=1,FBDIM
           do ifs=1,nsend(3)
              s_buff3(ns_buff)=send_fb(ifs,ifb,3)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(s_buff3,scount,mpi_real,peabove,itag,
     +                    r_buff3,rcount,mpi_real,pebelow,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,FBDIM
          do ir=1,nrecv(3)
             recv_fb(ir,ifb,3)=r_buff3(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7703
        scount=nsend(3)*IFBDIM
        rcount=nrecv(3)*IFBDIM
        allocate(isbuff3(scount))
        allocate(ir_buff3(rcount))
        ns_buff=1
        do ifb=1,IFBDIM
           do ifs=1,nsend(3)
              isbuff3(ns_buff)=isend_fb(ifs,ifb,3)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(isbuff3,scount,mpi_integer,peabove,itag,
     +                    ir_buff3,rcount,mpi_integer,pebelow,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,IFBDIM
          do ir=1,nrecv(3)
             irecv_fb(ir,ifb,3)=ir_buff3(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7804
        scount=nsend(4)*FBDIM
        rcount=nrecv(4)*FBDIM
        allocate(s_buff4(scount))
        allocate(r_buff4(rcount))
        ns_buff=1
        do ifb=1,FBDIM
           do ifs=1,nsend(4)
              s_buff4(ns_buff)=send_fb(ifs,ifb,4)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(s_buff4,scount,mpi_real,pebelow,itag,
     +                    r_buff4,rcount,mpi_real,peabove,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,FBDIM
          do ir=1,nrecv(4)
             recv_fb(ir,ifb,4)=r_buff4(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo
c
        itag=(timestep*10000)+7704
        scount=nsend(4)*IFBDIM
        rcount=nrecv(4)*IFBDIM
        allocate(isbuff4(scount))
        allocate(ir_buff4(rcount))
        ns_buff=1
        do ifb=1,IFBDIM
           do ifs=1,nsend(4)
              isbuff4(ns_buff)=isend_fb(ifs,ifb,4)
              ns_buff=ns_buff+1
           enddo
        enddo
        call MPI_sendrecv(isbuff4,scount,mpi_integer,pebelow,itag,
     +                    ir_buff4,rcount,mpi_integer,peabove,itag,
     +                    mpi_comm_world,status,ierr)
       nr_buff=1
       do ifb=1,IFBDIM
          do ir=1,nrecv(4)
             irecv_fb(ir,ifb,4)=ir_buff4(nr_buff)
             nr_buff=nr_buff+1
          enddo
       enddo

!         do ifb=1,8
!            call MPI_reduce(nrecv(ifb),nrecv_oa,1,mpi_integer,mpi_sum
!     +                      ,0,mpi_comm_world,ierr)
!            call MPI_reduce(nsend(ifb),nsend_oa,1,mpi_integer,mpi_sum
!     +                      ,0,mpi_comm_world,ierr)
!            if (mpi_rank.EQ.0 .AND. nsend_oa.GT.0) then
!c               print*,"nrecv_overall(",ifb,")=",nrecv_oa
!               print*,"nsend_overall(",ifb,")=",nsend_oa
!            endif
!         enddo
c
c udpate brand,ibrand & fb_imm,ifb_imm
c
      deallocate(s_buff1,r_buff1,isbuff1,ir_buff1)
      deallocate(s_buff2,r_buff2,isbuff2,ir_buff2)
      deallocate(s_buff3,r_buff3,isbuff3,ir_buff3)
      deallocate(s_buff4,r_buff4,isbuff4,ir_buff4)
      deallocate(s_buff5,r_buff5,isbuff5,ir_buff5)
      deallocate(s_buff6,r_buff6,isbuff6,ir_buff6)
      deallocate(s_buff7,r_buff7,isbuff7,ir_buff7)
      deallocate(s_buff8,r_buff8,isbuff8,ir_buff8) 
!36    continue 

      if (nbmax_imm1+imm_new.GT.0) then
 
      allocate(upd_fb_imm(nbmax_imm1+imm_new,FBDIM))
      allocate(iupd_fb_imm(nbmax_imm1+imm_new,IFBDIM))

!cc fb_imm,ifb_imm
!c    
!37    continue
c      print*,"nbmax_imm1,mpi_rank",nbmax_imm1,mpi_rank
      if (nbmax_imm1.NE.0) then
      do imf=1,nbmax_imm1
          upd_fb_imm(imf,1:FBDIM) = fb_imm(imf,1:FBDIM) 
         iupd_fb_imm(imf,1:IFBDIM)=ifb_imm(imf,1:IFBDIM)
         iupd_fb_imm(imf,4)=mpi_rank
      enddo
      endif
cc
      do imf=1,imm_new
          upd_fb_imm(nbmax_imm1+imf,1:FBDIM) = fb_imm_new(imf,1:FBDIM)
         iupd_fb_imm(nbmax_imm1+imf,1:IFBDIM) =ifb_imm_new(imf,1:IFBDIM)
         iupd_fb_imm(nbmax_imm1+imf,4) = mpi_rank
      enddo

      if (nbmax_imm1.GT.0) then
         deallocate(fb_imm)
         deallocate(ifb_imm)
      endif

      if (nbmax_imm1+imm_new.GT.0) then
         allocate(fb_imm(nbmax_imm1+imm_new,FBDIM))
         allocate(ifb_imm(nbmax_imm1+imm_new,IFBDIM))
         fb_imm=upd_fb_imm
         ifb_imm=iupd_fb_imm
         deallocate(upd_fb_imm)
         deallocate(iupd_fb_imm)
      endif
 
      endif ! nbmax_imm1+imm_new.GT.0        
cc
cc brand, ibrand 
      nbmax_new=nbmax2+nrecv_total
      nbmax_new_imm=nbmax_imm1+imm_new

      if(nbmax_new.GT.0) then

      allocate(upd_fb(nbmax_new,FBDIM))
      allocate(iupd_fb(nbmax_new,IFBDIM))
      do ib=1,nbmax2
          upd_fb(ib,1:FBDIM) = brand(ib,1:FBDIM)
         iupd_fb(ib,1:IFBDIM)=ibrand(ib,1:IFBDIM)
         iupd_fb(ib,4) = mpi_rank
      enddo
c 
      irecv_tl=0
      do ir=1,8
         irecv=0
         do ib=1,nrecv(ir)
            irecv_tl=irecv_tl+1
            upd_fb(nbmax2+irecv_tl,1:FBDIM)=recv_fb(ib,1:FBDIM,ir)
            iupd_fb(nbmax2+irecv_tl,1:IFBDIM)=irecv_fb(ib,1:IFBDIM,ir)
            iupd_fb(nbmax2+irecv_tl,2)=1
            iupd_fb(nbmax2+irecv_tl,4)=mpi_rank
         enddo
      enddo
      if (irecv_tl.NE.nrecv_total) then
      print*,"counted nrecv_tot;",irecv_tl,"but nrec_total:",nrecv_total
      stop
      endif 
c
39    continue
      if(nbmax_o.GT.0) then
         deallocate(brand)
         deallocate(ibrand)
      endif

      allocate(brand(nbmax_new,FBDIM))
      allocate(ibrand(nbmax_new,IFBDIM))
          brand=upd_fb
         ibrand=iupd_fb
      deallocate(upd_fb)
      deallocate(iupd_fb)
c   
      else
         if(nbmax_o.GT.0) then
              deallocate(brand)
              deallocate(ibrand)
          endif
      endif

40    continue 

      deallocate(fb_imm_new)
      deallocate(ifb_imm_new)
      deallocate(send_fb)
      deallocate(isend_fb)
      deallocate(recv_fb)
      deallocate(irecv_fb)

cc update nbrand, nfb_imm_
      do k=1,l
!         do j=jscan_start,jscan_end
!            do i=iscan_start,iscan_end
          do j=1,mp
             do i=1,np
               nbrand(i,j,k)=0
!               nfb_imm(i,j,k)=0
            enddo
         enddo
      enddo
 
      do nb=1,nbmax_new
        bx=brand(nb,1) 
        by=brand(nb,2)
        bz=brand(nb,3)
           ia=nint(bx/dx+0.5*n+0.5)
           ja=nint(by/dy+0.5*m+0.5)
        ib=ia-((npos-1)*np)
        jb=ja-((mpos-1)*mp)

      if(ib.LT.0 .or. ib.GT.np) then
         print*,"[updd fb] ib=",ib
         print*,"x,y,z",bx,by,bz
         print*,"ijb,ija,mpi",ib,jb,ia,ja,mpi_rank
      endif
      if(jb.LT.0 .or. jb.GT.mp) then
         print*,"[updd fb] jb=",jb
         print*,"x,y,z",bx,by,bz
         print*,"ijb,ija,mpi",ib,jb,ia,ja,mpi_rank
      endif

        ik=1
        do while(bz.GE.cl_hgt(ib,jb,ik))
           ik=ik+1
        enddo
        kb=ik

      nbrand(ib,jb,kb)=nbrand(ib,jb,kb)+1
      enddo

!      do nb=1,nbmax_new_imm
!        bx=fb_imm(nb,1)
!        by=fb_imm(nb,2)
!        bz=fb_imm(nb,3)
!c            ia=int(bx/dx+0.5*n+0.5)
!c            ja=int(by/dy+0.5*m+0.5)
!           ia=nint(bx/dx+0.5*n+0.5)
!           ja=nint(by/dy+0.5*m+0.5)
!        ib=ia-((npos-1)*np)
!        jb=ja-((mpos-1)*mp)
!ccc it is possible that change mpi_rank on landing. 
!      if(ib.LT.0) ib = 0 
!      if(ib.GT.np) ib= np  
!      if(jb.LT.0) jb=0 
!      if(jb.GT.mp) jb = mp
!
!        ik=1
!        do while(bz.GE.cl_hgt(ib,jb,ik))
!           ik=ik+1
!        enddo
!        kb=ik
!
!56    continue
!c      nfb_imm(ia,ja,kb)=nfb_imm(ia,ja,kb)+1
!      enddo
      nbmax_o=nbmax_new
      nbmax_imm=nbmax_imm1+imm_new
66    continue
      deallocate(status)
      return
      end subroutine fb_update

c234567890123456789012345678901234567890123456789012345678901234567890**

      subroutine fb_ignite(n_crt,n_ign)

      use metryic
      use gridsetup
      use xvo
      use msga
      use turba
      use fireteca
      use firebranda
      !use ignite

      Implicit None
      !include 'mpif.h'

      integer nbmax_im,n_crt,n_ign
      integer i,j,k
      real temp_cri

      n_crt=0
      n_ign=0
      do k=1,l
        do j=1,mp
           do i=1,np
             if(nfb_imm(i,j,k).GT.2*nb_per_cell) then
                 n_crt = n_crt+1
!                print*,'nfb_imm',i,j,k,nfb_imm(i,j,k),temps(1,i,j,k)
                 if((temps(i,j,k).LT.1000 !targettemp(1) !targettemp=1000
     +      .and.   rhof(i,j,k).GT.fden_limit) )then

                   temp_cri=((nfb_imm(i,j,k)+0.1)/nb_per_cell-2.0)*350
     +                      +300
                   if(temp_cri.GT.1000) temp_cri=1000
                   if(temps(i,j,k).LT.temp_cri) then
                      temps(i,j,k)=temps(i,j,k)+
     +                ((0.1*350  ! ramprate(1)          !ramprate(1)=350 degree/sec
     +                *(((nfb_imm(i,j,k)+0.1)/nb_per_cell)-2))*0.5
     +                *dt)
                      sies(i,j,k)=cpsolid(i,j,k)*temps(i,j,k)
                      n_ign= n_ign+1
!                      print*,'new temps',temps(1,i,j,k),temps(2,i,j,k)
                    endif ! temps < temp_cri

                endif ! temps < tagettemp, which is a cap.

             endif ! nfb_imm>2*nb_per_cell 
           enddo
        enddo
      enddo


      return
      end subroutine fb_ignite
c234567890123456789012345678901234567890123456789012345678901234567890**


c234567890123456789012345678901234567890123456789012345678901234567890**

      subroutine fb_writio(wiunit1,wiunit2,it3)
      use metryic
      use gridsetup
      use msga
      use fireteca
      use firebranda  

      Implicit None

      !include 'mpif.h'
c
      character base*5 /'fbvo.'/, fname1*7
      character basei*5 /'imfb.'/,fnamei1*7
      character base_t*5 /'fbtx.'/,fname_t1*7
      character basei_t*5 /'imtx.'/,fnamei_t1*7
      character basen*5 /'nbpc.'/,fnamen1*7 
      character fname2*8,fnamei2*8,fname_t2*8,fnamei_t2*8
      character fname3*9,fnamei3*9,fname_t3*9,fnamei_t3*9
      character fname4*10,fnamei4*10,fname_t4*10,fnamei_t4*10
      character fname5*11,fnamei5*11,fname_t5*11,fnamei_t5*11
      character fname6*12,fnamei6*12,fname_t6*12,fnamei_t6*12
      character fnamen2*8
      character fnamen3*9
      character fnamen4*10
      character fnamen5*11
      character fnamen6*12
 
      integer,allocatable::nbmaxs(:),nbmaxDIM(:),inbmaxDIM(:)
      integer,allocatable::nb_imms(:),nb_immDIM(:),inb_immDIM(:)
      integer,allocatable::disp(:),idp(:),disp_im(:),idp_im(:)

      real,allocatable::sbf_fb(:),rbf_fb(:),sbf_im(:),rbf_im(:)
      integer,allocatable::isbf_fb(:),irbf_fb(:),isbf_im(:),irbf_im(:)

      real,allocatable::g_brand(:,:),g_fb_imm(:,:)
      integer,allocatable::g_ibrand(:,:),g_ifb_imm(:,:)
  
      integer nbmax_tl,nb_imm_tl,it3,wiunit1,wiunit2
c      integer data_fb,data_im,idata_fb,idata_im,ncnt,icnt
      integer nbmaxD,nb_immD,inbmaxD,inb_immD
      integer ipr,n1,n2,n3,ierr
      !logical :: new

      nbmaxD=nbmax*FBDIM
      inbmaxD=nbmax*IFBDIM
      nb_immD=nb_imm*FBDIM
      inb_immD=nb_imm*IFBDIM
      if (mpi_rank.EQ.0) then
          allocate(nbmaxs(0:nproc-1))
          allocate(nbmaxDIM(0:nproc-1))
          allocate(inbmaxDIM(0:nproc-1))
          allocate(disp(0:nproc-1))
          allocate(idp(0:nproc-1))
          disp(0)=0
          idp(0)=0
         if(.NOT. nfb_io) then
          allocate(nb_imms(0:nproc-1))
          allocate(nb_immDIM(0:nproc-1))
          allocate(inb_immDIM(0:nproc-1))
          allocate(disp_im(0:nproc-1))
          allocate(idp_im(0:nproc-1))
          disp_im(0)=0
          idp_im(0)=0
         endif
      endif
      call mpi_gather(nbmax,1,mpi_integer,nbmaxs,1,mpi_integer,
     +                0,mpi_comm_world,ierror)
      if(.NOT. nfb_io) call mpi_gather(nb_imm,1,mpi_integer,nb_imms
     +               ,1,mpi_integer,
     +                0,mpi_comm_world,ierror)
      if(mpi_rank.EQ.0) then
         do ipr=0,nproc-1
            nbmaxDIM(ipr)=nbmaxs(ipr)*FBDIM
            inbmaxDIM(ipr)=nbmaxs(ipr)*IFBDIM
            if(.NOT. nfb_io) then
               nb_immDIM(ipr)=nb_imms(ipr)*FBDIM
               inb_immDIM(ipr)=nb_imms(ipr)*IFBDIM
            endif
            if (ipr.GT.0) then
               disp(ipr)=disp(ipr-1)+nbmaxDIM(ipr-1)
               idp(ipr)=idp(ipr-1)+inbmaxDIM(ipr-1)
               if(.NOT. nfb_io) then
               disp_im(ipr)=disp_im(ipr-1)+nb_immDIM(ipr-1)
               idp_im(ipr)=idp_im(ipr-1)+inb_immDIM(ipr-1)
               endif
            endif
         enddo

         nbmax_tl=sum(nbmaxs)
         if(.NOT. nfb_io) nb_imm_tl=sum(nb_imms)

         allocate(g_brand(nbmax_tl+1,FBDIM))
         allocate(g_ibrand(nbmax_tl+1,IFBDIM))
         allocate(rbf_fb((nbmax_tl+1)*FBDIM))
         allocate(irbf_fb((nbmax_tl+1)*IFBDIM))
      endif
c
        allocate(sbf_fb(nbmax*FBDIM+1))
        allocate(isbf_fb(nbmax*IFBDIM+1))
c  
        do n1=1,nbmax
          do n2=1,FBDIM
             sbf_fb(n2+(n1-1)*FBDIM)=brand(n1,n2)
          enddo
          do n3=1,IFBDIM
             isbf_fb(n3+(n1-1)*IFBDIM)=ibrand(n1,n3)
          enddo
        enddo
       call mpi_gatherv(sbf_fb,nbmaxD,mpi_real,rbf_fb,
     +                   nbmaxDIM,disp,mpi_real,0,mpi_comm_world,ierr)
       call mpi_gatherv(isbf_fb,inbmaxD,mpi_integer,irbf_fb,
     +              inbmaxDIM,idp,mpi_integer,0,mpi_comm_world,ierr)
c       if (mpi_rank.EQ.0) print*,"gathering brand",nbmax_tl,it3
c
      if (mpi_rank.eq.0) then
         if (it3.LT.100) then
          write (fname1,'(a,i2.2)') base,it3
          if(tx_out) write (fname_t1,'(a,i2.2)') base_t,it3
          open (wiunit1,file=fname1,form='unformatted',status='unknown')
        if(tx_out) open (wiunit1+2,file=fname_t1,form='formatted',status='unknown')
         elseif (it3.LT.1000) then
          write (fname2,'(a,i3.3)') base,it3
          if(tx_out) write (fname_t2,'(a,i3.3)') base_t,it3
          open (wiunit1,file=fname2,form='unformatted',status='unknown')
        if(tx_out) open (wiunit1+2,file=fname_t2,form='formatted',status='unknown')
         elseif (it3.LT.10000) then
          write (fname3,'(a,i4.4)') base,it3
          if(tx_out) write (fname_t3,'(a,i4.4)') base_t,it3
          open (wiunit1,file=fname3,form='unformatted',status='unknown')
        if(tx_out) open (wiunit1+2,file=fname_t3,form='formatted',status='unknown')
         elseif (it3.LT.100000) then
          write (fname4,'(a,i5.5)') base,it3
          if(tx_out) write (fname_t4,'(a,i5.5)') base_t,it3
          open (wiunit1,file=fname4,form='unformatted',status='unknown')
        if(tx_out) open (wiunit1+2,file=fname_t4,form='formatted',status='unknown')
         elseif (it3.LT.1000000) then
          write (fname5,'(a,i6.6)') base,it3
          if(tx_out) write (fname_t5,'(a,i6.6)') base_t,it3
          open (wiunit1,file=fname5,form='unformatted',status='unknown')
        if(tx_out) open (wiunit1+2,file=fname_t5,form='formatted',status='unknown')
         elseif (it3.LT.10000000) then
          write (fname6,'(a,i7.7)') base,it3
          if(tx_out) write (fname_t6,'(a,i7.7)') base_t,it3
          open (wiunit1,file=fname6,form='unformatted',status='unknown')
        if(tx_out) open (wiunit1+2,file=fname_t6,form='formatted',status='unknown')
         endif
          if(tx_out) write (wiunit1+2,*) 'i',nbmax_tl,FBDIM,IFBDIM
          write (wiunit1) nbmax_tl,FBDIM,IFBDIM
        do n1=1,nbmax_tl
           do n2=1,FBDIM
                g_brand(n1,n2)=rbf_fb(n2+(n1-1)*FBDIM)
           enddo
           do n3=1,IFBDIM
                g_ibrand(n1,n3)=irbf_fb(n3+(n1-1)*IFBDIM)
           enddo
           if(tx_out) then
c     +     write (wiunit1+2,119) n1,g_brand(n1,1),g_brand(n1,2),
c     +      g_brand(n1,3),g_brand(n1,4),g_brand(n1,5),g_brand(n1,6),
c     +      g_brand(n1,7),g_brand(n1,8)*1000,g_brand(n1,9)*1000,
c     +      g_brand(n1,10),g_brand(n1,11),g_brand(n1,12),g_brand(n1,13),
c     +      g_brand(n1,14),g_brand(n1,15),g_brand(n1,16),g_brand(n1,17),
c     +      g_brand(n1,18)*1000,g_brand(n1,19)*1000,g_brand(n1,20),
c     +      g_ibrand(n1,1),g_ibrand(n1,2),g_ibrand(n1,3),g_ibrand(n1,4),
c     +      g_ibrand(n1,5),g_ibrand(n1,6)
       write (wiunit1+2,*) 'b',g_brand(n1,1),g_brand(n1,2),g_brand(n1,3)
            endif
c           write (wiunit1,118) n1,g_ibrand(n1,1),g_ibrand(n1,2),
c     +                            g_ibrand(n1,3),g_ibrand(n1,4)

        enddo 
c        print*, "writing brand: it=",it3,"nbmax=",nbmax_tl
      write (wiunit1) g_brand
      write (wiunit1) g_ibrand
        deallocate(g_brand) 
        deallocate(g_ibrand)
        deallocate(rbf_fb)
        deallocate(irbf_fb)
      endif
      deallocate(sbf_fb)
      deallocate(isbf_fb)
cc
cc gather and write g_fb_im
c   
      if(.NOT. nfb_io) then

      if(mpi_rank.EQ.0) then
          allocate(g_fb_imm(nb_imm_tl+1,FBDIM))
          allocate(g_ifb_imm(nb_imm_tl+1,IFBDIM))
          allocate(rbf_im((nb_imm_tl+1)*FBDIM))
          allocate(irbf_im((nb_imm_tl+1)*IFBDIM))
      endif
      allocate(sbf_im(nb_imm*FBDIM+1))
      allocate(isbf_im(nb_imm*IFBDIM+1))
c
      do n1=1,nb_imm
         do n2=1,FBDIM
            sbf_im(n2+(n1-1)*FBDIM)=fb_imm(n1,n2)
         enddo
         do n3=1,IFBDIM
            isbf_im(n3+(n1-1)*IFBDIM)=ifb_imm(n1,n3)
        enddo
      enddo
      call mpi_gatherv(sbf_im,nb_immD,mpi_real,rbf_im,
     +              nb_immDIM,disp_im,mpi_real,0,mpi_comm_world,ierr)
      call mpi_gatherv(isbf_im,inb_immD,mpi_integer,irbf_im,
     +        inb_immDIM,idp_im,mpi_integer,0,mpi_comm_world,ierr)
c
      if (mpi_rank.eq.0) then
c         print*,"gthering imms",nb_imm_tl,it3
        if (it3.LT.100) then
           write (fnamei1,'(a,i2.2)') basei,it3
           if(tx_out) write (fnamei_t1,'(a,i2.2)') basei_t,it3
         open (wiunit2,file=fnamei1,form='unformatted',status='unknown')
       if(tx_out) open (wiunit2+2,file=fnamei_t1,form='formatted',status='unknown')
        elseif (it3.LT.1000) then
           write (fnamei2,'(a,i3.3)') basei,it3
           if(tx_out) write (fnamei_t2,'(a,i3.3)') basei_t,it3
         open (wiunit2,file=fnamei2,form='unformatted',status='unknown')
       if(tx_out) open (wiunit2+2,file=fnamei_t2,form='formatted',status='unknown')
        elseif (it3.LT.10000) then
           write (fnamei3,'(a,i4.4)') basei,it3
           if(tx_out) write (fnamei_t3,'(a,i4.4)') basei_t,it3
         open (wiunit2,file=fnamei3,form='unformatted',status='unknown')
       if(tx_out) open (wiunit2+2,file=fnamei_t3,form='formatted',status='unknown')
        elseif (it3.LT.100000) then
           write (fnamei4,'(a,i5.5)') basei,it3
           if(tx_out) write (fnamei_t4,'(a,i5.5)') basei_t,it3
         open (wiunit2,file=fnamei4,form='unformatted',status='unknown')
       if(tx_out) open (wiunit2+2,file=fnamei_t4,form='formatted',status='unknown')
        elseif (it3.LT.1000000) then
           write (fnamei5,'(a,i6.6)') basei,it3
           if(tx_out) write (fnamei_t5,'(a,i6.6)') basei_t,it3
         open (wiunit2,file=fnamei5,form='unformatted',status='unknown')
       if(tx_out) open (wiunit2+2,file=fnamei_t5,form='formatted',status='unknown')
        elseif (it3.LT.10000000) then
           write (fnamei6,'(a,i7.7)') basei,it3
           if(tx_out) write (fnamei_t6,'(a,i7.7)') basei_t,it3
         open (wiunit2,file=fnamei6,form='unformatted',status='unknown')
       if(tx_out) open (wiunit2+2,file=fnamei_t6,form='formatted',status='unknown')
        endif
         write (wiunit2) nb_imm_tl,FBDIM,IFBDIM
         if(tx_out) write (wiunit2+2,*) 'i',nb_imm_tl,FBDIM,IFBDIM
        do n1=1,nb_imm_tl
           do n2=1,FBDIM
               g_fb_imm(n1,n2)=rbf_im(n2+(n1-1)*FBDIM)
           enddo
           do n3=1,IFBDIM
               g_ifb_imm(n1,n3)=irbf_im(n3+(n1-1)*IFBDIM)
           enddo

!        if((g_ifb_imm(n1,2).EQ.3).or.(g_ifb_imm(n1,2).EQ.5).or.(g_ifb_imm(n1,2).EQ.6).and.(tx_out) )then 
         if((g_ifb_imm(n1,2).EQ.5).and.(tx_out))then

       write (wiunit2+2,119) n1,g_fb_imm(n1,1),g_fb_imm(n1,2),
     +  g_fb_imm(n1,3),g_fb_imm(n1,4),g_fb_imm(n1,5),g_fb_imm(n1,6),
     +  g_fb_imm(n1,7),g_fb_imm(n1,8)*1000,g_fb_imm(n1,9)*1000,
     +  g_fb_imm(n1,10),g_fb_imm(n1,11),g_fb_imm(n1,12),g_fb_imm(n1,13),
     +  g_fb_imm(n1,14),g_fb_imm(n1,15),g_fb_imm(n1,16),g_fb_imm(n1,17),
     +  g_fb_imm(n1,18)*1000,g_fb_imm(n1,19)*1000,g_fb_imm(n1,20),
     +  g_ifb_imm(n1,1),g_ifb_imm(n1,2),g_ifb_imm(n1,3),g_ifb_imm(n1,4),
     +  g_ifb_imm(n1,5),g_ifb_imm(n1,6)
c         write (wiunit2+2,*) 'b',g_fb_imm(n1,1),g_fb_imm(n1,2),g_fb_imm(n1,3)
c        write (wiunit2,118) n1,g_ifb_imm(n1,1),g_ifb_imm(n1,2),
c     +                           g_ifb_imm(n1,3),g_ifb_imm(n1,4)
         endif
          enddo
         write (wiunit2) g_fb_imm
         write (wiunit2) g_ifb_imm
          deallocate(g_fb_imm)
          deallocate(g_ifb_imm)
          deallocate(rbf_im)
          deallocate(irbf_im)
         endif
        deallocate(sbf_im)
        deallocate(isbf_im)

       else !nfb_io==.true. 

       if (it3.LT.100) then
cJEC ! KOO open commented out.  041514 
         write (fnamen1,'(a,i2.2)') basen,it3
       elseif (it3.LT.1000) then
         write (fnamen2,'(a,i3.3)') basen,it3
       elseif (it3.LT.10000) then
         write (fnamen3,'(a,i4.4)') basen,it3
       elseif (it3.LT.100000) then
         write (fnamen4,'(a,i5.5)') basen,it3
       elseif (it3.LT.1000000) then
         write (fnamen5,'(a,i6.6)') basen,it3
       elseif (it3.LT.10000000) then
         write (fnamen6,'(a,i7.7)') basen,it3
cJEC
       endif
      !new=.true.
!! KOO 041514 MPI_IO 
!!!
!!! MPI_IO write - writeio3D_fb 

       if (it3.LT.100) then
!      call writeio3D_fb(trim(fnamen1))
!       call writeio3D_fb(fnamen1)
       elseif (it3.LT.1000) then
!      call writeio3D_fb(trim(fnamen2))
!       call writeio3D_fb(fnamen2)
       elseif (it3.LT.10000) then
!      call writeio3D_fb(trim(fnamen3))
!       call writeio3D_fb(fnamen3)
       elseif (it3.LT.100000) then
!      call writeio3D_fb(trim(fnamen4))
!       call writeio3D_fb(fnamen4)
       elseif (it3.LT.1000000) then
!      call writeio3D_fb(trim(fnamen5))
!       call writeio3D_fb(fnamen5)
      elseif (it3.LT.10000000) then
!      call writeio3D_fb(trim(fnamen6))
!       call writeio3D_fb(fnamen6)
       endif
!!!

      endif !nfb_io 
c
cc
       if (mpi_rank.EQ.0) then
         deallocate(nbmaxs)
         deallocate(nbmaxDIM)
         deallocate(inbmaxDIM)
         deallocate(disp)
         deallocate(idp)
         close(wiunit1)
         if(tx_out) close(wiunit1+2)
         close(wiunit2)

         if(.NOT. nfb_io) then
           deallocate(nb_imms)
           deallocate(nb_immDIM)
           deallocate(inb_immDIM)
           deallocate(disp_im)
           deallocate(idp_im)
           if(tx_out) close(wiunit2+2)
         endif

       endif
119    format(I10,F9.2,F9.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,
     +           F9.2,F9.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,
     +           I10,I8,I8,I8,I3,I3)
c
c118    format(I5,I10,I8,I8,I8)
c
      return
      end subroutine fb_writio

c234567890123456789012345678901234567890123456789012345678901234567890**


      subroutine fb_getwind(bx,by,bz,bu,bv,bw)

      use workavg
      use metryic
      use xvo
      use gridsetup
      use msga
      use firebranda


      Implicit None

      integer ia1,ja1,ib1,jb1,kb1,kb2,kb3,kb4,ik,i_vel
c      real bx_l,by_l
      real ddx,ddy,ddz1,ddz2,ddz3,ddz4,dzk1,dzk2,dzk3,dzk4
      real bx,by,bz,bu,bv,bw
      real vol_tt,vol(8)
      real vel_cel(8,3)
      integer isw

      isw=0
!      if(((bx.EQ.-604) .AND.((by.GT.23).AND.(by.LT.24))
!     +    .AND. ((bz.GT.104) .AND. (bz.LT.105)))) isw=1
           ia1=nint(bx/dx+0.5*n)
           ja1=nint(by/dy+0.5*m)

      ib1=ia1-((npos-1)*np)
      jb1=ja1-((mpos-1)*mp)
      kb1=l
      kb2=l
      kb3=l
      kb4=l

      if(ib1.LT.0 .or. ib1.GT.np) then
         print*,"[getwind]ib1=",ib1
         print*,"x,y,z",bx,by,bz
         print*,"ijb1,ija1,mpi",ib1,jb1,ia1,ja1,mpi_rank
      endif
      if(jb1.LT.0 .or. jb1.GT.mp) then
         print*,"[getwind]jb1=",jb1
         print*,"x,y,z",bx,by,bz
         print*,"ijb1,ija1,mpi",ib1,jb1,ia1,ja1,mpi_rank
      endif
!      if(isw.EQ.1) then
!         print*,"[getwind]ib1=",ib1
!         print*,"[getwind]jb1=",jb1
!         print*,"ijb1,ija1,mpi",ib1,jb1,ia1,ja1,mpi_rank
!      endif
!1501  continue
      ik=1
      do while(bz.GE.zposition(ib1,jb1,ik))
!         if(ik.GE.l+1) print*,ik,bz,zposition(ib1,jb1,ik),1501 
         ik=ik+1
      enddo
      kb1=ik-1  ! be carefull with -1....not comparing with cl_hgt

!1502  continue
      ik=1
      do while(bz.GE.zposition(ib1,jb1+1,ik))
!         if(ik.GE.l+1) print*,ik,bz,zposition(ib1,jb1+1,ik),1502
         ik=ik+1
      enddo
      kb2=ik-1

!1503  continue
      ik=1
      do while(bz.GE.zposition(ib1+1,jb1,ik))
!         if(ik.GE.l+1) print*,ik,bz,zposition(ib1+1,jb1,ik),1503
         ik=ik+1
      enddo
      kb3=ik-1

!1504  continue
      ik=1
      do while(bz.GE.zposition(ib1+1,jb1+1,ik))
!         if(ik.GE.l+1) print*,ik,bz,zposition(ib1+1,jb1+1,ik),1504
         ik=ik+1
      enddo
      kb4=ik-1

! NaN Check koo

!      if((ib1.eq.0 .or.ib1.eq.1) .and. (jb1.eq.0 .or. jb1.eq.1))
!     +   then
!       print*,"isw on,",ia1,ja1,kb1 
!       isw=1
!      endif 

      ddx=bx-(ia1-0.5*n-0.5)*dx
      ddy=by-(ja1-0.5*m-0.5)*dy

      ddz1=bz-zposition(ib1,jb1,kb1)
      dzk1=dzk(ib1,jb1,kb1)

      ddz2=bz-zposition(ib1,jb1+1,kb2)
      dzk2=dzk(ib1,jb1+1,kb2)

      ddz3=bz-zposition(ib1+1,jb1,kb3)
      dzk3=dzk(ib1+1,jb1,kb3)

      ddz4=bz-zposition(ib1+1,jb1+1,kb4)
      dzk4=dzk(ib1+1,jb1+1,kb4)

!      if(isw.EQ.1) then
!         print*,"bx,by,x,y",bx,by,x(ia1),y(ja1)
!         print*,"dx,dy,ddx,ddy",dx,dy,ddx,ddy
!         print*,"ddz",ddz1,ddz2,ddz3,ddz4
!         print*,"dzk",dzk1,dzk2,dzk3,dzk4
!      endif

      vol(1)=(dx-ddx)*(dy-ddy)*(dzk4-ddz4)
      vol(2)=(dx-ddx)*ddy*(dzk3-ddz3)
      vol(3)=ddx*(dy-ddy)*(dzk2-ddz2)
      vol(4)=ddx*ddy*(dzk1-ddz1)
      vol(5)=(dx-ddx)*(dy-ddy)*ddz4
      vol(6)=(dx-ddx)*ddy*ddz3
      vol(7)=ddx*(dy-ddy)*ddz2
      vol(8)=ddx*ddy*ddz1
      vol_tt=vol(1)+vol(2)+vol(3)+vol(4)+vol(5)+vol(6)+vol(7)+vol(8)
c  
      do i_vel=1,3
       if(kb1.GT.0.AND.kb1.LT.l) then
        vel_cel(1,i_vel)=xvb(ib1,jb1,kb1,i_vel)/xvb(ib1,jb1,kb1,nv)
        vel_cel(5,i_vel)=xvb(ib1,jb1,kb1+1,i_vel)/xvb(ib1,jb1,kb1+1,nv)
       elseif(kb1.EQ.0) then
        vel_cel(1,i_vel)=0.0
        vel_cel(5,i_vel)=xvb(ib1,jb1,kb1+1,i_vel)/xvb(ib1,jb1,kb1+1,nv)
       elseif(kb1.EQ.l) then
        vel_cel(1,i_vel)=xvb(ib1,jb1,kb1,i_vel)/xvb(ib1,jb1,kb1,nv)
        vel_cel(5,i_vel)=vel_cel(1,i_vel)
       endif

       if(kb2.GT.0.AND.kb2.LT.l) then
        vel_cel(2,i_vel)=xvb(ib1,jb1+1,kb2,i_vel)/xvb(ib1,jb1+1,kb2,nv)
        vel_cel(6,i_vel)
     +              =xvb(ib1,jb1+1,kb2+1,i_vel)/xvb(ib1,jb1+1,kb2+1,nv)
       elseif(kb2.EQ.0) then
        vel_cel(2,i_vel)=0.0
        vel_cel(6,i_vel)
     +              =xvb(ib1,jb1+1,kb2+1,i_vel)/xvb(ib1,jb1+1,kb2+1,nv)
       elseif(kb2.EQ.l) then
        vel_cel(2,i_vel)=xvb(ib1,jb1+1,kb2,i_vel)/xvb(ib1,jb1+1,kb2,nv)
        vel_cel(6,i_vel)=vel_cel(2,i_vel)
       endif

       if(kb3.GT.0.AND.kb3.LT.l) then
        vel_cel(3,i_vel)=xvb(ib1+1,jb1,kb3,i_vel)/xvb(ib1+1,jb1,kb3,nv)
        vel_cel(7,i_vel)
     +              =xvb(ib1+1,jb1,kb3+1,i_vel)/xvb(ib1+1,jb1,kb3+1,nv)
       elseif(kb3.EQ.0) then
        vel_cel(3,i_vel)=0.0
        vel_cel(7,i_vel)
     +              =xvb(ib1+1,jb1,kb3+1,i_vel)/xvb(ib1+1,jb1,kb3+1,nv)
       elseif(kb3.EQ.l) then
        vel_cel(3,i_vel)=xvb(ib1+1,jb1,kb3,i_vel)/xvb(ib1+1,jb1,kb3,nv)
        vel_cel(7,i_vel)=vel_cel(3,i_vel)
       endif


       if(kb4.GT.0.AND.kb4.LT.l) then
        vel_cel(4,i_vel)
     +              =xvb(ib1+1,jb1+1,kb4,i_vel)/xvb(ib1+1,jb1+1,kb4,nv)
        vel_cel(8,i_vel)
     +          =xvb(ib1+1,jb1+1,kb4+1,i_vel)/xvb(ib1+1,jb1+1,kb4+1,nv)
       elseif(kb4.EQ.0) then
        vel_cel(4,i_vel)=0.0
        vel_cel(8,i_vel)
     +          =xvb(ib1+1,jb1+1,kb4+1,i_vel)/xvb(ib1+1,jb1+1,kb4+1,nv)
       elseif(kb4.EQ.l) then
        vel_cel(4,i_vel)
     +              =xvb(ib1+1,jb1+1,kb4,i_vel)/xvb(ib1+1,jb1+1,kb4,nv)
        vel_cel(8,i_vel)=vel_cel(4,i_vel)
       endif
      enddo
!      if(isw.EQ.1)then
!       print*," in getwind",ib1,jb1,kb1,mpi_rank
!       print*,"vel1:",vel_cel(1,1),vel_cel(1,2),vel_cel(1,3)
!       print*,"vel2:",vel_cel(2,1),vel_cel(2,2),vel_cel(2,3)
!       print*,"vel3:",vel_cel(3,1),vel_cel(3,2),vel_cel(3,3)
!       print*,"vel4:",vel_cel(4,1),vel_cel(4,2),vel_cel(4,3)
!       print*,"vel5:",vel_cel(5,1),vel_cel(5,2),vel_cel(5,3)
!       print*,"vel6:",vel_cel(6,1),vel_cel(6,2),vel_cel(6,3)
!       print*,"vel7:",vel_cel(7,1),vel_cel(7,2),vel_cel(7,3)
!       print*,"vel8:",vel_cel(8,1),vel_cel(8,2),vel_cel(8,3)
!       print*,"vols1-4",vol(1),vol(2),vol(3),vol(4)
!       print*,"vols5-8",vol(5),vol(6),vol(7),vol(8)
!       print*,"total vol",vol_tt
!      endif
c      
         bu= ( vol(1)*vel_cel(1,1) + vol(2)*vel_cel(2,1) +
     +         vol(3)*vel_cel(3,1) + vol(4)*vel_cel(4,1) +
     +         vol(5)*vel_cel(5,1) + vol(6)*vel_cel(6,1) +
     +         vol(7)*vel_cel(7,1) + vol(8)*vel_cel(8,1))/vol_tt
c
         bv= ( vol(1)*vel_cel(1,2) + vol(2)*vel_cel(2,2) +
     +         vol(3)*vel_cel(3,2) + vol(4)*vel_cel(4,2) +
     +         vol(5)*vel_cel(5,2) + vol(6)*vel_cel(6,2) +
     +         vol(7)*vel_cel(7,2) + vol(8)*vel_cel(8,2))/vol_tt
c
         bw= ( vol(1)*vel_cel(1,3) + vol(2)*vel_cel(2,3) +
     +         vol(3)*vel_cel(3,3) + vol(4)*vel_cel(4,3) +
     +         vol(5)*vel_cel(5,3) + vol(6)*vel_cel(6,3) +
     +         vol(7)*vel_cel(7,3) + vol(8)*vel_cel(8,3))/vol_tt
!      if(isw.EQ.1) then
!          print*,"bu,bv,bw",bu,bv,bw
!          print*,"bx,by,bz",bx,by,bz
!          print*,"ia,ja,ib,jb,mpi_rank",ia1,ja1,ib1,jb1,mpi_rank
!          print*,"kbs:",kb1,kb2,kb3,kb4
!          print*,"vols1",vol(1),vol(2),vol(3),vol(4)
!          print*,"vols5",vol(5),vol(6),vol(7),vol(8)
!          print*,"vels1:",vel_cel(1,1),vel_cel(2,1),vel_cel(3,1),vel_cel(4,1)
!          print*,"vels5:",vel_cel(5,1),vel_cel(6,1),vel_cel(7,1),vel_cel(8,1)
!       print*,"ddzs:",ddz1,ddz2,ddz3,ddz4
!       print*,"dzks:",dzk1,dzk2,dzk3,dzk4
!       endif
       return
       end subroutine fb_getwind

c234567890123456789012345678901234567890123456789012345678901234567890**

      subroutine fb_property(ishape2,mo_f,dm_f,r_f,t_f,rho_f,
     +                        rho_g,w_abs)

      use metryic
      use gridsetup
      use xvo
      use msga
      use firebranda

      Implicit None

      real visair,fo,B_wood,gamma_wood,w_abs
      real mo_f,dm_f,r_f,t_f,rho_f,rho_g
      integer ishape2

cccc  properties
c
        visair=2.3e-5
        B_wood=1.2
        gamma_wood=0.5
c
c     
      if (ishape2.LE.5) then
           mo_f=3.141592*t_f*(r_f**2)*rho_f
      else
           mo_f=4*3.141592*(r_f**3)*rho_f/3
      endif

       if (mod(ishape2,3).EQ.0) then
          dm_f = 0
c   
ccc   Dr.Woychesse's model - disk dh/dt
ca
       elseif (ishape2.EQ.1) then
          fo = -0.39924
          if(w_ter) w_abs=sqrt(2*9.81*t_f*rho_f/(rho_g*Cd_dn))
          t_f=t_f+dt*fo*8*rho_g/rho_f*sqrt(w_abs*visair/r_f)/3
          dm_f=rho_f*3.141592*t_f*(r_f**2)-mo_f
c
ccc   Emmons' model I - disk dr/dt
c
       elseif (ishape2.EQ.2) then
          fo = -0.52055
          if(w_ter) w_abs=sqrt(2*9.81*t_f*rho_f/(rho_g*Cd_dn))
          r_f=r_f+dt*fo*rho_g/rho_f*sqrt(w_abs*visair/t_f)
          dm_f=rho_f*3.141592*t_f*(r_f**2)-mo_f
c
ccc   Emmons' model II - cylinder dh/dt
c
       elseif (ishape2.EQ.4) then
          fo = -0.52055
          if(w_ter) w_abs=sqrt(3.141592*9.81*r_f*rho_f/(rho_g*Cd_cn)) 
          t_f=t_f+dt*fo*1.113*rho_g/rho_f*sqrt(w_abs*visair/r_f)
          dm_f=rho_f*3.141592*t_f*(r_f**2)-mo_f
c
ccc   Emmons' model III - cylinder dr/dt`
c
      elseif (ishape2.EQ.5) then
          fo = -0.52055
          if(w_ter) w_abs=sqrt(3.141592*9.81*r_f*rho_f/(rho_g*Cd_cn))
          r_f=r_f+dt*fo*rho_g/rho_f*sqrt(w_abs*visair/(r_f*3.141592))
          dm_f=rho_f*3.141592*t_f*(r_f**2)-mo_f
c 
ccc  droplet model - sphere dr/dt
c
      elseif (ishape2.EQ.7) then
          r_f=r_f-dt*visair*(rho_g/rho_f)*log(1+B_wood)/r_f
          dm_f=rho_f*4*3.141592*(r_f**3)/3-mo_f
      endif
30    continue
      if(dm_f.GT.0) then
c        print*,"adding mass!"
c        print*,"ishape2,mo_f,dm_f,r_f,t_f,rho_f,rho_g,w_abs"
c        print*,ishape2,mo_f,dm_f,r_f,t_f,rho_f,rho_g,w_abs
      endif  
c      if(r_f.LT.0) r_f=0
      return
      end subroutine fb_property

!!!! © Copyright 2007 Los Alamos National Security, LLC All rights reserved, Authors: Eunmo Koo (koo_e@lanl.gov) 


