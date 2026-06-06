c2345678***************************************************
       module advectionScheme
       Implicit None

       save

       contains


!*****************************************************************************************!
!     this routine computes the advection on the large time steps and
!     the nv variables. It can use either donorcell or fct
!
!****************************************************************************************!
      subroutine advec(xv,il,iu,jl,ju,lls)
      use gridsetup
      use xve
      use metryic
      use msga
      use workavg

      Implicit None

      integer,intent(in) :: il,iu,jl,ju,lls

      real xv(il:iu,jl:ju,lls,nv)
      integer::kv

      !if (mpi_rank.eq.0) write(6,*) 'begin advec' 
      call setAdvectiveVelocities(uavg,vavg,oavg,gc1,gc2,gc3,il,iu,jl,ju,lls,iord) ! this is advvel in code. 
       
      !if (mpi_rank.eq.0) write(6,*) 'advective velocities computed' 
      if(iwallclock.EQ.1) call computeWallTime(wtime2,4,'advl')
      if(iwallclock.EQ.1) wtime2=MPI_Wtime() 

      if (impdataold==0) then
!!!FP092019 the new implementation of mpdata might be incorrect in case
!of topo 
        do kv=1,nv
          if (kv.ne.3) then 
            call mpdata(xv(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),h,il,iu,jl,ju,lls,0)
          else  ! specific bc for w for bottom (w(0)=-w(1))
            call mpdata(xv(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),h,il,iu,jl,ju,lls,1)
          endif
          !call updated(xv(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
        enddo
      else
        call mpdataold(xv,h,xe,il,iu,jl,ju,lls,nv)
      endif

      if(iwallclock.EQ.1) call computeWallTime(wtime2,11,'advc')
      !if (mpi_rank.eq.0) write(6,*) 'end advection'
      return
      end subroutine advec
!*****************************************************************************************!
      subroutine mpdataold(x,h,xe,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use msga
      use advo ! code2-where u1,u2,u3 are defined. 

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

!      real,dimension(il:iu+1,jl:ju,lls) :: u1        ! code2-they are defined as global arrays
!      real,dimension(il:iu,jl:ju+1,lls) :: u2
!      real,dimension(il:iu,jl:ju,lls+1) :: u3
      real,dimension(il:iu,jl:ju,lls,nvp) :: x,xe
      real,dimension(il:iu,jl:ju,lls) :: h
      real,allocatable:: v1mp(:,:,:,:)                 ! v123 is declared as 3d array  
      real,allocatable:: v2mp(:,:,:,:)                 ! changed as v123mp for code2
      real,allocatable:: v3mp(:,:,:,:)
      real,allocatable:: f1(:,:,:,:)
      real,allocatable:: f2(:,:,:,:)
      real,allocatable:: f3(:,:,:,:)
      real,allocatable:: f1o(:,:,:,:)
      real,allocatable:: f2o(:,:,:,:)
      real,allocatable:: f3o(:,:,:,:)
      real,allocatable:: cpmp(:,:,:,:)                 ! code2-they are defined as global arrays
      real,allocatable:: cnmp(:,:,:,:)
      real,allocatable:: mxo(:,:,:,:)
      real,allocatable:: mno(:,:,:,:)
      real,allocatable::  a(:,:,:,:)
      real,allocatable:: mx(:,:,:,:)
      real,allocatable:: mn(:,:,:,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: n3,n3m,kv,i,j,k,ip,im,jp,jm,kp,km,itrfct
      integer :: itr,ilft,illim,iulim,jllim,julim,ia,ja
      real :: ep,rhoin,rhoout,ain,aout,rmxuse,rmnuse,a1p,a2p,a1n,a2n
      real :: tmpp,tmpn,tmp,v1d,v2d,c1,c2

      n3=l+1
      n3m=n3-1
      ep=1.e-10

!      if(j3.eq.1) then
      allocate (v1mp(1-ih:np+ih+1,1-ih:mp+ih, l,nv))
      allocate (v2mp(1-ih:np+ih,1-ih:mp+ih+1, l,nv))
      allocate (v3mp(1-ih:np+ih,1-ih:mp+ih, l+1,nv))
      allocate (f1(1-ih:np+ih+1,1-ih:mp+ih, l,nv))
      allocate (f2(1-ih:np+ih,1-ih:mp+ih+1, l,nv))
      allocate (f3(1-ih:np+ih,1-ih:mp+ih, l+1,nv))
      allocate (f1o(1-ih:np+ih+1,1-ih:mp+ih, l,nv))
      allocate (f2o(1-ih:np+ih,1-ih:mp+ih+1, l,nv))
      allocate (f3o(1-ih:np+ih,1-ih:mp+ih, l+1,nv))
      allocate (cpmp(1-ih:np+ih, 1-ih:mp+ih, l,nv))
      allocate (cnmp(1-ih:np+ih, 1-ih:mp+ih, l,nv))
      allocate (mxo(1-ih:np+ih, 1-ih:mp+ih, l,nv))
      allocate (mno(1-ih:np+ih, 1-ih:mp+ih, l,nv))
      allocate (a(1-ih:np+ih, 1-ih:mp+ih, l,nv-1))
      allocate (mx(1-ih:np+ih, 1-ih:mp+ih, l,nv-1))
      allocate (mn(1-ih:np+ih, 1-ih:mp+ih, l,nv-1))
!      endif ! code2 j3 is removed 
c
c
c     transfer data from shared arrays to msg arrays
c
c updated (array,dim in x, dim in y, dim in z, ngp in x, ngp in y,
c0-east west and north-south, 1 all directions, 2 east-west only,
c3-north-south only

      do kv=1,nv
      call updated(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1,0)  ! code2-0 is added

      do k=1,n3m+1
       do j=1,mp
        do i=1,np
            v3mp(i,j,k,kv) = u3(i,j,k)
          enddo
        enddo
      end do

      do k=1,n3m
        do j=1,mp
          do i=1,np+1
            v1mp(i,j,k,kv) = u1(i,j,k)
          end do
        end do

        do i=1,np
          do j=1,mp+1
            v2mp(i,j,k,kv) = u2(i,j,k)
          end do
        end do
      end do
      enddo

      if(nonos.eq.1) then
       do kv=1,nv-1
        do k=1,l
         do j=1,mp
          do i=1,np
            a(i,j,k,kv)=x(i,j,k,kv)/x(i,j,k,nv)
          enddo ! i=1:np
         enddo ! j=1:mp
       enddo ! k=1:l
       call updated(a(1-ih,1-ih,1,kv),a(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1,0) ! code2-0 is added

         do k=1,n3m
            km=max0(k-1,1  )
            kp=min0(k+1,n3m)
            do j=1,mp
               if (botdedge.eq.1 .and. j.eq.1) then !add d rrl
                  jm = ibcy*(-1) + ibyo*1
               else
                  jm = j - 1
               end if
               if (topdedge.eq.1 .and. j.eq.mp) then !add d rrl
                  jp = ibcy*(mp+2) + ibyo*mp
               else
                  jp = j + 1
               end if
               do i=1,np
                  if (leftdedge.eq.1 .and. i.eq.1) then !add d rrl
                     im = ibcx*(-1) + ibxo*1
                  else
                     im = i - 1
                  end if
                  if (rightdedge.eq.1 .and. i.eq.np) then !add d rrl
                     ip = ibcx*(np+2) + ibxo*np
                  else
                     ip = i + 1
                  end if
           mx(i,j,k,kv)=amax1(a(im,j,k,kv),a(i,j,k,kv),a(ip,j,k,kv),
     .      a(i,jm,k,kv),a(i,jp,k,kv),a(i,j,kp,kv),a(i,j,km,kv))
           mn(i,j,k,kv)=amin1(a(im,j,k,kv),a(i,j,k,kv),a(ip,j,k,kv),
     .      a(i,jm,k,kv),a(i,jp,k,kv),a(i,j,kp,kv),a(i,j,km,kv))
               end do ! i=1:np
            end do ! j=1:mp
         end do ! k=1:n3m(l)
      enddo ! kv=1:nv-1 
      end if ! nonos==1

      if(nonosold.eq.1) then
         do kv=1,nv
         do k=1,n3m
            km=max0(k-1,1  )
            kp=min0(k+1,n3m)
            do j=1,mp
               if (botdedge.eq.1 .and. j.eq.1) then !add d rrl
                  jm = ibcy*(-1) + ibyo*1
               else
                  jm = j - 1
               end if
               if (topdedge.eq.1 .and. j.eq.mp) then !add d rrl
                  jp = ibcy*(mp+2) + ibyo*mp
               else
                  jp = j + 1
               end if
               do i=1,np
                  if (leftdedge.eq.1 .and. i.eq.1) then !add d rrl
                     im = ibcx*(-1) + ibxo*1
                  else
                     im = i - 1
                  end if
                  if (rightdedge.eq.1 .and. i.eq.np) then !add d rrl
                     ip = ibcx*(np+2) + ibxo*np
                  else
                     ip = i + 1
                  end if
           mxo(i,j,k,kv)=amax1(x(im,j,k,kv),x(i,j,k,kv),x(ip,j,k,kv),
     .             x(i,jm,k,kv),x(i,jp,k,kv),x(i,j,kp,kv),x(i,j,km,kv))
           mno(i,j,k,kv)=amin1(x(im,j,k,kv),x(i,j,k,kv),x(ip,j,k,kv),
     .             x(i,jm,k,kv),x(i,jp,k,kv),x(i,j,kp,kv),x(i,j,km,kv))
               end do ! i=1:np
            end do ! j=1:mp
         end do ! k=1:nm3(l)
      enddo ! kv=1:nv
      end if ! nonosold==1

      c1=1.
      c2=0.

                         do 30 itr=1,iord ! loop to LINE 1402  
      do kv=1,nv

      call updated(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1,0) ! code2-0 is added

       ilft=1+leftdedge*1                       !add d rrl
         do k=1,n3m
            do j=1,mp
               do i=ilft,np
                  f1(i,j,k,kv)=donor(c1*x(i-1,j,k,kv)+c2,
     .                               c1*x(i,j,k,kv)+c2,
     .                 v1mp(i,j,k,kv))
               end do
            end do
         end do

       if (rightdedge.eq.0) then                       !add d rrl
          call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
       else
          call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
       end if


      if (leftdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do j=1,mp
               f1(1 ,j,k,kv)= ibcx*f1(-1,j,k,kv)
     .              +ibxo*donor(c1*xe(1,j,k,kv)+c2,
     .                          c1*x(1,j,k,kv)+c2,
     .              v1mp(1,j,k,kv))
            end do
         end do
      end if

      if (rightdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do j=1,mp
               f1(np+1,j,k,kv)= ibcx*f1(np+3,j,k,kv)
     .  +ibxo*donor(c1*x(np,j,k,kv)+c2,
     .              c1*xe(np,j,k,kv)+c2,
     .         v1mp(np+1,j,k,kv))
            end do
         end do
      end if


      if (topdedge.eq.0) then                       !add d rrl
          call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
      else
          call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
      end if

      if (botdedge.eq.0) then                       !add d rrl
         do k=1,n3m
            do j=1,mp
               do i=1,np
                  f2(i,j,k,kv)=donor(c1*x(i,j-1,k,kv)+c2,
     .                               c1*x(i,j,k,kv)+c2,
     .                 v2mp(i,j,k,kv))
               end do
            end do
         end do
      end if
      if (botdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do j=2,mp
               do i=1,np
                  f2(i,j,k,kv)=donor(c1*x(i,j-1,k,kv)+c2,
     .                              c1*x(i,j,k,kv)+c2,
     .                 v2mp(i,j,k,kv))
               end do
            end do
         end do
      end if

      if (topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
      end if

      if (botdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do i=1,np
               f2(i,1 ,k,kv)= ibcy*f2(i,-1,k,kv)
     .              +ibyo*donor(c1*xe(i,1,k,kv)+c2,
     .                          c1*x(i,1,k,kv)+c2,
     .              v2mp(i,1,k,kv))

            end do
         end do
      end if
      if (topdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do i=1,np
               f2(i,mp+1,k,kv)= ibcy*f2(i,mp+3,k,kv)
     .              +ibyo*donor(c1*x(i,mp,k,kv)+c2,
     .                          c1*xe(i,mp,k,kv)+c2,
     .              v2mp(i,mp+1,k,kv))
            end do
         end do
      end if


      do 333 k=2,n3m
      do 333 j=1,mp
      do 333 i=1,np
  333 f3(i,j,k,kv)=donor(c1*x(i,j,k-1,kv)+c2,c1*x(i,j,k,kv)+c2,
     .v3mp(i,j,k,kv))

         if(ibctopbot.eq.0) then
         do j=1,mp
            do i=1,np
               f3(i,j, 1,kv)=-f3(i,j,2,kv)
               f3(i,j,n3,kv)=-f3(i,j,n3m,kv)
            enddo
         enddo
         else
         do j=1,mp
            do i=1,np
               f3(i,j, 1,kv)=0.
               f3(i,j,n3,kv)=0.
            enddo
         enddo
         endif


      if (rightdedge.eq.0) then                       !add d rrl
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
      else
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
      endif
      if(topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
      end if

      do k=1,n3m
      do j=1,mp
      do i=1,np
         x(i,j,k,kv)=x(i,j,k,kv)-( f1(i+1,j,k,kv)-f1(i,j,k,kv)
     .                            +f2(i,j+1,k,kv)-f2(i,j,k,kv)
     .                            +f3(i,j,k+1,kv)-f3(i,j,k,kv) )
     .  /h(i,j,k)
      end do
      end do
      end do

      enddo ! kv=1:nv


      if(itr.eq.iord) go to 6  !! Skip if iord.eq.1: start of mpdata when iord=2: escaping to LINE 1403
       c1=0.
       c2=1.
       do kv=1,nv
       do k=1,n3m
       do j=1,mp
       do i=1,np+1
         f1(i,j,k,kv)=v1mp(i,j,k,kv)
         v1mp(i,j,k,kv)=0.
       end do
       end do
       end do

       do k=1,n3m
       do j=1,mp+1
       do i=1,np
         f2(i,j,k,kv)=v2mp(i,j,k,kv)
         v2mp(i,j,k,kv)=0.
       end do
       end do
       end do
 
       do k=1,n3
       do j=1,mp
       do i=1,np
         f3(i,j,k,kv)=v3mp(i,j,k,kv)
         v3mp(i,j,k,kv)=0.
       end do
       end do
       end do

!!! compute antidiffusive velocities in x direction
      call updated(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1,0) ! code2-0 is added

      if (topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1,0) ! code2-0 is added
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1,0) ! code2-0 is added 
      end if
      call updated(f3(1-ih,1-ih,1,kv),f3(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1,0) ! code2-0 is added

      illim = 1 + 1*leftdedge                       !add d rrl
      iulim = np
      jllim = 1  + 1*botdedge                       !add d rrl
      julim = mp - 1*topdedge                       !add d rrl

      do k=2,n3-2
      do j=jllim,julim
      do i=illim,iulim
        v1mp(i,j,k,kv)=vdyf(x(i-1,j,k,kv),x(i,j,k,kv),f1(i,j,k,kv),
     *              .5*(h(i-1,j,k)+h(i,j,k)))
     *   +vcorr(f1(i,j,k,kv),
     *      f2(i-1,j,k,kv)+f2(i-1,j+1,k,kv)
     *     +f2(i,j+1,k,kv)+f2(i,j,k,kv),
     *      x(i-1,j-1,k,kv),x(i,j-1,k,kv),
     *      x(i-1,j+1,k,kv),x(i,j+1,k,kv),
     *    .5*(h(i-1,j,k)+h(i,j,k)))
     *   +vcorr(f1(i,j,k,kv),
     *      f3(i-1,j,k,kv)+f3(i-1,j,k+1,kv)
     *     +f3(i,j,k+1,kv)+f3(i,j,k,kv),
     *       x(i-1,j,k-1,kv),x(i,j,k-1,kv),
     *       x(i-1,j,k+1,kv),x(i,j,k+1,kv),
     *    .5*(h(i-1,j,k)+h(i,j,k)))
      end do
      end do
      end do

      if (rightdedge.eq.0) then                       !add d rrl
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0) ! code2-0 is added
      else
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0) ! code2-0 is added
      end if

      if(ibcy.eq.1) then

         if (botdedge.eq.1) then                       !add d rrl
            illim = 1  + 1*leftdedge                       !add d rrl
            iulim = np
            do k=2,n3-2
               do i=illim,iulim
                  v1mp(i,1,k,kv)=vdyf(x(i-1,1,k,kv),x(i,1,k,kv),
     *            f1(i,1,k,kv),
     *         .5*(h(i-1,1,k)+h(i,1,k)))
     *     +vcorr(f1(i,1,k,kv),
     *            f2(i-1,1,k,kv)+f2(i-1,2,k,kv)
     *           +f2(i,2,k,kv)+f2(i,1,k,kv),
     *             x(i-1,-1,k,kv),x(i,-1,k,kv),
     *             x(i-1,2,k,kv),x(i,2,k,kv),
     *         .5*(h(i-1,1,k)+h(i,1,k)))
     *     +vcorr(f1(i,1,k,kv),
     *            f3(i-1,1,k,kv)+f3(i-1,1,k+1,kv)
     *           +f3(i,1,k+1,kv)+f3(i,1,k,kv),
     *             x(i-1,1,k-1,kv),x(i,1,k-1,kv),
     *             x(i-1,1,k+1,kv),x(i,1,k+1,kv),
     *             .5*(h(i-1,1,k)+h(i,1,k)))
               enddo
            enddo
         end if

         if (rightdedge.eq.0) then                       !add d rrl
            call updated(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0) ! code2-0 is added
         else
            call updated(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0) ! code2-0 is added
         end if

         if (topdedge.eq.1) then                       !add d rrl
            illim = 1 + 1*leftdedge                       !add d rrl
            iulim = np
            do k=2,n3-2
               do i=illim,iulim
                  v1mp(i,mp,k,kv)=v1mp(i,mp+1,k,kv)
               enddo
            enddo
         end if
      end if

      if(idiv.eq.1) then
         illim = 1 + 1*leftdedge                       !add d rrl
         iulim = np
         jllim = 1  + (1-ibcy)*botdedge                       !add d rrl
         julim = mp + (-1+ibcy)*topdedge                       !add d rrl
      do 511 k=2,n3-2
      do 511 j=jllim,julim
      do 511 i=illim,iulim
      v1d=-vdiv1(f1(i-1,j,k,kv),f1(i,j,k,kv),f1(i+1,j,k,kv),
     *                 .5*(h(i-1,j,k)+h(i,j,k)))
     *  -vdiv2(f1(i,j,k,kv),f2(i-1,j+1,k,kv),
     *         f2(i,j+1,k,kv),f2(i-1,j,k,kv),
     *   f2(i,j,k,kv),   .5*(h(i-1,j,k)+h(i,j,k)))
     *  -vdiv2(f1(i,j,k,kv),f3(i-1,j,k+1,kv),
     *         f3(i,j,k+1,kv),f3(i-1,j,k,kv),
     *   f3(i,j,k,kv),   .5*(h(i-1,j,k)+h(i,j,k)))
  511 v1mp(i,j,k,kv)=v1mp(i,j,k,kv)+
     *     (pp(v1d)*x(i-1,j,k,kv)-pn(v1d)*x(i,j,k,kv))
      endif ! idiv.eq.1

!!! compute antidiffusive velocities in y direction
         illim = 1  + 1*leftdedge                       !add d rrl
         iulim = np - 1*rightdedge                       !add d rrl
         jllim = 1 + 1*botdedge                       !add d rrl
         julim = mp

      do k=2,n3-2
      do j=jllim,julim
      do i=illim,iulim
      v2mp(i,j,k,kv)=vdyf(x(i,j-1,k,kv),x(i,j,k,kv),f2(i,j,k,kv),
     *              .5*(h(i,j-1,k)+h(i,j,k)))
     * +vcorr(f2(i,j,k,kv),
     *      f1(i,j-1,k,kv)+f1(i,j,k,kv)
     *     +f1(i+1,j,k,kv)+f1(i+1,j-1,k,kv),
     *       x(i-1,j-1,k,kv),x(i-1,j,k,kv),
     *       x(i+1,j-1,k,kv),x(i+1,j,k,kv),
     *               .5*(h(i,j-1,k)+h(i,j,k)))
     * +vcorr(f2(i,j,k,kv),
     *      f3(i,j-1,k,kv)+f3(i,j,k,kv)
     *     +f3(i,j,k+1,kv)+f3(i,j-1,k+1,kv),
     *       x(i,j-1,k-1,kv),x(i,j,k-1,kv),
     *       x(i,j-1,k+1,kv),x(i,j,k+1,kv),
     *               .5*(h(i,j-1,k)+h(i,j,k)))
      end do
      end do
      end do

      if(ibcx.eq.1) then
         jllim = 1 + 1*botdedge                       !add d rrl
         julim = mp
         if (leftdedge.eq.1) then                       !add d rrl
            do k=2,n3-2
               do j=jllim,julim
                  v2mp(1,j,k,kv)=vdyf(x(1,j-1,k,kv),
     *                     x(1,j,k,kv),f2(1,j,k,kv),
     *                 .5*(h(1,j-1,k)+h(1,j,k)))
     *             +vcorr(f2(1,j,k,kv),
     *                    f1(1,j-1,k,kv)+f1(1,j,k,kv)
     *                   +f1(2,j,k,kv)+f1(2,j-1,k,kv),
     *                     x(-1,j-1,k,kv),x(-1,j,k,kv),
     *                     x(2,j-1,k,kv),x(2,j,k,kv),
     *             .5*(h(1,j-1,k)+h(1,j,k)))
     *             +vcorr(f2(1,j,k,kv),
     *                    f3(1,j-1,k,kv)+f3(1,j,k,kv)
     *                   +f3(1,j,k+1,kv)+f3(1,j-1,k+1,kv),
     *                     x(1,j-1,k-1,kv),x(1,j,k-1,kv),
     *                     x(1,j-1,k+1,kv),x(1,j,k+1,kv),
     *             .5*(h(1,j-1,k)+h(1,j,k)))
               end do
            end do
         end if

         if (topdedge.eq.0) then                       !add d rrl
            call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,2,0) ! code2-0 is added
         else
            call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,2,0) ! code2-0 is added
         end if

         if (rightdedge.eq.1) then                       !add d rrl
            do k=2,n3-2
               do j=jllim,julim
                  v2mp(np,j,k,kv)=v2mp(np+1,j,k,kv)
               end do
            end do
         end if
      end if


      if(idiv.eq.1) then
         illim = 1  + (1-ibcx)*leftdedge                       !add d rrl
         iulim = np + (-1+ibcx)*rightdedge                       !add d rrl
         jllim = 1 + 1*botdedge                       !add d rrl
         julim = mp
         do k=2,n3-2
            do j=jllim,julim
               do i=illim,iulim
                  v2d=-vdiv1(f2(i,j-1,k,kv),
     *                       f2(i,j,k,kv),f2(i,j+1,k,kv),
     *             .5*(h(i,j-1,k)+h(i,j,k)))
     *                -vdiv2(f2(i,j,k,kv),
     *                       f1(i+1,j-1,k,kv),f1(i+1,j,k,kv),
     *                       f1(i,j-1,k,kv),f1(i,j,k,kv),
     *                 .5*(h(i,j-1,k)+h(i,j,k)))
     *                -vdiv2(f2(i,j,k,kv),f3(i,j-1,k+1,kv),
     *                 f3(i,j,k+1,kv),f3(i,j-1,k,kv),
     *                 f3(i,j,k,kv),    .5*(h(i,j-1,k)+h(i,j,k)))
                  v2mp(i,j,k,kv)=v2mp(i,j,k,kv)+(pp(v2d)*x(i,j-1,k,kv)-
     *                                       pn(v2d)*x(i,j,k,kv))
               end do
            end do
         end do
      endif ! idiv.eq.1 

!!! compute antidiffusive velocities in z direction

      illim = 1  + 1*leftdedge                       !add d rrl
      iulim = np - 1*rightdedge                       !add d rrl
      jllim = 1  + 1*botdedge                       !add d rrl
      julim = mp - 1*topdedge                       !add d rrl


      do 53 k=2,n3m
      do 53 j=jllim,julim
      do 53 i=illim,iulim

   53 v3mp(i,j,k,kv)=vdyf(x(i,j,k-1,kv),x(i,j,k,kv),f3(i,j,k,kv),
     *                 .5*(h(i,j,k-1)+h(i,j,k)))
     * +vcorr(f3(i,j,k,kv),
     *        f1(i,j,k-1,kv)+f1(i,j,k,kv)
     *       +f1(i+1,j,k,kv)+f1(i+1,j,k-1,kv),
     *         x(i-1,j,k-1,kv),x(i-1,j,k,kv),
     *         x(i+1,j,k-1,kv),x(i+1,j,k,kv),
     *               .5*(h(i,j,k-1)+h(i,j,k)))
     * +vcorr(f3(i,j,k,kv),
     *        f2(i,j,k-1,kv)+f2(i,j+1,k-1,kv)
     *       +f2(i,j+1,k,kv)+f2(i,j,k,kv),
     *         x(i,j-1,k-1,kv),x(i,j-1,k,kv),
     *         x(i,j+1,k-1,kv),x(i,j+1,k,kv),
     *               .5*(h(i,j,k-1)+h(i,j,k)))
c     if(mpi_rank.eq.1.and.kv.eq.4) then 
c     print*,f3(1,5,4,4),'v3'
c     print*,v3(1,5,4,4),'v3'
c     endif

      if(ibcx.eq.1) then
         jllim = 1  + 1*botdedge                       !add d rrl
         julim = mp - 1*topdedge                       !add d rrl
      if (leftdedge.eq.1) then                       !add d rrl
         do k=2,n3m
            do j=jllim,julim
               v3mp(1,j,k,kv)=vdyf(x(1,j,k-1,kv),
     *                           x(1,j,k,kv),f3(1,j,k,kv),
     *             .5*(h(1,j,k-1)+h(1,j,k)))
     *             +vcorr(f3(1,j,k,kv),
     *                    f1(1,j,k-1,kv)+f1(1,j,k,kv)
     *                   +f1(2,j,k,kv)+f1(2,j,k-1,kv),
     *                     x(-1,j,k-1,kv),x(-1,j,k,kv),
     *                     x(2,j,k-1,kv),x(2,j,k,kv),
     *             .5*(h(1,j,k-1)+h(1,j,k)))
     *             +vcorr(f3(1,j,k,kv),
     *                    f2(1,j,k-1,kv)+f2(1,j+1,k-1,kv)
     *                   +f2(1,j+1,k,kv)+f2(1,j,k,kv),
     *                     x(1,j-1,k-1,kv),x(1,j-1,k,kv),
     *                     x(1,j+1,k-1,kv),x(1,j+1,k,kv),
     *             .5*(h(1,j,k-1)+h(1,j,k)))
            enddo
         enddo
      endif
      call updated(v3mp(1-ih,1-ih,1,kv),v3mp(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,2,0) ! code2-0 is added

      if (rightdedge.eq.1) then                       !add d rrl
         do k=2,n3m
            do j=jllim,julim
               v3mp(np,j,k,kv)=v3mp(np+1,j,k,kv)
c               ja=(mpos-1)*mp + j
c               v3(np,j,k)=vdyf(x(np+1,j,k-1),x(np+1,j,k),f3(np+1,j,k),
c     *             .5*(h(np+1,j,k-1)+h(np+1,j,k)))
c     *             +vcorr(f3(np+1,j,k),
c     *             f1(np+1,j,k-1)+f1(np+1,j,k)+f1(np+2,j,k)+
c     *             f1(np+2,j,k-1),
c     *             x(np-2,j,k-1),x(np-2,j,k),x(np+2,j,k-1),x(np+2,j,k),
c     *             .5*(h(np+1,j,k-1)+h(np+1,j,k)))
c     *             +vcorr(f3(np+1,j,k),
c     *             f2(np+1,j,k-1)+f2(np+1,j+1,k-1)+
c     *             f2(np+1,j+1,k)+f2(np+1,j,k),
c     *             x(np+1,j-1,k-1),x(np+1,j-1,k),x(np+1,j+1,k-1),
c     *             x(np+1,j+1,k),
c     *             .5*(h(np+1,j,k-1)+h(np+1,j,k)))
            enddo
         enddo
      endif
      end if

      if(ibcy.eq.1) then
         illim = 1  + 1*leftdedge                       !add d rrl
         iulim = np - 1*rightdedge                       !add d rrl
         if (botdedge.eq.1) then                       !add d rrl
            do k=2,n3m
               do i=illim,iulim

                  v3mp(i,1,k,kv)=vdyf(x(i,1,k-1,kv),
     *             x(i,1,k,kv),f3(i,1,k,kv),
     *             .5*(h(i,1,k-1)+h(i,1,k)))
     *             +vcorr(f3(i,1,k,kv),
     *                    f1(i,1,k-1,kv)+f1(i,1,k,kv)
     *                   +f1(i+1,1,k,kv)+f1(i+1,1,k-1,kv),
     *                     x(i-1,1,k-1,kv),x(i-1,1,k,kv),
     *                     x(i+1,1,k-1,kv),x(i+1,1,k,kv),
     *             .5*(h(i,1,k-1)+h(i,1,k)))
     *             +vcorr(f3(i,1,k,kv),
     *                    f2(i,1,k-1,kv)+f2(i,2,k-1,kv)
     *                   +f2(i,2,k,kv)+f2(i,1,k,kv),
     *                     x(i,-1,k-1,kv),x(i,-1,k,kv),
     *                     x(i,2,k-1,kv),x(i,2,k,kv),
     *             .5*(h(i,1,k-1)+h(i,1,k)))
               enddo
            enddo
         end if

      call updated(v3mp(1-ih,1-ih,1,kv),v3mp(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,3,0) ! code2-0 is added

         if (topdedge.eq.1) then                       !add d rrl
            do k=2,n3m
               do i=illim,iulim
                  v3mp(i,mp,k,kv)=v3mp(i,mp+1,k,kv)

c                  v3(i,mp,k)=vdyf(x(i,mp+1,k-1),x(i,mp+1,k),
c     *             f3(i,mp+1,k),
c     *             .5*(h(i,mp+1,k-1)+h(i,mp+1,k)))
c     *             +vcorr(f3(i,mp+1,k),
c     *             f1(i,mp+1,k-1)+f1(i,mp+1,k)+f1(i+1,mp+1,k)+
c     *             f1(i+1,mp+1,k-1),
c     *             x(i-1,mp+1,k-1),x(i-1,mp+1,k),x(i+1,mp+1,k-1),
c     *             x(i+1,mp+1,k),
c     *             .5*(h(i,mp+1,k-1)+h(i,mp+1,k)))
c     *             +vcorr(f3(i,mp+1,k),
c     *             f2(i,mp+2,k-1)+f2(i,mp+3,k-1)+f2(i,mp+3,k)+
c     *             f2(i,mp+2,k),
c     *             x(i,mp-1,k-1),x(i,mp-1,k),x(i,mp+2,k-1),x(i,mp+2,k),
c     *             .5*(h(i,mp+1,k-1)+h(i,mp+1,k)))
               enddo
            enddo
         end if
      if(ibcx.eq.1) then
         if (leftdedge.eq.1 .and. botdedge.eq.1) then !add d rrl
            do k=2,n3m
               v3mp(1,1,k,kv)=vdyf(x(1,1,k-1,kv),
     *          x(1,1,k,kv),f3(1,1,k,kv),
     *              .5*(h(1,1,k-1)+h(1,1,k)))
     *              +vcorr(f3(1,1,k,kv),
     *                     f1(1,1,k-1,kv)+f1(1,1,k,kv)
     *                    +f1(2,1,k,kv)+f1(2,1,k-1,kv),
     *                      x(-1,1,k-1,kv),x(-1,1,k,kv),
     *                      x(2,1,k-1,kv),x(2,1,k,kv),
     *              .5*(h(1,1,k-1)+h(1,1,k)))
     *              +vcorr(f3(1,1,k,kv),
     *                     f2(1,1,k-1,kv)+f2(1,2,k-1,kv)
     *                    +f2(1,2,k,kv)+f2(1,1,k,kv),
     *                      x(1,-1,k-1,kv),x(1,-1,k,kv),
     *                      x(1,2,k-1,kv),x(1,2,k,kv),
     *              .5*(h(1,1,k-1)+h(1,1,k)))
            end do
         end if

         call updated(v3mp(1-ih,1-ih,1,kv),v3mp(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1,0) ! code2-0 is added
         if (rightdedge.eq.1 .and. topdedge.eq.1) then  !add d rrl
            do k=2,n3m
               v3mp(np,mp,k,kv)=v3mp(np+1,mp+1,k,kv)
            enddo
         endif
         if (rightdedge.eq.1 .and. botdedge.eq.1) then  !add d rrl
            do k=2,n3m
               v3mp(np,1,k,kv)=v3mp(np+1,1,k,kv)
            enddo
         endif
         if (leftdedge.eq.1 .and. topdedge.eq.1) then   !add d rrl
            do k=2,n3m
               v3mp(1,mp,k,kv)=v3mp(1,mp+1,k,kv)
            enddo
         endif
      end if ! ibcx.eq.1
      end if ! ibcy.eq.1 

      if(idiv.eq.1) then
        illim = 1  + (1-ibcx)*leftdedge                      !add d rrl
        iulim = np + (-1+ibcx)*rightdedge                    !add d rrl
        jllim = 1  + (1-ibcy)*botdedge                       !add d rrl
        julim = mp + (-1+ibcy)*topdedge                      !add d rrl

        do 531 k=2,n3m
        do 531 j=jllim,julim
        do 531 i=illim,iulim
         ia=(npos-1)*np + i
         ja=(mpos-1)*mp + j

         v2d=-vdiv1(f3(i,j,k-1,kv),
     *              f3(i,j,k,kv),f3(i,j,k+1,kv),
     *        .5*(h(i,j,k-1)+h(i,j,k)))
     *        -vdiv2(f3(i,j,k,kv),
     *               f1(i+1,j,k-1,kv),f1(i+1,j,k,kv),
     *               f1(i,j,k-1,kv),
     *               f1(i,j,k,kv),  .5*(h(i,j,k-1)+h(i,j,k)))
     *        -vdiv2(f3(i,j,k,kv),
     *               f2(i,j+1,k-1,kv),f2(i,j+1,k,kv),
     *        f2(i,j,k-1,kv),
     *     f2(i,j,k,kv),    .5*(h(i,j,k-1)+h(i,j,k)))
  531 v3mp(i,j,k,kv)=v3mp(i,j,k,kv)
     *    +(pp(v2d)*x(i,j,k-1,kv)-pn(v2d)*x(i,j,k,kv))

      endif ! idiv==1
! HERE 
      if (rightdedge.eq.0) then                       !add d rrl
         call updated(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0)  ! code2-0 is added
      else
         call updated(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0)  ! code2-0 is added
      end if

       if(ibcx.eq.1) then
          if (leftdedge.eq.1) then                       !add d rrl
             do k=1,n3m
                do j=1,mp
                   v1mp(1,j,k,kv)=v1mp(-1,j,k,kv)
                end do
             end do
          end if
          if (rightdedge.eq.1) then                       !add d rrl
             do k=1,n3m
                do j=1,mp
                   v1mp(np+1,j,k,kv)=v1mp(np+3,j,k,kv)
                enddo
             enddo
          end if
       end if

      if (topdedge.eq.0) then                       !add d rrl
         call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0)  ! code2-0 is added
      else
         call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0)  ! code2-0 is added 
      end if

       if(ibcy.eq.1) then
          if (botdedge.eq.1) then                       !add d rrl
             do k=1,n3m
                do i=1,np
                   v2mp(i,1,k,kv)=v2mp(i,-1,k,kv)
                end do
             end do
          end if
          if (topdedge.eq.1) then                       !add d rrl
             do k=1,n3m
                do i=1,np
                   v2mp(i,mp+1,k,kv)=v2mp(i,mp+3,k,kv)
                enddo
             enddo
          end if
       end if
       enddo ! kv=1:nv

c     if(mpi_rank.eq.1) print*,v3(1,5,4,4),'v3',f3(1,5,4,4)
                  if(nonosold.eq.1) then  !                       non-osscilatory option

      do kv=1,nv
      do 401 k=1,n3m
       km=max0(k-1,1  )
       kp=min0(k+1,n3m)
      do 401 j=1,mp
         if (botdedge.eq.1 .and. j.eq.1) then   !add d rrl
            jm = ibcy*(-1) + ibyo*1
         else
            jm = j - 1
         end if
         if (topdedge.eq.1 .and. j.eq.mp) then  !add d rrl
            jp = ibcy*(mp+2) + ibyo*mp
         else
            jp = j + 1
         end if
      do 401 i=1,np
         if (leftdedge.eq.1 .and. i.eq.1) then !add d rrl
            im = ibcx*(-1) + ibxo*1
         else
            im = i - 1
         end if
         if (rightdedge.eq.1 .and. i.eq.np) then !add d rrl
            ip = ibcx*(np+2) + ibxo*np
         else
            ip = i + 1
         end if
         mxo(i,j,k,kv)=amax1(x(im,j,k,kv),x(i,j,k,kv),
     .                       x(ip,j,k,kv),mxo(i,j,k,kv),
     .                       x(i,jm,k,kv),x(i,jp,k,kv),
     .                       x(i,j,kp,kv),x(i,j,km,kv))
 401     mno(i,j,k,kv)=amin1(x(im,j,k,kv),x(i,j,k,kv),
     .                       x(ip,j,k,kv),mno(i,j,k,kv),
     .                       x(i,jm,k,kv),x(i,jp,k,kv),
     .                       x(i,j,kp,kv),x(i,j,km,kv))

      illim = 1
      iulim = np + 1*rightdedge                       !add d rrl
      jllim = 1
      julim = mp

      do 402 k=1,n3m
      do 402 j=jllim,julim
      do 402 i=illim,iulim
  402 f1(i,j,k,kv)=donor(c2,c2,v1mp(i,j,k,kv))

      illim = 1
      iulim = np
      jllim = 1
      julim = mp + 1*topdedge                       !add d rrl

      do 403 k=1,n3m
      do 403 j=jllim,julim
      do 403 i=illim,iulim
  403 f2(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))

      do 4033 k=1,n3
      do 4033 j=1,mp
      do 4033 i=1,np
 4033 f3(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))

      if (rightdedge.eq.0) then                       !add d rrl
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
      else
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
      end if

      if (topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
      end if

      do 444 k=1,n3m
      do 444 j=1,mp
      do 444 i=1,np
      cpmp(i,j,k,kv)=(mxo(i,j,k,kv)-x(i,j,k,kv))*h(i,j,k)/
     1( pn(f1(i+1,j,k,kv))+pp(f1(i,j,k,kv))
     2 +pn(f2(i,j+1,k,kv))+pp(f2(i,j,k,kv))
     3 +pn(f3(i,j,k+1,kv))+pp(f3(i,j,k,kv))+ep)

      cnmp(i,j,k,kv)=(x(i,j,k,kv)-mno(i,j,k,kv))*h(i,j,k)/
     1( pp(f1(i+1,j,k,kv))+pn(f1(i,j,k,kv))
     2 +pp(f2(i,j+1,k,kv))+pn(f2(i,j,k,kv))
     3 +pp(f3(i,j,k+1,kv))+pn(f3(i,j,k,kv))+ep)
 444  continue

      call updated(cpmp(1-ih,1-ih,1,kv),cpmp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0,0) ! code2-0 is added
      call updated(cnmp(1-ih,1-ih,1,kv),cnmp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0,0) ! code2-0 is added

      illim = 1 + 1*leftdedge                       !add d rrl
      do k=1,n3m
        do j=1,mp
          do i=illim,np
            v1mp(i,j,k,kv)= pp(v1mp(i,j,k,kv))
     *                *amin1(1.,cpmp(i,j,k,kv),cnmp(i-1,j,k,kv))
     *                 -pn(v1mp(i,j,k,kv))
     *                *amin1(1.,cpmp(i-1,j,k,kv),cnmp(i,j,k,kv))
          end do
        end do
      end do
      if (ibcx.eq.1) then
         if (rightdedge.eq.0) then                       !add d rrl
            call updated(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0) ! code2-0 is added
         else
            call updated(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0) ! code2-0 is added
         end if
         do k=1,n3m
            do j=1,mp
               if (leftdedge.eq.1) then                       !add d rrl
                  v1mp(1 ,j,k,kv)=v1mp(-1,j,k,kv)
               end if
               if (rightdedge.eq.1) then                       !add d rrl
                  v1mp(np+1,j,k,kv)=v1mp(np+3  ,j,k,kv)
               end if
            end do
         end do
      end if

      jllim = 1 + 1*botdedge                       !add d rrl
      do k=1,n3m
        do j=jllim,mp
          do i=1,np
            v2mp(i,j,k,kv)= pp(v2mp(i,j,k,kv))
     .                *amin1(1.,cpmp(i,j,k,kv),cnmp(i,j-1,k,kv))
     .                -pn(v2mp(i,j,k,kv))
     .                *amin1(1.,cpmp(i,j-1,k,kv),cnmp(i,j,k,kv))
          end do
        end do
      end do

      if (ibcy.eq.1) then
         if (topdedge.eq.0) then                       !add d rrl
            call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1,0) ! code2-0 is added
         else
            call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1,0) ! code2-0 is added
         end if

        do k=1,n3m
          do i=1,np
             if (botdedge.eq.1) then                       !add d rrl
                v2mp(i, 1,k,kv)=v2mp(i,-1,k,kv)
             end if
             if (topdedge.eq.1) then                       !add d rrl
                v2mp(i,mp+1,k,kv)=v2mp(i,mp+3,k,kv)
             end if
          end do
        end do
      end if

      do k=2,n3m
        do j=1,mp
          do i=1,np
            v3mp(i,j,k,kv)= pp(v3mp(i,j,k,kv))
     .             *amin1(1.,cpmp(i,j,k,kv),cnmp(i,j,k-1,kv))
     *                -pn(v3mp(i,j,k,kv))
     .             *amin1(1.,cpmp(i,j,k-1,kv),cnmp(i,j,k,kv))
          end do
        end do
      end do

      enddo ! kv=1:nv

                  endif  !! nonosold.eq.1 

      if(nonos.eq.1) then
       do 1000 itrfct=1,nfct !!nfct=1 is default

       if(itrfct.eq.1) then
        do kv=1,nv

         illim = 1
         iulim = np + 1*rightdedge                       !add d rrl
         jllim = 1
         julim = mp
         do 502 k=1,l
         do 502 j=jllim,julim
         do 502 i=illim,iulim
  502 f1(i,j,k,kv)=donor(c2,c2,v1mp(i,j,k,kv))

         illim = 1
         iulim = np
         jllim = 1
         julim = mp + 1*topdedge                       !add d rrl
         do 503 k=1,n3m
         do 503 j=jllim,julim
         do 503 i=illim,iulim
  503 f2(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))

         do 5033 k=1,n3
         do 5033 j=1,mp
         do 5033 i=1,np
 5033 f3(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))

         if (rightdedge.eq.0) then                       !add d rrl
          call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
         else
          call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
         end if

         if (topdedge.eq.0) then                       !add d rrl
          call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
         else
          call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
         end if
        enddo !kv=1:nv
      endif ! itrfct.eq.1

!!! correction coefficients for variables nv-1
      do kv=1,nv-1
      do k=1,l
      do j=1,mp
      do i=1,np
       if(abs(mx(i,j,k,kv)).lt.ep) mx(i,j,k,kv)=0.
       if(abs(mn(i,j,k,kv)).lt.ep) mn(i,j,k,kv)=0.
       if(k.eq.1) f3(i,j,k,kv)=-f3(i,j,k+1,kv)
       if(k.eq.1) f3(i,j,k,nv)=-f3(i,j,k+1,nv)
       if(k.eq.l) f3(i,j,k+1,kv)=-f3(i,j,k,kv)
       if(k.eq.l) f3(i,j,k+1,nv)=-f3(i,j,k,nv)
       rhoin=  (pn(f1(i+1,j,k,nv))+pp(f1(i,j,k,nv))
     .        +pn(f2(i,j+1,k,nv))+pp(f2(i,j,k,nv))
     .        +pn(f3(i,j,k+1,nv))+pp(f3(i,j,k,nv)))
       rhoout=-(pp(f1(i+1,j,k,nv))+pn(f1(i,j,k,nv))
     .        +pp(f2(i,j+1,k,nv))+pn(f2(i,j,k,nv))
     .        +pp(f3(i,j,k+1,nv))+pn(f3(i,j,k,nv)))
       ain=    (pn(f1(i+1,j,k,kv))+pp(f1(i,j,k,kv))
     .        +pn(f2(i,j+1,k,kv))+pp(f2(i,j,k,kv))
     .        +pn(f3(i,j,k+1,kv))+pp(f3(i,j,k,kv)))
       aout=  -(pp(f1(i+1,j,k,kv))+pn(f1(i,j,k,kv))
     .        +pp(f2(i,j+1,k,kv))+pn(f2(i,j,k,kv))
     .        +pp(f3(i,j,k+1,kv))+pn(f3(i,j,k,kv)))
          cpmp(i,j,k,kv)=
     . pp(mx(i,j,k,kv)*x(i,j,k,nv)-x(i,j,k,kv))*h(i,j,k)
     .       /(ain-pp(mx(i,j,k,kv))*rhoout
     .            +pn(mx(i,j,k,kv))*rhoin+ep)
         cnmp(i,j,k,kv)=
     . pp(x(i,j,k,kv)-mn(i,j,k,kv)*x(i,j,k,nv))*h(i,j,k)
     .     /(-aout+pp(mn(i,j,k,kv))*rhoin
     .            -pn(mn(i,j,k,kv))*rhoout+ep)
      enddo ! i=1:np
      enddo ! j=1:mp
      enddo ! k=1:l

      call updated(cpmp(1-ih,1-ih,1,kv),cpmp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0,0) ! code2-0 is added
      call updated(cnmp(1-ih,1-ih,1,kv),cnmp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0,0) ! code2-0 is added
      enddo ! kv=1:nv-1

!!! correction coeffiecients for rho
      do kv=1,nv-1
      do k=1,l
      do j=1,mp
      do i=1,np
c  JLW - Fix from Jon Reisner 9/2005
c     a1p=cpmp(i,j,k,kv)+amax1(0.,sign(1., mx(i,j,k,kv)))
c     a2p=cnmp(i,j,k,kv)+amax1(0.,sign(1.,-mn(i,j,k,kv)))
c     a1n=cnmp(i,j,k,kv)+amax1(0.,sign(1., mn(i,j,k,kv)))
c     a2n=cpmp(i,j,k,kv)+amax1(0.,sign(1.,-mx(i,j,k,kv)))
       rmxuse=-1.*mx(i,j,k,kv)
       if(mx(i,j,k,kv).eq.0.) rmxuse=0.
       rmnuse=-1.*mn(i,j,k,kv)
       if(mn(i,j,k,kv).eq.0.) rmnuse=0.
       a1p=cpmp(i,j,k,kv)+amax1(0.,sign(1., mx(i,j,k,kv)))
       a2p=cnmp(i,j,k,kv)+amax1(0.,sign(1.,rmnuse))
       a1n=cnmp(i,j,k,kv)+amax1(0.,sign(1., mn(i,j,k,kv)))
       a2n=cpmp(i,j,k,kv)+amax1(0.,sign(1.,rmxuse))
       tmpp=amin1(a1p,a2p)
       tmpn=amin1(a1n,a2n)
       if(kv.eq.1) cpmp(i,j,k,nv)=tmpp
       if(kv.eq.1) cnmp(i,j,k,nv)=tmpn
       cpmp(i,j,k,nv)=amin1(tmpp,cpmp(i,j,k,nv))
       cnmp(i,j,k,nv)=amin1(tmpn,cnmp(i,j,k,nv))
      enddo !i=1:np
      enddo !j=1:mp
      enddo !k=1:l
      enddo !kv=1:nv-1 

      call updated(cpmp(1-ih,1-ih,1,nv),cpmp(1-ih,1-ih,1,nv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0,0) ! code2-0 is added
      call updated(cnmp(1-ih,1-ih,1,nv),cnmp(1-ih,1-ih,1,nv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0,0) ! code2-0 is added

c     if(mpi_rank.eq.1) then
c     print*,v1(1,5,3,4),'v1'
c     print*,v1(2,5,3,4),'v1'
c     print*,v2(1,5,3,4),'v2'
c     print*,v2(1,6,3,4),'v2'
c     print*,v3(1,5,3,4),'v3'
c     print*,v3(1,5,4,4),'v3'
c     print*, pp(v1(1,5,3,4))
c    *        *amin1(1.,cpmp(1,5,3,4),cnmp(0,5,3,4))
c    *                 -pn(v1(1,5,3,4))
c    *        *amin1(1.,cpmp(0,5,3,4),cnmp(1,5,3,4))
c     endif

      do kv=1,nv
      do k=1,l
        do j=1,mp
          do i=1+leftdedge,np                        !add d rrl
          v1mp(i,j,k,kv)= pp(v1mp(i,j,k,kv))
     *        *amin1(1.,cpmp(i,j,k,kv),cnmp(i-1,j,k,kv))
     *                 -pn(v1mp(i,j,k,kv))
     *        *amin1(1.,cpmp(i-1,j,k,kv),cnmp(i,j,k,kv))
          end do
        end do
      end do

      if (ibcx.eq.1) then
         if (rightdedge.eq.0) then                       !add d rrl
            call updated(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0) ! code2-0 is added
         else
            call updated(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1,0) ! code2-0 is added 
         end if
         do k=1,n3m
            do j=1,mp
               if (leftdedge.eq.1) then                       !add d rrl
                  v1mp(1 ,j,k,kv)=v1mp(-1,j,k,kv)
               end if                       !add d rrl
               if (rightdedge.eq.1) then                       !add d rrl
                  v1mp(np+1,j,k,kv)=v1mp(np+3  ,j,k,kv)
               end if
            end do
         end do
      end if

      do k=1,l
        do j=1+botdedge,mp                        !add d rrl
          do i=1,np
            v2mp(i,j,k,kv)= pp(v2mp(i,j,k,kv))
     *             *amin1(1.,cpmp(i,j,k,kv),cnmp(i,j-1,k,kv))
     *                   -pn(v2mp(i,j,k,kv))
     *             *amin1(1.,cpmp(i,j-1,k,kv),cnmp(i,j,k,kv))
          end do
        end do
      end do

      if (ibcy.eq.1) then
         if (topdedge.eq.0) then                       !add d rrl
            call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1,0) ! code2-0 is added
         else
            call updated(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1,0) ! code2-0 is added
         end if

        do k=1,n3m
          do i=1,np
             if (botdedge.eq.1) then                       !add d rrl
                v2mp(i, 1,k,kv)=v2mp(i,-1,k,kv)
             end if
             if (topdedge.eq.1) then                       !add d rrl
                v2mp(i,mp+1,k,kv)=v2mp(i,mp+3,k,kv)
             end if
          end do
        end do
      end if

      do k=2,l
        do j=1,mp
          do i=1,np
            v3mp(i,j,k,kv)= pp(v3mp(i,j,k,kv))
     *             *amin1(1.,cpmp(i,j,k,kv),cnmp(i,j,k-1,kv))
     *                   -pn(v3mp(i,j,k,kv))
     *             *amin1(1.,cpmp(i,j,k-1,kv),cnmp(i,j,k,kv))
          end do
        end do
      end do
      enddo
c     if(mpi_rank.eq.1) then
c     print*,v1(1,5,3,4),'v1'
c     print*,v1(2,5,3,4),'v1'
c     print*,cn(1,5,3,4),'cn'
c     print*,cn(0,5,3,4),'cn'
c     print*,cp(1,5,3,4),'cp'
c     print*,cp(0,5,3,4),'cp'
c     print*,v2(1,5,3,4),'v2'
c     print*,v2(1,6,3,4),'v2'
c     print*,v3(1,5,3,4),'v3'
c     print*,v3(1,5,4,4),'v3'
c     endif

      if(itrfct.lt.nfct) then
       do kv=1,nv
       do 602 k=1,l
       do 602 j=1,mp
       do 602 i=1,np+1
       tmp=f1(i,j,k,kv)
       f1o(i,j,k,kv)=donor(c2,c2,v1mp(i,j,k,kv))
  602 f1(i,j,k,kv)=tmp-f1o(i,j,k,kv)

       do 603 k=1,l
       do 603 j=1,mp+1
       do 603 i=1,np
       tmp=f2(i,j,k,kv)
       f2o(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))
  603 f2(i,j,k,kv)=tmp-f2o(i,j,k,kv)

       do 6033 k=1,n3
       do 6033 j=1,mp
       do 6033 i=1,np
       tmp=f3(i,j,k,kv)
       f3o(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))
 6033 f3(i,j,k,kv)=tmp-f3o(i,j,k,kv)

       if(rightdedge.eq.0) then                       !add d rrl
        call updated(f1o(1-ih,1-ih,1,kv),f1o(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
       else
        call updated(f1o(1-ih,1-ih,1,kv),f1o(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2,0) ! code2-0 is added
       endif
 
       if(botdedge.eq.0) then                       !add d rrl
        call updated(f2o(1-ih,1-ih,1,kv),f2o(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
       else
        call updated(f2o(1-ih,1-ih,1,kv),f2o(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3,0) ! code2-0 is added
       endif

       do k=1,l
       do j=1,mp
       do i=1,np
      x(i,j,k,kv)=x(i,j,k,kv)-( f1o(i+1,j,k,kv)-f1o(i,j,k,kv)
     .                         +f2o(i,j+1,k,kv)-f2o(i,j,k,kv)
     .                         +f3o(i,j,k+1,kv)-f3o(i,j,k,kv) )
     .                          /h(i,j,k)
       enddo !i
       enddo !j
       enddo !k 

       enddo !kv=1:nv
      endif ! itrfct < nfct

1000  continue ! itrfct=1:nfct
      endif ! nonos.eq.1 

   30                      continue   ! loop from LINE 268 - do 30 itr=1,iord  
    6 continue                        ! goto from LINE 455 - if(itr.eq.irod) goto 6

      do kv=1,nv
       call updated(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1,0) ! code2-0 is added
      enddo

      deallocate (v1mp) ! code2 - was v1
      deallocate (v2mp) ! code2 - was v2
      deallocate (v3mp) ! code2 - was v3
      deallocate (f1)
      deallocate (f2)
      deallocate (f3)
      deallocate (f1o)
      deallocate (f2o)
      deallocate (f3o)
      deallocate (cpmp) ! code2 - was cp 
      deallocate (cnmp) ! code2 - was cn 
      deallocate (mxo)
      deallocate (mno)
      deallocate (mx)
      deallocate (mn)
      deallocate (a)
2500  continue

      return
      end subroutine mpdataold

!*****************************************************************************************!
      
!****************************************************************************************
! high order advection scheme for large time step, based on mpdata -
! this is FP's version, which doesn't have nonos, nonosold option
!****************************************************************************************
      subroutine mpdata(xv,xe,h,il,iu,jl,ju,lls,isW) ! xe deleted from 2nd argument
      use gridsetup
      use advo
      use msga
      Implicit None
      integer,intent(in) :: il,iu,jl,ju,lls
      real,dimension(il:iu,jl:ju,lls) :: xv,xe
      real,dimension(il:iu,jl:ju,lls) :: h
      integer :: i,j,k, kp1,km1    !,k0
      real :: ha, fluxIn, fluxOut
      real :: eps,vda
      integer:: isW  !if isW=1 bc is xv(0)=-xv(1), otherwhise xv(0)=xv(1)
      real:: plusOrMinus1

      eps = 1e-15

      call updated(xv,xe,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
! compute donor cell      
      do k=1,l
       do j=1,mp
        do i=1,np+1
         fd1(i,j,k)=donor(xv(i-1,j,k),xv(i,j,k),u1(i,j,k))
        enddo
       enddo
      enddo
      do k=1,l
        do j=1,mp+1
         do i=1,np
          fd2(i,j,k)=donor(xv(i,j-1,k),xv(i,j,k),u2(i,j,k))
         enddo
        enddo
      enddo
      do k=2,l
       do j=1,mp
        do i=1,np
         fd3(i,j,k)=donor(xv(i,j,k-1),xv(i,j,k),u3(i,j,k))
        enddo
       enddo
      enddo
      if(ibctopbot.eq.0) then
       stop
      else
       do j=1,mp
        do i=1,np
          fd3(i,j, 1)=0.
          fd3(i,j,l+1)=0.
        enddo
       enddo
      endif

!     MPDATA method
      if (iord.eq.2) then
! compute pmx,pmn
      do k=1,l
       km1 = max(k-1,1)
       kp1 = min(k+1,l)
       do j=1,mp
         do i=1,np
           pmx(i,j,k)=max(xv(i,j,k),xv(i-1,j,k),
     +         xv(i+1,j,k),xv(i,j-1,k),xv(i,j+1,k),
     +         xv(i,j,km1),xv(i,j,kp1))
           pmn(i,j,k)=min(xv(i,j,k),xv(i-1,j,k),
     +         xv(i+1,j,k),xv(i,j-1,k),xv(i,j+1,k),
     +         xv(i,j,km1),xv(i,j,kp1))
         enddo
       enddo
      enddo
! compute temporary values based on donor cell

      do k=1,l
       do j=1,mp
        do i=1,np
         xv(i,j,k)=xv(i,j,k)-(
     &                 fd1(i+1,j,k)-fd1(i,j,k)
     &                +fd2(i,j+1,k)-fd2(i,j,k)
     &                +fd3(i,j,k+1)-fd3(i,j,k))/h(i,j,k)
        enddo
       enddo
      enddo 

      call updated(xv,xe,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)

! update of pmx and pmn with new xv from donor cell
      do k=1,l
        km1 = max(k-1,1)
        kp1 = min(k+1,l)
        do j=1,mp
         do i=1,np
           pmx(i,j,k)=max(xv(i,j,k),xv(i-1,j,k),
     +         xv(i+1,j,k),xv(i,j-1,k),xv(i,j+1,k),
     +         xv(i,j,km1),xv(i,j,kp1),pmx(i,j,k))
           pmn(i,j,k)=min(xv(i,j,k),xv(i-1,j,k),
     +         xv(i+1,j,k),xv(i,j-1,k),xv(i,j+1,k),
     +         xv(i,j,km1),xv(i,j,kp1),pmn(i,j,k))
         enddo
        enddo
      enddo
! computation of antidiffusive velocities     
      ! the updated of u1, u2, and u3 is done in setAdvectiveVelocities
      !call updated(u1,u1,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1,1)
      !call updated(u2,u2,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1,2)
      !call updated(u3,u3,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1,0)
      ! x velocities
      do k=1,l    
       if (k.eq.1) then
        km1 = 1
        plusOrMinus1 = - 2* isW + 1  !-1 when isW=1, 1 otherwhise
       else
        km1=k-1
        plusOrMinus1 = 1
       endif
       ! FIXME: xv(0) is not xv(1) when topo (see getSubGroundvalue in
       ! field update)
       kp1 = min(k+1,l)
       do j=1,mp
        do i=1,np+1
         ha = 0.5*(h(i-1,j,k)+h(i,j,k))
         v1(i,j,k)=vdyf(xv(i-1,j,k),xv(i,j,k),u1(i,j,k),ha)
     + -vcorr_fp(u1(i,j,k),u2(i-1,j,k)+u2(i-1,j+1,k)+u2(i,j+1,k)+u2(i,j,k), ! vcorr is changed to vcorr_fp to avoid confusion
     +        xv(i-1,j-1,k),xv(i,j-1,k),xv(i-1,j+1,k),xv(i,j+1,k),ha)
     + -vcorr_fp(u1(i,j,k),u3(i-1,j,k)+u3(i-1,j,k+1)+u3(i,j,k+1)+u3(i,j,k),
     +      plusOrMinus1*xv(i-1,j,km1),plusOrMinus1*xv(i,j,km1),
     +         xv(i-1,j,kp1),xv(i,j,kp1),ha)
        enddo
       enddo
      enddo
      if (idiv.eq.1) then
      do k=1,l
       kp1 = min(k+1,l)
       do j=1,mp
        do i=1,np+1
         ha = 0.5*(h(i-1,j,k)+h(i,j,k))
         vda = -vdiv1(u1(i-1,j,k),u1(i,j,k),u1(i+1,j,k),ha)
     +-vdiv2(u1(i,j,k),u2(i-1,j+1,k),u2(i,j+1,k),u2(i-1,j,k),u2(i,j,k),ha)
     +-vdiv2(u1(i,j,k),u3(i-1,j,kp1),u3(i,j,kp1),u3(i-1,j,k),u3(i,j,k),ha)
         v1(i,j,k)=v1(i,j,k)+donor(xv(i-1,j,k),xv(i,j,k),vda)
        enddo
       enddo
      enddo
      endif
      ! y velocities
      do k=1,l    
       if (k.eq.1) then
        km1 = 1
        plusOrMinus1 = - 2* isW + 1  !-1 when isW=1, 1 otherwhise
       else
        km1=k-1
        plusOrMinus1 = 1
       endif
       !km1 = max(k-1,1)
       kp1 = min(k+1,l)
       do j=1,mp+1
        do i=1,np
         ha = 0.5*(h(i,j-1,k)+h(i,j,k))
         v2(i,j,k)=vdyf(xv(i,j-1,k),xv(i,j,k),u2(i,j,k),ha)
     + -vcorr_fp(u2(i,j,k),u1(i,j-1,k)+u1(i+1,j-1,k)+u1(i+1,j,k)+u1(i,j,k),
     +        xv(i-1,j-1,k),xv(i-1,j,k),xv(i+1,j-1,k),xv(i+1,j,k),ha)
     + -vcorr_fp(u2(i,j,k),u3(i,j-1,k)+u3(i,j,k)+u3(i,j,k+1)+u3(i,j-1,k+1),
     +      plusOrMinus1 * xv(i,j-1,km1),plusOrMinus1 * xv(i,j,km1)
     +        ,xv(i,j-1,kp1),xv(i,j,kp1),ha)
        enddo
       enddo
      enddo
      if (idiv.eq.1) then
       do k=1,l
       kp1 = min(k+1,l)
       do j=1,mp+1
        do i=1,np
         ha = 0.5*(h(i,j-1,k)+h(i,j,k))
         vda = -vdiv1(u2(i,j-1,k),u2(i,j,k),u2(i,j+1,k),ha)
     +-vdiv2(u2(i,j,k),u1(i+1,j-1,k),u1(i+1,j,k),u1(i,j-1,k),u1(i,j,k),ha)
     +-vdiv2(u2(i,j,k),u3(i,j-1,kp1),u3(i,j,kp1),u3(i,j-1,k),u3(i,j,k),ha)
         v2(i,j,k)=v2(i,j,k)+donor(xv(i,j-1,k),xv(i,j,k),vda)
        enddo
       enddo
      enddo
      endif
      ! z velocities
      do k=2,l    
       do j=1,mp
        do i=1,np
         ha = 0.5*(h(i,j,k-1)+h(i,j,k))
         v3(i,j,k)=vdyf(xv(i,j,k-1),xv(i,j,k),u3(i,j,k),ha)
     + -vcorr_fp(u3(i,j,k),u1(i,j,k-1)+u1(i,j,k)+u1(i+1,j,k)+u1(i+1,j,k-1),
     +     xv(i-1,j,k-1),xv(i-1,j,k),xv(i+1,j,k-1),xv(i+1,j,k),ha)
     + -vcorr_fp(u3(i,j,k),u2(i,j,k-1)+u2(i,j+1,k-1)+u2(i,j+1,k)+u2(i,j,k),
     +     xv(i,j-1,k-1),xv(i,j-1,k),xv(i,j+1,k-1),xv(i,j+1,k),ha)
        enddo
       enddo
      enddo
      if (idiv.eq.1) then
      do k=2,l
       do j=1,mp
        do i=1,np
         ha = 0.5*(h(i,j,k-1)+h(i,j,k))
         vda = -vdiv1(u3(i,j,k-1),u3(i,j,k),u3(i,j,k+1),ha)
     +-vdiv2(u3(i,j,k),u2(i,j+1,k-1),u2(i,j+1,k),u2(i,j,k-1),u2(i,j,k),ha)
     +-vdiv2(u3(i,j,k),u1(i+1,j,k-1),u1(i+1,j,k),u1(i,j,k-1),u1(i,j,k),ha)
         v3(i,j,k)=v3(i,j,k)+donor(xv(i,j,k-1),xv(i,j,k),vda)
        enddo
       enddo
      enddo
      endif
      !v3 in k=1 and k=l+1 is never used
      !do j=1,mp
      !do i=1, np
      !  v3(i,j,1)=0.0
      !  v3(i,j,l+1) =0.0
      !enddo
      !enddo

      ! computation of temporary fluxes for oscilation limiter
      do k=1,l
       do j=1,mp
        do i=1,np+1
         fd1(i,j,k)=donor1(v1(i,j,k))
        enddo
       enddo
      enddo
      do k=1,l
        do j=1,mp+1
         do i=1,np
          fd2(i,j,k)=donor1(v2(i,j,k))
         enddo
        enddo
      enddo
      do k=2,l
       do j=1,mp
        do i=1,np
         fd3(i,j,k)=donor1(v3(i,j,k))
        enddo
       enddo
      enddo
      ! limiter coefficients
      do k=1,l
       do j=1,mp
        do i=1,np
         fluxIn=pp(fd1(i,j,k)) + pn(fd1(i+1,j,k))
     +          + pp(fd2(i,j,k)) + pn(fd2(i,j+1,k))
     +          + pp(fd3(i,j,k)) + pn(fd3(i,j,k+1))
         fluxOut=pp(fd1(i+1,j,k)) + pn(fd1(i,j,k))
     +         + pp(fd2(i,j+1,k)) + pn(fd2(i,j,k))
     +         + pp(fd3(i,j,k+1)) + pn(fd3(i,j,k))
         cp(i,j,k) = (pmx(i,j,k)-xv(i,j,k))*h(i,j,k)/
     +                (fluxIn+eps)
         cn(i,j,k) = (xv(i,j,k)-pmn(i,j,k))*h(i,j,k)/
     +                (fluxOut+eps)
        enddo
       enddo
      enddo
      call updated(cp,cp,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
      call updated(cn,cn,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
      do k=1,l
       do j=1,mp
        do i=1,np+1
         v1(i,j,k)=donor(min(1.0,cp(i,j,k),cn(i-1,j,k)),
     .        min(1.0,cp(i-1,j,k),cn(i,j,k)),v1(i,j,k))
        enddo
       enddo
      enddo
      do k=1,l
        do j=1,mp+1
         do i=1,np
         v2(i,j,k)=donor(min(1.0,cp(i,j,k),cn(i,j-1,k)),
     .       min(1.0,cp(i,j-1,k),cn(i,j,k)),v2(i,j,k))
         enddo
        enddo
      enddo
      do k=2,l
       do j=1,mp
        do i=1,np
         v3(i,j,k)=donor(min(1.0,cp(i,j,k),cn(i,j,k-1)),
     .      min(1.0,cp(i,j,k-1),cn(i,j,k)),v3(i,j,k))
        enddo
       enddo
      enddo

! donor cell on final antidiffusive velocities
      do k=1,l
       do j=1,mp
        do i=1,np+1
         fd1(i,j,k)=donor1(v1(i,j,k))
        enddo
       enddo
      enddo
      do k=1,l
        do j=1,mp+1
         do i=1,np
          fd2(i,j,k)=donor1(v2(i,j,k))
         enddo
        enddo
      enddo
      do k=2,l
       do j=1,mp
        do i=1,np
         fd3(i,j,k)=donor1(v3(i,j,k))
        enddo
       enddo
      enddo
      endif !iord.eq.2, end mpdata method

! update of xv      
      do k=1,l
       do j=1,mp
        do i=1,np
         xv(i,j,k)=xv(i,j,k)-(
     &                 fd1(i+1,j,k)-fd1(i,j,k)
     &                +fd2(i,j+1,k)-fd2(i,j,k)
     &                +fd3(i,j,k+1)-fd3(i,j,k))/h(i,j,k)
        enddo
       enddo
      enddo 

      return
      end subroutine mpdata

!*******************************************************************************************!
! from u, v, o divided by gi and cell centered, this routine set u1,u2,
! u3 (module advo), the contravariant face coordinates time dt/dx * 1/gi required for
! advective schemes
! when ihighorder.eq.2, u1,u2 and u3 are updated (for antidiffusive vel
! in mpdata...
!*******************************************************************************************!

      subroutine setAdvectiveVelocities(u,v,o,gcx,gcy,gcz,il,iu,jl,ju,lls,ihighorder)
      use metryic
      use gridsetup
      use msga
      use advo
      use xve

! u (uavg from advec) is uavg_code, u1 is uab in code.
! v (vavg from advec) is vavg_code, u2 is vab in code.
! o (oavg from advec) is ovag_code, u3 is oab in code   

      Implicit None

      integer,intent(in) :: il,iu,jl,ju,lls, ihighorder
      integer :: i,j,k
      real :: gcx,gcy,gcz
      real u(il:iu,jl:ju,lls),v(il:iu,jl:ju,lls),o(il:iu,jl:ju,lls)
 
      ! compute contravariant cell face in the x direction
      call updated(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,2,0)
      do k=1,l
       do j=1,mp
        do i=1,np+1
         u1(i,j,k)=0.5*(u(i,j,k)+u(i-1,j,k))*gcx           
        enddo
       enddo
      enddo
      if(leftdedge.eq.1) then
       do k=1,l
        do j=1,mp
          ! 0 gradient between xe, cause xe(0 is not defined)
          u1(1,j,k)=xe(1,j,k,1)/xe(1,j,k,nv)*gcx/gi(i,j,k)
        enddo
       enddo
      endif
      if(rightdedge.eq.1) then
       if (ibclatopen.eq.0) then
         do k=1,l
           do j=1,mp
             u1(np+1,j,k)=xe(np,j,k,1)/xe(np,j,k,nv)*gcx/gi(i,j,k)
             ! here we assume 0 gradient between xe
           enddo
         enddo
       else if (ibclatopen.eq.1) then
          do k=1,l
           do j=1,mp
             u1(np+1,j,k)=max(u1(np,j,k),0.0)
           enddo
         enddo
       else if (ibclatopen.eq.2) then
          do k=1,l
           do j=1,mp
             u1(np+1,j,k)=u1(np,j,k)
           enddo
         enddo
       else
          write(6,*) 'wrong ubclapoption'
          stop
       endif 
      endif


      ! compute contravariant cell face in the y direction
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,3,0)
      do k=1,l
        do j=1,mp+1
         do i=1,np
          u2(i,j,k)=0.5*(v(i,j,k)+v(i,j-1,k))*gcy
         enddo
        enddo
       enddo
       if(botdedge.eq.1)then
        if (ibclatopen.eq.0) then
        do k=1,l
         do i=1,np
           ! here we assume 0 gradient between xe
           u2(i,1,k)=xe(i,1,k,2)/xe(i,1,k,nv)*gcy/gi(i,j,k)
         enddo
        enddo
        else if (ibclatopen.eq.1) then
        do k=1,l
         do i=1,np
           u2(i,1,k)=min(0.0,u2(i,2,k))
         enddo
        enddo
        else if (ibclatopen.eq.2) then
        do k=1,l
         do i=1,np
           u2(i,1,k)=u2(i,2,k)
         enddo
        enddo
        else 
          write(6,*) 'wrong ubclapoption'
          stop
        endif
        
       endif
       if(topdedge.eq.1)then
        if (ibclatopen.eq.0) then
        do k=1,l
         do i=1,np
           ! here we assume 0 gradient between xe
          u2(i,mp+1,k)=xe(i,mp,k,2)/xe(i,mp,k,nv)*gcy/gi(i,j,k)
         enddo
        enddo
        else if (ibclatopen.eq.1) then
        do k=1,l
         do i=1,np
           u2(i,mp+1,k)=max(0.0,u2(i,mp,k))
         enddo
        enddo
        else if (ibclatopen.eq.2) then
        do k=1,l
         do i=1,np
           u2(i,mp+1,k)=u2(i,mp,k)
         enddo
        enddo
        else 
          write(6,*) 'wrong ubclapoption'
          stop
        endif
       endif

      ! compute contravariant cell face in the z direction
      do k=2,l
       do j=1,mp
        do i=1,np
         u3(i,j,k)=0.5*(o(i,j,k)+o(i,j,k-1))*gcz
        enddo
       enddo
      enddo
      if(ibctopbot.eq.0) then
       do j=1,mp
        do i=1,np
         !u3(i,j,1)=xe(i,j,1,3)/xe(i,mp,1,5)*gcz/gi(i,j,k)
         !u3(i,j,l+1)=xe(i,j,l,3)/xe(i,mp,l,5)*gcz/gi(i,j,k)
         !FIXME : - or +? probably +
         u3(i,j,1)=-u3(i,j,2)
         u3(i,j,l+1)=-u3(i,j,l)
         !u3(i,j,1)=u3(i,j,2)
         !u3(i,j,l+1)=u3(i,j,l)
        enddo
       enddo
      else
       do j=1,mp
        do i=1,np
         u3(i,j,1)=0.0
         u3(i,j,l+1)=0.0
        enddo
       enddo
      endif !if(ibctopbot.eq.0)

      ! the following updated are required for antidiffusive vel in
      ! mpdata
      if (ihighorder.eq.2) then
        call updated(u1,u1,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1,1)
        call updated(u2,u2,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1,2)
        call updated(u3,u3,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1,0)
      endif

      end subroutine setAdvectiveVelocities
      
!*****************************************************************************************!
      real function pp(y)
      Implicit None
      real :: y

      pp = amax1(0.0,y)

      end function pp
      !*************************************
      real function pn(y)
      Implicit None
      real :: y

      pn = -amin1(0.0,y)

      end function pn
      !*************************************
      real function donor(y1,y2,a10)
      Implicit None
      real :: y1,y2,a10
      ! donor=y1 for positive velocity

      donor = pp(a10)*y1-pn(a10)*y2

      end function donor
      !*************************************
      real function donor1(a10)
      Implicit None
      real :: a10
      ! donor=1 for positive velocity

      donor1 = pp(a10)-pn(a10)

      end function donor1
      !*************************************
      real function vdiv1(a1,a2,a3,r)
      Implicit None
      real :: a1,a2,a3,r

      vdiv1 = 0.25*a2*(a3-a1)/r

      end function vdiv1
      !*************************************
      real function vdiv2(aa,b1,b2,b3,b4,r)
      Implicit None
      real :: aa,b1,b2,b3,b4,r

      vdiv2 = 0.25*aa*(b1+b2-b3-b4)/r

      end function vdiv2
      !*************************************
      real function rat2(z1,z2)
      Implicit None
      real :: z1,z2

      rat2 = (z2-z1)*0.5

      end function rat2
      !*************************************
      real function rat4(z0,z1,z2,z3)
      Implicit None
      real :: z0,z1,z2,z3

      rat4 = (z3+z2-z1-z0)*0.25

      end function rat4
      !*************************************
      !*************************************
      real function vdyf(x1,x2,aa,r)
      Implicit None
      real :: x1,x2,aa,r

      vdyf = (abs(aa)-aa**2/r)*rat2(x1,x2)

      end function vdyf
      !*************************************
      real function vcorr(aa,b,y0,y1,y2,y3,r)
      Implicit None
      real :: aa,b,y0,y1,y2,y3,r

      vcorr = -0.125*aa*b/r*rat4(y0,y1,y2,y3)

      end function vcorr
      !*************************************
      !*************************************
      real function vcorr_fp(aa,b,y0,y1,y2,y3,r)
      Implicit None
      real :: aa,b,y0,y1,y2,y3,r

      vcorr_fp = 0.125*aa*b/r*rat4(y0,y1,y2,y3)

      end function vcorr_fp
!*****************************************************************************************!
      end module advectionScheme
!*****************************************************************************************!
