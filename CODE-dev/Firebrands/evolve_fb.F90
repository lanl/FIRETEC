!-----------------------------------------------------------------------
! evolveFB consumes and moves firebrands
!-----------------------------------------------------------------------
subroutine evolveFB(u,v,w,rhog,tempg,ilh,iuh,jlh,juh,klh,lls)
  use gridlist_variables, only : wp,n,m,l,dx,dy,dz,ibc_left,ibc_right, &
    ibc_top,ibc_bot,dt
  use firebrand_gridlist_variables, only : shapeFB,dsizeFB
  use firebrand_general_variables, only : Cd_dsk,Cd_cyl,Cd_sph, &
    numBrands,xvbrand,nvbrand,ix,iy,iz,iuvel,ivvel,iwvel,irad,ihght,irho
  use gridsetup, only : np,mp
  use msga_variables, only : npos,mpos
  use metric_variables, only : zs,zedge
  use constants, only : pi,g,tolerance
  use firebrand_general_variables, only : xvbrand,numBrands,nvbrand
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,klh,lls
  real(wp),intent(in) :: u(ilh:iuh,jlh:juh,klh:lls)
  real(wp),intent(in) :: v(ilh:iuh,jlh:juh,klh:lls)
  real(wp),intent(in) :: w(ilh:iuh,jlh:juh,klh:lls)
  real(wp),intent(in) :: rhog(ilh:iuh,jlh:juh,lls)
  real(wp),intent(in) :: tempg(ilh:iuh,jlh:juh,lls)

  integer :: i,j,k,it
  integer :: numBrand,count
  real(wp) :: xvu,xvv,xvw
  real(wp) :: w_bx,w_by,w_bz,w_abs=0._wp,w_absz!,w_bh
  real(wp) :: nu,fo,Ac=0._wp
  real(wp) :: A_ek,massFB
  real(wp),external :: zcart
  real(wp),allocatable :: xvtmpbrand(:,:)
  logical :: retain(numBrands)

  ! Executable Code
  retain = .true.
  numBrand=numBrands
  do it=1,numBrands
    ! Find grid coordinates of firebrand
    i=floor(xvbrand(it,ix)/dx)-(npos-1)*np+1
    j=floor(xvbrand(it,iy)/dy)-(mpos-1)*mp+1
    do k=1,lls
      if(zcart(zedge(k+1),i,j).gt.xvbrand(it,iz)) exit
    enddo
    call GetWind(xvbrand(it,ix),xvbrand(it,iy), &
      xvbrand(it,iz),i,j,k,u,v,w,ilh,iuh,jlh,juh,klh,lls,xvu,xvv,xvw)
    
    ! Velocity differences between particle and field
    w_bx=(xvbrand(it,iuvel)-xvu)**2
    w_by=(xvbrand(it,ivvel)-xvv)**2
    w_bz=(xvbrand(it,iwvel)-xvw)**2
    w_absz=sqrt(w_bx+w_by+w_bz)
    !w_bh  =max(0.00001,sqrt(w_bx+w_by))

    ! Consume firebrands through burning
    nu=(1.7616E-5_wp*(tempg(i,j,k)/273.16_wp)**1.5_wp*(273.16_wp+110.4_wp)/ &
      (tempg(i,j,k)+110.4_wp))/rhog(i,j,k) ! Kinematic air viscosity (m2/s)
    if(dsizeFB.eq.1)then
      if(shapeFB.eq.1)then ! Emmons' model - disk dr/dt
        fo=0.52055_wp
        Ac=pi*xvbrand(it,irad)**2*Cd_dsk
        w_abs=sqrt(2*g*xvbrand(it,ihght)*xvbrand(it,irho)/ &
          (rhog(i,j,k)*Cd_dsk))
        xvbrand(it,irad)=xvbrand(it,irad)-fo*dt* &
          rhog(i,j,k)/xvbrand(it,irho)*sqrt(w_abs*nu/xvbrand(it,ihght))
      elseif(shapeFB.eq.2)then ! Emmons' model - cylinder dr/dt
        fo=0.52055_wp
        Ac=2*xvbrand(it,irad)*xvbrand(it,ihght)*Cd_cyl
        w_abs=sqrt(pi*g*xvbrand(it,irad)*xvbrand(it,irho)/ &
          (rhog(i,j,k)*Cd_cyl))
        xvbrand(it,irad)=xvbrand(it,irad)-fo*dt* &
          rhog(i,j,k)/xvbrand(it,irho)* &
          sqrt(w_abs*nu/(xvbrand(it,irad)*pi))
      elseif(shapeFB.eq.3)then ! droplet model - sphere dr/dt
        w_abs=sqrt(w_bx+w_by+w_bz)
        Ac=pi*xvbrand(it,irad)**2*Cd_sph
        xvbrand(it,irad)=xvbrand(it,irad)-dt*nu* &
          rhog(i,j,k)/xvbrand(it,irho)* &
          log(2.2_wp)/xvbrand(it,irad)
      endif 
    elseif(dsizeFB.eq.2)then
      if(shapeFB.eq.1)then ! Dr.Woychesse's model - disk dh/dt
        fo=0.39924_wp
        Ac=pi*xvbrand(it,irad)**2*Cd_dsk
        w_abs=sqrt(2*g*xvbrand(it,ihght)*xvbrand(it,irho)/ &
          (rhog(i,j,k)*Cd_dsk))
        xvbrand(it,ihght)=xvbrand(it,ihght)-8.0_wp/3.0_wp*fo*dt* &
          rhog(i,j,k)/xvbrand(it,irho)* &
          sqrt(w_abs*nu/xvbrand(it,irad))
      elseif(shapeFB.eq.2)then ! Emmons' model - cylinder dh/dt
        fo=0.52055_wp
        Ac=2*xvbrand(it,ihght)*xvbrand(it,irad)*Cd_cyl
        w_abs=sqrt(pi*g*xvbrand(it,irad)*xvbrand(it,irho)/ &
          (rhog(i,j,k)*Cd_cyl))
        xvbrand(it,ihght)=xvbrand(it,ihght)-1.113_wp*fo*dt* &
          rhog(i,j,k)/xvbrand(it,irho)*sqrt(w_abs*nu/xvbrand(it,irad))
      endif 
    endif
    
    ! Identify burned out brands
    if(xvbrand(it,ihght).lt.tolerance.or.xvbrand(it,irad).lt.tolerance &
      .or.xvbrand(it,irho).lt.tolerance)then ! Burnout
      retain(it)=.false.
      numBrand=numBrand-1
      cycle
    endif

    ! Determine new change in firebrand velocity
    select case(shapeFB)
      case(1,2)
        massFB=xvbrand(it,irho)*pi*xvbrand(it,irad)**2*xvbrand(it,ihght)
      case(3)
        massFB=xvbrand(it,irho)*4.0_wp/3.0_wp*pi*xvbrand(it,irad)**3
      case default
        massFB=tolerance
    end select
    A_ek=0.5_wp*Ac*rhog(i,j,k)/massFB

    xvbrand(it,iuvel)=xvbrand(it,iuvel)+ &
      (xvu-xvbrand(it,iuvel))*A_ek*w_abs*dt!*(sin(ak)-cos(ak)*w_bz/w_bh)
    xvbrand(it,ivvel)=xvbrand(it,ivvel)+ &
      (xvv-xvbrand(it,ivvel))*A_ek*w_abs*dt!*(sin(ak)-cos(ak)*w_bz/w_bh)
    xvbrand(it,iwvel)=xvbrand(it,iwvel)-g*dt+ &
      (xvw-xvbrand(it,iwvel))*A_ek*w_absz*dt!*(sin(ak)-cos(ak)*w_bz/w_bh)
   
    ! Move firebrands with given velocities
    xvbrand(it,ix)=xvbrand(it,ix)+xvbrand(it,iuvel)*dt
    xvbrand(it,iy)=xvbrand(it,iy)+xvbrand(it,ivvel)*dt
    xvbrand(it,iz)=xvbrand(it,iz)+xvbrand(it,iwvel)*dt
    
    ! Identify brands no longer in flight
    if(xvbrand(it,iz).le.zs(i,j).or.xvbrand(it,iz).ge.l*dz)then ! Landed on ground or lofted beyond vertical domain
      retain(it)=.false.
      numBrand=numBrand-1
      cycle
    endif

    ! Check left and right boundaries
    if(xvbrand(it,ix).lt.0)then ! Outside x domain on left
      if(ibc_left.eq.0)then
        xvbrand(it,ix)=n*dx+xvbrand(it,ix)
      elseif(ibc_left.eq.5) then
        xvbrand(it,ix)= -xvbrand(it,ix)
        xvbrand(it,iuvel)=-xvbrand(it,iuvel)
      else
        retain(it)=.false.
        numBrand=numBrand-1
        cycle
      endif
    else if(xvbrand(it,ix).gt.(n*dx))then ! Outside x domain on right
      if(ibc_right.eq.0)then
        xvbrand(it,ix)=xvbrand(it,ix)-n*dx
      elseif(ibc_right.eq.5) then
        xvbrand(it,ix)=n*dx-(xvbrand(it,ix)-(n*dx))
        xvbrand(it,iuvel)=-xvbrand(it,iuvel)
      else
        retain(it)=.false.
        numBrand=numBrand-1
        cycle
      endif
    endif
    
    ! Check top and bottom boundaries
    if(xvbrand(it,iy).lt.0)then ! Outside y domain on bottom
      if(ibc_bot.eq.0)then
        xvbrand(it,iy)=m*dy+xvbrand(it,iy)
      elseif(ibc_bot.eq.5) then
        xvbrand(it,iy)=-xvbrand(it,iy)
        xvbrand(it,ivvel)=-xvbrand(it,ivvel)
      else
        retain(it)=.false.
        numBrand=numBrand-1
        cycle
      endif
    else if(xvbrand(it,iy).gt.m*dy)then ! Outside y domain on top
      if(ibc_top.eq.0)then
        xvbrand(it,iy)=xvbrand(it,iy)-m*dy
      elseif(ibc_top.eq.5) then
        xvbrand(it,iy)=m*dy-(xvbrand(it,iy)-(m*dy))
        xvbrand(it,ivvel)=-xvbrand(it,ivvel)
      else
        retain(it)=.false.
        numBrand=numBrand-1
        cycle
      endif
    endif
  enddo

  ! Remove irrelevant brands
  allocate(xvtmpbrand(numBrand,nvbrand))
  count=1
  do it=1,numBrands
    if(retain(it))then
      xvtmpbrand(count,:)=xvbrand(it,:)
      count=count+1
    endif
  enddo
  call move_alloc(xvtmpbrand,xvbrand)
  numBrands=numBrand

  call updateFB

end subroutine evolveFB

!-----------------------------------------------------------------------
! updateFB updates the location and distribution of firebrands across
! processors
!-----------------------------------------------------------------------
subroutine updateFB
  use gridlist_variables, only : wp,dx,dy,nprocx
  use gridsetup, only : np,mp
  use firebrand_general_variables, only : xvbrand,numBrands,nvbrand, &
    ix,iy
  use msga_variables, only : numprocs,mpi_integer,mpi_max, &
    mpi_comm_world,ierror
  use msga_variables, only : mpi_wp
  use mpi
  Implicit None

  ! Local Variables
  integer :: it,iproc
  integer :: npos,mpos,mpi_rank_fb
  integer :: mxbrand,nsize
  integer :: sndBrands(numprocs),rcvBrands(numprocs)
  real(wp),allocatable :: xvsndBrands(:,:,:),xvrcvBrands(:,:,:)

  ! Executable Code
  ! Count and communicate the number of brands sent and received 
  ! between all the processors
  sndBrands=0
  do it=1,numBrands
    npos=floor(xvbrand(it,ix)/dx/np)
    mpos=floor(xvbrand(it,iy)/dy/mp)
    mpi_rank_fb=mpos*nprocx+npos
    sndBrands(mpi_rank_fb+1)=sndBrands(mpi_rank_fb+1)+1
  enddo
  call mpi_alltoall(sndBrands,1,mpi_integer,rcvBrands,1,mpi_integer, &
    mpi_comm_world,ierror)
  call mpi_allreduce(maxval(sndBrands),mxbrand,1,mpi_integer, &
    mpi_max,mpi_comm_world,ierror)

  ! Communicate brand data between processors
  allocate(xvsndBrands(mxbrand,nvbrand,numprocs))
  sndBrands=0
  do it=1,numBrands
    npos=floor(xvbrand(it,ix)/dx/np)
    mpos=floor(xvbrand(it,iy)/dy/mp)
    mpi_rank_fb=mpos*nprocx+npos
    sndBrands(mpi_rank_fb+1)=sndBrands(mpi_rank_fb+1)+1
    xvsndBrands(sndBrands(mpi_rank_fb+1),:,mpi_rank_fb+1)= &
      xvbrand(it,:)
  enddo
  allocate(xvrcvBrands(mxbrand,nvbrand,numprocs))
  nsize=mxbrand*nvbrand
  call mpi_alltoall(xvsndBrands,nsize,mpi_wp,xvrcvBrands,nsize,mpi_wp, &
    mpi_comm_world,ierror)
  deallocate(xvsndBrands)

  ! Reallocate brands for moved brands
  deallocate(xvbrand)
  numBrands=sum(rcvBrands)
  allocate(xvbrand(numBrands,nvbrand))
  it=1
  do iproc=1,numprocs
    if(rcvBrands(iproc).gt.0) then
      xvbrand(it:it+rcvBrands(iproc)-1,:)= &
        xvrcvBrands(1:rcvBrands(iproc),:,iproc)
      it=it+rcvBrands(iproc)
    endif
  enddo
  deallocate(xvrcvBrands)

end subroutine updateFB
