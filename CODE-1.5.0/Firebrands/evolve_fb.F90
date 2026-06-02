!-----------------------------------------------------------------------
! evolveFB consumes and moves firebrands
!-----------------------------------------------------------------------
subroutine evolveFB(u,v,w,rhog,tempg,ilh,iuh,jlh,juh,klh, &
    il,iu,jl,ju,lls)
  use gridlist_variables, only : prec,n,m,l,dx,dy,dz,ibcx,ibcy
  use firebrand_gridlist_variables, only : shapeFB,dsizeFB
  use firebrand_general_variables, only : Cd_dsk,Cd_cyl,Cd_sph, &
    numBrands,xvbrand,nvbrand,ix,iy,iz,iuvel,ivvel,iwvel,irad,ihght,irho
  use gridsetup, only : np,mp,dt
  use msga_variables, only : npos,mpos
  use metric_variables, only : zedge
  use zcart_function, only : zcart
  use constants, only : pi,g
  use firebrand_general_variables, only : xvbrand,numBrands,nvbrand
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,klh
  integer,intent(in) :: il,iu,jl,ju,lls
  real(prec),intent(in) :: u(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: v(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: w(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: rhog(il:iu,jl:ju,lls)
  real(prec),intent(in) :: tempg(il:iu,jl:ju,lls)

  integer :: i,j,k,it
  integer :: numBrand,count
  real(prec) :: xvu,xvv,xvw
  real(prec) :: w_bx,w_by,w_bz,w_abs,w_absz!,w_bh
  real(prec) :: nu,fo,Ac
  real(prec) :: A_ek
  real(prec),allocatable :: xvtmpbrand(:,:)
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
    nu=(1.7616E-5*(tempg(i,j,k)/273.16)**1.5*(273.16+110.4)/ &
      (tempg(i,j,k)+110.4))/rhog(i,j,k) ! Kinematic air viscosity (m2/s)
    if(dsizeFB.eq.1)then
      if(shapeFB.eq.1)then ! Emmons' model - disk dr/dt
        fo=0.52055
        Ac=pi*xvbrand(it,irad)**2*Cd_dsk
        w_abs=sqrt(2*g*xvbrand(it,ihght)*xvbrand(it,irho)/ &
          (rhog(i,j,k)*Cd_dsk))
        xvbrand(it,irad)=xvbrand(it,irad)-fo*dt* &
          rhog(i,j,k)/xvbrand(it,irho)*sqrt(w_abs*nu/xvbrand(it,ihght))
      elseif(shapeFB.eq.2)then ! Emmons' model - cylinder dr/dt
        fo=0.52055
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
          log(2.2)/xvbrand(it,irad)
      endif 
    elseif(dsizeFB.eq.2)then
      if(shapeFB.eq.1)then ! Dr.Woychesse's model - disk dh/dt
        fo=0.39924
        Ac=pi*xvbrand(it,irad)**2*Cd_sph
        w_abs=sqrt(2*g*xvbrand(it,ihght)*xvbrand(it,irho)/ &
          (rhog(i,j,k)*Cd_dsk))
        xvbrand(it,ihght)=xvbrand(it,ihght)-8./3.*fo*dt* &
          rhog(i,j,k)/xvbrand(it,irho)* &
          sqrt(w_abs*nu/xvbrand(it,irad))
      elseif(shapeFB.eq.2)then ! Emmons' model - cylinder dh/dt
        fo=0.52055
        Ac=2*xvbrand(it,ihght)*xvbrand(it,irad)*Cd_cyl
        w_abs=sqrt(pi*g*xvbrand(it,irad)*xvbrand(it,irho)/ &
          (rhog(i,j,k)*Cd_cyl))
        xvbrand(it,ihght)=xvbrand(it,ihght)-1.113*fo*dt* &
          rhog(i,j,k)/xvbrand(it,irho)*sqrt(w_abs*nu/xvbrand(it,irad))
      endif 
    endif

    ! Determine new change in firebrand velocity
    A_ek=0.5*Ac*rhog(i,j,k)

    xvbrand(it,iuvel)=xvbrand(it,iuvel)+xvu*A_ek*w_abs*dt!*(sin(ak)-cos(ak)*w_bz/w_bh)
    xvbrand(it,ivvel)=xvbrand(it,ivvel)+xvv*A_ek*w_abs*dt!*(sin(ak)-cos(ak)*w_bz/w_bh)
    xvbrand(it,iwvel)=xvbrand(it,iwvel)-g*dt+xvw*A_ek*w_absz*dt!*(sin(ak)-cos(ak)*w_bz/w_bh)

    ! Move firebrands with given velocities
    xvbrand(it,ix)=xvbrand(it,ix)+xvbrand(it,iuvel)*dt
    xvbrand(it,iy)=xvbrand(it,iy)+xvbrand(it,ivvel)*dt
    xvbrand(it,iz)=xvbrand(it,iz)+xvbrand(it,iwvel)*dt
    
    ! Identify irrelevant brands
    if(xvbrand(it,ihght).lt.0.000001.or.xvbrand(it,irad).lt.0.000001 &
      .or.xvbrand(it,irho).lt.0.0001)then ! Burnout
      retain(it)=.false.
      numBrand=numBrand-1
      cycle
    else if(xvbrand(it,iz).le.0.or.xvbrand(it,iz).ge.l*dz)then ! Landed on ground or lofted beyond vertical domain
      retain(it)=.false.
      numBrand=numBrand-1
    else if(xvbrand(it,ix).le.0)then ! Outside x domain
      if(ibcx.eq.0)then
        xvbrand(it,ix)=n*dx+xvbrand(it,ix)
      else
        retain(it)=.false.
        numBrand=numBrand-1
      endif
    else if(xvbrand(it,ix).ge.(n*dx))then ! Outside x domain
      if(ibcx.eq.0)then
        xvbrand(it,ix)=xvbrand(it,ix)-n*dx
      else
        retain(it)=.false.
        numBrand=numBrand-1
      endif
    else if(xvbrand(it,iy).le.0)then ! Outside y domain
      if(ibcy.eq.0)then
        xvbrand(it,iy)=m*dy+xvbrand(it,iy)
      else
        retain(it)=.false.
        numBrand=numBrand-1
      endif
    else if(xvbrand(it,iy).ge.m*dy)then ! Outside y domain
      if(ibcy.eq.0)then
        xvbrand(it,iy)=xvbrand(it,iy)-m*dy
      else
        retain(it)=.false.
        numBrand=numBrand-1
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
  deallocate(xvbrand)
  allocate(xvbrand(numBrand,nvbrand))
  numBrands=numBrand
  xvbrand=xvtmpbrand
  deallocate(xvtmpbrand)

  call updateFB

end subroutine evolveFB

!-----------------------------------------------------------------------
! updateFB updates the location and distribution of firebrands across
! processors
!-----------------------------------------------------------------------
subroutine updateFB
  use gridlist_variables, only : prec,dx,dy,nprocx
  use gridsetup, only : np,mp
  use firebrand_general_variables, only : xvbrand,numBrands,nvbrand, &
    ix,iy
  use msga_variables, only : numprocs,mpi_integer,mpi_max, &
    mpi_comm_world,ierror
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  Implicit None

  ! Local Variables
  integer :: it,iproc
  integer :: npos,mpos,mpi_rank_fb
  integer :: mxbrand,nsize
  integer :: sndBrands(numprocs),rcvBrands(numprocs) 
  real(prec),allocatable :: xvsndBrands(:,:,:),xvrcvBrands(:,:,:)

  ! Executable Code
  ! Count and communicate the number of brands sent and received 
  ! between all the processors
  sndBrands=0
  do it=1,numBrands
    npos=floor(xvbrand(it,ix)/dx/np)
    mpos=floor(xvbrand(it,iy)/dy/mp)
    mpi_rank_fb=mpos*nprocx+npos
    if(mpi_rank_fb.gt.3) &
      print*,'ERROR0',npos,mpos,xvbrand(it,ix),xvbrand(it,iy)
    sndBrands(mpi_rank_fb+1)=sndBrands(mpi_rank_fb+1)+1
  enddo
  do iproc=1,numprocs
    call mpi_scatter(sndBrands,1,mpi_integer,rcvBrands(iproc), &
      1,mpi_integer,iproc-1,mpi_comm_world,ierror)
  enddo
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
  do iproc=1,numprocs
#ifdef DBL_PREC
    call mpi_scatter(xvsndBrands,nsize,mpi_double_precision, &
      xvrcvBrands(:,:,iproc),nsize,mpi_double_precision,iproc-1, &
      mpi_comm_world,ierror)
#else
    call mpi_scatter(xvsndBrands,nsize,mpi_real,xvrcvBrands(:,:,iproc),&
      nsize,mpi_real,iproc-1,mpi_comm_world,ierror)
#endif
  enddo
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
