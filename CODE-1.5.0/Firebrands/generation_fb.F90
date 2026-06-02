!-----------------------------------------------------------------------
! launchFB generates new firebrands
!-----------------------------------------------------------------------
subroutine launchFB(u,v,w,rhog,ilh,iuh,jlh,juh,klh,lls, &
    rhof,temps,il,iu,jl,ju,lfuel,nfuel)
  use gridlist_variables, only : prec,dx,dy
  use firebrand_gridlist_variables, only : tempHot,shapeFB, &
    radius,height,rhoFB
  use firebrand_general_variables, only : Cd_dsk,Cd_cyl,Cd_sph, &
    numBrands,xvbrand,nvbrand,ix,iy,iz,iuvel,ivvel,iwvel,irad,ihght, &
    itemp,irho,itheta,iphi
  use gridsetup, only : np,mp
  use msga_variables, only : npos,mpos
  use metric_variables, only : zedge
  use zcart_function, only : zcart
  use constants, only : pi,g
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,klh,lls
  integer,intent(in) :: il,iu,jl,ju,lfuel,nfuel
  real(prec),intent(in) :: u(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: v(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: w(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: rhog(il:iu,jl:ju,lls)
  real(prec),intent(in) :: rhof(nfuel,il:iu,jl:ju,lfuel)
  real(prec),intent(in) :: temps(nfuel,il:iu,jl:ju,lfuel)

  integer :: i,j,k,ia,ja,it
  integer :: newBrands
  logical :: newCheck(il:iu,jl:ju,lfuel)
  real(prec) :: force_drag,force_grav
  real(prec) :: randi,randj,randk
  real(prec) :: tempBrand(numBrands,nvbrand)
  real(prec),allocatable :: newBrand(:,:)

  ! Executable Code
  newBrands=0
  newCheck=.false.
  do k=1,lfuel
    do j=jl,ju
      do i=il,iu
        ! Check criteria for consideration for launching
        if(any(temps(:,i,j,k).gt.tempHot) &
          .and.any(rhof(:,i,j,k).gt.0.00001) &
          .and.(w(i,j,k).gt.0)) then

          ! Calculate drag and gravitational forces ! TODO Oriented with wind
          select case(shapeFB)
            case(1) ! DISK
              force_drag=0.5*Cd_dsk*(2*radius*height) &
                *rhog(i,j,k)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
              force_grav=g*rhoFB*(pi*radius**2*height)
            case(2) ! CYLINDER
              force_drag=0.5*Cd_cyl*(2*pi*radius**2) &
                *rhog(i,j,k)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
              force_grav=rhoFB*g*(pi*radius**2*height)
            case(3) ! SPHERE
              force_drag=0.5*Cd_sph*(2*pi*radius**2) &
                *rhog(i,j,k)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
              force_grav=rhoFB*g*4./3.*pi*radius**3.
          end select
          
          ! Launch Firebrand (TODO Assumes a specific height and/or radius)
          if(force_drag.gt.force_grav)then
            newBrands=newBrands+1
            newCheck(i,j,k)=.true.
          endif
        endif
      enddo
    enddo
  enddo

  allocate(newBrand(newBrands,nvbrand))
  it=1
  do k=1,lfuel
    do j=jl,ju
      do i=il,iu
        if(newCheck(i,j,k))then
          ! Cartesian coordinates of new firebrand
          ia=(npos-1)*np+i
          ja=(mpos-1)*mp+j
          call random_number(randi)
          newBrand(it,ix)=(ia-1)*dx+randi*dx
          call random_number(randj)
          newBrand(it,iy)=(ja-1)*dy+randj*dy
          call random_number(randk)
          newBrand(it,iz)=zcart(zedge(k),i,j)+ &
            randk*(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))

          ! Velocities of new firebrand (averaged between no movement and wind velocity field)
          call GetWind(newBrand(it,ix),newBrand(it,iy), &
            newBrand(it,iz),i,j,k,u,v,w,ilh,iuh,jlh,juh,klh,lls, &
            newBrand(it,iuvel),newBrand(it,ivvel),newBrand(it,iwvel))

          ! Characteristic size of new firebrand
          newBrand(it,irad) = radius
          newBrand(it,ihght) = height

          ! Temperature, density, and angle of alignment of new firebrand
          newBrand(it,itemp) = maxval(temps(:,i,j,k))
          newBrand(it,irho) = rhoFB
          newBrand(it,itheta)= 0. ! TODO use this
          newBrand(it,iphi)= 0. ! TODO use this
          it=it+1
        endif
      enddo
    enddo
  enddo

  tempBrand=xvbrand
  deallocate(xvbrand)
  allocate(xvbrand(numBrands+newBrands,nvbrand))
  xvbrand(1:numBrands,:)=tempBrand
  xvbrand(numBrands+1:numBrands+newBrands,:)=newBrand
  numBrands=numBrands+newBrands
  deallocate(newBrand)

end subroutine launchFB 
