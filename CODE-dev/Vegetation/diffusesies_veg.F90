!-----------------------------------------------------------------------
! diffusesies spreads creeping fire laterally through a diffusion 
! of temperature technique (modified from original sies-based diffusion)
!-----------------------------------------------------------------------
subroutine diffusesies(xvrhoFuel,xvrhoWater,frhof,fsies, &
    il,iu,jl,ju,lls,nfuel,ih)
  use gridlist_variables, only : wp
  use gridsetup, only : dxi,dyi
  use fuel_variables, only : hf,temps
  use metric_variables, only : gi
  Implicit None

  ! Input/Output Arguments
  integer,intent(in) :: il,iu,jl,ju,lls,nfuel,ih
  real(wp),intent(inout) :: frhof(nfuel,il:iu,jl:ju,lls)
  real(wp),intent(inout) :: fsies(nfuel,il:iu,jl:ju,lls)
  real(wp),intent(in) :: xvrhoFuel(nfuel,il:iu,jl:ju,lls)
  real(wp),intent(in) :: xvrhoWater(nfuel,il:iu,jl:ju,lls)
  
  ! Local Variables
  integer :: ift,i,j,k
  real(wp) :: coef,rhos
  real(wp),allocatable :: tempsdiff(:,:,:)
  real(wp),allocatable :: hxc(:,:,:),hyc(:,:,:)

  ! Executable Code
  allocate(tempsdiff(il:iu,jl:ju,lls)); tempsdiff=0._wp

  allocate(hxc(il:iu+1,jl:ju,lls)); hxc=0._wp
  allocate(hyc(il:iu,jl:ju+1,lls)); hyc=0._wp
   
  do ift=1,nfuel
    call update(temps(ift,:,:,:),iu,ju,lls, &
      il-ih,iu+ih,jl-ih,ju+ih,1,0,0)

    ! compute contravariant component x-flux at (i-1/2,j,k) hxc (sqrt(g) burried in...
    ! top bc of phi is hardcoded kp1=l when k=l
    do k=1,lls
      do j=jl,ju
        do i=il,iu+1
          if(temps(ift,i,j,k).LT.800) then                                               ! The orginal version of the code checked for when the gas temperature was less than 800 instead of solid temp
            coef=1.25_wp/dxi                                                                ! W/(m*K) thermal conductivity by matching the heat flux from convection to the heat flux from conductivity in this modeling approach we get the coefficient by rearranging h(dT) = (d/dx)(k*dT)
          else
            coef=0.0_wp
          endif 
          hxc(i,j,k)=coef*(temps(ift,i,j,k)-temps(ift,i-1,j,k))*dxi                   ! units [J*s^-1*m^-2]
        enddo
      enddo
    enddo

    ! compute contravariant component y-flux at (i,j-1/2,k) hyc (sqrt(g) burried in...)
    ! top bc of phi is hardcoded kp1=l when k=l
    do k=1,lls
      do j=jl,ju+1
        do i=il,iu  
          if(temps(ift,i,j,k).LT.800) then
            coef=1.25_wp/dyi                                                              ! W/(m*K) thermal conductivity ... the value of 1.25 was chosen so that it is effectively the same coefficient (0.001) as what was used in the original diffusesies which was diffusing sies
          else
            coef=0.0_wp
          endif
          hyc(i,j,k)=coef*(temps(ift,i,j,k)-temps(ift,i,j-1,k))*dyi                   ! units [J*s^-1*m^-2]
        enddo
      enddo
    enddo
      
    ! compute Laplacian term
    do k=1,lls
      do j=jl,ju
        do i=il,iu
          tempsdiff(i,j,k)=tempsdiff(i,j,k) &                                         ! units [J*s^-1*m^-3]
            +max(0._wp,gi(i,j,k)*((hxc(i+1,j,k)-hxc(i,j,k))*dxi &                        ! originally there was a random 2 in front of max that was removed and combined into the coef term
            +(hyc(i,j+1,k)-hyc(i,j,k))*dyi))
        enddo
      enddo
    enddo
  enddo ! ift

  do k=1,lls
    do j=jl,ju
      do i=il,iu
        do ift=1,nfuel
          if(xvrhoFuel(ift,i,j,k).gt.1.e-6_wp)then
            rhos=xvrhoFuel(ift,i,j,k)+xvrhoWater(ift,i,j,k)
            frhof(ift,i,j,k)=frhof(ift,i,j,k)-tempsdiff(i,j,k)/hf                     ! tempsdiff/hf units ... [J*s^-1*m^-3][J-1*kg].     = [kg*m^-3*s^-1]
            fsies(ift,i,j,k)=fsies(ift,i,j,k)+tempsdiff(i,j,k)/rhos                   ! tempsdiff/hf units ... [J*s^-1*m^-3][m^3*kg^-1]   = [J*kg^-1*s^-1]
          endif
        enddo
      enddo
    enddo
  enddo

  deallocate(hxc,hyc,tempsdiff)
      
end subroutine diffusesies
