!-----------------------------------------------------------------------
! diffusesies spreads creeping fire laterally through a diffusion 
! of ignitions technique. Not complete diffusion in that no 
! energy is moved between cells, but the effect of diffusion is 
! computed and then positive heat spread is treated like an 
! ignition consuming fuel to supply heat in adjacent cells, hence
! maintaining the mass and energy balance
!-----------------------------------------------------------------------
subroutine diffusesies(xvsies,xvrhoFuel,xvrhoWater,frhof,fsies,tempg, &
    il,iu,jl,ju,lls,nfuel,ih)
  use gridlist_variables, only : prec
  use gridsetup, only : dxi,dyi
  use fuel_variables, only : hf
  use metric_variables, only : gi
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nfuel,ih
  real(prec),intent(inout) :: frhof(nfuel,il:iu,jl:ju,lls)
  real(prec),intent(inout) :: fsies(nfuel,il:iu,jl:ju,lls)
  real(prec),intent(in) :: xvsies(nfuel,il:iu,jl:ju,lls)
  real(prec),intent(in) :: xvrhoFuel(nfuel,il:iu,jl:ju,lls)
  real(prec),intent(in) :: xvrhoWater(nfuel,il:iu,jl:ju,lls)
  real(prec),intent(in) :: tempg(il-ih:iu+ih,jl-ih:ju+ih,lls)
  
  integer :: ift,i,j,k
  real(prec) :: coef,rhos
  real(prec),allocatable :: sies(:,:,:,:),siesdiff(:,:,:)
  real(prec),allocatable :: hxc(:,:,:),hyc(:,:,:)

  ! Executable Code
  allocate(sies(nfuel,il-ih:iu+ih,jl-ih:ju+ih,lls)); sies=0.0
  sies(:,il:iu,jl:ju,:)=xvsies
  allocate(siesdiff(il:iu,jl:ju,lls)); siesdiff=0.

  allocate(hxc(il:iu+1,jl:ju,lls)); hxc=0.0
  allocate(hyc(il:iu,jl:ju+1,lls)); hyc=0.0
   
  do ift=1,nfuel
    call update(sies(ift,:,:,:),sies(ift,:,:,:),iu,ju,lls, &
      il-ih,iu+ih,jl-ih,ju+ih,1,0)

    ! compute contravariant component x-flux at (i-1/2,j,k) hxc (sqrt(g) burried in...
    ! top bc of phi is hardcoded kp1=l when k=l
    do k=1,lls
      do j=jl,ju
        do i=il,iu+1
          if(tempg(i,j,k).LT.800) then 
            coef=0.001
          else
            coef=0.0
          endif 
          hxc(i,j,k)=coef*(sies(ift,i,j,k)-sies(ift,i-1,j,k))*dxi
        enddo
      enddo
    enddo

    ! compute contravariant component y-flux at (i,j-1/2,k) hyc (sqrt(g) burried in...)
    ! top bc of phi is hardcoded kp1=l when k=l
    do k=1,lls
      do j=jl,ju+1
        do i=il,iu  
          if(tempg(i,j,k).LT.800) then
            coef=0.001
          else
            coef=0.0
          endif
          hyc(i,j,k)=coef*(sies(ift,i,j,k)-sies(ift,i,j-1,k))*dyi
        enddo
      enddo
    enddo
      
    ! compute Laplacian term
    do k=1,lls
      do j=jl,ju
        do i=il,iu
          siesdiff(i,j,k)=siesdiff(i,j,k) &
            +2.*max(0.,gi(i,j,k)*((hxc(i+1,j,k)-hxc(i,j,k))*dxi &
            +(hyc(i,j+1,k)-hyc(i,j,k))*dyi))
        enddo
      enddo
    enddo
  enddo ! ift

  do k=1,lls
    do j=jl,ju
      do i=il,iu
        do ift=1,nfuel
          if(xvrhoFuel(ift,i,j,k).gt.1.d-6)then
            rhos=xvrhoFuel(ift,i,j,k)+xvrhoWater(ift,i,j,k)
            frhof(ift,i,j,k)=frhof(ift,i,j,k)-siesdiff(i,j,k)*rhos/hf
            fsies(ift,i,j,k)=fsies(ift,i,j,k)+siesdiff(i,j,k)
          endif
        enddo
      enddo
    enddo
  enddo

  deallocate(hxc,hyc,sies,siesdiff)
      
end subroutine diffusesies
