!----------------------------------------------------------------
! diffusesies spreads creeping fire laterally through a diffusion 
! of ignitions technique. Not complete diffusion in that no 
! energy is moved between cells, but the effect of diffusion is 
! computed and then positive heat spread is treated like an 
! ignition consuming fuel to supply heat in adjacent cells, hence
! maintaining the mass and energy balance
!----------------------------------------------------------------
subroutine diffusesies 
  use gridlist_variables, only : ih,dx,dy,nfuel
  use gridsetup, only : np,mp
  use thermo_variables, only : tempg
  use metric_variables, only : gi
  use fuel_variables, only : lfuel,sies,siesdiff
  Implicit None

  ! Local Variables
  integer :: ift,i,j,k
  real :: coef
  real,allocatable :: hxc(:,:,:),hyc(:,:,:)

  ! Executable Code
  allocate(hxc(1-ih:np+ih,1-ih:mp+ih,lfuel))
  allocate(hyc(1-ih:np+ih,1-ih:mp+ih,lfuel))
  
  do ift=1,nfuel
    call update(sies(ift,:,:,:),sies(ift,:,:,:),np,mp,lfuel, &
      1-ih,np+ih,1-ih,mp+ih,1,0)

! compute contravariant component x-flux at (i-1/2,j,k) hxc (sqrt(g) burried in...
! top bc of phi is hardcoded kp1=l when k=l
    do k=1,lfuel
      do j=1,mp
        do i=1,np+1
          if(tempg(i,j,k).LT.800) then 
            coef=0.001
          else
            coef=0.0
          endif 
          hxc(i,j,k)=coef*(sies(ift,i,j,k)-sies(ift,i-1,j,k))/dx
        enddo
      enddo
    enddo

! compute contravariant component y-flux at (i,j-1/2,k) hyc (sqrt(g) burried in...)
! top bc of phi is hardcoded kp1=l when k=l
    do k=1,lfuel
      do i=1,np
        do j=1,mp+1  
          if(tempg(i,j,k).LT.800) then
            coef=0.001
          else
            coef=0.0
          endif
          hyc(i,j,k)=coef*(sies(ift,i,j,k)-sies(ift,i,j-1,k))/dy
        enddo
      enddo
    enddo
      
! compute Laplacian term
    siesdiff(ift,:,:,:)=0.0
    do k=1,lfuel
      do j=1,mp
        do i=1,np
          siesdiff(ift,i,j,k)=2.*max(0.,gi(i,j,k)* &
            ((hxc(i+1,j,k)-hxc(i,j,k))/dx+(hyc(i,j+1,k)-hyc(i,j,k))/dy))
        enddo
      enddo
    enddo
  enddo ! ift

  deallocate (hxc,hyc)
      
end subroutine diffusesies
