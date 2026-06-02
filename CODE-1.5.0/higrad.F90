!-donorcell: donorcell advection scheme for individual variable
!-gradmoa: pressure force for moa

!-----------------------------------------------------------------------
! compute the donorcell advection scheme on field xf base on
! advective velocities u1,u2,u3 (face contravariant *dt/dx*1/gi)
! this scheme is used in the inner loop of moa and in force.f
! (forcing advection in compress)
!-----------------------------------------------------------------------
subroutine donorcell(xf,il,iu,jl,ju,lls)
  use gridlist_variables, only : l,prec
  use forcings, only : u1,u2,u3
  use metric_variables, only : gi
  use gridsetup, only : np,mp
  use advection_functions, only : donor
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real(prec),intent(inout) :: xf(il:iu,jl:ju,lls)

  integer :: i,j,k
  !real(prec),external :: donor
  real(prec),allocatable :: fd(:,:,:,:)

  ! Executable Code
  allocate(fd(3,np+1,mp+1,l+1)); fd=0.
  do k=1,l
    do j=1,mp
      do i=1,np+1
        fd(1,i,j,k)=donor(xf(i-1,j,k),xf(i,j,k),u1(i,j,k))
      enddo
    enddo
  enddo
  do k=1,l
    do j=1,mp+1
      do i=1,np
        fd(2,i,j,k)=donor(xf(i,j-1,k),xf(i,j,k),u2(i,j,k))
      enddo
    enddo
  enddo
  do k=2,l
    do j=1,mp
      do i=1,np
        fd(3,i,j,k)=donor(xf(i,j,k-1),xf(i,j,k),u3(i,j,k))
      enddo
    enddo
  enddo

  do k=1,l
    do j=1,mp
      do i=1,np
        xf(i,j,k)=xf(i,j,k)-gi(i,j,k)*(fd(1,i+1,j,k)-fd(1,i,j,k) &
          +fd(2,i,j+1,k)-fd(2,i,j,k)+fd(3,i,j,k+1)-fd(3,i,j,k))
      enddo
    enddo
  enddo
  deallocate(fd)
    
end subroutine donorcell

!-----------------------------------------------------------------------
! gradmoa computes the pressure force on the small time steps and 
! add those terms to the overall forcing terms
! NB : the lower bc on pr assumes that current is pr(0)=pr(1)
! and not that pr(1/2)=pr(1)....
!-----------------------------------------------------------------------
subroutine gradmoa()
  use gridlist_variables, only : l,ih,prec
  use forcings, only : forceSI_xv
  use xvall, only : xvtmp,iuvel,ivvel,iwvel,itemp
  use metric_variables, only : gmul,gi,c13,c23
  use gridsetup, only : np,mp,dxi,dyi,dzi
  use msga_variables, only : leftdedge,rightdedge,botdedge,topdedge
  use thermo_variables, only : pre,cp_gas,cv_gas,mw_gas
  use constants, only : Rgas,Pref
  Implicit None

  ! Local Variables
  real(prec) :: g13,g23,g33
  integer :: i,j,k
  real*8 :: pr_double,rg_over_prrcp_gas
  real(prec),allocatable:: pr(:,:,:),px(:, :,:),py(:, :,:),pz(:, :,:)

  ! Executable Code
  allocate (pr(1-ih:np+ih,1-ih:mp+ih,l)); pr = 0.0
  allocate (px(1-ih:np+ih,1-ih:mp+ih,l)); px = 0.0
  allocate (py(1-ih:np+ih,1-ih:mp+ih,l)); py = 0.0
  allocate (pz(1-ih:np+ih,1-ih:mp+ih,l)); pz = 0.0

  do k=1,l
    do j=1,mp
      do i=1,np
        rg_over_prrcp_gas=Rgas/mw_gas(i,j,k)* &
          Pref**(-Rgas/cp_gas(i,j,k)/mw_gas(i,j,k))
        pr_double=(xvtmp(i,j,k,itemp)* &
          rg_over_prrcp_gas)**(cp_gas(i,j,k)/cv_gas(i,j,k))-pre(i,j,k) 
#ifdef DBL_PREC
        pr(i,j,k)=pr_double
#else
        pr(i,j,k)=real(pr_double)
#endif
      enddo
    enddo
  enddo
  call update(pr,pr,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
  
  do i=1+leftdedge,np-rightdedge                       
    do k=1,l
      do j=1,mp
        px(i,j,k)=0.5*(pr(i+1,j,k)-pr(i-1,j,k))*dxi
      enddo
    enddo
  enddo
  if(leftdedge.eq.1) then
    do k=1,l
      do j=1,mp
        px(1,j,k)=(pr(2,j,k)-pr(1,j,k))*dxi
      enddo
    enddo
  endif
  if(rightdedge.eq.1) then                     
    do k=1,l
      do j=1,mp
        px(np,j,k)=(pr(np,j,k)-pr(np-1,j,k))*dxi
      enddo
    enddo
  endif

  do k=1,l
    do i=1,np
      do j=1+botdedge,mp-topdedge                  
        py(i,j,k)=0.5*(pr(i,j+1,k)-pr(i,j-1,k))*dyi
      enddo
    enddo
  enddo
  if(botdedge.eq.1) then                     
    do k=1,l
      do i=1,np
        py(i,1,k)=(pr(i,2,k)-pr(i,1,k))*dyi
      enddo
    enddo
  endif
  if(topdedge.eq.1) then
    do k=1,l
      do i=1,np
        py(i,mp,k)=(pr(i,mp,k)-pr(i,mp-1,k))*dyi
      enddo
    enddo
  endif

  do k=2,l-1
    do j=1,mp
      do i=1,np
        pz(i,j,k)=0.5*(pr(i,j,k+1)-pr(i,j,k-1))*dzi
      enddo
    enddo
  enddo
  
  do j=1,mp
    do i=1,np
      pz(i,j,1)= 0.5*(pr(i,j,2)-pr(i,j,1))*dzi
      pz(i,j,l)= 0.5*(pr(i,j,l)-pr(i,j,l-1))*dzi
    enddo
  enddo
  
  do k=1,l
    do j=1,mp
      do i=1,np
        g13=c13(i,j)*gmul(k)
        g23=c23(i,j)*gmul(k)
        g33= gi(i,j,k)

        forceSI_xv(i,j,k,iuvel)=forceSI_xv(i,j,k,iuvel) &
          -(px(i,j,k)+g13*pz(i,j,k))
        forceSI_xv(i,j,k,ivvel)=forceSI_xv(i,j,k,ivvel) &
          -(py(i,j,k)+g23*pz(i,j,k))
        forceSI_xv(i,j,k,iwvel)=forceSI_xv(i,j,k,iwvel) &
          -(g33*pz(i,j,k))
      enddo
    enddo
  enddo

  deallocate (pr,px,py,pz)
end subroutine gradmoa
