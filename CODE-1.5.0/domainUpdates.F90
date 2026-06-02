module SubGround_function
  use gridlist_variables, only : prec

  contains
  real(prec) function getSubGroundValue(u,uim1,uip1,ujm1,ujp1,g13,g23,g33)
    use gridlist_variables, only : dx,dy,dz
    Implicit None

    ! Local variables
    real(prec) u,uim1,uip1,ujm1,ujp1,g13,g23,g33
    real(prec) ux1,uy1 ! 2dx - u derivatives in cell 1 on the model grid
    real(prec) uzground ! gradient at the wall on the model grid

    ! Executable Code
    ux1=0.5*(uip1-uim1)/dx
    uy1=0.5*(ujp1-ujm1)/dy
    uzground=-(g13*ux1+g23*uy1)/g33 !uz : grad(u) 's normal to ground component is 0
    getSubGroundValue=u-dz*uzground

  end function getSubGroundValue
end module SubGround_function

!-----------------------------------------------------------------------
! small_domain_updates updates all environmental arrays and variables to
! be used across this large time step 
!-----------------------------------------------------------------------
subroutine implicit_domain_updates
  use forcings, only : forceSI_xv,forceSI_xvfuel
  use gridlist_variables, only : l,ih,nfuel,ifire
  use gridsetup, only : np,mp
  use xvall, only : xvtmp,xvfueltmp,xvrho,nv,irho,itemp, &
    iEmitStart,iEmitStop,irhof,irhow,isies
  use fuel_variables, only : temps,lfuel
  use thermo_variables, only : pr,tempg,cv_gas,cp_gas,mw_gas
  use emission_gridlist_variables, only : nEmit
  use emission_general_variables, only : cpEmit,mwEmit
  Implicit None

  ! Local Variables
  integer :: i,j,k,kv,ift

  ! Executable Code
  forceSI_xv=0.
  forceSI_xvfuel=0.
  
  do i=1,np
    do j=1,mp
      do k=1,l
        call updateGasThermProps( &
          xvtmp(i,j,k,iEmitStart:iEmitStop)/xvtmp(i,j,k,irho), &
          cpEmit,mwEmit,nEmit,xvtmp(i,j,k,itemp)/xvtmp(i,j,k,irho), &
          xvtmp(i,j,k,irho),cp_gas(i,j,k),cv_gas(i,j,k),mw_gas(i,j,k), &
          pr(i,j,k),tempg(i,j,k))
        do kv=1,nv
          xvrho(i,j,k,kv)=xvtmp(i,j,k,kv)/xvtmp(i,j,k,irho)
        enddo
      enddo
    enddo
  enddo
  do kv=1,nv
    call update(xvrho(:,:,1:l,kv),xvrho(:,:,1:l,kv),np,mp,l, &
      1-ih,np+ih,1-ih,mp+ih,1,0)
    call update(tempg,tempg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
  enddo
  
  if(ifire.eq.1)then
    do i=1,np
      do j=1,mp
        do k=1,lfuel
          do ift=1,nfuel
            call updateSolidThermProps(xvfueltmp(ift,i,j,k,irhof), &
              xvfueltmp(ift,i,j,k,irhow),xvfueltmp(ift,i,j,k,isies), &
              temps(ift,i,j,k))
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine implicit_domain_updates

!-----------------------------------------------------------------------
! large_domain_updates updates all environmental arrays and variables to
! be used across this large time step 
!-----------------------------------------------------------------------
subroutine explicit_domain_updates
  use gridlist_variables, only : l,ilspgf,frqlspgf,ixevariation,ih, &
    iheatsource,nfuel,iturb,ifire
  use gridsetup, only : np,mp,ittot
  use metric_variables, only : gmul,c13,c23,gi
  use forcings, only : forceLE_xv,forceSE_xv,forceLE_xvfuel, &
    forceSE_xvfuel
  use xvall, only : xvfuel,xv,xe,iuvel,ivvel,iwvel,itemp,irho, &
    iEmitStart,iEmitStop, &
    nv,irhof,irhow,xvrho,xvfuel,isies
  use fuel_variables, only : lfuel,temps
  use thermo_variables, only : pr,tempg,cp_gas,cv_gas,mw_gas
  use emission_gridlist_variables, only : nEmit
  use emission_general_variables, only : cpEmit,mwEmit
  use SubGround_function
  Implicit None

  ! Local Variables
  integer :: i,j,k,ift,kv
  real(prec):: g13, g23, g33 ! metric tensor 
  real(prec)::uground,vground,wground ! cartesian at ground level     

  ! Executable Code
  forceSE_xv=0.
  forceLE_xv=0.
  forceSE_xvfuel=0.
  forceLE_xvfuel=0.

  do i=1,np
    do j=1,mp
      do k=1,l
        call updateGasThermProps( &
          xv(i,j,k,iEmitStart:iEmitStop)/xv(i,j,k,irho), &
          cpEmit,mwEmit,nEmit,xv(i,j,k,itemp)/xv(i,j,k,irho), &
          xv(i,j,k,irho),cp_gas(i,j,k),cv_gas(i,j,k),mw_gas(i,j,k), &
          pr(i,j,k),tempg(i,j,k))
        do kv=1,nv
          xvrho(i,j,k,kv)=xv(i,j,k,kv)/xv(i,j,k,irho)
        enddo
      enddo
    enddo
  enddo
  do kv=1,nv
    call update(xvrho(:,:,1:l,kv),xvrho(:,:,1:l,kv),np,mp,l, &
      1-ih,np+ih,1-ih,mp+ih,1,0)
    call update(tempg,tempg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
  enddo
  k=1 !bottombc
  do j=1-ih+1,mp+ih-1
    do i=1-ih+1,np+ih-1
      g13 = c13(i,j)*gmul(1)
      g23 = c23(i,j)*gmul(1)
      g33 = g13**2+g23**2+gi(i,j,1)**2
      do kv=1,nv
        xvrho(i,j,0,kv)=getSubGroundValue(xvrho(i,j,k,kv), &
          xvrho(i-1,j,k,kv),xvrho(i+1,j,k,kv),xvrho(i,j-1,k,kv), &
          xvrho(i,j+1,k,kv),g13,g23,g33)
      enddo
      uground=0.5*(xvrho(i,j,1,iuvel)+xvrho(i,j,0,iuvel))
      vground=0.5*(xvrho(i,j,1,ivvel)+xvrho(i,j,0,ivvel))
      wground = - (g13 * uground + g23 * vground)/gi(i,j,1)
      xvrho(i,j,0,iwvel) = 2*wground-xvrho(i,j,1,iwvel)
    enddo
  enddo

  if(ifire.eq.1)then
    do i=1,np
      do j=1,mp
        do k=1,lfuel
          do ift=1,nfuel
            call updateSolidThermProps(xvfuel(ift,i,j,k,irhof), &
              xvfuel(ift,i,j,k,irhow),xvfuel(ift,i,j,k,isies), &
              temps(ift,i,j,k))
          enddo
        enddo
      enddo
    enddo
  endif
  if(iturb.eq.2) call turbFieldUpdate()
  if(ixevariation.ge.1) call xevariation(ixevariation,ih, &
    xe(1:np,1:mp,:,:),1,np,1,mp,l,nv)
  if(ilspgf.ge.1.and.mod(ittot,frqlspgf).eq.0) call largeScalePGF
  if(iheatsource.gt.0) call updateHeatSourcePosition(iheatsource, &
    xvfuel(:,:,:,:,irhof),nfuel,1,np,1,mp,lfuel)

end subroutine explicit_domain_updates
