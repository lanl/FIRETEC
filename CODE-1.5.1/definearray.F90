!-----------------------------------------------------------------------
! Allocate global arrays
!-----------------------------------------------------------------------
subroutine definearray
  use gridlist_variables, only : nts,ih,iturb,ilspgf,l,n,m,wp
  use forcings, only : forceSE_xv,forceSI_xv,forceLE_xv,xvLim, &
    uavg,vavg,wavg,oavg,wt,wtf,u1,u2,u3
  use xvall, only : xv,xe,xvtmp,relaxxv,nv,xvrho,xvrhotmp
  use metric_variables, only : x,y,z,zedge,zs,gmul,c13,c23,h,gi
  use thermo_variables, only : tempg,pr,pre,cp_gas,cv_gas,mw_gas
  use lspgf_variables, only : sintheta,sinthetaf,flspgf
  use gridsetup, only : np,mp
  Implicit None

  ! Executable Code
  call xvpointers
  
  allocate(xe(np,mp,l,nv)); xe=0._wp
  allocate(xv(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(xvLim(2,nv)); xvLim(1,:)=-1.e20_wp; xvLim(2,:)=1.e20_wp
  allocate(xvrho(1-ih:np+ih,1-ih:mp+ih,0:l,nv-1))
  allocate(xvtmp(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(xvrhotmp(1-ih:np+ih,1-ih:mp+ih,0:l,nv-1))
  allocate(wt(nts+1)); wt=1./nts
  allocate(wtf(nts)); wtf=1./nts
  allocate(uavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(vavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(wavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(oavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(forceSE_xv(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(forceSI_xv(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(forceLE_xv(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(c13(1-ih:np+ih,1-ih:mp+ih))
  allocate(c23(1-ih:np+ih,1-ih:mp+ih))
  allocate(gi(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(gmul(l))
  allocate(h(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(zs(1-ih:np+ih,1-ih:mp+ih)); zs=0.0_wp
  allocate(x(1-ih:n+ih))
  allocate(y(1-ih:m+ih))
  allocate(z(l))
  allocate(zedge(l+1))
  allocate(relaxxv(np,mp,l,4)); relaxxv=0.0_wp
  if(iturb.gt.0) call defineTurbArray
  allocate(pr(np,mp,l))
  allocate(mw_gas(np,mp,l))
  allocate(cp_gas(np,mp,l))
  allocate(cv_gas(np,mp,l))
  allocate(tempg(1-ih:np+ih,1-ih:mp+ih,l))
  if(ilspgf.ge.1) then
    allocate(sintheta(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(sinthetaf(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(flspgf(np,mp))
  endif
  allocate(pre(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(u1(np+1,mp,l))
  allocate(u2(np,mp+1,l))
  allocate(u3(np,mp,l+1))

end subroutine definearray

!-----------------------------------------------------------------------
! xvpointers sets up an array of indexes in the xv array for a more
! dynamic capability in the transporting array
!-----------------------------------------------------------------------
subroutine xvpointers
  use gridlist_variables, only : iturb,iemissions
  use vegetation_gridlist_variables, only : ifire,inonlocal
  use radiation_gridlist_variables, only : irad,isolar
  use emission_gridlist_variables, only : nEmit,nAero,nMAero
  use emission_general_variables, only : spEmit,spAero
  use xvall
  Implicit None

  ! Local Variables
  integer :: i,j
  integer :: indx

  ! Executable Code
  ! Transported array pointers
  iuvel = 1
  ivvel = 2
  iwvel = 3
  itheta = 4
  nv = 4
  if(iturb.gt.0)then
    nv = nv+1
    ika = nv
    if(iturb.gt.1)then
      nv = nv+1
      ikb = nv
    endif
  endif
  if(iemissions.gt.0)then
    iEmitStart = nv+1
    if(nEmit.gt.0)then
      iEmitStop  = nv+nEmit
      nv = nv+nEmit
      do i=1,nEmit
        if(spEmit(i).eq.'O2')  iO2  = iEmitStart-1+i
        if(spEmit(i).eq.'H2O') iH2O = iEmitStart-1+i
      enddo
    endif
    if(nAero.gt.0)then
      iEmitStop  = nv+nAero*nMAero
      nv = nv + nAero*nMAero
    endif
  endif
  if(inonlocal.eq.1)then
    nv = nv+1
    imixfrac = nv
  endif
  nv = nv+1
  irho = nv
  
  allocate(xv_list(nv))
  xv_list(irho) = 'rho'
  if(iemissions.gt.0)then
    if(nAero.gt.0)then
      do i=1,nAero
        do j=1,nMAero
          indx=index(spAero(i),' ')
          write(xv_list(iEmitStart+nEmit+nMAero*(i-1)+j-1),"(A,I0)") &
            spAero(i)(1:indx),j
        enddo
      enddo
    endif
    if(nEmit.gt.0)then
      do i=1,nEmit
        xv_list(iEmitStart-1+i) = spEmit(i)
      enddo
    endif
  endif
  xv_list(imixfrac) = 'mixfrac'
  xv_list(ikb) = 'kb'
  xv_list(ika) = 'ka'
  xv_list(itheta) = 'theta'
  xv_list(iwvel) = 'w'
  xv_list(ivvel) = 'v'
  xv_list(iuvel) = 'u'

  ! Non-advected solid array pointers
  irhof = 1
  irhow = 2
  nvfuel = 2
  if(isolar.eq.1.or.ifire.eq.1.or.irad.gt.0)then
    isies = nvfuel+1
    nvfuel = nvfuel+1
  endif
  if(ifire.eq.1)then
    ipsiw = nvfuel+1
    nvfuel = nvfuel+1
  endif

  allocate(xvfuel_list(nvfuel))
  xvfuel_list(ipsiw) = 'psiwmax'
  xvfuel_list(isies) = 'sies'
  xvfuel_list(irhow) = 'rhoWater'
  xvfuel_list(irhof) = 'rhoFuel'

end subroutine xvpointers
