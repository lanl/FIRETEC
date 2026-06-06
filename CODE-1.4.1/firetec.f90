!----------------------------------------------------------------
! firetec subroutine resolves burning, heat transfer between
! solid and gas phase, and radiation heat transport
!----------------------------------------------------------------
subroutine firetec
  use gridlist_variables, only : ntp,ih,irad,idiffsies,irhovapor,nfuel, &
    icallrad,l,ifire
  use xvall, only : iO2,ivapor,xv,irho,ika,ikb
  use thermo_variables, only : tempg
  use fuel_variables, only : lfuel,temps,min_rhoFuel,rhoFuel,siesDiff,fcorr,rhoWater,actualFuelDepth
  use gridsetup, only : np,mp,ittot
  use msga_variables, only : mpi_rank,ierror
  use turba
  use workavg
  Implicit None

  ! Local Variables
  integer :: i,j,k,ift,its
  real :: rneteng,ff_sum,fw_sum
  real :: rhoFuelold
  character(len=10) :: iftname
  
  ! Executable Code
  do ift=1,nfuel
    write(iftname,"(i0)") ift
    call rmaxmin(temps(ift,1:np,1:mp,:),'temps_'//iftname, &
      1,np,1,mp,lfuel,1)
  enddo
  call rmaxmin(tempg,'tempg',1-ih,np+ih,1-ih,mp+ih,l,1)

  ! Radiation
  if(mod(ittot,icallrad).eq.0)then
    if(irad.eq.1)then ! Diffusive Radiation Scheme
      if(mpi_rank.eq.0) print*,'Diffusive Radiation Not implemented'
      call mpi_finalize(ierror)
      STOP
    elseif(irad.eq.2)then ! Monte-Carlo Radiation Scheme
      call firerad_MC
    elseif(irad.eq.3)then ! Radiation Sink
      call radiationSink
    endif
  endif
      
  if(idiffsies.eq.1)then
    call diffusesies
    call rmaxmin(siesdiff,'siesdiff',1,np,1,mp,lfuel,1)
  endif
  
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        if (sum(rhoFuel(:,i,j,k)).gt.min_rhoFuel) then
          ! begin firetec small timestepping loop
          do its=1,ntp
            ff_sum=0.0;fw_sum=0.0;rneteng=0.0
            rhoFuelold=sum(rhoFuel(:,i,j,k))
            call updateGasThermProps(ifire*xv(i,j,k,iO2)/xv(i,j,k,irho), &
              irhovapor*xv(i,j,k,ivapor)/xv(i,j,k,irho))
            call fuelDragCorrection(fcorr,nfuel,actualFuelDepth(:,i,j),rhoFuel(:,i,j,k),rhoWater(:,i,j,k),k)
            do ift=1,nfuel
              if (rhoFuel(ift,i,j,k).gt.min_rhoFuel/nfuel) then
                call convection(ift,i,j,k)
                call fuel(ift,i,j,k,rneteng,ff_sum,fw_sum,rhoFuelold)
              endif
            enddo
            call gas_dtp(i,j,k,rneteng,ff_sum,fw_sum)
          enddo ! loop on its
        else ! rhoFuel(i,j,k).le.minrhoFuel, loop on large time step
          call gas_ltp(i,j,k)
        endif
      enddo !i
    enddo !j
  enddo !k
  
  ! LOOP ON LARGE TIME STEP ELSEWHERE
  do k=1+lfuel,l
    do j=1,mp
      do i=1,np
        call gas_ltp(i,j,k)
      enddo !i
    enddo !j
  enddo !k

  ! Ignite new cells
  call ignite

end subroutine firetec

!------------------------------------------------------------------
! gas_dtp compute source terms for gas equations on small time step   
!------------------------------------------------------------------
subroutine gas_dtp(i,j,k,rneteng,ff_sum,fw_sum)
  use gridlist_variables, only : irhovapor,cpwood,Water2WoodRatio
  use forcings, only : force
  use gridsetup, only : dtp
  use xvall, only : xv,nv,itemp,iO2,irho,ivapor
  use thermo_variables, only : pr,cp_gas,rg_over_prrcp_gas, &
    cp_over_cv_gas,rg_over_cp_gas,cpwater
  use fuel_variables, only : rnfuel,rno,convht,tcrit,twvap
  use radiation_variables, only : firad
  use turba
  use workavg
  Implicit none

  ! Local Variables
  integer,intent(in) ::i,j,k
  real,intent(in) :: rneteng
  real,intent(in) :: ff_sum,fw_sum
  integer :: kv
  real :: gammaterm,rwatergainht,rnetmass,energy
  real :: fi,frho,fox,frhovapor
  real :: cv_air=717.
  real :: cv_vapor=1996.

  ! Executable Code
  do kv=1,3 ! Winds
    xv(i,j,k,kv)=xv(i,j,k,kv)+0.5*force(i,j,k,kv)*dtp
  enddo
  do kv=4,nv-1 ! Everything else
    xv(i,j,k,kv)=max(1.e-8,xv(i,j,k,kv)+0.5*force(i,j,k,kv)*dtp)
  enddo
    
  frho=(ff_sum*rnfuel+fw_sum)*dtp ! Mass released to gas phase
  fox=-ff_sum*rno*dtp ! Mass of oxygen change
  if(irhovapor.eq.1) then
    frhovapor=(fw_sum+ff_sum*rnfuel*Water2WoodRatio)*dtp
    energy=energy+0.5*force(i,j,k,ivapor)*(cv_vapor-cv_air)  ! vapor diffusion term
  endif
    
  ! Computation of theta eq rhs
  gammaterm=-sum(convht(:,i,j,k))+firad(i,j,k)  !-qxt
  rwatergainht=fw_sum*cpwater*twvap
  rnetmass=ff_sum*rnfuel*tcrit*cpwood+rwatergainht
  energy=gammaterm+rneteng+rnetmass
  pr(i,j,k)=(xv(i,j,k,itemp)*rg_over_prrcp_gas)**cp_over_cv_gas
  fi=energy*dtp/cp_gas*(1.e5/pr(i,j,k))**rg_over_cp_gas

  xv(i,j,k,itemp)=xv(i,j,k,itemp)+fi
  xv(i,j,k,iO2)=xv(i,j,k,iO2)+fox
  xv(i,j,k,irho)=xv(i,j,k,irho)+frho
  if(irhovapor.eq.1) xv(i,j,k,ivapor)=xv(i,j,k,ivapor)+frhovapor

end subroutine gas_dtp

!----------------------------------------------------------------
! gas_ltp compute source terms for gas equations on small time step   
!----------------------------------------------------------------
subroutine gas_ltp(i,j,k)
  use forcings, only : force
  use gridsetup, only : dt
  use xvall, only : xv,nv,itemp
  use thermo_variables, only : pr,rg_over_prrcp_gas,cp_over_cv_gas, &
    rg_over_cp_gas,cp_gas
  use radiation_variables, only : firad
  use turba
  use workavg
  Implicit none
  
  ! Local Variables
  integer,intent(in) :: i,j,k
  integer :: kv
  real :: energy,capqterm

  ! Executable Code
  energy=firad(i,j,k)
  pr(i,j,k)=(xv(i,j,k,itemp)*rg_over_prrcp_gas)**cp_over_cv_gas
  capqterm=energy*(1.e5/pr(i,j,k))**rg_over_cp_gas/cp_gas !rrl
  xv(i,j,k,itemp)=xv(i,j,k,itemp)+2.*capqterm*dt
  
  do kv=1,3 ! Winds
    xv(i,j,k,kv)=xv(i,j,k,kv)+0.5*force(i,j,k,kv)*dt
  enddo
  do kv=4,nv-1 ! Everything else
    xv(i,j,k,kv)=max(1.e-8,xv(i,j,k,kv)+0.5*force(i,j,k,kv)*dt)
  enddo

end subroutine gas_ltp

!----------------------------------------------------------------
! ignite computes ignition terms on large time step   
!----------------------------------------------------------------
subroutine ignite
  use gridlist_variables, only : cpwood,tambient,nfuel
  use gridsetup, only : np,mp,time
  use msga_variables, only : npos,mpos
  use ign_variables, only : ignLocation,ignTime,nIgn,targetTemp, &
    rampRate
  use fuel_variables, only : temps,sies,rhoWater
  Implicit None

  ! Local Variables
  integer :: iign,i,j,k,ift

  ! Executable Code
  if(time.le.maxval(igntime)+(targetTemp-tambient)/rampRate)then
    do iign=1,nIgn
      if(time.ge.ignTime(iign).and. &
        time.le.ignTime(iign)+(targetTemp-tambient)/rampRate.and. &
        npos.eq.ignLocation(iign,1)/np+1.and. &
        mpos.eq.ignLocation(iign,2)/mp+1)then
        i=mod(ignLocation(iign,1),np)+1
        j=mod(ignLocation(iign,2),mp)+1
        k=ignLocation(iign,3)
        do ift=1,nfuel
          temps(ift,i,j,k)=max(temps(ift,i,j,k),min((time-ignTime(iign)) &
            *rampRate+tambient,targetTemp))
          sies(ift,i,j,1)=cpwood*temps(ift,i,j,k)
          rhoWater(ift,i,j,k)=0.0
        enddo
      endif
    enddo
  endif

end subroutine ignite
