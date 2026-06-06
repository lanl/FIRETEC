!----------------------------------------------------------------
! dragm routine calculates changes and evolution of wind due to 
! drag effects
!----------------------------------------------------------------
subroutine dragm(xv,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : nfuel,rhoMicro
  use forcings, only : force
  use xvall, only: xvrho,iuvel,ivvel,iwvel,ika,ikb,irho
  use fuel_variables, only : lfuel,rhoFuel,min_rhoFuel,sizeScale, &
    actualFuelDepth
  use metric_variables, only : zedge,zs,z
  use gridsetup, only : np,mp,dt
  use turb_variables, only : cd
  use constants, only : pi
  use turba
  use workavg
  use metric_variables_old
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real,intent(inout) :: xv(il:iu,jl:ju,lls,nvp)

  integer :: ift,i,j,k,kv
  real :: sp,av
  real :: zla2,zla3,rinterp
  real :: dragxt,dragyt,dragzt
  real :: uFuel,vFuel,wFuel
  real :: bendh,bendx,bendy,bendz
  real,external :: zcart

  ! Executable Code
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        if(sum(rhoFuel(:,i,j,k)).gt.min_rhoFuel) then
          dragxt=0.
          dragyt=0.
          dragzt=0.
          if(k.eq.1) then
            do ift=1,nfuel
              ! Fuel velocities, interprets velocity at fuel top based on 2nd and 3rd cell velocities
              zla3=zcart(z(3),i,j)-zs(i,j)
              zla2=zcart(z(2),i,j)-zs(i,j)
              rinterp=(actualFuelDepth(ift,i,j)-zla2)/(zla3-zla2)
              uFuel=(xvrho(i,j,2,iuvel)+rinterp*(xvrho(i,j,3,iuvel)-xvrho(i,j,2,iuvel)))/2.
              vFuel=(xvrho(i,j,2,ivvel)+rinterp*(xvrho(i,j,3,ivvel)-xvrho(i,j,2,ivvel)))/2.
              wFuel=(xvrho(i,j,2,iwvel)+rinterp*(xvrho(i,j,3,iwvel)-xvrho(i,j,2,iwvel)))/2.
              sp=sqrt(uFuel**2.+vFuel**2.+wFuel**2.)
              
              ! Bending grass correlations
              bendh=cos(atan((uFuel**2.+vFuel**2.)/(3.+wFuel**2.)))
              bendx=1.-(1.-bendh)*(xvrho(i,j,k,iuvel)**2./ &
                (xvrho(i,j,k,iuvel)**2.+xvrho(i,j,k,ivvel)**2.+1e-6))
              bendy=1.-(1.-bendh)*(xvrho(i,j,k,ivvel)**2./ &
                (xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+1e-6))
              bendz=sin(atan((xvrho(i,j,k,iuvel)**2.+xvrho(i,j,k,ivvel)**2.)/ &
                (3.+xvrho(i,j,k,iwvel)**2.)))
              
              av=2./pi/sizeScale(ift,i,j,k)*rhoFuel(ift,i,j,k)/rhoMicro !correct for cylinder if ss is radius
              dragxt=dragxt+.5*cd*av*sp*bendx/sqrt(2.)
              dragyt=dragyt+.5*cd*av*sp*bendy/sqrt(2.)
              dragzt=dragzt+.5*cd*av*sp*(1./sqrt(2.)+(1.-1./sqrt(2.))*bendz)
            enddo
          else
            sp=sqrt(xvrho(i,j,k,iuvel)**2.+xvrho(i,j,k,ivvel)**2.+ &
              xvrho(i,j,k,iwvel)**2.)
            do ift=1,nfuel
              av=2./pi/sizeScale(ift,i,j,k)*rhoFuel(ift,i,j,k)/rhoMicro !correct for cylinder if ss is radius
              
              dragxt=dragxt+.5*cd*av*sp
              dragyt=dragyt+.5*cd*av*sp
              dragzt=dragzt+.5*cd*av*sp
            enddo
          endif
          force(i,j,k,iuvel)=force(i,j,k,iuvel)-2.*dragxt*xv(i,j,k,iuvel)*dt
          force(i,j,k,ivvel)=force(i,j,k,ivvel)-2.*dragyt*xv(i,j,k,ivvel)*dt
          force(i,j,k,iwvel)=force(i,j,k,iwvel)-2.*dragzt*xv(i,j,k,iwvel)*dt
        endif !rhoFuel.gt.rhoFuel_min
      enddo !i
    enddo !j
  enddo !k
 
end subroutine dragm

!----------------------------------------------------------------
! dragtk routine calculates changes and evolution of tka/tkb due 
! to drag effects
!----------------------------------------------------------------
subroutine dragtk(fk,tke,stk,tkeu,il,iu,jl,ju,ll)
  use gridlist_variables, only : l,nfuel,rhoMicro
  use xvall, only : xvrho,iuvel,ivvel,iwvel,xv,ika,ikb,irho
  use fuel_variables, only : sizeScale,rhoFuel,lfuel
  use gridsetup, only : np,mp,dt
  use turb_variables, only : cd,rtke_abc
  use constants, only : pi
  use turba
  Implicit None

  ! Local variables
  integer,intent(in) :: il,iu,jl,ju,ll
  real,intent(in) :: tke(il:iu,jl:ju,ll)
  real,intent(in) :: stk(il:iu,jl:ju,ll)
  real,intent(in) :: tkeu(il:iu,jl:ju,ll)
  real,intent(inout) :: fk(il:iu,jl:ju,ll)
  integer :: i,j,k,ift
  real :: sp,sqrtk
  real :: av,drag

  ! Executable Code
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        sp=sqrt(xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+xvrho(i,j,k,iwvel)**2)
        sqrtk=sqrt(rtke_abc(i,j,k)/xv(i,j,k,irho))
        drag=-sqrtk*tke(i,j,k)/stk(i,j,k)
        do ift=1,nfuel
          av=2./sizescale(ift,i,j,k)*rhoFuel(ift,i,j,k)/rhoMicro/pi
          drag = drag+av*sp*(0.25*tkeu(i,j,k)-tke(i,j,k))
        enddo
        fk(i,j,k)=fk(i,j,k)+2.*drag*dt
      enddo
    enddo
  enddo
  do k=lfuel+1,l
    do j=1,mp
      do i=1,np
        sqrtk=sqrt(rtke_abc(i,j,k)/xv(i,j,k,irho))
        drag=-sqrtk*tke(i,j,k)/stk(i,j,k)
        fk(i,j,k)=fk(i,j,k)+2.*drag*dt
      enddo
    enddo
  enddo
end subroutine dragtk

!----------------------------------------------------------------
! fuelDragCorrection calculates a correction term for wind velocities
! and turbulent kinetic energy in the fuel layers
!----------------------------------------------------------------
subroutine fuelDragCorrection(fcorr,nfuel,afd,rhoFuel,rhoWater,iflag)
  Implicit None

  ! Local Variables
  integer,intent(in) :: nfuel,iflag
  real,intent(in) :: afd(nfuel),rhoFuel(nfuel),rhoWater(nfuel)
  real,intent(inout) :: fcorr(nfuel,2)

  integer :: it,ift
  real :: totRhos,totVels,height
  integer,allocatable :: itfindex(:)
  real,allocatable :: afdtmp(:),velocityFuel(:)

  ! Executable Code
  if(iflag.ne.1) then ! Not bottom cell no need for correction
    fcorr(:,:)=1.
    return
  endif

  ! Sort individual fuel zones by afd
  allocate(afdtmp(nfuel),itfindex(nfuel),velocityFuel(nfuel))
  afdtmp=afd(:)
  do it=1,nfuel
    itfindex(minloc(afdtmp,dim=1))=it
    afdtmp(minloc(afdtmp,dim=1))=maxval(afd)+1
  enddo

  ! Total Fuel Mass Layers
  totRhos=0.
  do ift=1,nfuel
    totRhos=totRhos+(rhoFuel(ift)+rhoWater(ift))*afd(ift)
  enddo

  ! Velocity Fuel Layer Correction
  do ift=1,nfuel
    velocityFuel(ift)=0.
    do it=1,nfuel
      velocityFuel(ift)=velocityFuel(ift)+(rhoFuel(ift)+rhoWater(ift))*min(afd(it),afd(ift))/totRhos
    enddo
    fcorr(ift,1)=velocityFuel(ift)/2.
  enddo

  ! Energy Fuel Layer Correction
  totVels=0.
  do ift=1,nfuel
    totVels=totVels+velocityFuel(ift)**2.*afd(ift)
  enddo
  do ift=1,nfuel
    fcorr(ift,2)=velocityFuel(ift)**2.*afd(ift)/totVels
  enddo
  deallocate(afdtmp,itfindex,velocityFuel)

end subroutine fuelDragCorrection
