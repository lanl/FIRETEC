!-----------------------------------------------------------------------
! ignite computes ignition terms on large time step   
!-----------------------------------------------------------------------
subroutine ignite(force_sies,xvfuel,il,iu,jl,ju,lls,nfuel,nvfuel)
  use gridlist_variables, only : cpwood,tambient,prec
  use gridsetup, only : np,mp,time,dt
  use msga_variables, only : npos,mpos
  use ign_variables, only : ignLocation,ignTime,nIgn,targetTemp, &
    rampRate
  use fuel_variables, only : temps
  use thermo_variables, only : cpwater
  use xvall, only : isies,irhof,irhow
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nfuel,nvfuel
  real(prec),intent(in) :: xvfuel(nfuel,il:iu,jl:ju,lls,nvfuel)
  real(prec),intent(inout) :: force_sies(nfuel,il:iu,jl:ju,lls)

  integer :: iign,i,j,k,ift
  real(prec) :: cpsolid

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
          temps(ift,i,j,k)=max(temps(ift,i,j,k), &
            min((time-ignTime(iign))*rampRate+tambient,targetTemp))
          cpsolid=(xvfuel(ift,i,j,k,irhof)*cpwood &
            +xvfuel(ift,i,j,k,irhow)*cpwater)/(xvfuel(ift,i,j,k,irhof) &
            +xvfuel(ift,i,j,k,irhow))
          force_sies(ift,i,j,k)=force_sies(ift,i,j,k) &
            +(cpsolid*temps(ift,i,j,k)-xvfuel(ift,i,j,k,isies))/dt
        enddo
      endif
    enddo
  endif

end subroutine ignite
