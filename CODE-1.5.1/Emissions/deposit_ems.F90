!-----------------------------------------------------------------------
! deposit calculates deposition forcing of aerosols onto solid surfaces,
! ground and vegetation assuming Brownian motion of the particles and 
! some collision efficiency parameter.
!-----------------------------------------------------------------------
subroutine deposit(forcingAero,xvAero,Temp,rhoFuel,ilh,iuh,jlh,juh, &
    il,iu,jl,ju,lls,lfuel,nfuel)
  use gridlist_variables, only : wp
  use emission_gridlist_variables, only : nAero,nMAero
  use emission_general_variables, only : LagrangeInterp
  use constants, only : kB,pi
  use gridsetup, only : np,mp
  use vegetation_gridlist_variables, only : rhoMicro
  use fuel_variables, only : sizeScale
  use metric_variables, only : zedge
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,lls,il,iu,jl,ju,lfuel,nfuel
  real(wp),intent(inout) :: forcingAero(ilh:iuh,jlh:juh,lls,nAero*nMAero)
  real(wp),intent(in) :: xvAero(ilh:iuh,jlh:juh,lls,nAero*nMAero)
  real(wp),intent(in) :: Temp(ilh:iuh,jlh:juh,lfuel)
  real(wp),intent(in) :: rhoFuel(nfuel,il:iu,jl:ju,lfuel)

  integer :: i,j,k,r,it,ift
  real(wp) :: zeta=0.95_wp ! Water aerosols (Stulov, Murashkevich, and Fuchs, Journal of Aerosol Science (1978))
  real(wp) :: av
  real(wp) :: iM
  real(wp),external :: zcart

  ! Executable Code
  ! Aerosol deposition on vegetation
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        do it=1,nAero
          do r=1,nMAero
            iM=LagrangeInterp(nMAero,r-1._wp/2._wp, &
              xvAero(i,j,k,nMAero*(it-1)+1:nMAero*it))
            av=0._wp
            do ift=1,nfuel
              av=av+2._wp*(rhoFuel(ift,i,j,k)/rhoMicro)/ &
                sizeScale(ift,i,j,k)
            enddo
            forcingAero(i,j,k,nMAero*(it-1)+r)= &
              forcingAero(i,j,k,nMAero*(it-1)+r)- &
              zeta*av/4._wp*sqrt(8._wp*kB*Temp(i,j,k)/pi)*iM
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Aerosol deposition on the ground
  do j=1,mp
    do i=1,np
      do it=1,nAero
        do r=1,nMAero
          iM=LagrangeInterp(nMAero,r-1._wp/2._wp, &
            xvAero(i,j,1,nMAero*(it-1)+1:nMAero*it))
          av=1._wp/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))
          forcingAero(i,j,1,nMAero*(it-1)+r)= &
            forcingAero(i,j,1,nMAero*(it-1)+r)- &
            zeta*av/4._wp*sqrt(8._wp*kB*Temp(i,j,1)/pi)*iM
        enddo
      enddo
    enddo
  enddo

end subroutine deposit
