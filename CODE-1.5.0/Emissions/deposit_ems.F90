!-----------------------------------------------------------------------
! deposit calculates deposition forcing of aerosols onto solid surfaces,
! ground and vegetation assuming Brownian motion of the particles and 
! some collision efficiency parameter.
!-----------------------------------------------------------------------
subroutine deposit(forcingAero,xvAero,Temp,rhoFuel,sizeScale,rhoMicro, &
    il,iu,jl,ju,lls,nfuel)
  use gridlist_variables, only : prec
  use emission_gridlist_variables, only : nAero,nMAero
  use constants, only : kB,pi
  use metric_variables, only : zedge
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nfuel
  real(prec),intent(inout) :: forcingAero(il:iu,jl:ju,lls,nAero*nMAero)
  real(prec),intent(in) :: rhoMicro
  real(prec),intent(in) :: xvAero(il:iu,jl:ju,lls,nAero*nMAero)
  real(prec),intent(in) :: Temp(il:iu,jl:ju,lls)
  real(prec),intent(in) :: rhoFuel(nfuel,il:iu,jl:ju,lls)
  real(prec),intent(in) :: sizeScale(nfuel,il:iu,jl:ju,lls)

  integer :: i,j,k,r,it,ift
  real(prec) :: zeta=0.95 ! Water aerosols (Stulov, Murashkevich, and Fuchs, Journal of Aerosol Science (1978))
  real(prec) :: av
  real(prec) :: iM
  real(prec),external :: LagrangeInterp,zcart

  ! Executable Code
  ! Aerosol deposition on vegetation
  do i=il,iu
    do j=jl,ju
      do k=1,lls
        do it=1,nAero
          do r=1,nMAero
            iM=LagrangeInterp(nMAero,r-1./2., &
              xvAero(i,j,k,nMAero*(it-1)+1:nMAero*it))
            av=0.
            do ift=1,nfuel
              av=av+2.*(rhoFuel(ift,i,j,k)/rhoMicro)/ &
                sizeScale(ift,i,j,k)
            enddo
            forcingAero(i,j,k,nMAero*(it-1)+r)= &
              forcingAero(i,j,k,nMAero*(it-1)+r)- &
              zeta*av/4.*sqrt(8.*kB*Temp(i,j,k)/pi)*iM
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Aerosol deposition on the ground
  do i=il,iu
    do j=jl,ju
      do it=1,nAero
        do r=1,nMAero
          iM=LagrangeInterp(nMAero,r-1./2., &
            xvAero(i,j,1,nMAero*(it-1)+1:nMAero*it))
          av=1/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))
          forcingAero(i,j,1,nMAero*(it-1)+r)= &
            forcingAero(i,j,1,nMAero*(it-1)+r)- &
            zeta*av/4.*sqrt(8.*kB*Temp(i,j,1)/pi)*iM
        enddo
      enddo
    enddo
  enddo

end subroutine deposit
