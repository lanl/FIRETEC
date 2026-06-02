!-----------------------------------------------------------------------
! stokesDrag calculates settling forcing of aerosols using Stokes
! law drag for spherical particles in a fluid
!-----------------------------------------------------------------------
subroutine stokesDrag(forcingAero,xvAero,Temp,rho,il,iu,jl,ju,lls)
  use gridlist_variables, only : prec
  use emission_gridlist_variables, only : nAero,nMAero
  use emission_general_variables, only : rhoAero
  use constants, only : Rgas,g,pi
  use metric_variables, only : zedge,zs
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real(prec),intent(inout) :: forcingAero(il:iu,jl:ju,lls,nAero*nMAero)
  real(prec),intent(in) :: xvAero(il:iu,jl:ju,lls,nAero*nMAero)
  real(prec),intent(in) :: Temp(il:iu,jl:ju,lls)
  real(prec),intent(in) :: rho(il:iu,jl:ju,lls)

  integer :: i,j,k,r,it
  real(prec) :: mw_Air=28.96
  real(prec) :: mu,w_set
  real(prec) :: iM,iMp,iMm
  real(prec),external :: LagrangeInterp,zcart

  ! Executable Code
  do i=il,iu
    do j=jl,ju
      do k=2,lls-1
        mu=sqrt(8*Rgas*Temp(i,j,k)/(pi*mw_Air))
        do it=1,nAero
          w_set=2*(rhoAero(it)-rho(i,j,k))*g/(9.*mu)* &
            (3./(4.*pi*rhoAero(it)))**(2./3.) ! Non-moment settling velocity (m/s kg^2/3) 
          do r=1,nMAero
            iMm=LagrangeInterp(nMAero,r+2./3., &
              xvAero(i,j,k-1,nMAero*(it-1)+1:nMAero*it))
            iM=LagrangeInterp(nMAero,r+2./3., &
              xvAero(i,j,k,nMAero*(it-1)+1:nMAero*it))
            iMp=LagrangeInterp(nMAero,r+2./3., &
              xvAero(i,j,k+1,nMAero*(it-1)+1:nMAero*it))
            forcingAero(i,j,k,nMAero*(it-1)+r)= &
              forcingAero(i,j,k,nMAero*(it-1)+r)+w_set* &
              ((iMp-iM)/(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))- &
               (iM-iMm)/(zcart(zedge(k),i,j)-zcart(zedge(k-1),i,j)))
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Ground deposition from settling
  do i=il,iu
    do j=jl,ju
      mu=sqrt(8*Rgas*Temp(i,j,1)/(pi*mw_Air))
      do it=1,nAero
        w_set=2*(rhoAero(it)-rho(i,j,1))*g/(9.*mu)* &
          (3./(4.*pi*rhoAero(it)))**(2./3.) ! Non-moment settling velocity (m/s kg^2/3) 
        do r=1,nMAero
          iM=LagrangeInterp(nMAero,r+2./3., &
            xvAero(i,j,1,nMAero*(it-1)+1:nMAero*it))
          iMp=LagrangeInterp(nMAero,r+2./3., &
            xvAero(i,j,2,nMAero*(it-1)+1:nMAero*it))
          forcingAero(i,j,1,nMAero*(it-1)+r)= &
            forcingAero(i,j,1,nMAero*(it-1)+r)+w_set* &
            ((iMp-iM)/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))- &
              iM/(zcart(zedge(1),i,j)-zs(i,j)))
        enddo
      enddo
    enddo
  enddo


end subroutine stokesDrag
