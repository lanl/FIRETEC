!-----------------------------------------------------------------------
! stokesDrag calculates settling forcing of aerosols using Stokes
! law drag for spherical particles in a fluid
!-----------------------------------------------------------------------
subroutine stokesDrag(forcingAero,xvAero,Temp,rho,il,iu,jl,ju,lls)
  use gridlist_variables, only : wp,l
  use emission_gridlist_variables, only : nAero,nMAero
  use emission_general_variables, only : rhoAero,LagrangeInterp
  use constants, only : Rgas,g,pi
  use gridsetup, only : np,mp
  use metric_variables, only : zedge,zs
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real(wp),intent(inout) :: forcingAero(il:iu,jl:ju,lls,nAero*nMAero)
  real(wp),intent(in) :: xvAero(il:iu,jl:ju,lls,nAero*nMAero)
  real(wp),intent(in) :: Temp(il:iu,jl:ju,lls)
  real(wp),intent(in) :: rho(il:iu,jl:ju,lls)

  integer :: i,j,k,r,it
  real(wp) :: mw_Air=28.96_wp
  real(wp) :: mu,w_set
  real(wp) :: iM,iMp,iMm
  real(wp),external :: zcart

  ! Executable Code
  do k=2,l-1
    do j=1,mp
      do i=1,np
        mu=sqrt(8._wp*Rgas*Temp(i,j,k)/(pi*mw_Air))
        do it=1,nAero
          w_set=2._wp*(rhoAero(it)-rho(i,j,k))*g/(9._wp*mu)* &
            (3._wp/(4._wp*pi*rhoAero(it)))**(2._wp/3._wp) ! Non-moment settling velocity (m/s kg^2/3) 
          do r=1,nMAero
            iMm=LagrangeInterp(nMAero,r+2._wp/3._wp, &
              xvAero(i,j,k-1,nMAero*(it-1)+1:nMAero*it))
            iM=LagrangeInterp(nMAero,r+2._wp/3._wp, &
              xvAero(i,j,k,nMAero*(it-1)+1:nMAero*it))
            iMp=LagrangeInterp(nMAero,r+2._wp/3._wp, &
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
  do j=1,mp
    do i=1,np
      mu=sqrt(8._wp*Rgas*Temp(i,j,1)/(pi*mw_Air))
      do it=1,nAero
        w_set=2._wp*(rhoAero(it)-rho(i,j,1))*g/(9._wp*mu)* &
          (3._wp/(4._wp*pi*rhoAero(it)))**(2._wp/3._wp) ! Non-moment settling velocity (m/s kg^2/3) 
        do r=1,nMAero
          iM=LagrangeInterp(nMAero,r+2._wp/3._wp, &
            xvAero(i,j,1,nMAero*(it-1)+1:nMAero*it))
          iMp=LagrangeInterp(nMAero,r+2._wp/3._wp, &
            xvAero(i,j,2,nMAero*(it-1)+1:nMAero*it))
          forcingAero(i,j,1,nMAero*(it-1)+r)= &
            forcingAero(i,j,1,nMAero*(it-1)+r)+w_set* &
            ((iMp-iM)/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))- &
              iM/(zcart(zedge(2),i,j)-zs(i,j)))
        enddo
      enddo
    enddo
  enddo
end subroutine stokesDrag
