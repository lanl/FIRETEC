!-----------------------------------------------------------------------
! dragtk routine calculates changes and evolution of tka/tkb due 
! to drag effects
!-----------------------------------------------------------------------
subroutine dragtk(fk,tke,stk,tkeu,il,iu,jl,ju,ll)
  use gridlist_variables, only : l,wp
  use vegetation_gridlist_variables, only : nfuel,rhoMicro
  use xvall, only : xvrho,iuvel,ivvel,iwvel,xv,irho,xvfuel,irhof
  use fuel_variables, only : sizeScale,lfuel
  use gridsetup, only : np,mp
  use linn_turb_variables, only : rtke_abc
  use constants, only : pi
  Implicit None

  ! Local variables
  integer,intent(in) :: il,iu,jl,ju,ll
  real(wp),intent(in) :: tke(il:iu,jl:ju,ll)
  real(wp),intent(in) :: stk(il:iu,jl:ju,ll)
  real(wp),intent(in) :: tkeu(il:iu,jl:ju,ll)
  real(wp),intent(inout) :: fk(il:iu,jl:ju,ll)
  integer :: i,j,k,ift
  real(wp) :: sp,sqrtk
  real(wp) :: av

  ! Executable Code

  ! Turbulence baseline dissipation
  do k=1,l
    do j=1,mp
      do i=1,np
        sqrtk=sqrt(rtke_abc(i,j,k)/xv(i,j,k,irho))
        fk(i,j,k)=fk(i,j,k)-sqrtk*tke(i,j,k)/stk(i,j,k)
      enddo
    enddo
  enddo

  ! Drag cascading of turbulent scales
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        sp=sqrt(xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+ &
          xvrho(i,j,k,iwvel)**2)
        do ift=1,nfuel
          av=2._wp/sizescale(ift,i,j,k)*xvfuel(ift,i,j,k,irhof) &
            /rhoMicro/pi
          fk(i,j,k)=fk(i,j,k)+av*sp*(0.25_wp*tkeu(i,j,k)-tke(i,j,k))
        enddo
      enddo
    enddo
  enddo
end subroutine dragtk
