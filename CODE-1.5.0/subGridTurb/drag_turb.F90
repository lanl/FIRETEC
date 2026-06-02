!-----------------------------------------------------------------------
! dragtk routine calculates changes and evolution of tka/tkb due 
! to drag effects
!-----------------------------------------------------------------------
subroutine dragtk(fk,tke,stk,tkeu,il,iu,jl,ju,ll)
  use gridlist_variables, only : l,nfuel,rhoMicro,prec
  use xvall, only : xvrho,iuvel,ivvel,iwvel,xv,irho,xvfuel,irhof
  use fuel_variables, only : sizeScale,lfuel
  use gridsetup, only : np,mp
  use linn_turb_variables, only : rtke_abc
  use constants, only : pi
  Implicit None

  ! Local variables
  integer,intent(in) :: il,iu,jl,ju,ll
  real(prec),intent(in) :: tke(il:iu,jl:ju,ll)
  real(prec),intent(in) :: stk(il:iu,jl:ju,ll)
  real(prec),intent(in) :: tkeu(il:iu,jl:ju,ll)
  real(prec),intent(inout) :: fk(il:iu,jl:ju,ll)
  integer :: i,j,k,ift
  real(prec) :: sp,sqrtk
  real(prec) :: av,drag

  ! Executable Code
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        sp=sqrt(xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+ &
          xvrho(i,j,k,iwvel)**2)
        sqrtk=sqrt(rtke_abc(i,j,k)/xv(i,j,k,irho))
        drag=-sqrtk*tke(i,j,k)/stk(i,j,k)
        do ift=1,nfuel
          av=2./sizescale(ift,i,j,k)*xvfuel(ift,i,j,k,irhof) &
            /rhoMicro/pi
          drag = drag+av*sp*(0.25*tkeu(i,j,k)-tke(i,j,k))
        enddo
        fk(i,j,k)=fk(i,j,k)+drag
      enddo
    enddo
  enddo
  do k=lfuel+1,l
    do j=1,mp
      do i=1,np
        sqrtk=sqrt(rtke_abc(i,j,k)/xv(i,j,k,irho))
        drag=-sqrtk*tke(i,j,k)/stk(i,j,k)
        fk(i,j,k)=fk(i,j,k)+drag
      enddo
    enddo
  enddo
end subroutine dragtk
