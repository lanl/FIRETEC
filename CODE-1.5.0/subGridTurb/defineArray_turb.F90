!-----------------------------------------------------------------------
! Allocate subgrid turbulence arrays
!-----------------------------------------------------------------------
subroutine defineTurbArray
  use gridlist_variables, only : l,ih,iturb
  use gridsetup, only : np,mp
  use linn_turb_variables, only : sa,sb,saxy,saz,sqrtG_Kxy,sqrtG_Kz, &
    sqrtG_KG33,K_axy,K_az,K_b,rtke_abc
  Implicit None

  ! Executable Code
  if(iturb.eq.2)then
    allocate(sa(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(saxy(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(saz(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(sb(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(rtke_abc(1-ih:np+ih,1-ih:mp+ih,0:l)); rtke_abc=0.
    allocate(K_axy(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(K_az(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(K_b(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(sqrtG_Kxy(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(sqrtG_Kz(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(sqrtG_KG33(1-ih:np+ih,1-ih:mp+ih,l))
  endif

end subroutine defineTurbArray
