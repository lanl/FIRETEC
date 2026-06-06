!----------------------------------------------------------------
! rinitturb initializes all turbulence and drag arrays and 
! variables for use throughout the simulation
!----------------------------------------------------------------
subroutine rinitturb
  use gridlist_variables, only : isa,dx,dy,ih,l,iwindfieldin, &
    iturb,icfmeflag
  use gridsetup, only : np,mp
  use turb_variables, only : saxy,saz,sa,sb
  use metric_variables, only : z,zedge
  use xvall, only : xe,iuvel,ivvel,iwvel,ika,ikb,irho
  Implicit None

  ! Local Variables
  integer :: i,j,k
  real :: sp2
  real,external :: zcart

  ! Executable Code
  do k=1,l
    do j=1,mp
      do i=1,np
        if(isa.eq.1)then ! saxy=saz=cst
          saxy(i,j,k)=2.
          saz(i,j,k)=saxy(i,j,k)
        elseif(isa.eq.2)then ! saxy=saz depend on mesh
          saxy(i,j,k)=(dx*dy*(zcart(zedge(k+1),i,j)- &
            zcart(zedge(k),i,j)))**(1./3.)
          saz(i,j,k)=saxy(i,j,k)
        elseif(isa.eq.3)then ! saxy and saz differ
          saxy(i,j,k)=(dx*dy)**0.5
          saz(i,j,k)=zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)
        elseif(isa.eq.4)then ! similar to isa=2 with a correction in the lower part of the domain
          saxy(i,j,k)=(dx*dy*(zcart(zedge(k+1),i,j)- &
            zcart(zedge(k),i,j)))**(1./3.)
          saxy(i,j,k)=saxy(i,j,k)* &
            (1.+(saxy(i,j,k)/(0.4/0.23*zcart(z(k),i,j)))**2.)**(-0.5)
          saz(i,j,k)=saxy(i,j,k)
        elseif(isa.eq.5)then ! similar to isa=3 with a correction in the lower part of the domain
          saxy(i,j,k)=(dx*dy)**0.5* &
            (1.+(saxy(i,j,k)/(0.4/0.23*zcart(z(k),i,j)))**2.)**(-0.5)
          saz(i,j,k)=zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)* &
            (1.+(saxy(i,j,k)/(0.4/0.23*zcart(z(k),i,j)))**2.)**(-0.5)
        endif
        sa(i,j,k)=(saxy(i,j,k)**2.*saz(i,j,k))**(1./3.)
        sb(i,j,k)=0.25

        if(iwindfieldin.eq.0.or.icfmeflag.eq.0)then
          sp2=xe(i,j,k,iuvel)**2.+xe(i,j,k,ivvel)**2.+xe(i,j,k,iwvel)**2.
          if(iturb.gt.0)then
            xe(i,j,k,ika)=max(min(sp2*0.05,xe(i,j,k,irho)*2.0),xe(i,j,k,irho)*0.05)
          endif
          if(iturb.gt.1)then
            xe(i,j,k,ikb)=max(min(sp2*0.01,xe(i,j,k,irho)*0.5),xe(i,j,k,irho)*0.01)
          endif
        endif
      enddo
    enddo
  enddo

  if(iturb.gt.0)then
    call update(xe(:,:,:,ika),xe(:,:,:,ika),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
  endif
  if(iturb.gt.1)then
    call update(xe(:,:,:,ikb),xe(:,:,:,ikb),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
  endif

end subroutine rinitturb
