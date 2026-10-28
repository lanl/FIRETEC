!-----------------------------------------------------------------------
! rinitturb initializes all turbulence and drag arrays and 
! variables for use throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitSubTurb
  use gridlist_variables, only : dx,dy,l,iturb,wp
  use turb_gridlist_variables, only : isa
  use forcings, only : xvLim
  use gridsetup, only : np,mp
  use linn_turb_variables, only : saxy,saz,sa,sb
  use metric_variables, only : z,zedge,zs
  use xvall, only : xe,iuvel,ivvel,iwvel,ika,ikb,irho
  Implicit None

  ! Local Variables
  integer :: i,j,k
  real(wp) :: sp2
  real(wp),external :: zcart

  ! Executable Code
  select case(iturb)
    case(1)
      print*,'Smorgarinsky not implemented'
      stop
    case(2)
      do k=1,l
        do j=1,mp
          do i=1,np
            if(isa.eq.1)then ! saxy=saz=cst
              saxy(i,j,k)=2._wp
              saz(i,j,k)=saxy(i,j,k)
            elseif(isa.eq.2)then ! saxy=saz depend on mesh
              saxy(i,j,k)=(dx*dy*(zcart(zedge(k+1),i,j)- &
                zcart(zedge(k),i,j)))**(1._wp/3._wp)
              saz(i,j,k)=saxy(i,j,k)
            elseif(isa.eq.3)then ! saxy and saz differ
              saxy(i,j,k)=sqrt(dx*dy)
              saz(i,j,k)=zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)
            elseif(isa.eq.4)then ! similar to isa=2 with a correction in the lower part of the domain
              saxy(i,j,k)=(dx*dy*(zcart(zedge(k+1),i,j)- &
                zcart(zedge(k),i,j)))**(1._wp/3._wp)
              saxy(i,j,k)=saxy(i,j,k)/sqrt(1._wp+ &
                (saxy(i,j,k)/(0.4_wp/0.23_wp* &
                (zcart(z(k),i,j)-zs(i,j)))**2))
              saz(i,j,k)=saxy(i,j,k)
            elseif(isa.eq.5)then ! similar to isa=3 with a correction in the lower part of the domain
              saxy(i,j,k)=sqrt(dx*dy)
              saxy(i,j,k)=saxy(i,j,k)/sqrt(1._wp+ &
                (saxy(i,j,k)/(0.4_wp/0.23_wp* &
                (zcart(z(k),i,j)-zs(i,j)))**2))
              saz(i,j,k)=(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))/ &
                sqrt(1._wp+(saxy(i,j,k)/ &
                (0.4_wp/0.23_wp*(zcart(z(k),i,j)-zs(i,j)))**2))
            endif
            sa(i,j,k)=(saxy(i,j,k)**2._wp*saz(i,j,k))**(1._wp/3._wp)
            sb(i,j,k)=0.25_wp

            sp2=(xe(i,j,k,iuvel)**2+xe(i,j,k,ivvel)**2+ &
              xe(i,j,k,iwvel)**2)/xe(i,j,k,irho)**2
            if(iturb.eq.2)then
              xe(i,j,k,ika)=xe(i,j,k,irho)* &
                max(min(0.05_wp*sp2,2._wp),0.05_wp)
              xe(i,j,k,ikb)=xe(i,j,k,irho)* &
                max(min(0.01_wp*sp2,0.5_wp),0.01_wp)
            endif
          enddo
        enddo
      enddo
      xvLim(1,ika:ikb)=1.e-6_wp
      xvLim(2,ika:ikb)=80._wp
  end select

end subroutine rinitSubTurb
