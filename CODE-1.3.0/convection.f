ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c the subroutine convection() computes convht
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine convection(i,j,k)
      use turba
      use pres
      use fireteca
      use constants
      use gridsetup
      use xvo
      use xve
      use msga
!      use radiation,only:irad   !KOO irad in radiaion.mod

      Implicit None

      integer,intent(in) :: i,j,k
      real :: sstemp,rktemp,rhovapor
      real :: sp,re,h,av
      !if (k.gt.lfuel) return
            sstemp=sizescale(i,j,k)
            if(irhovapor.eq.1) rhovapor=xvb(i,j,k,8)

c computation of sp (velocity in fuel + turbulence)
            rktemp=xvfuel(i,j,k,4)/xvb(i,j,k,nv)+
     +             1.2*(xvfuel(i,j,k,5)/xvb(i,j,k,nv))
            sp=sqrt(rktemp)+
     &         sqrt((xvfuel(i,j,k,1)/xvb(i,j,k,nv))**2+
     &              (xvfuel(i,j,k,2)/xvb(i,j,k,nv))**2+
     &              (xvfuel(i,j,k,3)/xvb(i,j,k,nv))**2)

c computation of h and av
            re=sstemp*sp/2.e-05
            !h=.5*0.683*re**0.466*thermcondair/sstemp
            h=.25*0.683*re**0.466*thermcondair/sstemp  !FP rrl
            av=2.*(rhof(i,j,k)/rhomicro(i,j,k))/sstemp
c computation of tempg
            pr(i,j,k)=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
            tempg(i,j,k)=xvb(i,j,k,4)/xvb(i,j,k,nv)
     +                      *(pr(i,j,k)*1.e-5)**(rg/cp)

c computation of convht
            convht(i,j,k)=h*av*(tempg(i,j,k)-temps(i,j,k))    !wss
! FIXME
!        if(mpi_rank.eq.11 .and. i.eq.1 .and. j.eq.20 .and. k.eq.2) then
!        print*,'**************************************'
!        print*,'convht',convht(i,j,k),'h',h
!        print*,'tempg',tempg(i,j,k),'temps',temps(i,j,k)
!        print*,'re',re,'sp',sp,'rktemp',rktemp
!        print*,'tempg',xvb(i,j,k,4),xvb(i,j,k,nv),pr(i,j,k)
!        endif            

      return
      end subroutine convection
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

