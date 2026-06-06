      subroutine rinitturb()
      use gridsetup
      use xvo
      use xve
      use turba
      use constants
      use metryic
      use msga

      Implicit None

      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,nn
      real :: sp2
      real,external :: zcart
      if (iinra.eq.0) then
          sc=0.1
      else if (iinra.eq.1) then
          !TODO:check if 0.8*0.85 is not better compared to experiments (FP)
          sc=.08
      endif
 
      do k=1,l
      do j=1,mp
      do i=1,np
        if (isa.eq.1) then    ! saxy=saz=cst
          saxy(i,j,k)=2.0
          saz(i,j,k)=saxy(i,j,k)
        else if (isa.eq.2) then  ! saxy=saz depend on mesh
           saxy(i,j,k)=
     +        (dx*dy*(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)))**0.333   
           saz(i,j,k)=saxy(i,j,k)
        else if (isa.eq.3) then  ! saxy and saz differ
           saxy(i,j,k)=(dx*dy)**0.5
           saz(i,j,k)= zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)
        else if (isa.eq.4) then ! similar to isa=2 with a correction in the lower part of the domain
           saxy(i,j,k)=
     +        (dx*dy*(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)))**0.333   
     +        *(1+(saxy(i,j,k)/(0.4/0.23*zcart(z(k),i,j)))**2.0)**(-0.5)
           saz(i,j,k)=saxy(i,j,k)
        else if (isa.eq.5) then ! similar to isa=3 with a correction in the lower part of the domain
           saxy(i,j,k)=(dx*dy)**0.5
     +        *(1+(saxy(i,j,k)/(0.4/0.23*zcart(z(k),i,j)))**2.0)**(-0.5)
           saz(i,j,k)= zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)
     +        *(1+(saxy(i,j,k)/(0.4/0.23*zcart(z(k),i,j)))**2.0)**(-0.5)
        end if
        sa(i,j,k)=(saxy(i,j,k)*saxy(i,j,k)*saz(i,j,k))**0.333 
        sb(i,j,k)=.25  
        if (iturb.eq.1) sb(i,j,k)=0.0
!        if (mpi_rank.eq.1.and.i.eq.1.and.j.eq.1) 
!     +    write(6,*) 'turb scale in cell',k,':',sa(i,j,k),sb(i,j,k)
        if (idrag.eq.1) then
           cd(i,j,k)=1.0
        else if (idrag.eq.2) then
           cd(i,j,k)=0.15
        endif

      if (iwindfieldin.eq.0.or.icfmeflag.eq.0) then
!        if (uswitch.eq.2.or.vswitch.eq.0) then
!         xe(i,j,k,5)=1.0*xe(i,j,k,nv) !0.03
!         xe(i,j,k,6)=0.003
!        else
        sp2=xe(i,j,k,1)**2+xe(i,j,k,2)**2+xe(i,j,k,3)**2
       xe(i,j,k,5)=max(min(sp2*0.05,xe(i,j,k,nv)*2.0),xe(i,j,k,nv)*0.05)
       xe(i,j,k,6)=max(min(sp2*0.01,xe(i,j,k,nv)*0.5),xe(i,j,k,nv)*0.01)
!        end if
      end if

      enddo
      enddo
      enddo
      do nn=5,6
       call updated(xe(1-ih,1-ih,1,nn),xe(1-ih,1-ih,1,nn),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      enddo
      call updated(saxy,saxy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(saz,saz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(sa,sa,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(sb,sb,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      !call updated(sc,sc,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(cd,cd,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      return
      end
