!!-----------------------------------------------------------------------
!! Linear interpolation function
!!-----------------------------------------------------------------------
!module lin_interp_function
!  use gridlist_variables, only : prec
!
!  contains
!  real(prec) function lin_interp(xi,xl,xh,yl,yh)
!    Implicit None
!    real(prec) :: xi,xl,xh,yl,yh
!
!    lin_interp=(xi-xl)/(xh-xl)*(yh-yl)+yl
!
!  end function lin_interp
!
!end module lin_interp_function

!-----------------------------------------------------------------------
! ixevariation takes into account any evolution of xe 
!-----------------------------------------------------------------------
subroutine xevariation(ixevariation,ih,xe,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : n,m,u0,uswitch,v0,vswitch,zu,dx,dy,prec
  use uvRamp_variables, only : uramp,uramptime,vramp,vramptime
  use windfieldio, only : windspeedupfactor,itinterp,ibcells,jbcells, &
    xvdataname,xvdataold,xvdatanew,nvwind,w2f_index,itabsold,itabsnew
  use sensor_variables, only : nSensors,SensorTimes,SensorX,SensorY, &
    SensorZ,SensorUvel,SensorVvel,senDataOld,senDatanew,senw,freqSensor
  use xvall, only : xvfuel,irhof,iuvel,ivvel,irho
  use metric_variables, only : z
  use gridsetup, only : ittot,np,mp,time
  use msga_variables, only : mpi_rank,npos,mpos
  use fuel_variables, only : lfuel
  use zcart_function
  !use lin_interp_function
  Implicit None

  ! Local Variables
  integer,intent(in) :: ixevariation,ih
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(prec),intent(inout) :: xe(il:iu,jl:ju,lls,nvp)

  !real,external :: zcart,lin_interp
  integer :: i,j,k,ia,ja
  integer :: kv
  real(prec) :: totw
  character(len=257) :: fxvdataname

  ! Executable Code
  select case(ixevariation)
    case (1)
      if(uswitch==1)then
        do k=1,lfuel
          do j=1,mp
            do i=1,np
              xe(i,j,k,iuvel)=min(uramp*xe(i,j,k,irho), &
                (time*uramp/uramptime+u0)*xe(i,j,k,irho)) &
                *max(0.,(1.-1.5*sum(xvfuel(:,i,j,k,irhof))))  &
                *(min(100.,zcart(z(k),1,1))/zu)**(1./7.)
            enddo
          enddo
        enddo
        do k=lfuel+1,lls
          do j=1,mp
            do i=1,np
              xe(i,j,k,iuvel)=min(uramp*xe(i,j,k,irho), &
                (time*uramp/uramptime+u0)*xe(i,j,k,irho)) &
                *(min(100.,zcart(z(k),1,1))/zu)**(1./7.)
            enddo
          enddo
        enddo
      endif
      if(vswitch==1)then
        do k=1,lfuel
          do j=il+ih,iu-ih
            do i=jl+ih,ju-ih
              xe(i,j,k,ivvel)=min(vramp*xe(i,j,k,irho), &
                (time*vramp/vramptime+v0)*xe(i,j,k,irho)) &
                *max(0.,(1.-1.5*sum(xvfuel(:,i,j,k,irhof)))) &
                *(min(100.,zcart(z(k),1,1))/zu)**(1./7.)
            enddo
          enddo
        enddo
        do k=lfuel+1,lls
          do j=il+ih,iu-ih
            do i=jl+ih,ju-ih
              xe(i,j,k,ivvel)=min(vramp*xe(i,j,k,irho), &
                (time*vramp/vramptime+v0)*xe(i,j,k,irho)) &
                *(min(100.,zcart(z(k),1,1))/zu)**(1./7.)
            enddo
          enddo
        enddo
      endif
    case (2)
      if(mod(ittot,itinterp).eq.0) then
        xvdataold = xvdatanew
        itabsold = itabsnew
        itabsnew=itabsold+itinterp
        call namefile(itabsnew/windspeedupfactor,xvdataname,fxvdataname)
        if (mpi_rank.eq.0) write(6,*) 'reading file ',fxvdataname
        call windfld_read(fxvdataname)
      endif
      ! here we update xe 
      do k=1,lls
        do j=1,mp
          do i=1,np
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            if(ia.le.ibcells.or.ia.gt.n-ibcells &
              .or.ja.le.jbcells.or.ja.gt.m-jbcells)then
              do kv=1,nvwind
#ifdef DBL_PREC
                xe(i,j,k,w2f_index(kv))= &
                  (dble(ittot)-dble(itabsnew))/ &
                  (dble(itabsold)-dble(itabsnew))* &
                  (xvdataold(i,j,k,kv)-xvdatanew(i,j,k,kv)) &
                  +xvdatanew(i,j,k,kv)
                  !lin_interp(dble(ittot),dble(itabsnew),dble(itabsold) &
                  !,xvdatanew(i,j,k,kv),xvdataold(i,j,k,kv))

                  !(xi-xl)/real(xh-xl)*(yh-yl)+yl
                  !lin_interp(xi,xl,xh,yl,yh)
#else
                xe(i,j,k,w2f_index(kv))= &
                  (real(ittot)-real(itabsnew))/ &
                  (real(itabsold)-real(itabsnew))* &
                  (xvdataold(i,j,k,kv)-xvdatanew(i,j,k,kv)) &
                  +xvdatanew(i,j,k,kv)

                  !lin_interp(real(ittot),real(itabsnew),real(itabsold) &
                  !,xvdatanew(i,j,k,kv),xvdataold(i,j,k,kv))
#endif
              enddo
            endif
          enddo
        enddo
      enddo
      if(mpi_rank.eq.0)then
        do kv=1,nvwind
          write(6,*)'told,tot,new',itabsold,ittot,itabsnew
          write(6,*)'xvbdata: ',w2f_index(kv),'...',xvdataold(1,1,1,kv), &
            xe(1,1,1,w2f_index(kv)),xvdatanew(1,1,1,kv)
        enddo
      endif
    case (3)
      do k=1,lfuel
        do j=1,mp
          do i=1,np
            totw=0.
            do kv=1,nSensors
              senw(kv)=sqrt((time-SensorTimes(kv))**2.+ &
                ((i-1)*dx-SensorX(kv))**2.+((j-1)*dy-SensorY(kv))**2.)
              totw=totw+senw(kv)
            enddo
            xe(i,j,k,iuvel)=0.
            xe(i,j,k,ivvel)=0.
            do kv=1,nSensors
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+SensorUvel(kv)* &
                max(0.,1.-1.5*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+SensorVvel(kv)* &
                max(0.,1.-1.5*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(kv)/totw
            enddo
            xe(i,j,k,iuvel)=xe(i,j,k,iuvel)*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=xe(i,j,k,ivvel)*xe(i,j,k,irho)
          enddo
        enddo
      enddo 
      do k=lfuel+1,lls
        do j=1,mp
          do i=1,np
            totw=0.
            do kv=1,nSensors
              senw(kv)=sqrt((time-SensorTimes(kv))**2.+ &
                ((i-1)*dx-SensorX(kv))**2.+((j-1)*dy-SensorY(kv))**2.)
              totw=totw+senw(kv)
            enddo
            xe(i,j,k,iuvel)=0.
            xe(i,j,k,ivvel)=0.
            do kv=1,nSensors
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+SensorUvel(kv)* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)*senw(kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+SensorVvel(kv)* &
                 (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)*senw(kv)/totw
            enddo
            xe(i,j,k,iuvel)=xe(i,j,k,iuvel)*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=xe(i,j,k,ivvel)*xe(i,j,k,irho)
          enddo
        enddo
      enddo
    case (4)
      if(mod(ittot,freqSensor).eq.0) then
        senDataOld = senDataNew
        do kv=1,nSensors
          read(4884+kv,*) senDataNew(kv,1),senDataNew(kv,2), &
            senDataNew(kv,3),senDataNew(kv,4)
          if(mpi_rank.eq.0) print*,'Reading Sensor',kv,'datum', &
            ittot+freqSensor
        enddo
      endif
      do k=1,lfuel
        do j=1,mp
          do i=1,np
            totw=0.
            do kv=1,nSensors
              senw(kv)=sqrt((time-senDataNew(kv,1))**2.+ &
                ((i-1)*dx-SensorX(kv))**2.+((j-1)*dy-SensorY(kv))**2.)
              senw(nSensors+kv)=sqrt((time-senDataOld(kv,1))**2.+ &
                ((i-1)*dx-SensorX(kv))**2.+((j-1)*dy-SensorY(kv))**2.)
              totw=totw+senw(kv)+senw(nSensors+kv)
            enddo
            xe(i,j,k,iuvel)=0.
            xe(i,j,k,ivvel)=0.
            do kv=1,nSensors
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+senDataNew(kv,2)* &
                max(0.,1.-1.5*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+senDataNew(kv,3)* &
                max(0.,1.-1.5*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(kv)/totw
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+senDataOld(kv,2)* &
                max(0.,1.-1.5*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(nSensors+kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+senDataOld(kv,3)* &
                max(0.,1.-1.5*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(nSensors+kv)/totw
            enddo
            xe(i,j,k,iuvel)=xe(i,j,k,iuvel)*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=xe(i,j,k,ivvel)*xe(i,j,k,irho)
          enddo
        enddo
      enddo 
      do k=lfuel+1,lls
        do j=1,mp
          do i=1,np
            totw=0.
            do kv=1,nSensors
              senw(kv)=sqrt((time-senDataNew(kv,1))**2.+ &
                ((i-1)*dx-SensorX(kv))**2.+((j-1)*dy-SensorY(kv))**2.)
              senw(nSensors+kv)=sqrt((time-senDataOld(kv,1))**2.+ &
                ((i-1)*dx-SensorX(kv))**2.+((j-1)*dy-SensorY(kv))**2.)
              totw=totw+senw(kv)+senw(nSensors+kv)
            enddo
            xe(i,j,k,iuvel)=0.
            xe(i,j,k,ivvel)=0.
            do kv=1,nSensors
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+senDataNew(kv,2)* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+senDataNew(kv,3)* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(kv)/totw
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+senDataOld(kv,2)* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(nSensors+kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+senDataOld(kv,3)* &
                (zcart(z(k),j,k)/SensorZ(kv))**(1./7.)* &
                senw(nSensors+kv)/totw
            enddo
            xe(i,j,k,iuvel)=xe(i,j,k,iuvel)*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=xe(i,j,k,ivvel)*xe(i,j,k,irho)
          enddo
        enddo
      enddo 
      if(mpi_rank.eq.0)then
        do kv=1,nSensors
          print*,'Sensor',kv
          print*,'senDataOld:',senDataOld(kv,:)
          print*,'senDataNew:',senDataNew(kv,:)
        enddo
      endif
  end select 
end subroutine xevariation
