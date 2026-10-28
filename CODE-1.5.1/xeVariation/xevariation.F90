!-----------------------------------------------------------------------
! ixevariation takes into account any evolution of xe 
!-----------------------------------------------------------------------
subroutine xevariation(ixevariation,xe,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : n,m,u0,uswitch,v0,vswitch,zu,dx,dy,wp
  use uvRamp_variables, only : utarget,uramptime,vtarget,vramptime
  use windfieldio, only : windspeedupfactor,itinterp,ibcells,jbcells, &
    xvdataname,xvdataold,xvdatanew,nvwind,w2f_index,itabsold,itabsnew
  use sensor_variables, only : nSensors,SensorTimes,SensorX,SensorY, &
    SensorZ,SensorUvel,SensorVvel,senDataOld,senDatanew,senw
  use xvall, only : xvfuel,irhof,iuvel,ivvel,irho
  use metric_variables, only : z,zs
  use gridsetup, only : ittot,np,mp,time
  use msga_variables, only : mpi_rank,npos,mpos
  use fuel_variables, only : lfuel
  use constants, only : tolerance
  Implicit None

  ! Local Variables
  integer,intent(in) :: ixevariation
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(wp),intent(inout) :: xe(il:iu,jl:ju,lls,nvp)

  real(wp),external :: zcart
  integer :: i,j,k,ia,ja
  integer :: kv
  real(wp) :: totw
  character(len=257) :: fxvdataname

  ! Executable Code
  select case(ixevariation)
    case (1)
      if(uswitch==1)then
        do k=1,lfuel
          do j=jl,ju
            do i=il,iu
              xe(i,j,k,iuvel)=xe(i,j,k,irho) &
                *(u0+(utarget-u0)*min(1.,time/uramptime)) &
                *max(0._wp,(1._wp-1.5_wp*sum(xvfuel(:,i,j,k,irhof)))) &
                *((zcart(z(k),i,j)-zs(i,j))/zu)**(1._wp/7._wp)
            enddo
          enddo
        enddo
        do k=lfuel+1,lls
          do j=jl,ju
            do i=il,iu
              xe(i,j,k,iuvel)=xe(i,j,k,irho) &
                *(u0+(utarget-u0)*min(1.,time/uramptime)) &
                *((zcart(z(k),i,j)-zs(i,j))/zu)**(1._wp/7._wp)
            enddo
          enddo
        enddo
      endif
      if(vswitch==1)then
        do k=1,lfuel
          do j=jl,ju
            do i=il,iu
              xe(i,j,k,ivvel)=xe(i,j,k,irho) &
                *(v0+(vtarget-v0)*min(1.,time/vramptime)) &
                *max(0._wp,(1._wp-1.5_wp*sum(xvfuel(:,i,j,k,irhof)))) &
                *((zcart(z(k),i,j)-zs(i,j))/zu)**(1._wp/7._wp)
            enddo
          enddo
        enddo
        do k=lfuel+1,lls
          do j=jl,ju
            do i=il,iu
              xe(i,j,k,ivvel)=xe(i,j,k,irho) &
                *(v0+(vtarget-v0)*min(1.,time/vramptime)) &
                *((zcart(z(k),i,j)-zs(i,j))/zu)**(1._wp/7._wp)
            enddo
          enddo
        enddo
      endif
    case (2)
      if(ittot.ge.itabsnew) then
        xvdataold = xvdatanew
        itabsold = itabsnew
        itabsnew=itabsold+itinterp
        call namefile(itabsnew/windspeedupfactor,xvdataname,fxvdataname)
        if (mpi_rank.eq.0) write(6,*) 'reading file ',fxvdataname
        call windfld_read(fxvdataname,xvdatanew)
      endif
      ! here we update xe 
      do k=1,lls
        do j=jl,ju
          ja=(mpos-1)*mp+j
          do i=il,iu
            ia=(npos-1)*np+i
            if(ia.le.ibcells.or.ia.gt.n-ibcells &
              .or.ja.le.jbcells.or.ja.gt.m-jbcells)then
              do kv=1,nvwind
                xe(i,j,k,w2f_index(kv))= &
                  (real(ittot,kind=wp)-real(itabsnew,kind=wp))/ &
                  (real(itabsold,kind=wp)-real(itabsnew,kind=wp))* &
                  (xvdataold(i,j,k,kv)-xvdatanew(i,j,k,kv)) &
                  +xvdatanew(i,j,k,kv)
              enddo
            endif
          enddo
        enddo
      enddo
      if(mpi_rank.eq.0)then
        do kv=1,nvwind
          write(6,*)'told,tot,new',itabsold,ittot,itabsnew
          write(6,*)'xvbdata: ',w2f_index(kv),'...', &
            xvdataold(1,1,1,kv),xe(1,1,1,w2f_index(kv)), &
            xvdatanew(1,1,1,kv)
        enddo
      endif
    case (3)
      do k=1,lfuel
        do j=jl,ju
          ja=(mpos-1)*mp+j
          do i=il,iu
            ia=(npos-1)*np+i
            totw=0._wp
            do kv=1,nSensors
              senw(kv)=1._wp/(tolerance+sqrt((time-SensorTimes(kv))**2 &
                +((ia-1)*dx-SensorX(kv))**2+((ja-1)*dy-SensorY(kv))**2))
              totw=totw+senw(kv)
            enddo
            xe(i,j,k,iuvel)=0._wp
            xe(i,j,k,ivvel)=0._wp
            do kv=1,nSensors
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+SensorUvel(kv)* &
                max(0._wp,1._wp-1.5_wp*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+SensorVvel(kv)* &
                max(0._wp,1._wp-1.5_wp*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(kv)/totw
            enddo
            xe(i,j,k,iuvel)=xe(i,j,k,iuvel)*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=xe(i,j,k,ivvel)*xe(i,j,k,irho)
          enddo
        enddo
      enddo 
      do k=lfuel+1,lls
        do j=jl,ju
          ja=(mpos-1)*mp+j
          do i=il,iu
            ia=(npos-1)*np+i
            totw=0._wp
            do kv=1,nSensors
              senw(kv)=1._wp/(tolerance+sqrt((time-SensorTimes(kv))**2 &
                +((ia-1)*dx-SensorX(kv))**2+((ja-1)*dy-SensorY(kv))**2))
              totw=totw+senw(kv)
            enddo
            xe(i,j,k,iuvel)=0._wp
            xe(i,j,k,ivvel)=0._wp
            do kv=1,nSensors
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+SensorUvel(kv)* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+SensorVvel(kv)* &
                 (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                 senw(kv)/totw
            enddo
            xe(i,j,k,iuvel)=xe(i,j,k,iuvel)*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=xe(i,j,k,ivvel)*xe(i,j,k,irho)
          enddo
        enddo
      enddo
    case (4)
      do while(time.gt.senDataNew(1,1))
        senDataOld = senDataNew
        do kv=1,nSensors
          read(4884+kv,*) senDataNew(kv,1),senDataNew(kv,2), &
            senDataNew(kv,3),senDataNew(kv,4)
          if(mpi_rank.eq.0) print*,'Reading Sensor',kv,'datum', &
            senDataNew(kv,1)
        enddo
      enddo
      do k=1,lfuel
        do j=jl,ju
          ja=(mpos-1)*mp+j
          do i=il,iu
            ia=(npos-1)*np+i
            totw=0._wp
            do kv=1,nSensors
              senw(kv)=sqrt((time-senDataNew(kv,1))**2+ &
                ((ia-1)*dx-SensorX(kv))**2+((ja-1)*dy-SensorY(kv))**2)
              senw(nSensors+kv)=sqrt((time-senDataOld(kv,1))**2+ &
                ((ia-1)*dx-SensorX(kv))**2+((ja-1)*dy-SensorY(kv))**2)
              totw=totw+senw(kv)+senw(nSensors+kv)
            enddo
            xe(i,j,k,iuvel)=0._wp
            xe(i,j,k,ivvel)=0._wp
            do kv=1,nSensors
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+senDataNew(kv,2)* &
                max(0._wp,1._wp-1.5_wp*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+senDataNew(kv,3)* &
                max(0._wp,1._wp-1.5_wp*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(kv)/totw
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+senDataOld(kv,2)* &
                max(0._wp,1._wp-1.5_wp*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(nSensors+kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+senDataOld(kv,3)* &
                max(0._wp,1._wp-1.5_wp*sum(xvfuel(:,i,j,k,irhof)))* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(nSensors+kv)/totw
            enddo
            xe(i,j,k,iuvel)=xe(i,j,k,iuvel)*xe(i,j,k,irho)
            xe(i,j,k,ivvel)=xe(i,j,k,ivvel)*xe(i,j,k,irho)
          enddo
        enddo
      enddo 
      do k=lfuel+1,lls
        do j=jl,ju
          ja=(mpos-1)*mp+j
          do i=il,iu
            ia=(npos-1)*np+i
            totw=0._wp
            do kv=1,nSensors
              senw(kv)=sqrt((time-senDataNew(kv,1))**2+ &
                ((ia-1)*dx-SensorX(kv))**2+((ja-1)*dy-SensorY(kv))**2)
              senw(nSensors+kv)=sqrt((time-senDataOld(kv,1))**2+ &
                ((ia-1)*dx-SensorX(kv))**2+((ja-1)*dy-SensorY(kv))**2)
              totw=totw+senw(kv)+senw(nSensors+kv)
            enddo
            xe(i,j,k,iuvel)=0._wp
            xe(i,j,k,ivvel)=0._wp
            do kv=1,nSensors
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+senDataNew(kv,2)* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+senDataNew(kv,3)* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(kv)/totw
              xe(i,j,k,iuvel)=xe(i,j,k,iuvel)+senDataOld(kv,2)* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
                senw(nSensors+kv)/totw
              xe(i,j,k,ivvel)=xe(i,j,k,ivvel)+senDataOld(kv,3)* &
                (zcart(z(k),i,j)/SensorZ(kv))**(1._wp/7._wp)* &
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
