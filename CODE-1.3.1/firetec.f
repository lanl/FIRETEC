      subroutine firetec(it)
      use metryic
      use turba
      use pres
      use gridsetup
      use fireteca
      use xvo
      use xve
      use radiation
      use msga
      use ignite
      use nonlocalfunction
      use workavg
      use constants

      Implicit None


      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: it
      integer :: i,j,k,its
      real::rneteng,ff,fw  ! defined in fuel.f

      if(inonlocal.eq.1)then
       call firetecnonlocal(it)
      else
      call rmaxmin1(temps,'temps',1-ih,np+ih,1-ih,mp+ih,l)
      call rmaxmin1(tempg,'tempg',1-ih,np+ih,1-ih,mp+ih,l)

     
      ! here we update source position AND fuel properties in the heat source zone
      if (iheatsource.ge.1.and.((time+restarttime).ge.startigntime)) 
     + call updateHeatSourcePosition()

!FP09/2019 we don't want this supension anymore, which makes no sense in
!case of multiple ignitions
! - no reaction when ignited. 
!      if(((time+restarttime).ge.startigntime).and.  !line ignition
!     &        ((time+restarttime).le.(startigntime+
!     &                         (targettemp-tambient)/ramprate)))then
!         ffparam=0.0
!      else
!         ffparam=1.0
!      endif 
 

      
      if(iwallclock.EQ.1) wtime2=MPI_Wtime()
! fp moved the frhovapor to turb...
!      if(irhovapor.eq.1)then
!       frhovaporb=0.0
!      endif

      if(iwallclock.EQ.1) call computeWallTime(wtime2,12,'init in ftec')
! fp moved the foxb to turb...
! eventually computes oxygen diffusion : dti is put in the nts loop
!      foxb=0.0
!      if(ilapdo.eq.1) 
!     +   call lapdos(xvb,xvb(1-ih,1-ih,1,7),foxb, 1-ih,np+ih,1-ih,mp+ih,l,nv)
!      if(iwallclock.EQ.1) call computeWallTime(wtime2,13,'lapdo')

! eventually computes rhovapor diffusion : dti is put in the nts loop     
!      if(irhovapor.eq.1) 
!     +   call lapdos(xvb,xvb(1-ih,1-ih,1,8),frhovaporb, 1-ih,np+ih,1-ih,mp+ih,l,nv)


!  call to radiation every icallrad large time step     
      if(iwallclock.EQ.1) wtime2=MPI_Wtime() !KOO-walltimer 
       if ((it.eq.1.and.irst.eq.0).or.(mod(it,icallrad).eq.0)) then !KOO
         if (irad.eq.1) call fire_radiation()
         if (irad.eq.2) call firerad_MC()
         if (irad.eq.3) call radiationSink()
      endif   !end if (irad.eq.1.and.(it.eq.1.or.mod(it,10).eq.0))
      if(iwallclock.EQ.1) call computeWallTime(wtime2,14,'radiation')
       ! LOOP ON SMALL TIME STEP IN FUEL, LARGE TIME STEP ELSEWHERE
       !in CODE3, KOO moved back the loop for small timestep including
       ! above lfuel. FP09/2019 put it back as it is faster
      do k=1,lfuel
!       do k=1,l
        do j=1,mp
         do i=1,np
         !if (rhof(i,j,k).gt.min_rhof.or.(iheatsource.ge.1)) then
          if (rhof(i,j,k).gt.min_rhof) then
c      begin firetec small timestepping loop
          do its=1,ntp
cc call to fuel computes xvfuel (required in convection, fuel)
         !call computexvfuel(i,j,k)
            ff=0.0;fw=0.0;rneteng=0.0
            if (iheatsource.eq.0.or.(iheatsource.eq.2.and.k.gt.1)) then 
cc call to convection to define convht        
                call convection(i,j,k)
cc ff and fw, psif, psiw psiwmax, frhof,frhowater thanks to the pdf function
              !if (rhof(i,j,k).gt.min_rhof) call fuel(i,j,k,rneteng,ff,fw)
                call fuel(i,j,k,rneteng,ff,fw)
c compute frho, fox,frhovapor
            endif ! if (iheatsource.eq.0.or.(iheatsource.eq.2.and.k.gt.1)) then
            call gas_dtp(i,j,k,rneteng,ff,fw)
            !heatSourceRHS modify fi, frho
          enddo ! loop on its
          else ! rhof(i,j,k).le.minrhof, loop on large time step
            call gas_ltp(i,j,k)
          endif
         enddo !i
        enddo !j
      enddo !k
      ! LOOP ON LARGE TIME STEP ELSEWHERE
      do k=1+lfuel,l
        do j=1,mp
         do i=1,np
            call gas_ltp(i,j,k)
         enddo !i
        enddo !j
      enddo !k
      call ignite_pattern()
      !reset forcing terms to zero 

      ! call rmaxmin(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)
      ! call rmaxmin1(firad,'firad',1-ih,np+ih,1-ih,mp+ih,l)
      ! call rmaxmin1(temps,'temps at end of FIRETEC',1-ih,np+ih,1-ih,mp+ih,l)
      ! call rmaxmin1(frhosiesrad,'frhosiesrad',1-ih,np+ih,1-ih,mp+ih,l)
      ! call rmaxmin1(firad,'firad',10,40,4,16,10)

c walltimer !KOO         
      if(iwallclock.EQ.1) call computeWallTime(wtime2,17,'updates')
                                                             
      endif  !end of nonlocal shunt
!
      !JAS adding updates to keep halos in sync with real cells 7/7/06 
!      do nn=1,nv
!       call updated(xvb(1-ih,1-ih,1,nn),xe(1-ih,1-ih,1,nn),
!     +              np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
!      enddo
!      call updated(pr,pr,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c gas_dtp compute source terms for gas equations on small time step   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gas_dtp(i,j,k,rneteng,ff,fw)
      use turba
      use turbb
      use workavg
      use constants
      use xvo
      use xve
       use gridsetup
       use fireteca
      use pres
      use ignite
      Implicit none
      integer,intent(in)::i,j,k
      real,intent(inout):: rneteng !ff(i,j,k)*hf*(1.-thetasolid)
      real,intent(inout)::ff,fw !defined in fuel.f
      real:: gammaterm,rwatergainht,rnetmass,
     +      rnetengwater,energy,capqterm
      real::fi,frho,fox,frhovapor
        frho=ff*dtp*rnfuel
        !dti is already in it (turb.f) fox=(-ff*rno+0.5*foxb(i,j,k)*dti)*dtp
        fox=(-ff*rno+0.5*foxb(i,j,k))*dtp
! FP :This test is incorrect, as foxb can be positive
!       if(fox.GT.0) then
!         print*,"negative fox dpt",ff
!         print*,i,j,k,mpi_rank
!          print*,'xvb123',xvb(i,j,k,1),xvb(i,j,k,2),xvb(i,j,k,3)
!          print*,'xvb456',xvb(i,j,k,4),xvb(i,j,k,5),xvb(i,j,k,6)
!          print*,'xvb78p',xvb(i,j,k,7),xvb(i,j,k,8),pr(i,j,k)
!          STOP
!        endif

        if(irhovapor.eq.1)
     +    frhovapor=(fw+0.5*frhovaporb(i,j,k))*dtp
        ! computation of theta eq rhs
          ! here the basic sink qxt is now in firad (mode irad=3)
        gammaterm=-convht(i,j,k)+firad(i,j,k)  !-qxt
        rwatergainht=fw*cpwater*twvap      !same as below
        rnetmass=ff*rnfuel*tcrit*cpwood+rwatergainht
        rnetengwater=fw*cvvapor*twvap
        if (irhovapor.eq.0) rnetengwater=0.0
        energy=gammaterm+rneteng+rnetmass+rnetengwater
        pr(i,j,k)=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
        capqterm=energy*(1.e5/pr(i,j,k))**(rg/cp)/cp !rrl
        fi=fib(i,j,k)*dtp+2.*capqterm*dtp
        !TODO: should change this: forcing on fuel should be done in fuel.f
        if (iheatsource.ge.1.and.((time+restarttime).ge.startigntime))
     + call heatSourceRHS(i,j,k,fi,frho,dtp)

        xvb(i,j,k,1)=xvb(i,j,k,1)+0.5*f1avg(i,j,k)*dtp
        xvb(i,j,k,2)=xvb(i,j,k,2)+0.5*f2avg(i,j,k)*dtp
        xvb(i,j,k,3)=xvb(i,j,k,3)+0.5*f3avg(i,j,k)*dtp
        xvb(i,j,k,4)=xvb(i,j,k,4)+0.5*fi
        xvb(i,j,k,5)=xvb(i,j,k,5)+0.5*fka(i,j,k)*dtp
        xvb(i,j,k,6)=xvb(i,j,k,6)+0.5*fkb(i,j,k)*dtp
        xvb(i,j,k,7)=xvb(i,j,k,7)+fox
        xvb(i,j,k,nv)=xvb(i,j,k,nv)+frho
        if(irhovapor.eq.1)
     +    xvb(i,j,k,8)=xvb(i,j,k,8)+frhovapor
      end subroutine gas_dtp


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c gas_ltp compute source terms for gas equations on small time step   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine gas_ltp(i,j,k)
      use constants
      use turba
      use turbb
      use workavg
      use pres
      use xvo
      use xve
       use gridsetup
       use fireteca
       use ignite
      Implicit none
      integer,intent(in)::i,j,k
      real:: energy,capqterm
      real::fi,frho,fox,frhovapor

        frho=0.0
        fox=0.5*foxb(i,j,k)*dt
        if(irhovapor.eq.1)
     +    frhovapor=0.5*frhovaporb(i,j,k)*dt
        energy=firad(i,j,k)
        pr(i,j,k)=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
        capqterm=energy*(1.e5/pr(i,j,k))**(rg/cp)/cp !rrl
        fi=fib(i,j,k)*dt+2.*capqterm*dt
        if (iheatsource.ge.1.and.((time+restarttime).ge.startigntime).and.k.le.lHeatSource)
     +         call heatSourceRHS(i,j,k,fi,frho,dt)

        xvb(i,j,k,1)=xvb(i,j,k,1)+0.5*f1avg(i,j,k)*dt
        xvb(i,j,k,2)=xvb(i,j,k,2)+0.5*f2avg(i,j,k)*dt
        xvb(i,j,k,3)=xvb(i,j,k,3)+0.5*f3avg(i,j,k)*dt
        xvb(i,j,k,4)=xvb(i,j,k,4)+0.5*fi
        xvb(i,j,k,5)=xvb(i,j,k,5)+0.5*fka(i,j,k)*dt
        xvb(i,j,k,6)=xvb(i,j,k,6)+0.5*fkb(i,j,k)*dt
        xvb(i,j,k,7)=xvb(i,j,k,7)+fox
        !xvb(i,j,k,nv)=xvb(i,j,k,nv)+frho
        if(irhovapor.eq.1)
     +    xvb(i,j,k,8)=xvb(i,j,k,8)+frhovapor

      end subroutine gas_ltp


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ignite_pattern compute ignition terms on large time step   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
      subroutine ignite_pattern       
      use constants
      use metryic
      use turba
      use pres
      use gridsetup
      use fireteca
      use msga
      use ignite
       Implicit none
      integer :: i,j,k
      integer::ia,ja

!  JLW 05/09/13   aerial ignition
      integer :: iaerial,ii,jj
!  JLW 05/09/13   atv/driptorch ignition
      integer :: iatv


      !real,external :: zcart
      real ::tempsold
      real :: rignitionperiod,travelrate,pretime
      real :: finterp,finterp1,finterp2,ignduration,xcell,ycell,x1,y1,x2,y2,dist,dist1,dist2,ps1,ps2
       
c       BELOW THIS LINE : IGNITION ONLY       
       ignduration=10.0   ! FP: parameter used for ignition type 7 only
       if(igntype==2.and.((time+restarttime).lt.endigntime))then     !terratorch style ignition
        rignitionperiod=endigntime-startigntime
        travelrate=sqrt((xfirelinelow-xfirelinehigh)**2
     &             +(yfirelinelow-yfirelinehigh)**2)/rignitionperiod
        pretime=flamedistance/travelrate
        if(pretime.gt.startigntime)then
         startigntime=pretime
         endigntime=pretime+(rignitionperiod)
        endif
        if((time+restarttime).ge.startigntime-pretime.and.
     &     (time+restarttime).le.endigntime+pretime) then
         finterp=((time+restarttime)-startigntime)/(endigntime-startigntime)
         do k=1,l
          do j=1,mp
           do i=1,np
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            x1=(float(ia)-.5)*dx
            y1=(float(ja)-.5)*dy
            x2=(xfirelinehigh-xfirelinelow)*finterp+xfirelinelow
            y2=(yfirelinehigh-yfirelinelow)*finterp+yfirelinelow
            if(y1.ge.yfirelinelow.and.y1.le.yfirelinehigh.and.k.lt.4.)then
             dist=sqrt((x2-x1)**2+(y2-y1)**2)
             tempsold=temps(i,j,k)
             temps(i,j,k)=max(temps(i,j,k),min(targettemp,
     &                          targettemp-targettemp/flamedistance*dist))
             if(temps(i,j,k).ne.tempsold)then
              cpsolid(i,j,k)=cpwood
              sies(i,j,k)=cpsolid(i,j,k)*temps(i,j,k)
              rhowater(i,j,k)=0.0
              rmoist(i,j,k)=0.0
              rhos(i,j,k)=rhof(i,j,k)
             endif
            endif
           enddo
          enddo
         enddo
        endif !end if (time.ge.startigntime-pretime.and.time.le.endigntime+pretime)

!correction JLW: test on irst added (here added by JLD 12/06/07
       elseif((igntype.eq.1).and.((time+restarttime).ge.startigntime).and.  ! line ignition
     &        ((time+restarttime).le.(startigntime+
     &                         (targettemp-tambient)/ramprate)))then
        do k=1,l
         do j=1,mp
          do i=1,np
           if(ifirestart(i,j,k).eq.1)then
            !if(rhof(i,j,k).GT.min_rhof .and. k.LT.3) then  ! KOO 
            !FP092019 put igniteVerticalExtent as a parameter in
            !gridlist
            if (igniteVerticalExtent.eq.0.or.k.le.igniteVerticalExtent) then
            if(rhof(i,j,k).GT.min_rhof) then   
            temps(i,j,k)=max(temps(i,j,k),min(((time+restarttime)-startigntime)
     &                      *(ramprate)+tambient,targettemp))
            sies(i,j,k)=cpsolid(i,j,k)*temps(i,j,k)
!            print *,'JLW ignite',i,j,k,
!     &          time+restarttime,temps(i,j,1),sies(i,j,1),rhos(i,j,1)
            endif ! rhof(i,j,k).GT.min_rhof
            endif
           endif
          enddo
         enddo
        enddo

! JLW 05/09/13 
       elseif(igntype==4) then                                ! aerial ignition
        if (time+restarttime.le.igntime(naerial)+
     &                      (targettemp-tambient)/ramprate) then
         do iaerial=1,naerial
           if ((time+restarttime.ge.igntime(iaerial)).and.
     &         (time+restarttime.le.igntime(iaerial)+
     &                            (targettemp-tambient)/ramprate)) then
             do j=1,mp
               do i=1,np
                 ia=(npos-1)*np+i
                 ja=(mpos-1)*mp+j
                 if(ia.eq.ignloc(iaerial,1).and.
     &              ja.eq.ignloc(iaerial,2)) then
                   do jj=j,j+1
                     do ii=i,i+1
                       temps(ii,jj,1)=max(temps(ii,jj,1),
     &                                 min(((time+restarttime)-
     &                           igntime(iaerial))*ramprate+tambient,
     &                                                    targettemp))
                       cpsolid(ii,jj,1)=cpwood
                       sies(ii,jj,1)=cpsolid(ii,jj,1)*temps(ii,jj,1)
                       rhowater(ii,jj,1)=0.0
                       rmoist(ii,jj,1)=0.0
                       rhos(ii,jj,1)=rhof(ii,jj,1)
                     enddo
                   enddo
                 endif
               enddo
             enddo
           endif
         enddo
        endif

! JLW 05/09/13  
       elseif(igntype==5) then                         ! atv/driptorch ignition
        if (time+restarttime.lt.atvtime(natv,2)) then
         do iatv=1,natv
          if ((time+restarttime.ge.atvtime(iatv,1)).and.
     &             (time+restarttime.le.atvtime(iatv,2))) then
            finterp=((time+restarttime)-atvtime(iatv,1))/
     &                      (atvtime(iatv,2)-atvtime(iatv,1))
            do j=1,mp
              do i=1,np
                ia=(npos-1)*np+i
                ja=(mpos-1)*mp+j
                x1=(float(ia)-.5)*dx
                y1=(float(ja)-.5)*dy
                x2=(atvstop(iatv,1)-atvstart(iatv,1))*finterp+
     &                                          atvstart(iatv,1)
                y2=(atvstop(iatv,2)-atvstart(iatv,2))*finterp+
     &                                          atvstart(iatv,2)
                dist=sqrt((x2-x1)**2+(y2-y1)**2)
                tempsold=temps(i,j,1)
                temps(i,j,1)=max(temps(i,j,1),min(targettemp,
     &                    targettemp-targettemp/flamedistance*dist))
                if(temps(i,j,1).ne.tempsold)then
c       print *,'atv ',time,atvtime(iatv,1),x1,y1,x2,y2,dist,
c    &                temps(i,j,1),tempsold,finterp,i,j,
c    &                mpi_rank,targettemp,flamedistance
                  cpsolid(i,j,1)=cpwood
                  sies(i,j,1)=cpsolid(i,j,1)*temps(i,j,1)
                  rhowater(i,j,1)=0.0
                  rmoist(i,j,1)=0.0
                  rhos(i,j,1)=rhof(i,j,1)
                endif
              enddo
            enddo
          endif
         enddo
        endif
        !JLD two line igntion pattern 26/09/2013
        elseif(igntype==6)then
        if ((time+restarttime).lt.endigntime1)then     !two lines ignition pattern
        ! line 1!
        rignitionperiod=endigntime1-startigntime1
        travelrate=sqrt((xlineor1-xlineend1)**2
     &             +(ylineor1-ylineend1)**2)/rignitionperiod
        pretime=flamedistance/travelrate
        if(pretime.gt.startigntime1)then
         startigntime1=pretime
         endigntime1=pretime+(rignitionperiod)
        endif
        if((time+restarttime).ge.startigntime1-pretime.and.
     &     (time+restarttime).le.endigntime1+pretime) then
         finterp=((time+restarttime)-startigntime1)/(endigntime1-startigntime1)
         do k=1,l
          do j=1,mp
           do i=1,np
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            x1=(float(ia)-.5)*dx
            y1=(float(ja)-.5)*dy
            x2=(xlineend1-xlineor1)*finterp+xlineor1
            y2=(ylineend1-ylineor1)*finterp+ylineor1
            if(y1.ge.min(ylineor1,ylineend1).and.y1.le.max(ylineor1,ylineend1).and.k.lt.2)then
             dist=sqrt((x2-x1)**2+(y2-y1)**2)
             tempsold=temps(i,j,k)
             temps(i,j,k)=max(temps(i,j,k),min(targettemp,
     & (targettemp-tambient)*(1.-dist/flamedistance)+tambient))
             if(temps(i,j,k).ne.tempsold)then
              cpsolid(i,j,k)=cpwood
              sies(i,j,k)=cpsolid(i,j,k)*temps(i,j,k)
              rhowater(i,j,k)=0.0
              rmoist(i,j,k)=0.0
              rhos(i,j,k)=rhof(i,j,k)
             endif
            endif
           enddo
          enddo
         enddo
        endif !end if (time.ge.startigntime1-pretime.and.time.le.endigntime1+pretime)
        endif ! end if time lt endigntime1
        if ((time+restarttime).lt.endigntime2)then     !two lines ignition pattern
        ! line 2!
        rignitionperiod=endigntime2-startigntime2
        travelrate=sqrt((xlineor2-xlineend2)**2
     &             +(ylineor2-ylineend2)**2)/rignitionperiod
        pretime=flamedistance/travelrate
        if(pretime.gt.startigntime2)then
         startigntime2=pretime
         endigntime2=pretime+(rignitionperiod)
        endif
        if((time+restarttime).ge.startigntime2-pretime.and.
     &     (time+restarttime).le.endigntime2+pretime) then
         finterp=((time+restarttime)-startigntime2)/(endigntime2-startigntime2)
         do k=1,l
          do j=1,mp
           do i=1,np
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            x1=(float(ia)-.5)*dx
            y1=(float(ja)-.5)*dy
            x2=(xlineend2-xlineor2)*finterp+xlineor2
            y2=(ylineend2-ylineor2)*finterp+ylineor2
            if(y1.ge.min(ylineor2,ylineend2).and.y1.le.max(ylineor2,ylineend2).and.k.lt.2)then 
            dist=sqrt((x2-x1)**2+(y2-y1)**2)
             tempsold=temps(i,j,k)
             
             temps(i,j,k)=max(temps(i,j,k),min(targettemp,
     & (targettemp-tambient)*(1.-dist/flamedistance)+tambient))
             if(temps(i,j,k).ne.tempsold)then
              cpsolid(i,j,k)=cpwood
              sies(i,j,k)=cpsolid(i,j,k)*temps(i,j,k)
              rhowater(i,j,k)=0.0
              rmoist(i,j,k)=0.0
              rhos(i,j,k)=rhof(i,j,k)
             endif
            endif
           enddo
          enddo
         enddo
        endif !end if (time.ge.startigntime1-pretime.and.time.le.endigntime1+pretime)
        endif !end if
       else if (igntype==7.and.((time+restarttime).ge.startigntime).and.((time+restarttime).lt.(endigntime+ignduration)))then
        !modification of terratorch 2 style ignition : pimont 2014. Now ignition has a
        !duration ('ignduration' seconds to match what we do with igntype=1)
        rignitionperiod=endigntime-startigntime
        travelrate=sqrt((xfirelinelow-xfirelinehigh)**2
     &             +(yfirelinelow-yfirelinehigh)**2)/rignitionperiod
        finterp2=min(((time+restarttime)-startigntime)/(endigntime-startigntime),1.0)
        finterp1=max(((time+restarttime)-ignduration-startigntime)/(endigntime-startigntime),0.0)
        !ignition line at current time)
        x1=(xfirelinehigh-xfirelinelow)*finterp1+xfirelinelow
        y1=(yfirelinehigh-yfirelinelow)*finterp1+yfirelinelow
        x2=(xfirelinehigh-xfirelinelow)*finterp2+xfirelinelow
        y2=(yfirelinehigh-yfirelinelow)*finterp2+yfirelinelow
       ! if (mpi_rank.eq.0)  write(6,*) 'ignition line',x1,y1,x2,y2
         do k=1,3
          do j=1,mp
           do i=1,np
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            xcell=(float(ia)-.5)*dx
            ycell=(float(ja)-.5)*dy
            ! distance between cell and p1(x1,y1)
            dist1=sqrt((xcell-x1)**2+(ycell-y1)**2)
            ! distance between cell and p2(x2,y2)
            dist2=sqrt((xcell-x2)**2+(ycell-y2)**2)
            ! scalar product P1cell.P1P2
            ps1=(xcell-x1)*(x2-x1)+(ycell-y1)*(y2-y1)
            ! scalar product P2cell.P2P1
            ps2=(xcell-x2)*(x1-x2)+(ycell-y2)*(y1-y2)
            ! distance between cell and ignition line
            dist=abs((xcell-x1)*(y2-y1)-(ycell-y1)*(x2-x1))/
     +           sqrt((x2-x1)**2+(y2-y1)**2)
            if ((ps1>=0.and.ps2>=0.and.dist<=flamedistance).or.dist1<=flamedistance.or.dist2<=flamedistance) then
            if (temps(i,j,k)<targettemp) then
        !     write(6,*) ia,ja,temps(i,j,k),dist
             temps(i,j,k)=temps(i,j,k)+(targettemp-tambient)*0.01
              cpsolid(i,j,k)=cpwood
              sies(i,j,k)=cpsolid(i,j,k)*temps(i,j,k)
              rhowater(i,j,k)=0.0
              rmoist(i,j,k)=0.0
              rhos(i,j,k)=rhof(i,j,k)
            endif
            endif
           enddo
          enddo
         enddo


       endif    !end if different ignition styles
      end subroutine ignite_pattern         

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!      subroutine computexvfuel(i,j,k)
!      Implicit none
!      integer::i,j,k


c                 if (actualfueldepth(i,j,k).lt.0.or.k.gt.1) then
c         if(actualfueldepth(i,j,k).lt.0.or.k.lt.1)then
c          xvfuel(i,j,k,1)=xvb(i,j,k,1)
c          xvfuel(i,j,k,2)=xvb(i,j,k,2)
c          xvfuel(i,j,k,3)=xvb(i,j,k,3)
c          xvfuel(i,j,k,4)=xvb(i,j,k,5)
c          xvfuel(i,j,k,5)=xvb(i,j,k,6)
c         elseif(k.lt.1)then
c          ztopcell=zcart(zedge(k+1),i,j)-zs(i,j)
c          urhoa=xvb(i,j,k+1,1)
c          xvfuel(i,j,k,1)=xvb(i,j,k,1)*
c     &    sqrt(actualfueldepth(i,j,k)/ztopcell)
c          vrhoa=xvb(i,j,k+1,2)
c          if(abs(xvb(i,j,k+1,2)).le.abs(xvb(i,j,k,2)).or.
c     &       xvb(i,j,k,2)/xvb(i,j,k+1,2).le.0.)then
c           vrhoa=xvb(i,j,k,2)
c          endif
c          xvfuel(i,j,k,2)=xvb(i,j,k,2)*sqrt(actualfueldepth(i,j,k)/ztopcell)
c          xvfuel(i,j,k,3)=xvb(i,j,k,3)*sqrt(actualfueldepth(i,j,k)/ztopcell)
c          speedrho2=(xvb(i,j,k,1)**2+xvb(i,j,k,2)**2+xvb(i,j,k,3)**2)
c          speedfuelrho2=(xvfuel(i,j,k,1)**2+xvfuel(i,j,k,2)**2+xvfuel(i,j,k,3)**2)
c          xvfuel(i,j,k,4)=xvb(i,j,k,5)*speedfuelrho2/speedrho2
c          xvfuel(i,j,k,5)=xvb(i,j,k,6)*speedfuelrho2/speedrho2
c         endif
!      return

!      end subroutine computexvfuel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine radiationSink()
      use gridsetup
      use xve
      use fireteca
      use constants 
      use  pres
      use xvo
        Implicit none
        integer::i,j,k 
         do k=1,l
         do j=1,mp
         do i=1,np
             ! TODO : the update of tambientarray should be done in rinit or xevariation only
             tambientarray(i,j,k)=xe(i,j,k,4)/xe(i,j,k,nv)
     &                           *(pre(i,j,k)*1.0e-5)**(rg/cp)
             pr(i,j,k)=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
             tempg(i,j,k)=xvb(i,j,k,4)/xvb(i,j,k,nv)*(pr(i,j,k)*1.e-5)**(rg/cp)
             firad(i,j,k)=-rke*sigma/sqrt(dx*dy)*(tempg(i,j,k)**4
     &          -tambientarray(i,j,k)**4)
         enddo
         enddo
         enddo
        return
        end subroutine radiationSink

!iheatsource>=1
!This subroutine update the position of a heat source in the scene according to its ros
! htros, but also the fuel properties (burn fuel in the first cell k=1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine updateHeatSourcePosition()
       use gridsetup
       use fireteca
       use msga
       use turba
       Implicit none
       
       Integer :: i,j
       real::xa,ya
       xhsmax=min(dx*n,xhsmax+dts*nts*hsros)
       xhsmin=max(0.0,xhsmin+dts*nts*hsros)
       ! we force the source to be greater than 1 cell
       xhsmin=min(xhsmax-dx,xhsmin)
       ! here we update value of rhof in the heat source
       do j=1,mp
         do i=1,np
         xa=dx*((npos-1)*np+i)
         ya=dy*((mpos-1)*mp+j)
         if (xa.ge.xhsmin.and.xa.le.xhsmax
     +            .and.ya.ge.yhsmin.and.ya.le.yhsmax)
     +        rhof(i,j,1) = 0.01  ! burn fuel (drag is very small)
         enddo
       enddo  
      end subroutine updateHeatSourcePosition
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       heatSourceRHS modify fi, frho so that
c       - the fuel does not burn in k=1
c       - there is the appropriate enery realise in the heatsource
c       - if iheatsource=1, the canopy can not burn, no convection
c       - if iheatsource=2, the canopy can not burn, convection
c       - if iheatsource=3, the canopy can burn
c       - if ihsmass=1 there is also a source of mass in k.eq.1
       subroutine heatSourceRHS(i,j,k,fi,frho,dths)
       use gridsetup
       use fireteca
       use msga
       use constants
       use turba
       use xvo
       use pres
       use metryic
       Implicit none
       Integer,intent(in) :: i,j,k
       real,intent(inout)::fi,frho
       real,intent(in)::dths !time step
       Integer :: ia,ja
       real:: xhslen,yhslen ! dimension on heat source within cell i,j,k (<=dx,dy)
       real:: energy,capqterm,energyfrombyram,hc,hchar
       real, external ::zcart
       if (k.le.lHeatSource) then
         ia=(npos-1)*np+i
         ja=(mpos-1)*mp+j
         xhslen=max(0.0,min (ia*dx,xhsmax)-max((ia-1)*dx,xhsmin))  
         yhslen=max(0.0,min (ja*dy,yhsmax)-max((ja-1)*dy,yhsmin))  
         ! the heat source is uniform between xhsmin and xhsmax
         ! hsint is kW/m, energy is W/m3
         energyfrombyram = hsint*yhslen*xhslen /(xhsmax-xhsmin)*1000.0
     +    /(dx*dy*(zcart(zedge(lHeatSource+1),i,j)-zcart(zedge(1),i,j)))
         ! computation of energy really realised in gas phase
         hc=hf/rnfuel
         hchar=32e6
         energy=energyfrombyram !* (hc-hchar*rhohydrothresh)/hc
         pr(i,j,k)=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
         capqterm=energy*(1.e5/pr(i,j,k))**(rg/cp)/cp 
!         if (capqterm.ge.0.001) write(6,*),ia,ja,xhslen,energyfrombyram
!         if (ia.eq.75.and.ja.eq.75) then
!           write(6,*) xhslen,yhslen,hsint,
!     +         xhsmax-xhsmin,energyfrombyram,energy,capqterm
!           capqterm=0.0
!         endif
         fi=fib(i,j,k)*dths+2.*capqterm*dths
         convht(i,j,k)=0.0 ! for the fuel
         if (ihsmass.eq.0) then
            frho=0.0
         else
           !mass source from byram:
           frho=energyfrombyram/(hf/rnfuel)*dths
         endif !ihsmass.eq.0 
           ! if (frho.ne.0) write(6,*) ia,ja,frho
        else ! above plume (k.gt.lHeatSource)
         if (iheatsource.eq.1) then ! the canopy does not burn, no convection
             convht(i,j,k)=0.0  ! for the fuel
             fi=fib(i,j,k)*dths
         else if(iheatsource.eq.2) then
             energy=-convht(i,j,k)
             pr(i,j,k)=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
             capqterm=energy*(1.e5/pr(i,j,k))**(rg/cp)/cp 
             fi=fib(i,j,k)*dths+2.*capqterm*dths
          endif ! if iheatsource.ge.3=>canopy can burn
        end if
      
       end subroutine heatSourceRHS
