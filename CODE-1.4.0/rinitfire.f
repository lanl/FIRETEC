      subroutine rinitfire()
      use metryic
      use gridsetup
      use fireteca
      use constants
      use turba
      use pres
      use xvi
      use xvo
      use xve
      use msga
      use ignite
      use nonlocal
      use io

      Implicit None

      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,ift,ia,ja,ifuelsegments,ifuelseg,ifueltype,iindex !kv
      real :: zl,zla,ztopcell,zbottomcell
      real :: actualfuelheight,actualfuelbottom
      real :: fueldistributionslope,actualgroundload
      real,external :: zcart 
      real :: xhsmintmp,xhsmaxtmp,yhsmintmp,yhsmaxtmp ! extension of heat source (iheatsource>=1)
      real ::rhovaporFrac
c init oxygen
      xe(:,:,:,7)=0.23*xe(:,:,:,nv)*(1.-specifichumidity)

c init hydrocarbons for nonlocal chemistry
      if (inonlocal==1)  xe(:,:,:,8)=1.d-6
c init water vapor
      if(irhovapor.eq.1) xe(:,:,:,8)=specifichumidity*xe(:,:,:,nv)
      
      do k=1,l
         do j=1,mp
            do i=1,np
              rhomicro(:,i,j,k)=rhomicrovalue
            enddo
         enddo
      enddo
      cpsolid=1.
      ! FP below this line things are not done in case of real restart (irst.eq.1) 
      if (irst.eq.0.or.irst.eq.2) then
      temps=0.
      rhof=0.
      rhos=0.0
      rhowater=0.
      if (irst.eq.0) then ! temps already defined if irst.eq.2
        do k=1,l
         do j=1,mp
            do i=1,np
         rhovaporFrac=real(irhovapor)*xe(i,j,k,8)/xe(i,j,k,nv)
         call updateGasThermalProperties(rhovaporFrac)
         temps(:,i,j,k)=xe(i,j,k,4)*(pre(i,j,k)*1.e-5)**(rg_over_cp_gas)
     +                  /xe(i,j,k,nv)
            enddo
         enddo
        enddo
      endif !irst.eq.0

      rmoist=0   ! wss change 12/18/00
      rhofinitial=0.0
      rhof=0.0

      if (ivegread.eq.1) then 
         if (mpi_rank.eq.0) then
           open (unit=92,file='treesrhof.dat',form='unformatted',
     +           status='old')
           open (unit=94,file='treesmoist.dat',form='unformatted',
     +           status='old')
           open (unit=95,file='treesfueldepth.dat',form='unformatted',
     +           status='old')
          open (unit=93,file='treesss.dat',form='unformatted',
     +           status='old')
         endif
          do ift=1,nfuel
            call readio(rhof(ift,:,:,:),92,l,1-ih,np+ih,1-ih,mp+ih)
            call readio(sizescale(ift,:,:,:),93,l,1-ih,np+ih,1-ih,mp+ih) !FP
            call readio(rmoist(ift,:,:,:),94,l,1-ih,np+ih,1-ih,mp+ih) 
            call readio(actualfueldepth(ift,:,:,:),95,l,1-ih,np+ih,1-ih,mp+ih) 

         do k=1,l
          do j=1,mp
           do i=1,np
                  if(sizescale(ift,i,j,k).eq.0.0) sizescale(ift,i,j,k)=ss 
            if(k.eq.1)then
                    if(rhof(ift,i,j,k).eq.0.0) actualfueldepth(ift,i,j,k)=.05
             endif
           enddo
          enddo
         enddo
          enddo
          call rmaxmin2(rhof,'rhof',1-ih,np+ih,1-ih,mp+ih,l)
         close (92)
          call rmaxmin2(sizescale,'szsc',1-ih,np+ih,1-ih,mp+ih,l)
         close (93)   !FP
          call rmaxmin2(rmoist,'rmst',1-ih,np+ih,1-ih,mp+ih,l)
         close (94)
          call rmaxmin2(actualfueldepth,'afdp',1-ih,np+ih,1-ih,mp+ih,l)
         close (95)
      endif ! end ivegread.eq.1
         
      if (ivegread.eq.0) then
       do k=1,l
        do j=1,mp
          do i=1,np
            if (ifuelinra.eq.0) then
                  sizescale(:,i,j,k)=ss
            zl=zcart(z(k),i,j)
            zla=zcart(z(k),i,j)-zs(i,j)
            ztopcell=zcart(zedge(k+1),i,j)-zs(i,j)
            zbottomcell=zcart(zedge(k),i,j)-zs(i,j)
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
          !Setting up a special case of vegetaion for JLW Las Vegas NM runs
              ifuelsegments=1
              do 11 ifuelseg=1,ifuelsegments
                if (ifuelseg.eq.1)    ifueltype=3
c                if (ifuelseg.eq.1.and.k.eq.1)    ifueltype=3
c                if (ifuelseg.eq.2.and.k.eq.2)    ifueltype=20

                if (ifueltype.eq.3) then
                  actualfuelheight=    0.7
                  actualfuelbottom= 0.0            ! inserting 0 for ground 
                  fueldistributionslope=0.         !pos numbers mean the load is higher at the top kg/m^3/m
                  actualgroundload= 1.0       !fuel load at the ground if a linear extrapolation was done
                  rhomicro(1,i,j,k)=500.              !kg/m^3
                elseif (ifueltype.eq.20) then
                  actualfuelheight=16.0
                  actualfuelbottom= 0.0            ! inserting 0 for ground 
                  fueldistributionslope=0.         !pos numbers mean the load is higher at the top kg/m^3/m
                  actualgroundload=0.25 !1.0   !fuel load at the ground if a linear extrapolation was done
                  rhomicro(1,i,j,k)=500.              !kg/m^3
                elseif (ifueltype.eq.0) then
                  actualfuelheight=.7
                  actualfuelbottom= 0.            ! inserting 0 for ground
                  fueldistributionslope=0.         !pos numbers mean the load is higher at the top kg/m^3/m
                  actualgroundload=.00001       !fuel load at the ground if a linear extrapolation was done
                  rhomicro(1,i,j,k)=500.              !kg/m^3

                endif
                if (actualfuelheight.lt.ztopcell.and.k.eq.1) then
                    actualfueldepth(1,i,j,k)=actualfuelheight
                endif
                    if (zbottomcell.lt.actualfuelheight.and.ztopcell.gt.actualfuelbottom ) then
                      rhof(ifuelseg,i,j,k)=(min(actualfuelheight,ztopcell)-amax1(actualfuelbottom,zbottomcell))/
     &                  (ztopcell-zbottomcell)*(actualgroundload+.5*fueldistributionslope*
     &                  (min(actualfuelheight,ztopcell)+amax1(actualfuelbottom,zbottomcell)))
c                 if (i.eq.5.and.j.eq.5) write (*,*)'k=',k,
c    &                 zbottomcell,ztopcell,' rhof=',rhof(i,j,k)
                  rmoist(1,i,j,k)=.05
                  if (k.gt.2) rmoist(1,i,j,k)=.80

                endif
 11           continue
                elseif (ifuelinra.eq.1) then
            call defineFuelFP(i,j,k)
           endif  ! end ifuelinra.eq.1
          enddo   !enddo i
        enddo    !enddo j
       enddo     !enddo k
      endif                         !  end of ivegread loop
       !FPAZ adds these lines to be sure that actualfueldepth was properly defined and that a zone with just gas exist in cell k=1
          do j=1,mp
           do i=1,np
            do ift=1,nfuel
            zla=zcart(z(1),i,j)-zs(i,j)
              if (actualfueldepth(ift,i,j,1).ge.zla.or.actualfueldepth(ift,i,j,1).le.0.) actualfueldepth(ift,i,j,1)=zla-0.01
            enddo
           enddo
          enddo
 
      !set an appropriate fueldepth based on vegetation distribution!
      call set_fueldepth
      if (mpi_rank.eq.0) write(6,*) 'fueldepth,lfuel = ',fueldepth,lfuel 
      do k=1,lfuel
        do j=1,mp
          do i=1,np
              do ift=1,nfuel
            !zl=zcart(z(k),i,j)
            !zla=zcart(z(k),i,j)-zs(i,j)
            !ztopcell=zcart(zedge(k+1),i,j)-zs(i,j)
            !zbottomcell=zcart(zedge(k),i,j)-zs(i,j)
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
      
            !if (zbottomcell.le.fueldepth)  then
                rhof(ift,i,j,k)=amax1(rhof(ift,i,j,k),min_rhof/nfuel)
                sizescale(ift,i,j,k)=amax1(sizescale(ift,i,j,k),ss)
                rhowater(ift,i,j,k)=rmoist(ift,i,j,k)*rhof(ift,i,j,k)
                rhos(ift,i,j,k)=rhof(ift,i,j,k)*(1.+rmoist(ift,i,j,k))
                cpsolid(ift,i,j,k)=(rhof(ift,i,j,k)*cpwood
     &             +cpwater*rmoist(ift,i,j,k)*rhof(ift,i,j,k))/rhos(ift,i,j,k)
                rhofinitial(ift,i,j,k)=rhof(ift,i,j,k)
              ! temps =tempg here
              rhovaporFrac=real(irhovapor)*xe(i,j,k,8)/xe(i,j,k,nv)
              call updateGasThermalProperties(rhovaporFrac)
              temps(ift,i,j,k)=xe(i,j,k,4)
     &             *(pre(i,j,k)*1.e-5)**(rg_over_cp_gas)
     &                  /xe(i,j,k,nv)
              enddo
              if(sum(rhos(:,i,j,k)).ge.10.0) then 
               xe(i,j,k,1)=0.0
               xe(i,j,k,2)=0.0
               xv(i,j,k,1)=0.0
               xv(i,j,k,2)=0.0                
               xvb(i,j,k,1)=0.0
               xvb(i,j,k,2)=0.0
               xe(i,j,k+1,1)=0.5*xe(i,j,k+1,1)
               xv(i,j,k+1,1)=0.5*xv(i,j,k+1,1)
               xvb(i,j,k+1,1)=0.5*xvb(i,j,k+1,1)
               xe(i,j,k+1,2)=0.5*xe(i,j,k+1,2)
               xv(i,j,k+1,2)=0.5*xv(i,j,k+1,2)
               xvb(i,j,k+1,2)=0.5*xvb(i,j,k+1,2)
              endif
           enddo
         enddo
       enddo
      
      if (irst.eq.0) then ! already defined if irst.eq.2
         xvb(:,:,:,7)=0.23*xe(:,:,:,nv)
         if (inonlocal==1) xvb(:,:,:,8)=1.d-6
         if(irhovapor.eq.1) xvb(:,:,:,8)=0.0
      endif

c END OF DEFINITION COMMON TO irst=0 and 2
      else if (irst.eq.1) then
        call set_fueldepth
        do k=1,lfuel
          do j=1,mp
            do i=1,np
              do ift=1,nfuel
              sizescale(ift,i,j,k)=amax1(sizescale(ift,i,j,k),ss)
              rmoist(ift,i,j,k)=rhowater(ift,i,j,k)/rhof(ift,i,j,k)
              !rhowater(i,j,k)=rmoist(i,j,k)*rhof(i,j,k)
              rhos(ift,i,j,k)=rhof(ift,i,j,k)*(1.+rmoist(ift,i,j,k))
              cpsolid(ift,i,j,k)=(rhof(ift,i,j,k)*cpwood
     &             +cpwater*rmoist(ift,i,j,k)*rhof(ift,i,j,k))/rhos(ift,i,j,k)
              enddo
            enddo
          enddo
        enddo
        endif !irst
       
c     BELOW THIS LINE IGNITION PATTERN: all irst      
      if (iheatsource.ge.1) then !ifp=1 or ijld=1
         xhsmintmp=dx*n
         xhsmaxtmp=0.0
         yhsmintmp=dy*m
         yhsmaxtmp=0.0
      endif
      do k=1,l
       do j=1,mp
         do i=1,np
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            do ift=1,nfuel
            !if (rhof(nfuel,i,j,k).gt.min_rhof)  then
               if(igntype.eq.1)then
                do iindex=1,numignpts 
                  if(ia.eq.ignpts(iindex,1).and.ja.eq.ignpts(iindex,2)) then
                 !  write(6,*) 'setting ifirestart to 1 at',ia,i,ja,j,k,mpi_rank 
                    if (iheatsource.eq.0.and.rhof(ift,i,j,k).gt.min_rhof/nfuel) then ! default when ifp and ijld=0
                   ifirestart(i,j,k)=1
                   rmoist(ift,i,j,k)=0.00
                   rhowater(ift,i,j,k)=rmoist(ift,i,j,k)*rhof(ift,i,j,k)
                   rhos(ift,i,j,k)=rhof(ift,i,j,k)+rhowater(ift,i,j,k)
                   cpsolid(ift,i,j,k)=(rhof(ift,i,j,k)*cpwood
     &                   +cpwater*rmoist(ift,i,j,k)*rhof(ift,i,j,k))/rhos(ift,i,j,k)
       else if (iheatsource.ge.1) then !only if ifp=1 or ijld=1
        ! definition of heatsource extension
                 xhsmintmp=min(xhsmintmp,(ia-1)*dx+hsros*max(0,ittot-1)*dt)
                 xhsmaxtmp=max(xhsmaxtmp,ia*dx)
                 yhsmintmp=min(yhsmintmp,(ja-1)*dy)
                 yhsmaxtmp=max(yhsmaxtmp,ja*dy)
               endif  ! endif iheatsource.ge.1
                 endif   ! endif ignnpts
                enddo
              endif
              sies(ift,i,j,k)=temps(ift,i,j,k)*cpsolid(ift,i,j,k)
            enddo
            !endif  !end if(zbottomcell.le.fueldepth)
          enddo
        enddo
      enddo
      if (iheatsource.ge.1) then !ifp=1 or ijld=1
        call mpi_allreduce(xhsmintmp,xhsmin,1,mpi_real,mpi_min,mpi_comm_world,ierror)
        call mpi_allreduce(xhsmaxtmp,xhsmax,1,mpi_real,mpi_max,mpi_comm_world,ierror)
        call mpi_allreduce(yhsmintmp,yhsmin,1,mpi_real,mpi_min,mpi_comm_world,ierror)
        call mpi_allreduce(yhsmaxtmp,yhsmax,1,mpi_real,mpi_max,mpi_comm_world,ierror)
       if (mpi_rank.eq.0) then
          write(6,*)'the plume source extension is :[',xhsmin,xhsmax,']'
          write(6,*)'by [',yhsmin,yhsmax,']'
       endif 
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c this subroutine is used by FP to define hard coded fuel (depending on fuelinranumber)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defineFuelFP(i,j,k)
      use metryic
      use gridsetup
      use turba
      use fireteca
      use msga
      Implicit none
      integer::i,j,k,ia,iamin
      real, external:: zcart
            if (fuelinranumber.eq.0.and.k.le.11) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*) 
     +          'typical maritime pine canopy (see Dupont et al 2011)'
            ! here mesh should be aa1=0.1, dz=15,rhomicro=700,cd=0.25
        actualfueldepth(nfuel,i,j,1)=0.5
        rmoist(nfuel,i,j,k)=1.0
            ia=(npos-1)*np+i
            iamin=101  ! edge
            iamin=0   ! homogeneous forest
        sizescale(nfuel,i,j,k)=2.0/4000.0
            if (k.eq.1) then
          rhof(nfuel,i,j,k)=0.35*0.6  !drag has already an additive factor of 0.1 in dragm2
            else if (ia.ge.iamin.and.k.le.4) then
          rhof(nfuel,i,j,k)=0.35*0.013
            else if (ia.ge.iamin.and.k.eq.5) then
          rhof(nfuel,i,j,k)=0.35*0.005
            else if (ia.ge.iamin.and.k.eq.6) then
          rhof(nfuel,i,j,k)=0.35*0.015
            else if (ia.ge.iamin.and.k.eq.7) then
          rhof(nfuel,i,j,k)=0.35*0.027
            else if (ia.ge.iamin.and.k.eq.8) then
          rhof(nfuel,i,j,k)=0.35*0.20
            else if (ia.ge.iamin.and.k.eq.9) then
          rhof(nfuel,i,j,k)=0.35*0.37
            else if (ia.ge.iamin.and.k.eq.10) then
          rhof(nfuel,i,j,k)=0.35*0.15
            end if
            else if (fuelinranumber.eq.1.and.k.le.9) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*) 
     +          'deciduous canopy (see Pimont et al 2009)'
            ! here mesh should be aa1=0.1, dz=15,rhomicro=700,cd=0.25
        actualfueldepth(nfuel,i,j,1)=0.5
        rmoist(nfuel,i,j,k)=1.0
        sizescale(nfuel,i,j,k)=2.0/4000.0
            if (k.eq.1) then
          rhof(nfuel,i,j,k)=0.35*0.0471  
            else if (k.eq.2) then
          rhof(nfuel,i,j,k)=0.35*0.0540
            else if (k.eq.3) then
          rhof(nfuel,i,j,k)=0.35*0.0688
            else if (k.eq.4) then
          rhof(nfuel,i,j,k)=0.35*0.0943
            else if (k.eq.5) then
          rhof(nfuel,i,j,k)=0.35*0.1323
            else if (k.eq.6) then
          rhof(nfuel,i,j,k)=0.35*0.1466
            else if (k.eq.7) then
          rhof(nfuel,i,j,k)=0.35*0.1468
            else if (k.eq.8) then
          rhof(nfuel,i,j,k)=0.35*0.1127
            else if (k.eq.9) then
          rhof(nfuel,i,j,k)=0.35*0.0115
            end if
            else if (fuelinranumber.eq.2.and.k.le.9) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)
     +          'wind tunnel edge (see Pimont et al 2009)'
            ! here mesh should be aa1=0.1, dz=15,rhomicro=700,cd=0.25,n=384,m=128
        actualfueldepth(nfuel,i,j,1)=0.5
        rmoist(nfuel,i,j,k)=1.0
        sizescale(nfuel,i,j,k)=2.0/4000.0
            ia=(npos-1)*np+i
            if (k.eq.1.or.(k.le.7.and.(ia.le.53.or.ia.ge.251))) then  !hmax=13.25
          rhof(nfuel,i,j,k)=0.35*0.0834
            end if
            else if (fuelinranumber.eq.11.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*) 
     +          'peypin d''aygues grasses'
            ! wind 3m/s at 6m, fire width 80m
            ! ros 0.57-0.8m/s
        actualfueldepth(nfuel,i,j,1)=0.25
        rmoist(nfuel,i,j,k)=.1
        sizescale(nfuel,i,j,k)=2.0/4000.0
        rhof(nfuel,i,j,1)=0.65/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))  !0.65kg/m2
            else if (fuelinranumber.eq.12.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)
     +          'typical litter (catchpole, viguiere...)'
            ! wind 0m/s at 6m, fire width some m 
            ! ros 0.03m/s
        actualfueldepth(nfuel,i,j,1)=0.15
        rmoist(nfuel,i,j,k)=.1
        sizescale(nfuel,i,j,k)=2.0/4000.0
        
            else if (fuelinranumber.eq.130.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)
     +          'dense shrub plot 3 (spain, vega)'
            !plot3, backfire tests (FIRE PARADOX) !
            !Be careful 80x 100 m plot with safety breaks and ignition
            !over 80 m '
            ! wind 1.5m/s at 2m;, oblique (selected at 45 deg  in the backfire
            ! study !, fire width 80 m 
            ! head fire ros 0.09m/s, backfire ros 0.026 m/s
        actualfueldepth(nfuel,i,j,1)=0.6
        rmoist(nfuel,i,j,k)=.63
        sizescale(nfuel,i,j,k)=2.0/4000.0
        rhof(nfuel,i,j,1)=2.4/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))  !2.4kg/m2
            else if (fuelinranumber.eq.135.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*) 
     +          'dense shrub plot1 (spain, vega)'
            !plot1, backfire tests (FIRE PARADOX) !
            !Be careful 60x 130 m plot with safety breaks and ignition
            !over 80 m '
            ! wind 4.2 m/s at 2m, fire width 60 m 
            ! ros 0.33 m/s 
        actualfueldepth(nfuel,i,j,1)=0.5
        rmoist(nfuel,i,j,k)=.61
        sizescale(nfuel,i,j,k)=2.0/4000.0
        rhof(nfuel,i,j,1)=1.7/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))  !1.7kg/m2
                else if (fuelinranumber.eq.137.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)
     +          'gorse edres 1  (spain, vega)'
            ! plot30 x 40 m plot with safety breaks
            ! wind 4.2 m/s at 2m, fire width 30 m
               ! wind 5.6 m /s at 6 m 
            ! ros 0.27 m/s 
            ! be careful speculate on fuel : 3 t/ha and 50% dead 
        actualfueldepth(nfuel,i,j,1)=1.25
        rmoist(nfuel,i,j,k)=.7
        sizescale(nfuel,i,j,k)=2.0/4000.0
        rhof(nfuel,i,j,1)=2./(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))  !2.0 kg/m2
            else if (fuelinranumber.eq.14.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*) 
     +          'kermes oak garrigue, trou du ras france'
            ! wind 5.6m/s at 6m, fire width 10 m 
            ! ros 0.1m/s
        actualfueldepth(nfuel,i,j,1)=0.4
        rmoist(nfuel,i,j,k)=.7
        sizescale(nfuel,i,j,k)=2.0/6000.0
        rhof(nfuel,i,j,1)=0.8/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))  !0.8kg/m2
            else if (fuelinranumber.eq.15.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*) 
     +          'australian grassland, cheney'
            ! wind 3 and 9m/s at 2 m, fire width 50 m 
            ! ros 0.7-0.8m/s for 3m/s, 1.8-2.7m/s for 9m/s
        actualfueldepth(nfuel,i,j,1)=0.4
        rmoist(nfuel,i,j,k)=.05
        sizescale(nfuel,i,j,k)=2.0/12000.0
        rhof(nfuel,i,j,1)=0.4/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))  !0.4kg/m2
            else if (fuelinranumber.eq.150.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)
     +          'australian grassland, original tests Rod'
            ! wind 3 and 9m/s at 2 m, fire width 50 m 
            ! ros 0.7-0.8m/s for 3m/s, 1.8-2.7m/s for 9m/s
        actualfueldepth(nfuel,i,j,1)=0.7
        rmoist(nfuel,i,j,k)=.05
        sizescale(nfuel,i,j,k)=2.0/4000.0
        rhof(nfuel,i,j,1)=0.7/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))  !0.7kg/m2
            else if (fuelinranumber.eq.16.and.k.le.2) then
        if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*) 
     +          'bare soil'
        actualfueldepth(nfuel,i,j,1)=0.05
        rmoist(nfuel,i,j,k)=.1
        sizescale(nfuel,i,j,k)=2.0/4000.0
        rhof(nfuel,i,j,1)=0.1/(zcart(zedge(2),i,j)-zcart(zedge(1),i,j))  !0.1kg/m2
            endif ! fuelinranumber

      end  ! suboutine defineFuelFP
