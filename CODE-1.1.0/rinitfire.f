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
      integer :: i,j,k,ia,ja,ifuelsegments,ifuelseg,ifueltype,iindex !kv
      real :: zl,zla,ztopcell,zbottomcell
      real :: actualfuelheight,actualfuelbottom
      real :: fueldistributionslope,actualgroundload
      real,external :: zcart 
      real :: xhsmintmp,xhsmaxtmp,yhsmintmp,yhsmaxtmp ! extension of heat source (iheatsource>=1)
c init oxygen
      xe(:,:,:,7)=0.21*xe(:,:,:,nv)

c init hydrocarbons for nonlocal chemistry
      if (inonlocal==1)  xe(:,:,:,8)=1.d-6
c init water vapor
      if(irhovapor.eq.1) xe(:,:,:,8)=0.0
      
      do k=1,l
         do j=1,mp
            do i=1,np
        rhomicro(i,j,k)=rhomicrovalue
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
         temps(i,j,k)=xe(i,j,k,4)*(pre(i,j,k)*1.e-5)**(rg/cp)
     +                  /xe(i,j,k,nv)
            enddo
         enddo
        enddo
      endif !irst.eq.0

      rmoist=0   ! wss change 12/18/00
      if (idirt.eq.1) rhodirt=0   ! rrl change 9/20/01
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
         call readio(rhof,92,l,1-ih,np+ih,1-ih,mp+ih)  

!         print *,' rhof ',(rhof(400,400,k),k=1,l)
   
        call readio(sizescale,93,l,1-ih,np+ih,1-ih,mp+ih) !FP
         call readio(rmoist,94,l,1-ih,np+ih,1-ih,mp+ih) 
         call readio(actualfueldepth,95,l,1-ih,np+ih,1-ih,mp+ih) 
         !print *,' rmoist ',(rmoist(5,5,k),k=1,l)

         ! default value of ss to prevent division by 0.

      call rmaxmin1(rhof,'rhof',1-ih,np+ih,1-ih,mp+ih,l)
      call rmaxmin1(sizescale,'szsc',1-ih,np+ih,1-ih,mp+ih,l)
      call rmaxmin1(rmoist,'rmst',1-ih,np+ih,1-ih,mp+ih,l)
      call rmaxmin1(actualfueldepth,'afdp',1-ih,np+ih,1-ih,mp+ih,l)  


         do k=1,l
          do j=1,mp
           do i=1,np
            if(sizescale(i,j,k).eq.0.0) sizescale(i,j,k)=ss 
            if(k.eq.1)then
              if(rhof(i,j,k).eq.0.0) actualfueldepth(i,j,k)=.05
             endif

           enddo
          enddo
         enddo
         close (92)
         close (93)   !FP
         close (94)
         close (95)

      endif ! end ivegread.eq.1
         

      if (ivegread.eq.0) then
       do k=1,l
        do j=1,mp
          do i=1,np
            if (ifuelinra.eq.0) then
            sizescale(i,j,k)=ss
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
                  rhomicro(i,j,k)=500.              !kg/m^3
                elseif (ifueltype.eq.20) then
                  actualfuelheight=16.0
                  actualfuelbottom= 0.0            ! inserting 0 for ground 
                  fueldistributionslope=0.         !pos numbers mean the load is higher at the top kg/m^3/m
                  actualgroundload=0.25 !1.0   !fuel load at the ground if a linear extrapolation was done
                  rhomicro(i,j,k)=500.              !kg/m^3
                elseif (ifueltype.eq.0.and.rhodirt(i,j,k).lt.10) then
                  actualfuelheight=.7
                  actualfuelbottom= 0.            ! inserting 0 for ground
                  fueldistributionslope=0.         !pos numbers mean the load is higher at the top kg/m^3/m
                  actualgroundload=.00001       !fuel load at the ground if a linear extrapolation was done
                  rhomicro(i,j,k)=500.              !kg/m^3

                endif
                if (actualfuelheight.lt.ztopcell.and.k.eq.1) then
                    actualfueldepth(i,j,k)=actualfuelheight
                endif
                if (zbottomcell.lt.actualfuelheight.and.
     &                   ztopcell.gt.actualfuelbottom ) then
                  rhof(i,j,k)=(min(actualfuelheight,ztopcell)
     &                 -amax1(actualfuelbottom,zbottomcell))/
     &                  (ztopcell-zbottomcell)
     &                 *(actualgroundload +
     &                 .5*fueldistributionslope*
     &                 (min(actualfuelheight,ztopcell)
     &                 +amax1(actualfuelbottom,zbottomcell)))
c                 if (i.eq.5.and.j.eq.5) write (*,*)'k=',k,
c    &                 zbottomcell,ztopcell,' rhof=',rhof(i,j,k)
                  rmoist(i,j,k)=.05
                  if (k.gt.2) rmoist(i,j,k)=.80

                endif
 11           continue
           else if (ifuelinra.eq.1) then
      if (ifp.ge.1)call defineFuelFP(i,j,k)
           endif  ! end ifuelinra.eq.1
          enddo   !enddo i
        enddo    !enddo j
       enddo     !enddo k
      endif                         !  end of ivegread loop
       !FPAZ adds these lines to be sure that actualfueldepth was properly defined and that a zone with just gas exist in cell k=1
          do j=1,mp
           do i=1,np
            zla=zcart(z(1),i,j)-zs(i,j)
            if (actualfueldepth(i,j,1).ge.zla.or.actualfueldepth(i,j,1).le.0.) actualfueldepth(i,j,1)=zla-0.01
           enddo
          enddo
 
      !set an appropriate fueldepth based on vegetation distribution!
      call set_fueldepth
      if (mpi_rank.eq.0) write(6,*) 'fueldepth,lfuel = ',fueldepth,lfuel 
      do k=1,lfuel
        do j=1,mp
          do i=1,np
            !zl=zcart(z(k),i,j)
            !zla=zcart(z(k),i,j)-zs(i,j)
            !ztopcell=zcart(zedge(k+1),i,j)-zs(i,j)
            !zbottomcell=zcart(zedge(k),i,j)-zs(i,j)
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
      
            !if (zbottomcell.le.fueldepth)  then
              
              rhof(i,j,k)=amax1(rhof(i,j,k),min_rhof)
              sizescale(i,j,k)=amax1(sizescale(i,j,k),ss)
              rhowater(i,j,k)=rmoist(i,j,k)*rhof(i,j,k)
              rhos(i,j,k)=rhof(i,j,k)*(1.+rmoist(i,j,k))
              if (idirt.eq.1) rhos(i,j,k)=rhodirt(i,j,k)+rhos(i,j,k)
              cpsolid(i,j,k)=(rhof(i,j,k)*cpwood
     &             +cpwater*rmoist(i,j,k)*rhof(i,j,k))/rhos(i,j,k)
              if (idirt.eq.1) cpsolid(i,j,k)=cpsolid(i,j,k)+
     +              rhodirt(i,j,k)*cpdirt/rhos(i,j,k)
              rhofinitial(i,j,k)=rhof(i,j,k)
              temps(i,j,k)=xe(i,j,k,4)*(pre(i,j,k)*1.e-5)**(rg/cp)
     &                  /xe(i,j,k,nv)
              if(rhos(i,j,k).ge.10.0) then 
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
         xvb(:,:,:,7)=0.21*xe(:,:,:,nv)
         if (inonlocal==1) xvb(:,:,:,8)=1.d-6
         if(irhovapor.eq.1) xvb(:,:,:,8)=0.0
      endif

c END OF DEFINITION COMMON TO irst=0 and 2
      else if (irst.eq.1) then
        call set_fueldepth
        do k=1,lfuel
          do j=1,mp
            do i=1,np
              sizescale(i,j,k)=amax1(sizescale(i,j,k),ss)
              rmoist(i,j,k)=rhowater(i,j,k)/rhof(i,j,k)
              !rhowater(i,j,k)=rmoist(i,j,k)*rhof(i,j,k)
              rhos(i,j,k)=rhof(i,j,k)*(1.+rmoist(i,j,k))
              if (idirt.eq.1) rhos(i,j,k)=rhos(i,j,k)+rhodirt(i,j,k)
              cpsolid(i,j,k)=(rhof(i,j,k)*cpwood
     &             +cpwater*rmoist(i,j,k)*rhof(i,j,k))/rhos(i,j,k)
              if (idirt.eq.1) cpsolid(i,j,k)=cpsolid(i,j,k)+
     +              rhodirt(i,j,k)*cpdirt/rhos(i,j,k)
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
              !if (rhof(i,j,k).gt.min_rhof)  then
               if(igntype.eq.1)then
                do iindex=1,numignpts 
                 if(ia.eq.ignpts(iindex,1).and.
     &              ja.eq.ignpts(iindex,2)) then
                 !  write(6,*) 'setting ifirestart to 1 at',ia,i,ja,j,k,mpi_rank 
       if (iheatsource.eq.0.and.rhof(i,j,k).gt.min_rhof) then ! default when ifp and ijld=0
                   ifirestart(i,j,k)=1
                   rmoist(i,j,k)=0.00
                   rhowater(i,j,k)=rmoist(i,j,k)*rhof(i,j,k)
                   rhos(i,j,k)=rhof(i,j,k)+rhowater(i,j,k)
                   cpsolid(i,j,k)=(rhof(i,j,k)*cpwood
     &                   +cpwater*rmoist(i,j,k)*rhof(i,j,k))   /rhos(i,j,k)
              if (idirt.eq.1) cpsolid(i,j,k)=cpsolid(i,j,k)+
     +              rhodirt(i,j,k)*cpdirt/rhos(i,j,k)
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
              sies(i,j,k)=temps(i,j,k)*cpsolid(i,j,k)
            
            !endif  !end if(zbottomcell.le.fueldepth)
          enddo
        enddo
      enddo
      if (iheatsource.ge.1) then !ifp=1 or ijld=1
        call mpi_allreduce(xhsmintmp,xhsmin,1,mpi_real,mpi_min,
     +                   mpi_comm_world,ierror)
        call mpi_allreduce(xhsmaxtmp,xhsmax,1,mpi_real,mpi_max,
     +                   mpi_comm_world,ierror)
        call mpi_allreduce(yhsmintmp,yhsmin,1,mpi_real,mpi_min,
     +                   mpi_comm_world,ierror)
        call mpi_allreduce(yhsmaxtmp,yhsmax,1,mpi_real,mpi_max,
     +                   mpi_comm_world,ierror)
       if (mpi_rank.eq.0) then
        write(6,*) 'the plume source extension is :[',xhsmin,xhsmax,']'
        write(6,*) '             by [',yhsmin,yhsmax,']'
       endif 
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c this subroutine is used by FP (ifp=1) to define hard coded fuel (depending on fuelnumber)
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
            if (fuelnumber.eq.0.and.k.le.11) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'typical maritime pine canopy (see Dupont et al 2011)'
            ! here mesh should be aa1=0.1, dz=15,rhomicro=700,cd=0.25
            actualfueldepth(i,j,1)=0.5
            rmoist(i,j,k)=1.0
            ia=(npos-1)*np+i
            iamin=101
            iamin=0
            sizescale(i,j,k)=2.0/4000.0
            if (k.eq.1) then
                rhof(i,j,k)=0.35*0.6  !drag has already an additive factor of 0.1 in dragm2
            else if (ia.ge.iamin.and.k.le.4) then
                rhof(i,j,k)=0.35*0.013
            else if (ia.ge.iamin.and.k.eq.5) then
                rhof(i,j,k)=0.35*0.005
            else if (ia.ge.iamin.and.k.eq.6) then
                rhof(i,j,k)=0.35*0.015
            else if (ia.ge.iamin.and.k.eq.7) then
                rhof(i,j,k)=0.35*0.027
            else if (ia.ge.iamin.and.k.eq.8) then
                rhof(i,j,k)=0.35*0.20
            else if (ia.ge.iamin.and.k.eq.9) then
                rhof(i,j,k)=0.35*0.37
            else if (ia.ge.iamin.and.k.eq.10) then
                rhof(i,j,k)=0.35*0.15
            end if
            else if (fuelnumber.eq.1.and.k.le.9) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'deciduous canopy (see Pimont et al 2009)'
            ! here mesh should be aa1=0.1, dz=15,rhomicro=700,cd=0.25
            actualfueldepth(i,j,1)=0.5
            rmoist(i,j,k)=1.0
            sizescale(i,j,k)=2.0/4000.0
            if (k.eq.1) then
                rhof(i,j,k)=0.35*0.0471  
            else if (k.eq.2) then
                rhof(i,j,k)=0.35*0.0540
            else if (k.eq.3) then
                rhof(i,j,k)=0.35*0.0688
            else if (k.eq.4) then
                rhof(i,j,k)=0.35*0.0943
            else if (k.eq.5) then
                rhof(i,j,k)=0.35*0.1323
            else if (k.eq.6) then
                rhof(i,j,k)=0.35*0.1466
            else if (k.eq.7) then
                rhof(i,j,k)=0.35*0.1468
            else if (k.eq.8) then
                rhof(i,j,k)=0.35*0.1127
            else if (k.eq.9) then
                rhof(i,j,k)=0.35*0.0115
            end if
            else if (fuelnumber.eq.2.and.k.le.9) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*)
     +          'wind tunnel edge (see Pimont et al 2009)'
            ! here mesh should be aa1=0.1, dz=15,rhomicro=700,cd=0.25,n=384,m=128
            actualfueldepth(i,j,1)=0.5
            rmoist(i,j,k)=1.0
            sizescale(i,j,k)=2.0/4000.0
            ia=(npos-1)*np+i
            if (k.eq.1.or.(k.le.7.and.(ia.le.53.or.ia.ge.251))) then  !hmax=13.25
                rhof(i,j,k)=0.35*0.0834
            end if

            else if (fuelnumber.eq.11.and.k.le.2) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
            else if (fuelnumber.eq.11.and.k.le.2) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'peypin d''aygues grasses'
            ! wind 3m/s at 6m, fire width 80m
            ! ros 0.57-0.8m/s
            actualfueldepth(i,j,1)=0.25
            rmoist(i,j,k)=.1
            sizescale(i,j,k)=2.0/4000.0
            rhof(i,j,1)=0.65/(zcart(zedge(2),i,j)- 
     +        zcart(zedge(1),i,j))  !0.65kg/m2
            else if (fuelnumber.eq.12.and.k.le.2) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'typical litter (catchpole, viguiere...)'
            ! wind 0m/s at 6m, fire width some m 
            ! ros 0.03m/s
            actualfueldepth(i,j,1)=0.15
            rmoist(i,j,k)=.1
            sizescale(i,j,k)=2.0/4000.0
            rhof(i,j,1)=1.0/(zcart(zedge(2),i,j)-
     +        zcart(zedge(1),i,j))  !1.0kg/m2
            else if (fuelnumber.eq.130.and.k.le.2) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'dense shrub (spain, vega)'
            ! wind 1.5m/s at 2m, fire width 50 m 
            ! ros 0.1m/s
            actualfueldepth(i,j,1)=0.5
            rmoist(i,j,k)=.6
            sizescale(i,j,k)=2.0/4000.0
            rhof(i,j,1)=2.4/(zcart(zedge(2),i,j)-
     +        zcart(zedge(1),i,j))  !2.4kg/m2
            else if (fuelnumber.eq.135.and.k.le.2) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'dense shrub (spain, vega)'
            ! wind 1.5m/s at 2m, fire width 50 m 
            ! ros 0.33m/s
            actualfueldepth(i,j,1)=0.5
            rmoist(i,j,k)=.6
            sizescale(i,j,k)=2.0/4000.0
            rhof(i,j,1)=1.8/(zcart(zedge(2),i,j)-
     +        zcart(zedge(1),i,j))  !1.8kg/m2
            else if (fuelnumber.eq.14.and.k.le.2) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'kermes oak garrigue, trou du ras france'
            ! wind 5.6m/s at 6m, fire width 10 m 
            ! ros 0.1m/s
            actualfueldepth(i,j,1)=0.4
            rmoist(i,j,k)=.7
            sizescale(i,j,k)=2.0/4000.0
            rhof(i,j,1)=0.8/(zcart(zedge(2),i,j)-
     +        zcart(zedge(1),i,j))  !0.8kg/m2
            else if (fuelnumber.eq.15.and.k.le.2) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'australian grassland, cheney'
            ! wind 3 and 6m/s at 2m, fire width 50 m 
            ! ros 0.7-0.8m/s for 3m/s, 1.8-2.7m/s for 6m/s
            actualfueldepth(i,j,1)=0.7
            rmoist(i,j,k)=.05
            sizescale(i,j,k)=2.0/4000.0
            rhof(i,j,1)=0.7/(zcart(zedge(2),i,j)-
     +        zcart(zedge(1),i,j))  !0.7kg/m2
            else if (fuelnumber.eq.16.and.k.le.2) then
            if (mpi_rank.eq.0.and.i.eq.1.and.j.eq.1.and.k.eq.1)
     +       write(6,*) 
     +          'bare soil'
            actualfueldepth(i,j,1)=0.05
            rmoist(i,j,k)=.1
            sizescale(i,j,k)=2.0/4000.0
            rhof(i,j,1)=0.1/(zcart(zedge(2),i,j)- 
     +        zcart(zedge(1),i,j))  !0.1kg/m2

            endif ! fuelnumber

      end  ! suboutine defineFuelFP
