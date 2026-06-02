!-----------------------------------------------------------------------
! rinitrad initializes all arrays and variables for radiation 
! used throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitrad
  use gridlist_variables, only : irad,l,irandseed,iseed, &
    dx,dy,nfuel,nprocx,nprocy,iradeastflux,prec
  use gridsetup, only : np,mp
  use metric_variables, only : zedge
  use radiation_variables, only : volumeAZ,kfindex, &
    zCartEdgeAZGlobal,nphotbatch
  use fuel_variables, only : actualFuelDepth
  use msga_variables, only : mpi_rank,mpi_comm_world,numprocs,ierror
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use zcart_function
  Implicit None

  ! Local Variables
  integer :: i,j,k,it
  integer :: ia,ja,iprocx,jprocy,iproc,nsize
  real(prec),allocatable :: afdtmp(:)
  real(prec),allocatable :: zCartEdgeAZLocal(:,:,:,:)
  real(prec),allocatable :: zCartEdgeAZTmp(:,:,:,:)

  ! Executable Code
  !Allocate Monte Carlo  method radiation arrays
  call defineRadArray
  if(irad.eq.0) return
  if(iradeastflux.eq.1) nphotbatch=3.e5

  ! initialize the random number generator for Monte Carlo
  if(irandseed.eq.1) call system_clock(iseed)
  call random_seed(iseed)
  if(mpi_rank.eq.0) print*,'iseed=',iseed
           
  allocate(zCartEdgeAZLocal(4,np,mp,0:l+nfuel))
  zCartEdgeAZLocal = 0.0
  ! AZ  define zcartedge from top of the first cell (k=1) (real mesh) to the top 
  do k=1,l
    do j=1,mp
      do i=1,np
        zCartEdgeAZLocal(1,i,j,k+nfuel)=0.25* &
          (zcart(zedge(k+1),i-1,j-1)+zcart(zedge(k+1),i-1,j)+ &
           zcart(zedge(k+1),i,j-1)+zcart(zedge(k+1),i,j))  !top of k cell, x-, y-
        zCartEdgeAZLocal(2,i,j,k+nfuel)=0.25* &
          (zcart(zedge(k+1),i+1,j-1)+zcart(zedge(k+1),i+1,j)+ &
           zcart(zedge(k+1),i,j-1)+zcart(zedge(k+1),i,j))  !top of k cell, x+, y-
        zCartEdgeAZLocal(3,i,j,k+nfuel)=0.25* &
          (zcart(zedge(k+1),i-1,j+1)+zcart(zedge(k+1),i-1,j)+ &
           zcart(zedge(k+1),i,j+1)+zcart(zedge(k+1),i,j))  !top of k cell, x-, y+
        zCartEdgeAZLocal(4,i,j,k+nfuel)=0.25* &
          (zcart(zedge(k+1),i+1,j+1)+zcart(zedge(k+1),i+1,j)+ &
           zcart(zedge(k+1),i,j+1)+zcart(zedge(k+1),i,j))  !top of k cell, x+, y+
      enddo !i
    enddo !j
  enddo !k
  ! Ground Level AZ
  do j=1,mp
    do i=1,np
      zCartEdgeAZLocal(1,i,j,0)=0.25*(zcart(zedge(1),i-1,j-1)+ &
        zcart(zedge(1),i-1,j)+zcart(zedge(1),i,j-1)+zcart(zedge(1),i,j))  !top of k cell, x-, y-
      zCartEdgeAZLocal(2,i,j,0)=0.25*(zcart(zedge(1),i+1,j-1)+ &
        zcart(zedge(1),i+1,j)+zcart(zedge(1),i,j-1)+zcart(zedge(1),i,j))  !top of k cell, x+, y-
      zCartEdgeAZLocal(3,i,j,0)=0.25*(zcart(zedge(1),i-1,j+1)+ &
        zcart(zedge(1),i-1,j)+zcart(zedge(1),i,j+1)+zcart(zedge(1),i,j))  !top of k cell, x-, y+
      zCartEdgeAZLocal(4,i,j,0)=0.25*(zcart(zedge(1),i+1,j+1)+ & 
        zcart(zedge(1),i+1,j)+zcart(zedge(1),i,j+1)+zcart(zedge(1),i,j))  !top of k cell, x+, y+
    enddo !i
  enddo !j

  ! sort individual fuel zones by afd
  allocate(afdtmp(nfuel)); afdtmp = 0.0
  do j=1,mp
    do i=1,np
      afdtmp=actualFuelDepth(:,i,j)
      do k=1,nfuel
        kfindex(i,j,minloc(afdtmp,dim=1))=k
        afdtmp(minloc(afdtmp,dim=1))=maxval(actualFuelDepth(:,i,j))+1
      enddo
    enddo
  enddo 
  deallocate(afdtmp)
        
  ! AZ for individual fuel zones in bottom cell
  do k=1,nfuel
    do j=1,mp
      do i=1,np
        do it=1,4
#ifdef DBL_PREC
          zCartEdgeAZLocal(it,i,j,k)=max(zCartEdgeAZLocal(it,i,j,k-1)+ &
            1.d-4,min(zCartEdgeAZLocal(it,i,j,nfuel+1)- &
            1.d-4*(1+nfuel-k),zCartEdgeAZLocal(it,i,j,0)+ &
            actualFuelDepth(kfindex(i,j,k),i,j)))
#else
          zCartEdgeAZLocal(it,i,j,k)=max(zCartEdgeAZLocal(it,i,j,k-1)+ &
            1.e-4,min(zCartEdgeAZLocal(it,i,j,nfuel+1)- &
            1.e-4*(1+nfuel-k),zCartEdgeAZLocal(it,i,j,0)+ &
            actualFuelDepth(kfindex(i,j,k),i,j)))
#endif
        enddo !it 
      enddo !i
    enddo !j
  enddo !k

  ! AZ  volume
  do k=1,l+nfuel
    do j=1,mp
      do i=1,np
        volumeAZ(i,j,k)=dx*dy*.25*(sum(zCartEdgeAZLocal(:,i,j,k))- &
          sum(zCartEdgeAZLocal(:,i,j,k-1)))
      enddo !i
    enddo !j
  enddo !k

  ! this section to build zCartEdgeAZGlobal (n,m,l+nfuel) can probably be improved   
  allocate(zCartEdgeAZTmp(np,mp,0:l+nfuel,numprocs))
  zCartEdgeAZTmp = 0.0
  do it=1,4
    nsize=np*mp*(l+nfuel+1)
#ifdef DBL_PREC
    call mpi_allgather(zCartEdgeAZLocal(it,:,:,:),nsize, &
      mpi_double_precision,zCartEdgeAZTmp,nsize, &
      mpi_double_precision,mpi_comm_world,ierror)
#else
    call mpi_allgather(zCartEdgeAZLocal(it,:,:,:),nsize,mpi_real, &
      zCartEdgeAZTmp,nsize,mpi_real,mpi_comm_world,ierror)
#endif
    do iprocx=1,nprocx
      do jprocy=1,nprocy
        iproc=1+(iprocx-1)+(jprocy-1)*nprocx
        do k=0,l+nfuel
          do j=1,mp
            do i=1,np
              ia=(iprocx-1)*np+i
              ja=(jprocy-1)*mp+j
              zCartEdgeAZGlobal(it,ia,ja,k)=zCartEdgeAZTmp(i,j,k,iproc)
            enddo !i
          enddo !j
        enddo !k
      enddo !jprocy
    enddo !iprocx
  enddo !it
  deallocate(zCartEdgeAZLocal,zCartEdgeAZTmp)

end subroutine rinitrad

!-----------------------------------------------------------------------
! radiationSink estimates a heat loss due to radiation based on 
! comparing local temperatures against environmental temperatures
!-----------------------------------------------------------------------
subroutine radiationSink(force,il,iu,jl,ju,lls)
  use gridlist_variables, only : dx,dy,prec
  use thermo_variables, only : pre,pr,tempg,cp_gas,cv_gas
  use constants, only : sigma,pref
  use xvall, only : xe,itemp,irho
  Implicit none
  
  ! Local Variables  
  integer,intent(in) :: il,iu,jl,ju,lls
  real(prec),intent(inout) :: force(il:iu,jl:ju,lls)
  integer :: i,j,k 
  real(prec) :: tenv,firad

  ! Executable Code
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        tenv=xe(i,j,k,itemp)/xe(i,j,k,irho) &
          *(pre(i,j,k)/pref)**(1.-cv_gas(i,j,k)/cp_gas(i,j,k))
        firad=-2.*sigma/sqrt(dx*dy)*(tempg(i,j,k)**4-tenv**4)
        force(i,j,k)=force(i,j,k)+firad/cp_gas(i,j,k)* &
          (pref/pr(i,j,k))**(1.-cv_gas(i,j,k)/cp_gas(i,j,k))
      enddo
    enddo
  enddo
        
end subroutine radiationSink
        
!-----------------------------------------------------------------------
! firerad_MC initializes and dictates the running of the Monte-
! Carlo radiation scheme
!-----------------------------------------------------------------------
subroutine firerad_MC(force,il,iu,jl,ju,lls)
  use gridlist_variables, only : nfuel,icallrad,prec
  use gridsetup, only : ittot
  use fuel_variables, only : lfuel
  use radiation_variables, only : fsiesrad,Etot,Emin,firad
  use thermo_variables, only : pr,cp_gas,cv_gas
  use constants, only : pref
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real(prec),intent(inout) :: force(il:iu,jl:ju,lls)
  integer :: i,j,k,ift
  character(len=10) :: iftname

  ! Executable Code
  if(mod(ittot,icallrad).eq.0)then
    call papv_defs ! define papvtotAZ,papvgas,papvtot
    call source_defs ! define sourcegas, sourcesol, EAZ, Etot, and lmc
    if(Etot.ge.Emin) call montecarlo
    do ift=1,nfuel
      write(iftname,"(i0)") ift
      call rmaxmin_3D(fsiesrad(ift,:,:,:),'fsiesrad_'//iftname, &
        il,iu,jl,ju,lfuel,1)
    enddo
    call rmaxmin_3D(firad,'firad',il,iu,jl,ju,lls,1)
  endif
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        force(i,j,k)=force(i,j,k)+firad(i,j,k)/cp_gas(i,j,k)* &
          (pref/pr(i,j,k))**(1.-cv_gas(i,j,k)/cp_gas(i,j,k))
      enddo
    enddo
  enddo

end subroutine firerad_MC

!-----------------------------------------------------------------------
! papv_defs defines papv, papvgas, and papvtot for the photon 
! generation and shooting
!-----------------------------------------------------------------------
subroutine papv_defs
  use gridlist_variables, only : l,ih,nfuel,rhoMicro,isootmodel, &
    tambient,crad,nprocx,nprocy,prec
  use gridsetup, only : np,mp
  use msga_variables, only : numprocs,mpi_comm_world,ierror
#ifdef DBL_PREC  
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use radiation_variables, only : papvift,papvgas,kfindex, &
    papvtotAZ,volumeAZ
  use fuel_variables, only : temps,lfuel,min_rhoFuel,sizeScale
  use thermo_variables, only : tempg
  use xvall, only : xv,xvrho,iO2,irho,xvfuel,irhof
  Implicit None
    
  ! Local Variables
  integer :: i,j,k,kcell,ift
  integer :: ia,ja,iproc,iprocx,jprocy
  real(prec) :: rmaxsootcon=0.2 ! Maximum soot concentration
  real(prec) :: absemis=0.3
  real(prec) :: aemit
  real(prec) :: fvsoot
  real(prec),allocatable :: papvlocal(:,:,:),papvtotAZtmp(:,:,:,:)

  ! Executable Code
  do ift=1,nfuel
    call update(temps(ift,:,:,:),temps(ift,:,:,:),np,mp,lfuel, &
      1-ih,np+ih,1-ih,mp+ih,1,0)
  enddo
  call update(tempg,tempg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
  papvift  = 0.0
  papvgas  = 0.0
  papvtotAZ= 0.0
  allocate(papvlocal(np,mp,l+nfuel)); papvlocal = 0.0
  if(isootmodel.eq.1)then
    ! HERE FP DEFINE A NEW MODEL FOR ABSORPTION COEFFICIENT:
    ! soot volume fraction is defined as correlated to mixture fraction, which is a linear fonction of tempg and O2 deficit
    do k=1,l+nfuel
      if(k.le.nfuel)then
        kcell=1
      else
        kcell=k-nfuel
      endif
      do j=1,mp
        do i=1,np
#ifdef DBL_PREC
          fvsoot=0.899d-6*(max(3.05d-5*(tempg(i,j,kcell)-tambient) &
            +0.05*(1.-xv(i,j,kcell,iO2)/(0.233*xv(i,j,kcell,irho))) &
            -1.48d-4,0.0))**0.671 
#else
          fvsoot=0.899e-6*(max(3.05e-5*(tempg(i,j,kcell)-tambient) &
            +0.05*(1.-xv(i,j,kcell,iO2)/(0.233*xv(i,j,kcell,irho))) &
            -1.48e-4,0.0))**0.671 
#endif
          papvgas(i,j,k)=266.0*7.0*fvsoot*tempg(i,j,kcell)  
          papvlocal(i,j,k)=papvgas(i,j,k)
        enddo  !do i
      enddo  !do j
    enddo  !do k
  else
    do k=1,l+nfuel
      if(k.le.nfuel)then
        kcell=1
      else
        kcell=k-nfuel
      endif
      do j=1,mp
        do i=1,np
          !set absemis and aemit
          aemit=rmaxsootcon*sqrt(max(.233-xvrho(i,j,kcell,iO2),0.) &
            /.233)*xv(i,j,kcell,irho)
          papvgas(i,j,k)=crad*aemit*absemis*0.25
          papvlocal(i,j,k)=papvgas(i,j,k)
        enddo  !do i
      enddo  !do j
    enddo  !do k
  endif

  do k=1,nfuel
    do i=1,np
      do j=1,mp
        do ift=1,nfuel
          if(kfindex(i,j,ift).ge.k.and. &
            xvfuel(ift,i,j,1,irhof).gt.min_rhoFuel)then
            papvift(ift,i,j,k)=2./sizeScale(ift,i,j,1)* &
              xvfuel(ift,i,j,1,irhof)/(4.*rhoMicro)* &
              volumeAZ(i,j,k)/sum(volumeAZ(i,j,1:kfindex(i,j,ift)))
            papvlocal(i,j,k)=papvlocal(i,j,k)+papvift(ift,i,j,k)
          endif
        enddo
      enddo
    enddo
  enddo
  if(lfuel.gt.1)then
    do k=2,lfuel
      do j=1,mp
        do i=1,np
          do ift=1,nfuel
            if(xvfuel(ift,i,j,k,irhof).gt.min_rhoFuel)then
              papvift(ift,i,j,k+nfuel)=2./sizeScale(ift,i,j,k)* &
                xvfuel(ift,i,j,k,irhof)/(4.*rhoMicro)
              papvlocal(i,j,k+nfuel)=papvlocal(i,j,k+nfuel)+ &
                papvift(ift,i,j,k+nfuel)
            endif
          enddo ! do ift
        enddo  !do i
      enddo  !do j
    enddo  !do k
  endif

  allocate(papvtotAZtmp(np,mp,l+nfuel,numprocs))
  papvtotAZtmp = 0.0
#ifdef DBL_PREC
  call mpi_allgather(papvlocal,np*mp*(l+nfuel), &
    mpi_double_precision,papvtotAZtmp,np*mp*(l+nfuel), &
    mpi_double_precision,mpi_comm_world,ierror)
#else
  call mpi_allgather(papvlocal,np*mp*(l+nfuel),mpi_real, &
    papvtotAZtmp,np*mp*(l+nfuel),mpi_real,mpi_comm_world,ierror)
#endif
  do iprocx=1,nprocx
    do jprocy=1,nprocy
      iproc=iprocx+(jprocy-1)*nprocx
      do j=1,mp
        do i=1,np
          ia=(iprocx-1)*np+i
          ja=(jprocy-1)*mp+j
          papvtotAZ(ia,ja,:)=max(1.e-10,papvtotAZtmp(i,j,:,iproc))
        enddo
      enddo
    enddo
  enddo
  deallocate(papvlocal,papvtotAZtmp)
 
end subroutine papv_defs
      
!-----------------------------------------------------------------------
! source_defs defines sourcegas, sourcesol, EAZ, Etot, and
! lmc for photon generation and shooting
!-----------------------------------------------------------------------
subroutine source_defs()
  use gridlist_variables, only : n,m,l,tambient,nfuel,nprocx, &
    nprocy,rhoMicro,prec
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank,numprocs,mpi_sum,mpi_comm_world, &
    ierror
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use fuel_variables, only : lfuel,temps,sizeScale,min_rhoFuel
  use radiation_variables, only : sourcesol,sourcegas,papvift, &
    Etot,EAZ,volumeAZ,Emin,papvtotAZ,nphotbatch,lmc,papvgas
  use thermo_variables, only : tempg
  Implicit None 

  ! Local Variables
  integer :: i,j,k,ift,ia,ja
  integer :: iprocx,jprocy,iproc
  real(prec) :: csurr=.8  ! Surrounding weighting coeffient for gaussian subgrid temperature distribution
  real(prec) :: camb=.2  ! Local weighting coeffient for gaussian subgrid temperature distribution
  real(prec) :: nbrterm
  real(prec) :: facm=2.725 ! (3+2.449489743)/2
  real(prec) :: facp=0.275 ! (3-2.449489743)/2
  real(prec) :: localterm,wdth
  real(prec) :: t4bar
  real(prec) :: papvthresh
  real(prec) :: Etotlocal
  real(prec),allocatable :: EAZlocal(:,:,:),EAZtmp(:,:,:,:)

  ! Executable Code
  sourcesol=0.0
  sourcegas=0.0
  do k=1,nfuel
    do j=1,mp
      do i=1,np
        !define sourcesol in the bottom layer 
        do ift=1,nfuel
          if(temps(ift,i,j,1).gt.tambient+1)then
            !nbrterm accounts for neighbor influence on t4bar
            nbrterm=csurr/4.* &
              (abs(temps(ift,i,j,1)-temps(ift,i-1,j,1)) &
              +abs(temps(ift,i,j,1)-temps(ift,i+1,j,1)) &
              +abs(temps(ift,i,j,1)-temps(ift,i,j-1,1)) &
              +abs(temps(ift,i,j,1)-temps(ift,i,j+1,1)))
            if(lfuel.gt.1) &
              nbrterm=nbrterm*5./4.+csurr/5.* &
                (abs(temps(ift,i,j,1)-temps(ift,i,j,2)))
 
            !localterm accounts for local influence on t4bar
            localterm=camb*(temps(ift,i,j,1)-tambient)
            wdth=localterm+nbrterm
            t4bar=sqrt((temps(ift,i,j,1)**2+facp*wdth**2) &
              *(temps(ift,i,j,1)**2+facm*wdth**2))
#ifdef DBL_PREC
            sourcesol(ift,i,j,k)=4.0*papvift(ift,i,j,k) &
              *(t4bar-tambient**2)*(t4bar+tambient**2)*5.67051d-8
#else
            sourcesol(ift,i,j,k)=4.0*papvift(ift,i,j,k) &
              *(t4bar-tambient**2)*(t4bar+tambient**2)*5.67051e-8
#endif
          endif  !if(temps.gt.tambient)
        enddo ! loop on ift
      enddo ! do i
    enddo ! do j
  enddo ! do k
  if(lfuel.gt.1)then
    do k=2,lfuel
      do j=1,mp
        do i=1,np
          !define sourcesol above bottom layer
          do ift=1,nfuel
            if(temps(ift,i,j,k).gt.tambient+1)then
              !nbrterm accounts for neighbor influence on t4bar
              nbrterm=csurr/5.* &
                (abs(temps(ift,i,j,k)-temps(ift,i-1,j,k)) &
                +abs(temps(ift,i,j,k)-temps(ift,i+1,j,k)) &
                +abs(temps(ift,i,j,k)-temps(ift,i,j-1,k)) &
                +abs(temps(ift,i,j,k)-temps(ift,i,j+1,k)) &
                +abs(temps(ift,i,j,k)-temps(ift,i,j,k-1)))
              if(k.lt.lfuel) &
                nbrterm=nbrterm*6./5.+csurr/6.* &
                  (abs(temps(ift,i,j,k)-temps(ift,i,j,k+1)))
 
              !localterm accounts for local influence on t4bar
              localterm=camb*(temps(ift,i,j,k)-tambient)
              wdth=localterm+nbrterm
              t4bar=sqrt((temps(ift,i,j,k)**2+facp*wdth**2) &
                *(temps(ift,i,j,k)**2+facm*wdth**2))
#ifdef DBL_PREC
              sourcesol(ift,i,j,k+nfuel)=4.0*papvift(ift,i,j,k+nfuel) &
                *(t4bar-tambient**2)*(t4bar+tambient**2)*5.67051d-8
#else
              sourcesol(ift,i,j,k+nfuel)=4.0*papvift(ift,i,j,k+nfuel) &
                *(t4bar-tambient**2)*(t4bar+tambient**2)*5.67051e-8
#endif
            endif  !if(temps.gt.tambient)
          enddo ! loop on ift
        enddo ! do i
      enddo ! do j
    enddo ! do k
  endif

  do k=1,l-1
    do j=1,mp
      do i=1,np
        !define sourcegas
        if(tempg(i,j,k).gt.tambient+1)then
          !nbrterm accounts for neighbor influence on t4bar
          nbrterm=csurr/5.* &
            (abs(tempg(i,j,k)-tempg(i-1,j,k)) &
            +abs(tempg(i,j,k)-tempg(i+1,j,k)) &
            +abs(tempg(i,j,k)-tempg(i,j-1,k)) &
            +abs(tempg(i,j,k)-tempg(i,j+1,k)) &
            +abs(tempg(i,j,k)-tempg(i,j,k+1)))
          if(k.gt.1) then
            nbrterm=nbrterm*5./6.+csurr/6.* &
              (abs(tempg(i,j,k)-tempg(i,j,k-1)))
          endif
          
          ! localterm accounts for local influence on t4bar
          localterm=camb*(tempg(i,j,k)-tambient)
          wdth=localterm+nbrterm
          t4bar=sqrt((tempg(i,j,k)**2+facp*wdth**2) &
            *(tempg(i,j,k)**2+facm*wdth**2))
#ifdef DBL_PREC
          sourcegas(i,j,k)=4.0*papvgas(i,j,k+nfuel)* &
            (t4bar-tambient**2)*(t4bar+tambient**2)*5.67051d-8
#else
          sourcegas(i,j,k)=4.0*papvgas(i,j,k+nfuel)* &
            (t4bar-tambient**2)*(t4bar+tambient**2)*5.67051e-8
#endif
        endif  !if(tempg.gt.tambient)
      enddo !do i
    enddo !do j
  enddo !do k
      
  ! BELOW THAT LINE DEFINITION OF EAZ,EAZ2, Etot
  ! calculate the energy from sources in each cell on a per process subdomain
  allocate(EAZlocal(np,mp,l+nfuel)); EAZlocal = 0.0
  do k=2,l  ! FPAZ all real cell except the first one
    do j=1,mp
      do i=1,np
        EAZlocal(i,j,k+nfuel)=sourcegas(i,j,k)*volumeAZ(i,j,k+nfuel)
      enddo !do i
    enddo !do j
  enddo !do k
  if(lfuel.gt.1)then
    do k=2,lfuel
      do j=1,mp
        do i=1,np
          EAZlocal(i,j,k+nfuel)=EAZlocal(i,j,k+nfuel)+ &
            sum(sourcesol(:,i,j,k+nfuel))*volumeAZ(i,j,k+nfuel)
        enddo !do i
      enddo !do j
    enddo !do k
  endif
  do k=1,nfuel+1
    do j=1,mp
      do i=1,np
        EAZlocal(i,j,k)=(sourcegas(i,j,1)+ &
          sum(sourcesol(:,i,j,k)))*volumeAZ(i,j,k)
      enddo !do i
    enddo !do j
  enddo !do k

  ! gather the all subdomains to create full E array on all processes
  allocate(EAZtmp(np,mp,l+nfuel,numprocs)); EAZtmp = 0.0
#ifdef DBL_PREC
  call mpi_allgather(EAZlocal,np*mp*(l+nfuel),mpi_double_precision, &
    EAZtmp,np*mp*(l+nfuel),mpi_double_precision,mpi_comm_world,ierror)
#else
  call mpi_allgather(EAZlocal,np*mp*(l+nfuel),mpi_real, &
    EAZtmp,np*mp*(l+nfuel),mpi_real,mpi_comm_world,ierror)
#endif
  do iprocx=1,nprocx
    do jprocy=1,nprocy
      iproc=iprocx+(jprocy-1)*nprocx
      do j=1,mp
        do i=1,np
          ia=(iprocx-1)*np+i
          ja=(jprocy-1)*mp+j
          EAZ(ia,ja,:)=EAZtmp(i,j,:,iproc)
        enddo
      enddo
    enddo
  enddo
  Etotlocal=sum(EAZlocal(:,:,:))
  deallocate(EAZlocal,EAZtmp)
  
#ifdef DBL_PREC
  call mpi_reduce(Etotlocal,Etot,1,mpi_double_precision,mpi_sum,0, &
    mpi_comm_world,ierror)
#else
  call mpi_reduce(Etotlocal,Etot,1,mpi_real,mpi_sum,0,mpi_comm_world, &
    ierror)
#endif
  ! Workaround for a weird type mismatch bug in gfortran
#ifdef DBL_PREC
  call mpi_Bcast(Etot,1,mpi_double_precision,0,mpi_comm_world,ierror)
  papvthresh=2.*min_rhoFuel/(maxval(sizeScale)*rhoMicro)
#else
  call mpi_Bcast(Etot,1,mpi_real,0,mpi_comm_world,ierror)
  papvthresh=2.*min_rhoFuel/(maxval(sizeScale)*rhoMicro)
#endif
  if(Etot.LT.Emin)then
    if(mpi_rank.eq.0) print*,'Etot=',Etot,' is neglected'
  else
    lmc=5+nfuel
    do k=1,l+nfuel
      do j=1,m
        do i=1,n
          if(papvtotAZ(i,j,k).gt.papvthresh.or. &
            (EAZ(i,j,k)/Etot*nphotbatch).ge.1) then
            lmc = MAX(k,lmc)
          endif
        enddo !do i
      enddo !do j
    enddo !do k
    if(mpi_rank.eq.0) print*,'lmc',lmc 
  endif
      
end subroutine source_defs

!-----------------------------------------------------------------------
! main MC subroutine which resolves the photon dispersion scheme
! with a montecarlo resolution approach
!-----------------------------------------------------------------------
subroutine montecarlo
  use gridlist_variables, only : n,m,iradeastflux,prec
  use msga_variables, only : mpi_rank,mpi_integer,mpi_sum, &
    mpi_comm_world,ierror
  use radiation_variables, only : Etot,lmc,nphotbatch,nphotonAZ, &
    nphotonAZgather,nphotonEastFlux
  Implicit None

  ! Local Variables
  integer :: i,j,k
  integer :: nphottot
  integer :: iternum
  real(prec) :: maxpercabsdiff
  real(prec) :: Qphot
  real(prec) :: percabsAZold,percabsAZnew
  real(prec) :: accuracythreshold=5.e-5
  integer,allocatable :: nphotontmpEastFlux(:,:,:)

  ! Executable Code
  call allocate_photons

  if(iradeastflux.eq.1) accuracythreshold=1.e-5
  iternum=0
  nphottot=0
  maxpercabsdiff=1.
  percabsAZnew=0.
  do while(maxpercabsdiff.gt.accuracythreshold) ! loop for dynamic allocation of photons
    call shoot_photons
    nphottot=nphottot+nphotbatch
    iternum=iternum+1
    maxpercabsdiff=0.
    do k=1,lmc
      do j=1,m
        do i=1,n
          percabsAZold=percabsAZnew
          percabsAZnew=nphotonAZ(i,j,k)/nphottot
          maxpercabsdiff=max(maxpercabsdiff,percabsAZnew-percabsAZold)
        enddo
      enddo
    enddo
  enddo ! do while (maxpercabsdiff.gt.accuracythreshold) dynamic allocation
  
  Qphot=Etot/nphottot
  if(mpi_rank.eq.0) print*,'Energy of a photon=',Qphot,'Total energy=',&
    Etot,'Iterations=',iternum
  nphotonAZgather=0.
  call mpi_allreduce(nphotonAZ(:,:,1:lmc),nphotonAZgather(:,:,1:lmc), &
    lmc*n*m,mpi_integer,mpi_sum,mpi_comm_world,ierror)
  if(iradeastflux.eq.1)then
    nphotonEastFlux=0.
    allocate(nphotontmpEastFlux(n,m,lmc)); nphotontmpEastFlux = 0
    call mpi_allreduce(nphotonEastFlux(:,:,1:lmc),nphotontmpEastFlux, &
      lmc*n*m,mpi_integer,mpi_sum,mpi_comm_world,ierror)
    nphotonEastFlux=nphotontmpEastFlux
    deallocate(nphotontmpEastFlux)
  endif

  call updateRadArrays(Qphot)

end subroutine montecarlo

!-----------------------------------------------------------------------
! allocate_photons sets arrays and numbers associated with 
! photons used on each call of the radiation scheme
!-----------------------------------------------------------------------
subroutine allocate_photons
  use gridlist_variables, only : n,m,prec
  use msga_variables, only : numprocs,mpi_rank
  use radiation_variables, only : Etot,nphotbatch,nphotemisAZsplit, &
    lmc,EAZ,nphotonAZ
  Implicit None

  ! Local Variables
  integer :: i,j,k
  integer :: n2
  integer :: mypstrt,mypstp
  integer :: icpt,icptold
  real(prec) :: Qphotbatch
  integer,allocatable :: nperproc(:)
  real(prec),allocatable :: nphotemisAZ(:,:,:)

  ! Executable Code
  allocate(nperproc(0:numprocs-1),nphotemisAZ(n,m,lmc))
  nperproc = 0
  nphotemisAZ = 0.0

  Qphotbatch=Etot/nphotbatch
  nperproc(:)=nphotbatch/numprocs
  n2=MOD(nphotbatch,numprocs)
  nperproc(0:n2-1)=nperproc(0:n2-1)+n2
  
  ! Fill an integer array with the number of photons to be emitted in each cell
  do k=1,lmc
    do j=1,m
      do i=1,n
        nphotemisAZ(i,j,k)=EAZ(i,j,k)/Qphotbatch
      enddo
    enddo
  enddo
  if(mpi_rank.eq.0) print*,'There are',sum(nphotemisAZ), &
    ' photons in the batch'

  ! Define a start and stop point for the range of photons each process emits
  nphotemisAZsplit=0. ! For everything beyond lmc
  mypstrt=SUM(nperproc(0:mpi_rank-1))+1
  mypstp =mypstrt+nperproc(mpi_rank)-1
  icpt=0
  do k=1,lmc
    do j=1,m
      do i=1,n
        icptold=icpt
        icpt=icpt+int(nphotemisAZ(i,j,k))
        if(icptold.le.mypstrt)then
          if(icpt.ge.mypstp)then ! All photons go here
            nphotemisAZsplit(i,j,k)=mypstp-mypstrt+1
          else ! Only some photons go here
            nphotemisAZsplit(i,j,k)=icpt-mypstrt+1
          endif
        elseif(icptold.ge.mypstp)then
          if(icpt.ge.mypstp)then ! remainder of my photons go here
            nphotemisAZsplit(i,j,k)=mypstp-icptold
          else ! icpt is less than mystp so only some of my remaining photons go here
            nphotemisAZsplit(i,j,k)=icpt-icptold
          endif
        endif
      enddo
    enddo
  enddo
  deallocate(nperproc,nphotemisAZ)

  nphotonAZ=0.

end subroutine allocate_photons

!-----------------------------------------------------------------------
! this subroutine will shoot nphotemisAZsplit(ia,ja,k) and see 
! where they arrive (number added in nphotonAZ)
!-----------------------------------------------------------------------
subroutine shoot_photons
  use gridlist_variables, only : n,m,dx,dy,ibcx,ibcy,prec
  use constants, only : pi
  use radiation_variables, only : lmc,nphotemisAZsplit, &
    zCartEdgeAZGlobal,nphotonAZ,papvtotAZ
  Implicit None

  ! Local Variables
  integer :: i,j,k,iphot
  integer :: i1,j1,k1,kt
  real(prec) :: p1,p2,p3
  real(prec) :: Rtheta,Rphi,costheta,sintheta,phi
  real(prec) :: dir(3)
  real(prec) :: one=1.0
  real(prec) :: z1top,z2top,z3top,z4top,z1bot,z2bot,z3bot,z4bot
  real(prec) :: zbottomcell,ztopcell,zbottomreal,ztopreal
  real(prec) :: xf,yf,zf
  real(prec) :: lopt ! optical thickness
  real(prec) :: t,tx,ty,tz,tz1,tz2,delta
  real(prec) :: acoeftop,bcoeftop,ccoeftop,acoefbot,bcoefbot,ccoefbot
  !real(prec) :: zsurface

  ! Executable Code
  do k=1,lmc
    do j=1,m
      do i=1,n ! Loop for photons emitted by cell (i,j,k)
        do iphot=1,nphotemisAZsplit(i,j,k)
          ! Shoot direction
          call random_number(Rtheta)
          call random_number(Rphi)
          costheta=1.-2.*Rtheta
          sintheta=sqrt(1.-costheta**2.)
          phi=2.*pi*Rphi
          dir(1)=sintheta*cos(phi)
          dir(2)=sintheta*sin(phi)
          dir(3)=costheta
          !Determine initial photon starting positions (xf,yf,zf)
          z1top=zCartEdgeAZGlobal(1,i,j,k)
          z2top=zCartEdgeAZGlobal(2,i,j,k)
          z3top=zCartEdgeAZGlobal(3,i,j,k)
          z4top=zCartEdgeAZGlobal(4,i,j,k)
          z1bot=zCartEdgeAZGlobal(1,i,j,k-1)
          z2bot=zCartEdgeAZGlobal(2,i,j,k-1)
          z3bot=zCartEdgeAZGlobal(3,i,j,k-1)
          z4bot=zCartEdgeAZGlobal(4,i,j,k-1)
          zbottomcell=min(z1bot,z2bot,z3bot,z4bot)
          ztopcell=max(z1top,z2top,z3top,z4top)
          do
            call random_number(p1)
            call random_number(p2)
            call random_number(p3)
            xf=(i-1+p1)*dx
            yf=(j-1+p2)*dy
            zf=zbottomcell+p3*(ztopcell-zbottomcell)
            call zsurface(zbottomreal,p1,p2,z1bot,z2bot,z3bot,z4bot)
            call zsurface(ztopreal,p1,p2,z1top,z2top,z3top,z4top)
            if(zf.ge.zbottomreal.and.zf.le.ztopreal) exit
          enddo
          ! calculate which cell the photon arrives at
          i1=i
          j1=j
          k1=k
          kt=0
          call random_number(lopt)
          lopt=-log(lopt)
          do while(i1.ge.1.and.i1.le.n.and.j1.ge.1.and.j1.le.m.and. &
            k1.ge.1.and.k1.le.lmc.and.lopt.gt.0)
            if(dir(1).gt.0)then
              tx=(i1*dx-xf)/dir(1)
            else
              tx=(i1*dx-dx-xf)/dir(1)
            endif
            if(dir(2).gt.0)then
              ty=(j1*dy-yf)/dir(2)
            else
              ty=(j1*dy-dy-yf)/dir(2)
            endif
            tz=1e6
            tz1=0.
            tz2=0.
            z1top=zCartEdgeAZGlobal(1,i1,j1,k1)
            z2top=zCartEdgeAZGlobal(2,i1,j1,k1)
            z3top=zCartEdgeAZGlobal(3,i1,j1,k1)
            z4top=zCartEdgeAZGlobal(4,i1,j1,k1)
            z1bot=zCartEdgeAZGlobal(1,i1,j1,k1-1)
            z2bot=zCartEdgeAZGlobal(2,i1,j1,k1-1)
            z3bot=zCartEdgeAZGlobal(3,i1,j1,k1-1)
            z4bot=zCartEdgeAZGlobal(4,i1,j1,k1-1)
    
            ! Upward photon
            acoeftop=(z1top+z4top-z2top-z3top)*dir(1)*dir(2)/(dx*dy)
            bcoeftop=(-dir(3)+dir(2)/dy*(z3top-z1top)+dir(1)/dx* &
              (z2top-z1top)+(z1top+z4top-z2top-z3top)* &
              (dir(2)/dy*(xf/dx-(i1-1))+dir(1)/dx*(yf/dy-(j1-1))))
            ccoeftop=(z1top-zf+(z3top-z1top)*(yf/dy-(j1-1))+ &
              (z2top-z1top)*(xf/dx-(i1-1))+(z1top+z4top-z2top-z3top)* &
              (xf/dx-(i1-1))*(yf/dy-(j1-1)))
            if(kt.ge.0)then 
              if(abs(acoeftop).lt.1d-6)then    !a.eq.0
                if(ccoeftop*bcoeftop.lt.0) tz=-ccoeftop/bcoeftop  ! b.eq.0 correct intersection with tz>0
              else   ! the code never goes here in case of flat topo
                delta=bcoeftop**2-4*acoeftop*ccoeftop
                if(delta.gt.0.)then     !two solutions
                  tz1=-0.5*(bcoeftop-sqrt(delta))/acoeftop
                  tz2=-0.5*(bcoeftop+sqrt(delta))/acoeftop
                  if(tz1.gt.0.)then   ! choose the smallest positive root
                    tz=tz1 
                    if(tz2.gt.0.) tz=min(tz1,tz2)
                  elseif(tz2.gt.0.)then
                    tz=tz2
                  endif
                endif
              endif
            elseif(abs(acoeftop).gt.1d-6)then         ! look for potential second intersection
              ! here we use the fact that ccoef should be equal to 0 in this case (k1=k-1)
              if(bcoeftop*acoeftop.lt.0) tz=-bcoeftop/acoeftop  ! choose the second root when positive (first should be 0)
            endif

            ! Downward photon
            acoefbot=(z1bot+z4bot-z2bot-z3bot)*dir(1)*dir(2)/(dx*dy)
            bcoefbot=(-dir(3)+dir(2)/dy*(z3bot-z1bot)+dir(1)/dx* &
              (z2bot-z1bot)+(z1bot+z4bot-z2bot-z3bot)*(dir(2)/dy* &
              (xf/dx-(i1-1))+dir(1)/dx*(yf/dy-(j1-1))))
            ccoefbot=(z1bot-zf+(z3bot-z1bot)*(yf/dy-(j1-1))+ &
              (z2bot-z1bot)*(xf/dx-(i1-1))+(z1bot+z4bot-z2bot-z3bot)* &
              (xf/dx-(i1-1))*(yf/dy-(j1-1)))
         
            if(kt.le.0)then
              if(abs(acoeftop).lt.1d-6)then    !a.eq.0
                if(ccoefbot*bcoefbot.lt.0.) tz=-ccoefbot/bcoefbot ! correct intersection with tz>0
              else
                delta=bcoefbot**2-4*acoefbot*ccoefbot
                if(delta.gt.0.)then     !two solutions
                  tz1=-0.5*(bcoefbot-sqrt(delta))/acoefbot
                  tz2=-0.5*(bcoefbot+sqrt(delta))/acoefbot
                  if(tz1.gt.0.)then   ! choose the smallest positive root
                    tz=tz1
                    if(tz2.gt.0.) tz=min(tz1,tz2)
                  elseif(tz2.gt.0.)then
                    tz=tz2
                  endif
                endif
              endif
            elseif(abs(acoeftop).gt.1d-6)then         ! look for potential second intersection
              ! here we use the fact that ccoef should be equal to 0 in this case (suivvox(3)=-1)
              if(bcoefbot*acoefbot.lt.0) tz=-bcoefbot/acoefbot  ! choose the second root when positive (first should be 0)
            endif

            ! Intersection with cell face
            t=min(tx,ty,tz)
            lopt=lopt-t*papvtotAZ(i1,j1,k1)
            if(lopt.ge.0)then
              kt=0
              if(abs(t-tx).lt.1.e-8)then
                i1=i1+int(sign(one,dir(1)))
              elseif(abs(t-ty).lt.1.e-8)then
                j1=j1+int(sign(one,dir(2)))
              else
                k1=k1+int(sign(one,dir(3)))
                kt=int(sign(one,dir(3)))
              endif
           
              xf=xf+dir(1)*t
              yf=yf+dir(2)*t
              zf=zf+dir(3)*t
              
              ! Cyclic boundary conditions implementation
              if(ibcx.eq.0)then
                if(i1.lt.1)then
                  i1=i1+n
                  xf=xf+n*dx
                elseif(i1.gt.n)then
                  i1=i1-n
                  xf=xf-n*dx
                endif
              endif
              if(ibcy.eq.0)then
                if(j1.lt.1)then
                  j1=j1+m
                  yf=yf+m*dy
                elseif(j1.gt.m)then
                  j1=j1-m
                  yf=yf-m*dy
                endif
              endif
            endif ! lopt.gt.0
          enddo
          if(i1.ge.1.and.i1.le.n.and.j1.ge.1.and.j1.le.m.and. &
            k1.ge.1.and.k1.le.lmc) &
            nphotonAZ(i1,j1,k1)=nphotonAZ(i1,j1,k1)+1
        enddo ! do iphot
      enddo ! do i
    enddo ! do j
  enddo ! do k
      
end subroutine shoot_photons

!-----------------------------------------------------------------------
! this subroutine update values of firad and fsiesrad and 
! nphotonEastFlux according to nphotonAZtot and Qphot
!-----------------------------------------------------------------------
subroutine updateRadArrays(Qphot)
  use gridlist_variables, only : dy,nfuel,iradeastflux,prec
  use gridsetup, only : np,mp
  use msga_variables, only : npos,mpos
  use metric_variables, only : zedge
  use fuel_variables, only : lfuel
  use radiation_variables, only : lmc,papvtotAZ,firad, &
    sourcegas,volumeAZ,nphotonAZgather,papvgas,fsiesrad, &
    papvift,eastFlux,nphotonEastFlux,sourcesol
  use zcart_function, only : zcart
  Implicit None

  ! Local Variables
  real(prec),intent(in) :: Qphot

  integer :: i,j,k,kcell,ift,ia,ja
#ifdef DBL_PREC
  real(prec) :: papvthresh=8.d-6
#else
  real(prec) :: papvthresh=8.e-6
#endif
               
  ! Executable Code
  firad=-sourcegas
  fsiesrad=0.
  do k=1,lmc  !FPAZ
    if(k.le.nfuel)then
      kcell=1
    else
      kcell=k-nfuel
    endif
    do j=1,mp
      ja=(mpos-1)*mp+j
      do i=1,np
        ia=(npos-1)*np+i 
        if(papvtotAZ(ia,ja,k).gt.papvthresh)then
          firad(i,j,kcell)=firad(i,j,kcell)+Qphot* &
            nphotonAZgather(ia,ja,k)* &
            papvgas(i,j,k)/papvtotAZ(ia,ja,k)/volumeAZ(i,j,k)
          if(k.le.lfuel+nfuel)then
            do ift=1,nfuel
              fsiesrad(ift,i,j,kcell)=fsiesrad(ift,i,j,kcell)- &
                sourcesol(ift,i,j,k) &
                +Qphot*nphotonAZgather(ia,ja,k)*papvift(ift,i,j,k) &
                /papvtotAZ(ia,ja,k)/volumeAZ(i,j,k)
            enddo !do ift
          endif
        endif !papvtot2.gt.papvthresh
        if(iradeastflux.eq.1) &
          eastFlux(i,j,k)=Qphot*nphotonEastFlux(ia,ja,k+nfuel)/ &
          (dy*(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)))
      enddo !do i
    enddo !do j
  enddo !do k
       
end subroutine updateRadArrays

!-----------------------------------------------------------------------
! Utility Functions
!-----------------------------------------------------------------------

subroutine zsurface(zsurf,p1,p2,z1,z2,z3,z4)
use gridlist_variables, only : prec
  Implicit None

  real(prec),intent(inout) :: zsurf

  ! Local Variables
  real(prec),intent(in) :: p1,p2
  real(prec),intent(in) :: z1,z2,z3,z4
  
  ! Executable Code  
  zsurf=((z1+z4-z2-z3)*p2+z2-z1)*p1+(z3-z1)*p2+z1

end subroutine zsurface
