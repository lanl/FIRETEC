!-----------------------------------------------------------------------
! firerad_MC initializes and dictates the running of the Monte-
! Carlo radiation scheme
!-----------------------------------------------------------------------
subroutine firerad_MC(force,force_fuel,xvsolid,ilh,iuh,jlh,juh, &
    il,iu,jl,ju,lls,nfuel)
  use gridlist_variables, only : wp
  use radiation_gridlist_variables, only : icallrad
  use gridsetup, only : ittot
  use fuel_variables, only : lfuel
  use radiation_general_variables, only : fsiesrad,firad
  use radiation_montecarlo_variables, only : Etot,Emin
  use thermo_variables, only : pr,cp_gas,cv_gas
  use constants, only : pref
  use rmaxmin_mod, only : rmaxmin
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,il,iu,jl,ju,lls,nfuel
  real(wp),intent(in) :: xvsolid(nfuel,il:iu,jl:ju,lfuel)
  real(wp),intent(inout) :: force(ilh:iuh,jlh:juh,lls)
  real(wp),intent(inout) :: force_fuel(nfuel,il:iu,jl:ju,lfuel)

  integer :: i,j,k,ift
  character(len=10) :: iftname

  ! Executable Code
  if(mod(ittot,icallrad).eq.0)then
    call papv_defs ! define papvtotAZ,papvgas,papvtot
    call source_defs ! define sourcegas, sourcesol, EAZ, Etot, and lmc
    if(Etot.ge.Emin) then
      call montecarlo
    else
      return
    endif
    do ift=1,nfuel
      write(iftname,"(i0)") ift
      call rmaxmin(fsiesrad(ift,:,:,:),'fsiesrad_'//iftname)
    enddo
    call rmaxmin(firad,'firad')
  endif
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        force(i,j,k)=force(i,j,k)+firad(i,j,k)/cp_gas(i,j,k)* &
          (pref/pr(i,j,k))**(1._wp-cv_gas(i,j,k)/cp_gas(i,j,k))
      enddo
    enddo
  enddo
  do k=1,lfuel
    do j=jl,ju
      do i=il,iu
        do ift=1,nfuel
          force_fuel(ift,i,j,k)=force_fuel(ift,i,j,k)+ &
            fsiesrad(ift,i,j,k)/xvsolid(ift,i,j,k)
        enddo
      enddo
    enddo
  enddo
end subroutine firerad_MC

!-----------------------------------------------------------------------
! papv_defs defines papv, papvgas, and papvtot for the photon 
! generation and shooting
!-----------------------------------------------------------------------
subroutine papv_defs
#ifdef __GFORTRAN__
  use mpi, only : mpi_allgather
#endif
  use gridlist_variables, only : l,ih,tambient,nprocx,nprocy,wp
  use radiation_gridlist_variables, only : isootmodel,crad
  use vegetation_gridlist_variables, only : nfuel,rhoMicro
  use gridsetup, only : np,mp
  use msga_variables, only : numprocs,mpi_comm_world,ierror
  use msga_variables, only : mpi_wp
  use radiation_montecarlo_variables, only : papvift,papvgas, &
    papvtotAZ,volumeAZ
  use fuel_variables, only : temps,lfuel,min_rhoFuel,sizeScale, &
    kfrank
  use thermo_variables, only : tempg
  use xvall, only : xv,xvrho,iO2,irho,xvfuel,irhof
  Implicit None
    
  ! Local Variables
  integer :: i,j,k,kcell,ift
  integer :: ia,ja,iproc,iprocx,jprocy
  real(wp) :: rmaxsootcon=0.2_wp ! Maximum soot concentration
  real(wp) :: absemis=0.3_wp
  real(wp) :: aemit
  real(wp) :: fvsoot
  real(wp) :: o2frac
  real(wp),allocatable :: papvlocal(:,:,:),papvtotAZtmp(:,:,:,:)

  ! Executable Code
  do ift=1,nfuel
    call update(temps(ift,:,:,:),np,mp,lfuel, &
      1-ih,np+ih,1-ih,mp+ih,1,0,0)
  enddo
  call update(tempg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0,0)
  papvift  = 0._wp
  papvgas  = 0._wp
  papvtotAZ= 0._wp
  allocate(papvlocal(np,mp,l+nfuel)); papvlocal = 0._wp
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
          o2frac=0.233_wp
          if(iO2.gt.0) o2frac=xv(i,j,kcell,iO2)/xv(i,j,kcell,irho)
          fvsoot=0.899e-6_wp*(max(3.05e-5_wp* &
            (tempg(i,j,kcell)-tambient)+ &
            0.05_wp*(1._wp-o2frac/0.233_wp) &
            -1.48e-4_wp,0._wp))**0.671_wp
          papvgas(i,j,k)=266.0_wp*7._wp*fvsoot*tempg(i,j,kcell)  
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
          o2frac=0.233_wp
          if(iO2.gt.0) o2frac=xvrho(i,j,kcell,iO2)
          aemit=rmaxsootcon*xv(i,j,kcell,irho)* &
            sqrt(max(0.233_wp-o2frac,0._wp)/0.233_wp)
          papvgas(i,j,k)=crad*aemit*absemis*0.25_wp
          papvlocal(i,j,k)=papvgas(i,j,k)
        enddo  !do i
      enddo  !do j
    enddo  !do k
  endif

  do k=1,nfuel
    do j=1,mp
      do i=1,np
        do ift=1,nfuel
          if(kfrank(i,j,ift).ge.k.and. &
            xvfuel(ift,i,j,1,irhof).gt.min_rhoFuel)then
            papvift(ift,i,j,k)=2._wp/sizeScale(ift,i,j,1)* &
              xvfuel(ift,i,j,1,irhof)/(4._wp*rhoMicro)* &
              volumeAZ(i,j,k)/sum(volumeAZ(i,j,1:kfrank(i,j,ift)))
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
              papvift(ift,i,j,k+nfuel)=2._wp/sizeScale(ift,i,j,k)* &
                xvfuel(ift,i,j,k,irhof)/(4._wp*rhoMicro)
              papvlocal(i,j,k+nfuel)=papvlocal(i,j,k+nfuel)+ &
                papvift(ift,i,j,k+nfuel)
            endif
          enddo ! do ift
        enddo  !do i
      enddo  !do j
    enddo  !do k
  endif

  allocate(papvtotAZtmp(np,mp,l+nfuel,numprocs))
  papvtotAZtmp = 0.0_wp
  call mpi_allgather(papvlocal,np*mp*(l+nfuel),mpi_wp, &
    papvtotAZtmp,np*mp*(l+nfuel),mpi_wp,mpi_comm_world,ierror)
  do iprocx=1,nprocx
    do jprocy=1,nprocy
      iproc=iprocx+(jprocy-1)*nprocx
      do j=1,mp
        do i=1,np
          ia=(iprocx-1)*np+i
          ja=(jprocy-1)*mp+j
          papvtotAZ(ia,ja,:)=max(1.e-10_wp,papvtotAZtmp(i,j,:,iproc))
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
#ifdef __GFORTRAN__
  use mpi, only : mpi_allgather,mpi_allreduce
#endif
  use gridlist_variables, only : n,m,l,tambient,nprocx,nprocy,wp
  use vegetation_gridlist_variables, only : nfuel
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank,numprocs,mpi_comm_world,ierror
  use msga_variables, only : mpi_wp
  use fuel_variables, only : lfuel,temps
  use radiation_montecarlo_variables, only : sourcesol,sourcegas, &
    papvift,Etot,EAZ,volumeAZ,Emin,papvtotAZ,nphotbatch0,lmc,papvgas
  use radiation_general_variables, only : papvthresh
  use thermo_variables, only : tempg
  use constants, only : sigma
  Implicit None 

  ! Local Variables
  integer :: i,j,k,ift,ia,ja
  integer :: iprocx,jprocy,iproc
  real(wp) :: csurr=0.8_wp  ! Surrounding weighting coeffient for gaussian subgrid temperature distribution
  real(wp) :: camb=0.2_wp  ! Local weighting coeffient for gaussian subgrid temperature distribution
  real(wp) :: nbrterm
  real(wp) :: facm=2.725_wp ! (3+2.449489743)/2
  real(wp) :: facp=0.275_wp ! (3-2.449489743)/2
  real(wp) :: localterm,wdth
  real(wp) :: t4bar
  real(wp),allocatable :: EAZlocal(:,:,:),EAZtmp(:,:,:,:)

  ! Executable Code
  sourcesol=0.0_wp
  sourcegas=0.0_wp
  do k=1,nfuel
    do j=1,mp
      do i=1,np
        !define sourcesol in the bottom layer 
        do ift=1,nfuel
          if(temps(ift,i,j,1).gt.tambient+1)then
            !nbrterm accounts for neighbor influence on t4bar
            nbrterm=csurr/4._wp* &
              (abs(temps(ift,i,j,1)-temps(ift,i-1,j,1)) &
              +abs(temps(ift,i,j,1)-temps(ift,i+1,j,1)) &
              +abs(temps(ift,i,j,1)-temps(ift,i,j-1,1)) &
              +abs(temps(ift,i,j,1)-temps(ift,i,j+1,1)))
            if(lfuel.gt.1) &
              nbrterm=nbrterm*5._wp/4._wp+csurr/5._wp* &
                (abs(temps(ift,i,j,1)-temps(ift,i,j,2)))
 
            !localterm accounts for local influence on t4bar
            localterm=camb*(temps(ift,i,j,1)-tambient)
            wdth=localterm+nbrterm
            t4bar=sqrt((temps(ift,i,j,1)**2+facp*wdth**2) &
              *(temps(ift,i,j,1)**2+facm*wdth**2))
            sourcesol(ift,i,j,k)=4._wp*papvift(ift,i,j,k) &
              *(t4bar-tambient**2)*(t4bar+tambient**2)*sigma
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
              nbrterm=csurr/5._wp* &
                (abs(temps(ift,i,j,k)-temps(ift,i-1,j,k)) &
                +abs(temps(ift,i,j,k)-temps(ift,i+1,j,k)) &
                +abs(temps(ift,i,j,k)-temps(ift,i,j-1,k)) &
                +abs(temps(ift,i,j,k)-temps(ift,i,j+1,k)) &
                +abs(temps(ift,i,j,k)-temps(ift,i,j,k-1)))
              if(k.lt.lfuel) &
                nbrterm=nbrterm*6._wp/5._wp+csurr/6._wp* &
                  (abs(temps(ift,i,j,k)-temps(ift,i,j,k+1)))
 
              !localterm accounts for local influence on t4bar
              localterm=camb*(temps(ift,i,j,k)-tambient)
              wdth=localterm+nbrterm
              t4bar=sqrt((temps(ift,i,j,k)**2+facp*wdth**2) &
                *(temps(ift,i,j,k)**2+facm*wdth**2))
              sourcesol(ift,i,j,k+nfuel)=4._wp*papvift(ift,i,j,k+nfuel) &
                *(t4bar-tambient**2)*(t4bar+tambient**2)*sigma
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
          nbrterm=csurr/5.0_wp* &
            (abs(tempg(i,j,k)-tempg(i-1,j,k)) &
            +abs(tempg(i,j,k)-tempg(i+1,j,k)) &
            +abs(tempg(i,j,k)-tempg(i,j-1,k)) &
            +abs(tempg(i,j,k)-tempg(i,j+1,k)) &
            +abs(tempg(i,j,k)-tempg(i,j,k+1)))
          if(k.gt.1) then
            nbrterm=nbrterm*5._wp/6._wp+csurr/6._wp* &
              (abs(tempg(i,j,k)-tempg(i,j,k-1)))
          endif
          
          ! localterm accounts for local influence on t4bar
          localterm=camb*(tempg(i,j,k)-tambient)
          wdth=localterm+nbrterm
          t4bar=sqrt((tempg(i,j,k)**2+facp*wdth**2) &
            *(tempg(i,j,k)**2+facm*wdth**2))
          sourcegas(i,j,k)=4._wp*papvgas(i,j,k+nfuel)* &
            (t4bar-tambient**2)*(t4bar+tambient**2)*sigma
        endif  !if(tempg.gt.tambient)
      enddo !do i
    enddo !do j
  enddo !do k
      
  ! BELOW THAT LINE DEFINITION OF EAZ,EAZ2, Etot
  ! calculate the energy from sources in each cell on a per process subdomain
  allocate(EAZlocal(np,mp,l+nfuel))
  do k=2,l  ! gas energy source all real cell except the first one
    do j=1,mp
      do i=1,np
        EAZlocal(i,j,k+nfuel)=sourcegas(i,j,k)*volumeAZ(i,j,k+nfuel)
      enddo !do i
    enddo !do j
  enddo !do k
  if(lfuel.gt.1)then ! add in solid energy source all real cell except the first one
    do k=2,lfuel
      do j=1,mp
        do i=1,np
          EAZlocal(i,j,k+nfuel)=EAZlocal(i,j,k+nfuel)+ &
            sum(sourcesol(:,i,j,k+nfuel))*volumeAZ(i,j,k+nfuel)
        enddo !do i
      enddo !do j
    enddo !do k
  endif
  do k=1,nfuel+1 ! resolve the broken up first cell among it's different zones
    do j=1,mp
      do i=1,np
        EAZlocal(i,j,k)=(sourcegas(i,j,1)+ &
          sum(sourcesol(:,i,j,k)))*volumeAZ(i,j,k)
      enddo !do i
    enddo !do j
  enddo !do k

  ! gather the all subdomains to create full E array on all processes
  allocate(EAZtmp(np,mp,l+nfuel,numprocs)); EAZtmp = 0.0_wp
  call mpi_allgather(EAZlocal,np*mp*(l+nfuel),mpi_wp, &
    EAZtmp,np*mp*(l+nfuel),mpi_wp,mpi_comm_world,ierror)
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
  deallocate(EAZlocal,EAZtmp)
  
  Etot = sum(EAZ)
  if(Etot.LT.Emin)then
    if(mpi_rank.eq.0) print*,'Etot=',Etot,' is neglected'
  else
    lmc=5+nfuel
    do k=1,l+nfuel
      do j=1,m
        do i=1,n
          if(papvtotAZ(i,j,k).gt.papvthresh.or. &
            (EAZ(i,j,k)/Etot*nphotbatch0).ge.1) then
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
#ifdef __GFORTRAN__
  use mpi, only : mpi_allreduce
#endif
  use gridlist_variables, only : n,m,wp
  use msga_variables, only : mpi_rank,mpi_integer,mpi_sum,mpi_max, &
    mpi_wp,mpi_comm_world,ierror
  use radiation_montecarlo_variables, only : Etot,lmc,nphotbatch, &
    nphotonAZ,nphotonAZgather
  Implicit None

  ! Local Variables
  integer :: nphottot
  integer :: iternum
  real(wp) :: maxpercabsdiff,maxpercabsdiffLocal
  real(wp) :: Qphot
  real(wp) :: accuracythreshold=5.e-5_wp
  real(wp),allocatable :: percabsAZold(:,:,:),percabsAZnew(:,:,:)

  ! Executable Code
  call allocate_photons
  allocate(percabsAZold(n,m,lmc),percabsAZnew(n,m,lmc))

  iternum=0
  nphottot=0
  maxpercabsdiff=1._wp
  percabsAZnew=0._wp
  do while(maxpercabsdiff.gt.accuracythreshold) ! loop for dynamic allocation of photons
    call shoot_photons
    
    nphottot=nphottot+nphotbatch
    iternum=iternum+1
    percabsAZold=percabsAZnew
    percabsAZnew=real(nphotonAZ(:,:,1:lmc),kind=wp)/ &
      real(nphottot,kind=wp)
    maxpercabsdiffLocal=maxval(abs(percabsAZnew-percabsAZold))
    call mpi_allreduce(maxpercabsdiffLocal,maxpercabsdiff,1,mpi_wp, &
      mpi_max,mpi_comm_world,ierror)
  enddo ! do while (maxpercabsdiff.gt.accuracythreshold) dynamic allocation
  
  Qphot=Etot/nphottot
  if(mpi_rank.eq.0) print*,'Energy of a photon=',Qphot, &
    'Total energy=',Etot,'Iterations=',iternum
  nphotonAZgather=0
  call mpi_allreduce(nphotonAZ(:,:,1:lmc),nphotonAZgather(:,:,1:lmc), &
    lmc*n*m,mpi_integer,mpi_sum,mpi_comm_world,ierror)

  call updateRadArrays(Qphot)

end subroutine montecarlo

!-----------------------------------------------------------------------
! allocate_photons sets arrays and numbers associated with 
! photons used on each call of the radiation scheme
!-----------------------------------------------------------------------
subroutine allocate_photons
  use gridlist_variables, only : n,m,wp
  use msga_variables, only : numprocs,mpi_rank
  use radiation_montecarlo_variables, only : Etot,nphotbatch0, &
    nphotbatch,nphotemisAZsplit,lmc,EAZ,nphotonAZ
  Implicit None

  ! Local Variables
  integer :: i,j,k
  integer :: n2
  integer :: mypstrt,mypstp
  integer :: icpt,icptold
  integer :: overlapStart,overlapStop
  real(wp) :: Qphotbatch
  integer,allocatable :: nperproc(:)
  integer,allocatable :: nphotemisAZ(:,:,:)

  ! Executable Code
  allocate(nperproc(0:numprocs-1),nphotemisAZ(n,m,lmc))

  ! Fill an integer array with the number of photons to be emitted in each cell
  Qphotbatch=Etot/nphotbatch0
  do k=1,lmc
    do j=1,m
      do i=1,n
        nphotemisAZ(i,j,k)=int(EAZ(i,j,k)/Qphotbatch)
      enddo
    enddo
  enddo
  nphotbatch=sum(nphotemisAZ)
  nperproc(:)=int(nphotbatch/numprocs)
  n2=MOD(nphotbatch,numprocs)
  nperproc(0:n2-1)=nperproc(0:n2-1)+1

  ! Define a start and stop point for the range of photons each process emits
  nphotemisAZsplit=0 ! For everything beyond lmc
  mypstrt=SUM(nperproc(0:mpi_rank-1))
  mypstp =mypstrt+nperproc(mpi_rank)
  icpt=0
  k_loop: do k=1,lmc
    do j=1,m
      do i=1,n
        icptold=icpt
        icpt   =icpt+nphotemisAZ(i,j,k)
        overlapStart=max(icptold,mypstrt)
        overlapStop =min(icpt,mypstp)
        if(overlapStop.gt.overlapStart) then
          nphotemisAZsplit(i,j,k)=overlapStop-overlapStart
          if(overlapStop.eq.mypstp)then
            exit k_loop
          endif
        endif
      enddo
    enddo
  enddo k_loop
  deallocate(nperproc,nphotemisAZ)

  nphotonAZ=0

end subroutine allocate_photons

!-----------------------------------------------------------------------
! this subroutine will shoot nphotemisAZsplit(ia,ja,k) and see 
! where they arrive (number added in nphotonAZ)
!-----------------------------------------------------------------------
subroutine shoot_photons
  use gridlist_variables, only : n,m,dx,dy,ibc_left,ibc_right, &
    ibc_top,ibc_bot,wp
  use constants, only : pi,tolerance
  use radiation_montecarlo_variables, only : lmc,nphotemisAZsplit, &
    zCartEdgeAZGlobal,nphotonAZ,papvtotAZ
  use gridlist_variables, only : wp
  Implicit None

  ! Local Variables
  integer :: i,j,k,iphot
  integer :: i1,j1,k1,kt
  real(wp) :: p1,p2,p3
  real(wp) :: Rtheta,Rphi,costheta,sintheta,phi
  real(wp) :: dir(3)
  real(wp) :: one=1.0_wp
  real(wp) :: z1top,z2top,z3top,z4top,z1bot,z2bot,z3bot,z4bot
  real(wp) :: zbottomcell,ztopcell,zbottomreal,ztopreal
  real(wp) :: xf,yf,zf
  real(wp) :: lopt ! optical thickness
  real(wp) :: t,tx,ty,tz,tz1,tz2,delta
  real(wp) :: acoeftop,bcoeftop,ccoeftop,acoefbot,bcoefbot,ccoefbot

  ! Executable Code
  do k=1,lmc
    do j=1,m
      do i=1,n ! Loop for photons emitted by cell (i,j,k)
        do iphot=1,nphotemisAZsplit(i,j,k)
          ! Shoot direction
          call random_number(Rtheta)
          call random_number(Rphi)
          costheta=1._wp-2._wp*Rtheta
          sintheta=sqrt(1.0_wp-costheta**2.0_wp)
          phi=2._wp*pi*Rphi
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
          lopt=-log((1.-lopt))
          do while(i1.ge.1.and.i1.le.n.and.j1.ge.1.and.j1.le.m.and. &
            k1.ge.1.and.k1.le.lmc.and.lopt.gt.0)
            if(dir(1).gt.0)then
              tx=(i1*dx-xf)/dir(1)
            else
              tx=(i1*dx-dx-xf)/min(-tolerance,dir(1))
            endif
            if(dir(2).gt.0)then
              ty=(j1*dy-yf)/dir(2)
            else
              ty=(j1*dy-dy-yf)/min(-tolerance,dir(2))
            endif
            tz=1e6
            tz1=0.0_wp
            tz2=0.0_wp
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
              if(abs(acoeftop).lt.tolerance)then  !a.eq.0
                if(ccoeftop*bcoeftop.lt.0) tz=-ccoeftop/bcoeftop  ! b.eq.0 correct intersection with tz>0
              else   ! the code never goes here in case of flat topo
                delta=bcoeftop**2-4*acoeftop*ccoeftop
                if(delta.gt.0.0_wp)then     !two solutions
                  tz1=-0.5_wp*(bcoeftop-sqrt(delta))/acoeftop
                  tz2=-0.5_wp*(bcoeftop+sqrt(delta))/acoeftop
                  if(tz1.gt.0.0_wp)then   ! choose the smallest positive root
                    tz=tz1 
                    if(tz2.gt.0.0_wp) tz=min(tz1,tz2)
                  elseif(tz2.gt.0.0_wp)then
                    tz=tz2
                  endif
                endif
              endif
            elseif(abs(acoeftop).gt.tolerance)then  ! look for potential second intersection
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
              if(abs(acoefbot).lt.tolerance)then    !a.eq.0
                if(ccoefbot*bcoefbot.lt.0.0_wp) tz=-ccoefbot/bcoefbot ! correct intersection with tz>0
              else
                delta=bcoefbot**2-4*acoefbot*ccoefbot
                if(delta.gt.0.0_wp)then     !two solutions
                  tz1=-0.5_wp*(bcoefbot-sqrt(delta))/acoefbot
                  tz2=-0.5_wp*(bcoefbot+sqrt(delta))/acoefbot
                  if(tz1.gt.0.0_wp)then   ! choose the smallest positive root
                    tz=tz1
                    if(tz2.gt.0.0_wp) tz=min(tz1,tz2)
                  elseif(tz2.gt.0.0_wp)then
                    tz=tz2
                  endif
                endif
              endif
            elseif(abs(acoefbot).gt.tolerance)then   ! look for potential second intersection
              ! here we use the fact that ccoef should be equal to 0 in this case (suivvox(3)=-1)
              if(bcoefbot*acoefbot.lt.0) tz=-bcoefbot/acoefbot  ! choose the second root when positive (first should be 0)
            endif

            ! Intersection with cell face
            t=min(tx,ty,tz)
            lopt=lopt-t*papvtotAZ(i1,j1,k1)
            if(lopt.ge.0)then
              kt=0
              if(abs(t-tx).lt.tolerance)then 
                i1=i1+int(sign(one,dir(1)))
              elseif(abs(t-ty).lt.tolerance)then
                j1=j1+int(sign(one,dir(2)))
              else
                k1=k1+int(sign(one,dir(3)))
                kt=int(sign(one,dir(3)))
              endif
           
              xf=xf+dir(1)*t
              yf=yf+dir(2)*t
              zf=zf+dir(3)*t
              
              ! Cyclic boundary conditions implementation
              if(ibc_left.eq.0.and.i1.lt.1) then
                i1=i1+n
                xf=xf+n*dx
              endif
              if(ibc_right.eq.0.and.i1.gt.n) then
                i1=i1-n
                xf=xf-n*dx
              endif
              if(ibc_bot.eq.0.and.j1.lt.1)then
                j1=j1+m
                yf=yf+m*dy
              endif
              if(ibc_top.eq.0.and.j1.gt.m) then
                j1=j1-m
                yf=yf-m*dy
              endif
            endif ! lopt.gt.0
            ! Mirror / symmetric boundary conditions
            if(ibc_left.eq.5.and.xf.le.tolerance) dir(1) = -dir(1)
            if(ibc_right.eq.5.and.xf.ge.(n*dx-tolerance)) dir(1) = -dir(1)
            if(ibc_bot.eq.5.and.yf.le.tolerance) dir(2) = -dir(2)
            if(ibc_top.eq.5.and.yf.ge.(m*dy-tolerance)) dir(2) = -dir(2)
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
! this subroutine update values of firad and fsiesrad according to 
! nphotonAZtot and Qphot
!-----------------------------------------------------------------------
subroutine updateRadArrays(Qphot)
  use gridlist_variables, only : wp,dx,dy
  use vegetation_gridlist_variables, only : nfuel
  use gridsetup, only : np,mp
  use msga_variables, only : npos,mpos
  use fuel_variables, only : lfuel
  use radiation_general_variables, only : firad,fsiesrad,papvthresh
  use radiation_montecarlo_variables, only : lmc,papvtotAZ,sourcegas, &
    nphotonAZgather,papvgas,papvift,sourcesol
  use metric_variables, only : zedge
  Implicit None

  ! Local Variables
  real(wp),intent(in) :: Qphot

  integer :: i,j,k,kcell,ift,ia,ja
  real(wp) :: vCell
  real(wp),external :: zcart
               
  ! Executable Code
  firad=-sourcegas
  fsiesrad=0._wp
  do k=1,lmc  !FPAZ
    if(k.le.nfuel+1)then
      kcell=1
    else
      kcell=k-nfuel
    endif
    do j=1,mp
      ja=(mpos-1)*mp+j
      do i=1,np
        ia=(npos-1)*np+i 
        vCell=dx*dy*(zcart(zedge(kcell+1),i,j)-zcart(zedge(kcell),i,j)) ! Cell volume (m3)
        if(papvtotAZ(ia,ja,k).gt.papvthresh)then
          firad(i,j,kcell)=firad(i,j,kcell)+Qphot* &
            nphotonAZgather(ia,ja,k)* &
            papvgas(i,j,k)/papvtotAZ(ia,ja,k)/vCell
          if(k.le.lfuel+nfuel)then
            do ift=1,nfuel
              fsiesrad(ift,i,j,kcell)=fsiesrad(ift,i,j,kcell)- &
                sourcesol(ift,i,j,k) &
                +Qphot*nphotonAZgather(ia,ja,k)*papvift(ift,i,j,k) &
                /papvtotAZ(ia,ja,k)/vCell
            enddo !do ift
          endif
        endif !papvtot2.gt.papvthresh
      enddo !do i
    enddo !do j
  enddo !do k
       
end subroutine updateRadArrays

!-----------------------------------------------------------------------
! Utility Functions
!-----------------------------------------------------------------------
subroutine zsurface(zsurf,p1,p2,z1,z2,z3,z4)
use gridlist_variables, only : wp
  Implicit None

  real(wp),intent(inout) :: zsurf

  ! Local Variables
  real(wp),intent(in) :: p1,p2
  real(wp),intent(in) :: z1,z2,z3,z4
  
  ! Executable Code  
  zsurf=((z1+z4-z2-z3)*p2+z2-z1)*p1+(z3-z1)*p2+z1

end subroutine zsurface
