!-----------------------------------------------------------------------
! rinitLargeScalePGF sets up the large scale pressure gradient force 
! used with cyclic boundary conditions to maintain the wind's 
! driving force
!-----------------------------------------------------------------------
subroutine rinitLargeScalePGF
  use gridlist_variables, only : n,m,l,zgroundref,zab,izlspgf, &
    zu,irst,prec
  use lspgf_variables, only : iwindx,nMassFlux,massFluxTemp, &
    massFlux,targMassFlux,sintheta,sinthetaf,flspgf
  use msga_variables, only : mpi_rank,mpi_sum,ierror,mpi_comm_world
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use metric_variables, only : z,zs,zedge
  use gridsetup, only : np,mp
  use constants, only : pi,fcor3
  use xvall, only : xe,iuvel,ivvel
  use zcart_function
  Implicit none
  
  ! Local Variables  
  integer :: i,j,k
  real(prec) :: cosg,sing
  real(prec) :: intsintheta,intsinthetal,zla,gam,dzcell
  real(prec) :: z1,z2,z12,zcoef
  real(prec) :: targMassFluxSum
  !real(prec),external :: zcart

  ! Executable Code 
  ! compute wind direction
  cosg=xe(1,1,l,iuvel)/sqrt(xe(1,1,l,iuvel)**2+xe(1,1,l,ivvel)**2)
  sing=xe(1,1,l,ivvel)/sqrt(xe(1,1,l,iuvel)**2+xe(1,1,l,ivvel)**2)
     
  ! definition of massFlow integrals
  if(abs(xe(1,1,l,iuvel)*m).ge.abs(xe(1,1,l,ivvel)*n))then ! massFlow computed along y axis (as if wind aligned with xaxis)
    iwindx=1
    nMassFlux=m
    if(mpi_rank.eq.0) print*,'mass flux computed along yaxis over ' &
      ,nMassFlux,' cells'
  else  ! massFlux computed along x axis (as if wind aligned with y axis)
    iwindx=0
    nMassFlux = n
    if(mpi_rank.eq.0) print*, 'mass flux computed along xaxis over ' &
      ,nMassFlux,' cells'
  endif
  allocate(massFluxTemp(nMassFlux)); massFluxTemp=0.
  allocate(massFlux(nMassFlux)); massFlux=0.

  intsinthetal = 0.0 ! integrated value of sintheta over the vertical 
  ! definition of sintheta (vertical profile of the lspgf extrapolated from ekmann spiral) 
  do k=1,l
    do j=1,mp
      do i=1,np
        zla=zcart(z(k),i,j)-zgroundref !FP took zgroundref as a reference for pressure gradient
         ! LSPGF IS NOT COMPATIBLE WITH ROTATED GRAVITY BECAUSE W IS ) AT THE TOP
         ! when gravity is rotated, zla vary with slopeangle and slope azimuth:
        gam = pi/zab*zla
        sintheta(i,j,k)=0.5*exp(-gam)*sin(gam)/ &
          sqrt(1.0+exp(-2.0*gam)-2.0*exp(-gam)*cos(gam))
        dzcell=zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)
        if(izlspgf.eq.0) intsinthetal=intsinthetal+ &
          sintheta(i,j,k)*dzcell/(n*m)
      enddo
    enddo
  enddo
  if(izlspgf.eq.0)then  ! mean of integral over z of sintheta*dz
#ifdef DBL_PREC
    call mpi_allreduce(intsinthetal,intsintheta,1,mpi_double_precision, &
      mpi_sum,mpi_comm_world,ierror)
#else
    call mpi_allreduce(intsinthetal,intsintheta,1,mpi_real, &
      mpi_sum,mpi_comm_world,ierror)
#endif
  else !izlspgf.eq.1  : intsintheta=sintheta(zu)
    gam=pi/zab*(zu-zgroundref)
    intsintheta=0.5*exp(-gam)*sin(gam)/ &
      sqrt(1.+exp(-2.*gam)-2.*exp(-gam)*cos(gam))
  endif

  ! definition of flspgf and sinthetaf
  if(irst.eq.0) then ! definition of flspgf (n*m) array
    flspgf=8.*fcor3  ! initial value 
    do k=1,l
      do j=1,mp
        do i=1,np
          sinthetaf(i,j,k)=sintheta(i,j,k)*flspgf(i,j) 
        enddo
      enddo
    enddo
  else ! compute flspgf from sinthetaf in case of restart (from io.f)
    do j=1,mp
      do i=1,np
        flspgf(i,j) = sinthetaf(i,j,1)/sintheta(i,j,1)
      enddo
    enddo
  endif ! irst.eq.0
  if(mpi_rank.eq.0) print*,'ilspgf: flspgfini',flspgf(1,1)

  ! test of value of zu (should be between bottom and top of the
  ! domain for izlspgf=1
  if(izlspgf.eq.1) then
    if(zu<=zcart(zedge(2),i,j)-zs(i,j).or. &
      (zu>=zcart(zedge(l),i,j)-zs(i,j))) then
      print*,'izlspgf=1 incompatible with zu that should be &
        &above cell 1 and below cell l-1'
      stop
    endif
  endif
  ! computation of inital massFluxTemp on all sub domains (reduced with a mpi sum in array
  ! massFluxTarget later). Something similar is done in largeScalePGF to compute current value of
  ! massFlux (to compare it to the target define here).
  do k=1,l
    do j=1,mp
      do i=1,np
        z1=zcart(zedge(k),i,j)-zs(i,j)
        z2=zcart(zedge(k+1),i,j)-zs(i,j)
        if (izlspgf.eq.0) then ! massFlux on the whole domain
          call addToMassFluxTemp(i,j,k,1,z2-z1)
        else ! massFlux at a given height (zu) 
          ! compute the contributon of cell (i,j,k)  at height zu
          ! z1 is cell center of cell k-1 (or domain bottom)
          if(k.eq.1) then
            z1=0.0
          else
            z1=zcart(z(k-1),i,j)-zs(i,j)
          endif
          z2=zcart(z(min(k+1,l)),i,j)-zs(i,j)
          if(z1<=zu.and.zu<=z2)then
          ! cell i,j,k should contribute to massFlux
            z12 = zcart(z(k),i,j)-zs(i,j)
            ! zcoef is the weight of the cell to massFlux
            if(zu>=z12)then ! zu is between z12 and z2
              zcoef=(z2-zu)/(z2-z12)
            else ! zu is between z1 and z12
              zcoef=(zu-z1)/(z12-z1)
            endif
            call addToMassFluxTemp(i,j,k,1,zcoef)
          endif
        endif  ! izlspgf.eq.1
      enddo
    enddo
  enddo
  targMassFluxSum=sum(massFluxTemp)
#ifdef DBL_PREC
  call mpi_allreduce(targMassFluxSum,targMassFlux,1,mpi_double_precision, &
    mpi_sum,mpi_comm_world,ierror)
#else
  call mpi_allreduce(targMassFluxSum,targMassFlux,1,mpi_real, &
    mpi_sum,mpi_comm_world,ierror)
#endif
  if(mpi_rank.eq.0) print*,'total massflux initial:',targMassFlux
       
end subroutine rinitLargeScalePGF

!----------------------------------------------------------------
! addMToMassFluxTemp computes the contribution of cell i,j,k to mass flow in direction
! (cosg, sing), where g means geostrophic wind. Mass flow is computed in an array which is 
! perpendicular to wind direction (see definition of massFlow)
!
! isIni is a selector for xe (1) and xv (0)
! zfactor can either be dz (m) or a fraction of cell height (m/m)
!----------------------------------------------------------------
subroutine addToMassFluxTemp(i,j,k,isIni,zfactor)
  use gridlist_variables, only : l,dx,dy,prec
  use lspgf_variables, only : iwindx,cosg,sing,nMassFlux, &
    massFluxTemp
  use gridsetup, only : np,mp
  use msga_variables, only : npos,mpos
  use xvall, only : xe,xv,iuvel,ivvel
  Implicit None
  ! Local Variables
  integer,intent(in) :: i,j,k 
  integer,intent(in) :: isIni ! if 1, use xe, else use xv
  real(prec),intent(in) ::zfactor  !generally a dz or a cell fraction...
  integer :: ia,ja
  integer :: j0f
  real(prec) :: j0,dj0f  ! index in massflux array
  real(prec) :: rhovel1,rhovel2,rhovel
  ! iwindx=1: wind blowing mostly on x axis, rhovel1 is in cell i,j-1,k
  ! iwindx=0: wind blowing mostly on y axis, rhovel1 is in cell i-1,j,k

  ! Executable Code
  if(isIni==1)then ! initial wind vel
    rhovel1=xe(i-1+iwindx,j-iwindx,k,iuvel)*cosg+ &
      xe(i-1+iwindx,j-iwindx,k,ivvel)*sing  
    rhovel2=xe(i,j,k,iuvel)*cosg+xe(i,j,k,ivvel)*sing  
  else ! current windvel   
    rhovel1=xv(i-1+iwindx,j-iwindx,k,iuvel)*cosg+ &
      xv(i-1+iwindx,j-iwindx,k,ivvel)*sing
    rhovel2=xv(i,j,k,iuvel)*cosg+xv(i,j,k,ivvel)*sing 
  endif 
  ia=(npos-1)*np+i
  ja=(mpos-1)*mp+j
        
  ! j0f is the index in the massFluxTemp  with j0f = floor(j0)
  if(iwindx==1)then ! mass flux is computed for each j (as if wind aligned with x axis)
    j0=ja-xe(1,1,l,ivvel)/xe(1,1,l,iuvel)*(ia-1)
  else ! (iwindx==0) then ! mass flux is computed for each i (as if wind aligned with y axis)
    j0=ia-xe(1,1,l,iuvel)/xe(1,1,l,ivvel)*(ja-1)
  endif
        
  j0f=floor(j0)
  dj0f=j0-j0f  ! contribution of cell i-1+iwindx,j-iwindx,k  (0 if wind is aligned with x or y axis)
               ! j0f is the index on a cyclic array due to cyclic bc (length is nMassFlux)
  if(j0f.lt.1)then
    j0f=j0f+nMassFlux
  elseif(j0f.gt.nMassFlux)then
    j0f=j0f-nMassFlux
  endif 
  rhovel=rhovel1*dj0f+rhovel2*(1.-dj0f)
  massFluxTemp(j0f)=massFluxTemp(j0f)+zfactor*dx*dy*rhovel

end subroutine addToMassFluxTemp

!-----------------------------------------------------------------------
! this subroutine called by compress can increase or decrease the 
! pressure gradient to ensure a convergence of u
!-----------------------------------------------------------------------
subroutine largeScalePGF
  use gridlist_variables, only : l,izlspgf,zu,frqlspgf,ilspgf, &
    dx,dy,n,m,nprocx,nprocy,prec
  use constants, only : pi
  use gridsetup, only : np,mp,ittot,dt
  use metric_variables, only : z,zs,zedge
  use msga_variables, only : mpi_rank,mpi_sum,mpi_comm_world,ierror, &
    npos,mpos,numprocs
#ifdef DBL_PREC
  use msga_variables, only : mpi_double_precision
#else
  use msga_variables, only : mpi_real
#endif
  use lspgf_variables, only : totMassFlux,massFluxTemp,massFLux, &
    nMassFlux,flspgf,intsintheta,iwindx,targMassFlux,sinthetaf, &
    sintheta
  use xvall, only : xe,iuvel,ivvel
  use zcart_function
  Implicit None

  ! Implicit None
  integer :: i,j,k
  integer :: ia,ja,iprocx,jprocy,iproc
  integer :: ihalffootprint,jhalffootprint,ii,jj,ifp2,jfp
  integer :: j0f,j0fp1
  real(prec) :: j0,dj0f  ! index in massflux array
  real(prec) :: oldTotMassFlux
  real(prec) :: z1,z2,z12,zcoef
  real(prec) :: deltaf
  real(prec) :: currentFlux
  real(prec) :: normcoeff !normalisation coefficient for gaussian filtering
  real(prec)::sdi2, sdj2 !square standard deviation for gaussian filtering (ilspgf=3)
  !real(prec),external :: zcart
  real(prec),allocatable :: flspgftmp(:,:,:) ! for filtering flspgf in ilspgf.eq.3
  real(prec),allocatable :: flspgf2(:,:) ! for filtering flspgf in ilspgf.eq.3

  ! Executable Code
  oldTotMassFlux=totMassFlux
  ! here we compute temporary mass flux in wind direction for a subdomain
  massFluxTemp=0.0
  do k=1,l
    do j=1,mp
      do i=1,np
        if(izlspgf.eq.0)then ! massFlux computed in the whole domain
          z1=zcart(zedge(k),i,j)-zs(i,j)
          z2=zcart(zedge(k+1),i,j)-zs(i,j)
          call addToMassFluxTemp(i,j,k,0,z2-z1)
        else ! massFlux computed at ref height zu compute the contributon of cell (i,j,k) at height zu z1 is cell center of cell k-1 (or domain bottom)
          z1=zcart(z(max(k-1,1)),i,j)-zs(i,j)
          z2=zcart(z(min(k+1,l)),i,j)-zs(i,j)
          if(z1.le.zu.and.zu.le.z2)then 
            z12 = zcart(z(k),i,j)-zs(i,j)
            if(zu.ge.z12)then ! zu between z12 and z2
              zcoef=(z2-zu)/(z2-z12) ! zcoef is the weight of the cell to massFlux
            else ! zu between z1 and z12
              zcoef=(zu-z1)/(z12-z1)
            endif 
            call addToMassFluxTemp(i,j,k,0,zcoef)
          endif
        endif  ! izlspgf.eq.1
      enddo
    enddo
  enddo
! reduction to all domains
  do i=1,nMassFlux
#ifdef DBL_PREC
    call mpi_allreduce(massFluxTemp(i),massFlux(i),1,mpi_double_precision, &
      mpi_sum,mpi_comm_world,ierror)
#else
    call mpi_allreduce(massFluxTemp(i),massFlux(i),1,mpi_real, &
      mpi_sum,mpi_comm_world,ierror)
#endif
  enddo

! sum of the mass flux
  totMassFlux=0.0
  do j=1,nMassFlux
    totMassFlux=totMassFlux+massFlux(j)
  enddo
  if (ittot.le.frqlspgf) then ! we wait for turbulence development bef update
    if (mpi_rank.eq.0) then 
      print*,'ilspgf:current and target massflux ',totMassFlux,targMassFlux
      print*,'   no update of the forcing yet'
    endif
  else ! update of flspgf
    if(mpi_rank.eq.0)then 
      print*,'ilspgf:old,current and target massflux ',oldTotMassFlux,totMassFlux,targMassFlux
      print*,'   current flspgf is ',flspgf(1,1)
    endif  
    if(ilspgf.eq.1)then 
      deltaf=(targMassFlux-(totMassFlux+(totMassFlux-oldTotMassFlux))) &
        /(dx*dy*n*m*sqrt(xe(1,1,1,iuvel)**2.+xe(1,1,1,ivvel)**2.)*frqlspgf*dt*intsintheta)
      if(mpi_rank.eq.0) print*,'new deltaf, 1/10 flspgf:',deltaf,0.1*flspgf(1,1)
      flspgf=flspgf+deltaf
      if(totMassFlux.le.targMassFlux)then
        if(mpi_rank.eq.0) print*,'pressure gradient increased'
      else
        if(mpi_rank.eq.0) print*,'pressure gradient decreased'
      endif
    elseif(ilspgf.ge.2)then
      do j=1,mp
        do i=1,np
          ia=(npos-1)*np+i
          ja=(mpos-1)*mp+j
          ! here we compute index corresponding to i,j in the massFlux array
          if(iwindx==1)then ! mass flux is computed for each j (as if wind aligned with x axis)
            j0=ja-xe(1,1,1,ivvel)/xe(1,1,1,iuvel)*(ia-1)
          else ! (iwindx==0) then ! mass flux is computed for each i (as if wind aligned with y axis)
            j0=ia-xe(1,1,1,iuvel)/xe(1,1,1,ivvel)*(ja-1)
          endif
          ! interpolation between j0f and j0fp1
          j0f=floor (j0)
          dj0f=j0-j0f
          if(j0f.lt.1)then
            j0f=j0f+nMassFlux
          elseif(j0f.gt.nMassFlux)then
            j0f=j0f-nMassFlux
          endif
          j0fp1=j0f+1
          if(j0fp1.gt.nMassFlux)then
            j0fp1=j0fp1-nMassFlux
          endif
          currentFlux=massFlux(j0f)*(1-dj0f)+massFlux(j0fp1)*dj0f
          deltaf=(targMassFlux-nMassFlux*currentFlux) &
            /(dx*dy*n*m*sqrt(xe(1,1,1,iuvel)**2.+xe(1,1,1,ivvel)**2.)*frqlspgf*dt*intsintheta)
          flspgf(i,j)=flspgf(i,j)+deltaf
        enddo
      enddo 
    endif 
    if(ilspgf.eq.3)then ! gaussian filtering of flspgf array
      allocate(flspgftmp(np,mp,numprocs),flspgf2(n,m))
      flspgftmp = 0.0
      flspgf2 = 0.0
      do iproc=1,numprocs
#ifdef DBL_PREC
        call mpi_gather(flspgf,np*mp,mpi_double_precision,flspgftmp, &
          np*mp,mpi_double_precision,iproc-1,mpi_comm_world,ierror)
#else
        call mpi_gather(flspgf,np*mp,mpi_real,flspgftmp,np*mp,mpi_real, &
          iproc-1,mpi_comm_world,ierror)
#endif
      enddo
      do iprocx=1,nprocx
        do jprocy=1,nprocy
          iproc=iprocx+(jprocy-1)*nprocx
          do j=1,mp
            do i=1,np
              ia=(iprocx-1)*np+i
              ja=(jprocy-1)*mp+j
              flspgf2(ia,ja)=flspgftmp(i,j,iproc)
            enddo
          enddo
        enddo
      enddo
      ! footprint of filter n/5
      ihalffootprint=floor(n/5*0.5)
      ! footprint of filter m/5
      jhalffootprint=floor(m/5*0.5)
      ! standard deviation square:
      sdi2=0.1*ihalffootprint**2.0
      sdj2=0.1*jhalffootprint**2.0
      normcoeff=0.0
      do ifp2=-ihalffootprint,ihalffootprint
        do jfp=-jhalffootprint,jhalffootprint
          normcoeff=normcoeff &
#ifdef DBL_PREC
            +1./sqrt(2.*pi*sdi2)*exp(-dble(ifp2*ifp2)/(2.*sdi2)) &
            *1./sqrt(2.*pi*sdj2)*exp(-dble(jfp*jfp)/(2.*sdj2))
#else
            +1./sqrt(2.*pi*sdi2)*exp(-real(ifp2*ifp2)/(2.*sdi2)) &
            *1./sqrt(2.*pi*sdj2)*exp(-real(jfp*jfp)/(2.*sdj2))
#endif
        enddo
      enddo
      normcoeff=1./normcoeff
      do j=1,mp
        do i=1,np
          flspgf(i,j)=0.0
          ia=(npos-1)*np+i
          ja=(mpos-1)*mp+j
          do ifp2=-ihalffootprint,ihalffootprint
            do jfp=-jhalffootprint,jhalffootprint
              ii=ia+ifp2
              if(ii.lt.1) ii=ii+n
              if(ii.gt.n) ii=ii-n
              jj=ja+jfp
              if(jj.lt.1) jj=jj+m
              if(jj.gt.m) jj=jj-m
              flspgf(i,j)=flspgf(i,j)+flspgf2(ii,jj) &
#ifdef DBL_PREC
                *normcoeff/sqrt(2.*pi*sdi2)*exp(-dble(ifp2*ifp2)/(2.*sdi2)) &
                *1./sqrt(2.*pi*sdj2)*exp(-dble(jfp*jfp)/(2.*sdj2))
#else
                *normcoeff/sqrt(2.*pi*sdi2)*exp(-real(ifp2*ifp2)/(2.*sdi2)) &
                *1./sqrt(2.*pi*sdj2)*exp(-real(jfp*jfp)/(2.*sdj2))
#endif
            enddo
          enddo
        enddo
      enddo
      deallocate(flspgftmp,flspgf2)
    endif ! ilspgf.eq.3 (end filtering flspgf)
  endif !it.ge.frqlspgf
      
  do k=1,l
    do j=1,mp
      do i=1,np
        sinthetaf(i,j,k)=sintheta(i,j,k)*flspgf(i,j)
      enddo
    enddo
  enddo 
       
end subroutine largeScalePGF
