cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c this subroutine compute the contribution of cell i,j,k to mass flow in direction
c (cosg, sing), where g means geostrophic wind. Mass flow is computed in an array which is 
c perpendicular to wind direction (see definition of massFlow)

c  isIni is a selector for xe (1) and xv (0)
c zfactor can either be dz (m) or a fraction of cell height (m/m)
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine addToMassFluxTemp(i, j, k, isIni, zfactor)
        use xve
        use xvo
        use gridsetup
        use msga
        use lspgf
        Implicit None
        integer:: i, j, k 
        integer:: isIni ! if 1, use xe, else use xv
        real ::zfactor  !generally a dz or a cell fraction...
        integer :: ia, ja
        real :: j0  ! index in massflux array before floor
        integer:: j0f ! index in massflux array
        real:: dj0f
        real :: rhovel1,rhovel2,rhovel
        ! iwindx=1: wind blowing mostly on x axis, rhovel1 is in cell i,j-1,k
        ! iwindx=0: wind blowing mostly on y axis, rhovel1 is in cell i-1,j,k

        if (isIni==1) then ! initial wind vel
           rhovel1 = xe(i-1+iwindx,j-iwindx,k,1) * cosg + 
     +                xe(i-1+iwindx,j-iwindx,k,2) * sing  
           rhovel2 = xe(i,j,k,1) * cosg + xe(i,j,k,2) * sing  
        else ! current windvel   
           rhovel1 = xvb(i-1+iwindx,j-iwindx,k,1) * cosg +
     +                 xvb(i-1+iwindx,j-iwindx,k,2) * sing   
           rhovel2 = xvb(i,j,k,1) * cosg + xvb(i,j,k,2) * sing 
        endif 
!       write (*,*) rhovel1,rhovel2,i,j,iwindx,rhoug,rhovg        
        ia=(npos-1)*np + i
        ja=(mpos-1)*mp + j
        
        ! j0f is the index in the massFluxTemp  with j0f = floor(j0)
        if (iwindx==1) then ! mass flux is computed for each j (as if wind aligned with x axis)
          j0 = ja - rhovg/rhoug * (ia - 1)
        else ! (iwindx==0) then ! mass flux is computed for each i (as if wind aligned with y axis)
          j0 = ia - rhoug/rhovg * (ja - 1)
        endif
        
        j0f = floor (j0)
        dj0f =j0 - j0f  ! contribution of cell i-1+iwindx,j-iwindx,k  (0
                        ! if wind is aligned with x or y axis)
        ! j0f is the index on a cyclic array due to cyclic bc (length is
        ! nMassFlux)
        if (j0f < 1) then
           j0f = j0f + nMassFlux
        else if (j0f > nMassFlux) then
           j0f = j0f - nMassFlux
        endif 
        !write (*,*)  ia,ja,j0,j0f,dj0f       
 
        rhovel = rhovel1 * dj0f + rhovel2 * (1.0 - dj0f)
        massFluxTemp(j0f) = massFluxTemp(j0f) + zfactor * dx * dy * rhovel

        return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     fluxAlongWindDir compute for any i,j in a given subdomain, the appropriate
c     massFlux from the massFlux array
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function fluxAlongWindDir(i, j)
        use xve
        use xvo
        use gridsetup
        use msga
        use lspgf
        Implicit None
        integer:: i, j
        integer :: ia, ja
        real :: j0  ! index in massflux array before floor
        real:: j0f, j0fp1 ! index in massflux array
        real:: dj0f

        ia=(npos-1)*np + i
        ja=(mpos-1)*mp + j
        ! here we compute index corresponding to i,j in the massFlux array
        if (iwindx==1) then ! mass flux is computed for each j (as if wind aligned with x axis)
           j0 = ja - rhovg/rhoug * (ia - 1)
        else ! (iwindx==0) then ! mass flux is computed for each i (as if wind aligned with y axis)
           j0 = ia - rhoug/rhovg * (ja - 1)
        endif
        ! interpolation between j0f and j0fp1
        j0f = floor (j0)
        dj0f =j0 - j0f
        if (j0f < 1) then
          j0f = j0f + nMassFlux
        else if (j0f > nMassFlux) then
          j0f = j0f - nMassFlux
        endif
        j0fp1 = j0f + 1
        if (j0fp1 > nMassFlux) then
          j0fp1 = j0fp1 - nMassFlux
        endif
        fluxAlongWindDir = massFlux(j0f) * (1-dj0f) + massFlux(j0fp1) * dj0f
        return
      end function fluxAlongWindDir





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  this subroutine called by compress can increase or decrease the pressure
c gradient to ensure a convergence of u
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine largeScalePGF(it)
       use xve 
       use pres
       use gridsetup
       use msga
       use xvo
       use metryic
       use lspgf
       Implicit None
       integer, intent(in):: it
       integer :: i,j,k,ia,ja,ierr,ihalffootprint,jhalffootprint,ii,jj,ifp2,jfp
       real,allocatable::flspgftmp(:,:,:) ! for filtering flspgf in ilspgf.eq.3
       real,allocatable::flspgf2(:,:,:) ! for filtering flspgf in ilspgf.eq.3
       real::normcoeff !normalisation coefficient for gaussian filtering
       real,external::zcart, fluxAlongWindDir
       real::sdi2, sdj2 !square standard deviation for gaussian filtering (ilspgf=3)
       real:: currentFlux
       real:: z1,z2,z12,zcoef
       real :: deltaf
       real:: oldTotMassFlux
       integer :: ilspgfstart
       oldTotMassFlux = totMassFlux
c    here we compute temporary mass flux in wind direction for a subdomain
       massFluxTemp = 0.0
       do k=1,l
        do j=1,mp
         do i=1,np
            z1=zcart(zedge(k),i,j)-zs(i,j)
            z2=zcart(zedge(k+1),i,j)-zs(i,j)
            if (izlspgf.eq.0) then ! massFlux computed in the whole domain
              call addToMassFluxTemp(i, j, k, 0, z2-z1)
            else ! massFlux computed at ref height zu 
              ! compute the contributon of cell (i,j,k)  at height zu
              ! z1 is cell center of cell k-1 (or domain bottom)
              if (k.eq.1) then
                z1=0.0
              else
                z1=zcart(z(k-1),i,j)-zs(i,j)
              endif
              z2=zcart(z(min(k+1,l)),i,j)-zs(i,j)
              if (z1 <= zu .and. zu <= z2) then 
              ! cell i,j,k should contribute to massFlux
                z12 = zcart(z(k),i,j)-zs(i,j)
                ! zcoef is the weight of the cell to massFlux
                if (zu>=z12) then ! zu between z12 and z2
                   zcoef = (z2 - zu) / (z2 - z12)
                else ! zu between z1 and z12
                   zcoef = (zu - z1) / (z12 - z1)
                endif 
                call addToMassFluxTemp(i, j, k, 0, zcoef)
               endif
            end if  ! izlspgf.eq.1
         enddo
        enddo
       enddo
c reduction to all domains
       call mpi_allreduce(massFluxTemp,massFlux,nMassFlux,mpi_real,
     +                  mpi_sum,mpi_comm_world,ierr)

c sum of the mass flux
       totMassFlux=0.0
       do j=1,nMassFlux
         totMassFlux = totMassFlux + massFlux(j)
       enddo
       ilspgfstart=frqlspgf
       !if (mpi_rank.eq.0) then 
       !  write(6,*) 'it, ilspgfstart ',it,ilspgfstart
       !endif
       if (it.le.ilspgfstart) then ! we wait for turbulence development bef update
       if (mpi_rank.eq.0) then 
          write(6,*) 'ilspgf:current and target massflux ',totMassFlux,targMassFlux
          write(6,*) '   no update of the forcing yet'
       endif
       else ! update of flspgf
       if (mpi_rank.eq.0) then 
         write(6,*) 'ilspgf:old,current and target massflux ',oldTotMassFlux,totMassFlux,targMassFlux
          write(6,*) '   current flspgf is ',flspgf(1,1)
       endif  
       if (ilspgf.eq.1) then 
         deltaf = (targMassFlux-(totMassFlux+(totMassFlux-oldTotMassFlux)))
     +            /(dx*dy*n*m*sqrt(rhoug**2.+rhovg**2.)*tau*intsintheta)
         flspgf = flspgf+deltaf
         if (totMassFlux.le.targMassFlux) then
           if (mpi_rank.eq.0) write(6,*) 'pressure gradient increased'
         else
           if (mpi_rank.eq.0) write(6,*) 'pressure gradient decreased'
         endif      
       elseif (ilspgf.ge.2) then
         do j=1,mp
         do i=1,np
            currentFlux = fluxAlongWindDir(i, j)
            deltaf = (targMassFlux-nMassFlux*currentFlux)
     +            /(dx*dy*n*m*sqrt(rhoug**2.+rhovg**2.)*tau*intsintheta)
            flspgf(i,j) = flspgf(i,j) + deltaf
            !if (currentFlux.le.targMassFlux/nMassFlux) then
            !  flspgf(i,j)=1.1*flspgf(i,j)
            !else
            !  flspgf(i,j)=0.9*flspgf(i,j)
            !endif
        enddo
       enddo 
      endif 
      if (ilspgf.eq.3) then 
c     ! gaussian filtering of flspgf array
      allocate (flspgftmp(np,mp,1))  ! copy of flspgf on 1:np,1:mp
      do j=1,mp
        do i=1,np
           flspgftmp(i,j,1)=flspgf(i,j)
        enddo
      enddo
      allocate (flspgf2(n,m,1))
      call allgather3d(flspgftmp,flspgf2,1)
      do j=1,mp
        do i=1,np
          ia = (npos-1)*np + i
          ja = (mpos-1)*mp + j
           if (flspgftmp(i,j,1).ne.flspgf2(ia,ja,1))then 
         write(6,*) 'allgather3d pb:rank=',mpi_rank,i,j,ia,ja,flspgftmp(i,j,1),flspgf2(ia,ja,1)
        endif
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
        normcoeff=normcoeff+
     +     1./sqrt(2.0*3.14159027*sdi2)*exp(-real(ifp2*ifp2)/(2.0*sdi2))*
     +     1./sqrt(2.0*3.14159027*sdj2)*exp(-real(jfp*jfp)/(2.0*sdj2))
      enddo
      enddo
      normcoeff=1.0/normcoeff
      !write(6,*) 'normcoeff=',normcoeff
       do j=1,mp
         do i=1,np
!          if (k.eq.1.and.i.eq.1)
!     +     write(6,*) 'j,sintheta',j,sintheta(i,j,k)
          flspgf(i,j)=0.0
          ia = (npos-1)*np + i
          ja = (mpos-1)*mp + j
          do ifp2=-ihalffootprint,ihalffootprint
          do jfp=-jhalffootprint,jhalffootprint
            ii=ia+ifp2
            if (ii.lt.1) ii=ii+n
            if (ii.gt.n) ii=ii-n
            jj=ja+jfp
            if (jj.lt.1) jj=jj+m
            if (jj.gt.m) jj=jj-m
            flspgf(i,j)=flspgf(i,j)+flspgf2(ii,jj,1)*
     +       normcoeff/sqrt(2.0*3.14159027*sdi2)*exp(-real(ifp2*ifp2)/(2.0*sdi2))
     +       *1.0/sqrt(2.0*3.14159027*sdj2)*exp(-real(jfp*jfp)/(2.0*sdj2))
          enddo
          enddo
!          if (k.eq.1.and.i.eq.1)
!     +     write(6,*) 'j,sinthetaf',j,sintheta(i,j,k)
          enddo
c         if (mpi_rank.eq.0)
c     +    write(6,*) ja,xe1j(ja),xe1j2(ja),xe1j2f(ja)
         enddo
        endif ! ilspgf.eq.3 (end filtering flspgf)
      endif !it.ge.itlspgfstart
      do k=1,l
       do j=1,mp
        do i=1,np
          sinthetaf(i,j,k)=sintheta(i,j,k)*flspgf(i,j)
        enddo
       enddo
      enddo 
       return
       end
