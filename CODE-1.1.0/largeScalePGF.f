ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  this subroutine called by compress can increase or decrease the pressure
c gradient to ensure a convergence of u
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine largeScalePGF()
       use xve 
       use pres
       use gridsetup
       use msga
       use xvo
      use metryic
       Implicit None
       integer :: i,j,k,ierr,ia,ja,jhalffootprint,jj,jfp
       real ::xe1j2(m), xe2i2(n)  !mass flux arrays
       real,allocatable::sinthetatmp(:,:,:) ! for filtering sintheta in ilspgf.eq.4
       real,allocatable::sintheta2(:,:,:) ! for filtering sintheta in ilspgf.eq.4
       real::normcoeff !normalisation coefficient for gaussian filtering
       real ::xe1j2temp(m), xe2i2temp(n) !for gathering data
       real ::xe1tot2,xe2tot2,dzcell
       real,external::zcart
       real::sd2 !square standard deviation for gaussian filtering (ilspgf=4)
       real::coswind,sinwind ! cos and sin of geostrophic wind direction 
       real::fluxi,fluxc ! initial and current flux in wind direction 
c    here we compute temporary mass flux though boundary for largeScalePGF
        xe1j2=0.0
        xe1j2temp=0.0
c      if (leftedge.eq.1) then
        do k=1,l
         do j=1,mp
         do i=1,np
            dzcell=zcart(zedge(k+1),i,j)-
     +         zcart(zedge(k),i,j)
            ja=(mpos-1)*mp + j
            xe1j2temp(ja)=xe1j2temp(ja)+
     +          xvb(i,j,k,1)*dzcell
         enddo
       enddo
       enddo
c      endif
        call mpi_allreduce(xe1j2temp,xe1j2,m,mpi_real,
     +                  mpi_sum,mpi_comm_world,ierr)
         xe1tot2=0.0
         do ja=1,m
            xe1tot2=xe1tot2+xe1j2(ja)
         enddo
        xe2i2=0.0 
        xe2i2temp=0.0 
      ! if (botedge.eq.1) then
        do k=1,l
         do j=1,mp
         do i=1,np
            ia=(npos-1)*np + i
            dzcell=zcart(zedge(k+1),i,j)-
     +         zcart(zedge(k),i,j)
          xe2i2temp(ia)=xe2i2temp(ia)+xvb(i,j,k,2)
     +         *dzcell
         enddo
         enddo
       enddo
      !endif
        call mpi_allreduce(xe2i2temp,xe2i2,m,mpi_real,
     +                  mpi_sum,mpi_comm_world,ierr)
         xe2tot2=0.0
         do ia=1,n
            xe2tot2=xe2tot2+xe2i2(ia)
         enddo

        if (mpi_rank.eq.0) then 
       write(6,*) 'ilspgf:total massflux to west and south boundary',xe1tot2,xe2tot2
       write(6,*) 'initial massflux to west and south boundary',xe1tot,xe2tot
       write(6,*) 'ilspgf: current sintheta (k=1)',sintheta(1,1,1)

       endif  
       if (ilspgf.eq.2) then 
c        rrl 's modification        
c           sintheta=sintheta*xe1tot/xe1tot2
          coswind=xe(1,1,l,1)/sqrt(xe(1,1,l,1)**2+xe(1,1,l,2)**2)
          sinwind=xe(1,1,l,2)/sqrt(xe(1,1,l,1)**2+xe(1,1,l,2)**2)
           fluxc=xe1tot2*coswind+xe2tot2*sinwind
           fluxi=xe1tot*coswind+xe2tot*sinwind
           if (mpi_rank.eq.0) write(6,*) 'initial and current fluxes:',fluxi,fluxc
          if (fluxc.le.fluxi) then
           sintheta=1.1*sintheta
           !sintheta=1.5*sintheta
           if (mpi_rank.eq.0) write(6,*) 'pressure gradient increased'
         else
           sintheta=0.9*sintheta
           !sintheta=0.67*sintheta
           if (mpi_rank.eq.0) write(6,*) 'pressure gradient decreased'
         endif      
       elseif (ilspgf.ge.3) then
         
        do k=1,l
         do j=1,mp
         do i=1,np
            ja=(mpos-1)*mp + j
c        rrl 's modification        
c             sintheta(i,j,k)=xe1j(ja)*sintheta(i,j,k)/xe1j2(ja)
            if (xe1j2(ja).le.xe1j(ja)) then 
             sintheta(i,j,k)=1.1*sintheta(i,j,k)
            else
             sintheta(i,j,k)=0.9*sintheta(i,j,k)
            endif
         enddo
        enddo
       enddo 
      endif    
      if (ilspgf.eq.4) then 
c     ! gaussian filtering :
      allocate (sinthetatmp(np,mp,l))  ! copy of sintheta on 1:np,1:mp
      do k=1,l
       do j=1,mp
        do i=1,np
          sinthetatmp(i,j,k)=sintheta(i,j,k)
        enddo
       enddo
      enddo
      allocate (sintheta2(n,m,l))
      call allgather3d(sinthetatmp,sintheta2,l)
      ! footprint of filter m/5
      jhalffootprint=floor(m/5*0.5)
      ! standard deviation square:
      sd2=0.1*jhalffootprint**2.0
      normcoeff=0.0
      do jfp=-jhalffootprint,jhalffootprint
        normcoeff=normcoeff+
     +     1./sqrt(2.0*3.14159027*sd2)*exp(-real(jfp*jfp)/(2.0*sd2))
      enddo
      normcoeff=1.0/normcoeff
      !write(6,*) 'normcoeff=',normcoeff
      do k=1,l
       do j=1,mp
         do i=1,np
!          if (k.eq.1.and.i.eq.1)
!     +     write(6,*) 'j,sintheta',j,sintheta(i,j,k)
          sintheta(i,j,k)=0.0
          ja = (mpos-1)*mp + j
          do jfp=-jhalffootprint,jhalffootprint
            jj=ja+jfp
            if (jj.lt.1) jj=jj+m
            if (jj.gt.m) jj=jj-m
            sintheta(i,j,k)=sintheta(i,j,k)+sintheta2(i,jj,k)*
     +   normcoeff/sqrt(2.0*3.14159027*sd2)*exp(-real(jfp*jfp)/(2.0*sd2))
          enddo
!          if (k.eq.1.and.i.eq.1)
!     +     write(6,*) 'j,sinthetaf',j,sintheta(i,j,k)
          enddo
c         if (mpi_rank.eq.0)
c     +    write(6,*) ja,xe1j(ja),xe1j2(ja),xe1j2f(ja)
         enddo
         enddo
        endif
       return
       end
