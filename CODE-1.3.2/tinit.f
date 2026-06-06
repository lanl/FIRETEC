      subroutine tinit()
      use metryic
      use gridsetup
      use relax
      use xve
      use constants
      use msga

      Implicit None

      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,ia,ja,inv
      real :: zl,towi,t1
      real,external :: zcart
     
      real,allocatable:: relb(:)
      allocate (relb(nr+1))
      relb=0.
      if (irlx.ge.1.and.irly.ge.1) then
        do i=1,nr+1  
           if (irlx.ne.irly) then
             write(6,*) 'irlx and irly are not compatible'
             stop
           endif
           if (irlx.eq.1.or.irly.eq.1) then
              relb(i)=0.05*exp((-1.0*real(i))/(2.0*(real(nr-1)**2))*
     &             (real(i)-1.0/real(i))**2)
           else if (irlx.eq.2.or.irly.eq.2) then
        ! LATERAL DAMPERS WAS REFORMULATED FOLLOWING ARPS RECOMMANDATION
        ! here the formulation is such as tow=10 correspond to the recommended 
        ! value cbcdmp between dt/20 and dt/10 (here 0.15 dt) including the fact that 
        ! boundary is called twice in a time step
        ! NB : nr should should be between 5 and 7
             relb(i)=1.5*0.25/(tow*(1.0+(i/(0.5*nr-1.0))**2))
           else  if (irlx.ge.3.or.irly.ge.3) then
            relb(i)=EXP(-8.0*(real(i)-1.0)/real(nr-1))
           endif
           relb(1)=1.0 !for mass conservation
           relb(2)=1.0  ! FP : this is to account for the fact that xe assume no gradient
            !on the lateral side
         enddo
      endif ! end of definition of relb
      towi=1./tow
      relaxxv=0.
      do 22 k=1,l
      do 22 j=1,mp
      do 22 i=1,np
       zl=zcart(z(k),i,j)
       
      if (zl.ge.zab.and.ibctopopen.eq.0) then  ! FP upper dampers
       t1=amax1(0.,zl-zab)
      do 23 inv=1,4  
        if (iab.eq.1)   !initial formulation
     &     relaxxv(i,j,k,inv)=towi*towi*t1/((zb-zab)+1.e-10)
        if (iab.eq.2)  !JAS new exponential top absorber
     &     relaxxv(i,j,k,inv)=0.05*exp(-t1/(-(zb-zab)/log(towi*towi)))
        if (iab.ge.3.and.inv.le.3)  !FP
        ! here the formulation is such as tow=10 correspond to the
        ! recommended value of cfrdmp=1/(20*dt) in arps (including the
        ! fact that boundary is called twice over a time step
        !relaxxv(i,j,k,inv)=towi*0.25*(1.0-cos(pi*(zl-zab)/(zb-zab)))
        !correction to integrate the calls in higrad.f
     &   relaxxv(i,j,k,inv)=0.66*towi*0.25*(1.0-cos(pi*(zl-zab)/(zb-zab)))
   23 continue
      else if (zl.ge.zab.and.ibctopopen.eq.1) then
        ! in case of open top bc we relax theta only
        relaxxv(i,j,k,3)=0.66*towi*0.25*(1.0-cos(pi*(zl-zab)/(zb-zab)))
      endif
      if (iab.eq.4.and.zl.ge.zabt) then  !theta is relaxed lower
         relaxxv(i,j,k,4)=0.66*towi*0.25*(1.0-cos(pi*(zl-zabt)/(zb-zabt)))
      end if

      if (irlx.ge.1) then
       ia = (npos-1)*np + i
        if (ia.le.nr) then
         relaxxv(i,j,k,1)=max(relaxxv(i,j,k,1),relb(ia))
         relaxxv(i,j,k,2)=max(relaxxv(i,j,k,2),relb(ia))
         relaxxv(i,j,k,3)=max(relaxxv(i,j,k,3),relb(ia))
         relaxxv(i,j,k,4)=max(relaxxv(i,j,k,4),relb(ia))
        else if ((ia.gt.n-nr).and.(ibclatopen.eq.0)) then
         ! relaxation on the east bc is done only if ibclatopen is 0
         relaxxv(i,j,k,1)=max(relaxxv(i,j,k,1),relb(n-ia+1))
         relaxxv(i,j,k,2)=max(relaxxv(i,j,k,2),relb(n-ia+1))
         relaxxv(i,j,k,3)=max(relaxxv(i,j,k,3),relb(n-ia+1))
         relaxxv(i,j,k,4)=max(relaxxv(i,j,k,4),relb(n-ia+1))
        end if
       end if
       !if ((irly.eq.1).and.(ibclatopen.eq.0)) then
       if (irly.ge.1) then
         ! relaxation on the north and south bc is done only if ibclatopen is 0
        ja = (mpos-1)*mp + j
       if (ja.le.nr) then
        relaxxv(i,j,k,1)=max(relaxxv(i,j,k,1),relb(ja))
        relaxxv(i,j,k,2)=max(relaxxv(i,j,k,2),relb(ja))
        relaxxv(i,j,k,3)=max(relaxxv(i,j,k,3),relb(ja))
        relaxxv(i,j,k,4)=max(relaxxv(i,j,k,4),relb(ja))
       else if (ja.gt.m-nr) then
        relaxxv(i,j,k,1)=max(relaxxv(i,j,k,1),relb(m-ja+1))
        relaxxv(i,j,k,2)=max(relaxxv(i,j,k,2),relb(m-ja+1))
        relaxxv(i,j,k,3)=max(relaxxv(i,j,k,3),relb(m-ja+1))
        relaxxv(i,j,k,4)=max(relaxxv(i,j,k,4),relb(m-ja+1))
       end if
      end if
   22 continue
      deallocate (relb)
      return
      end
