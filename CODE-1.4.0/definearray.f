      subroutine definearray
      use gridsetup
      use constants
      use bc
      use workavg
      use forcinner
      use metryic
      use relax
      use radiation
      use turba
      use turbb
      use fireteca
      use xvi
      use xvo
      use xve
      use advo
      use pres
      use lspgf
      use weights
      use nonlocal
      use xvbin
      use updatedFields
      Implicit None

      allocate (uavg(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (vavg(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (wavg(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (oavg(1-ih:np+ih,1-ih:mp+ih,l))
      ! for f1avg, f2avg and f3avg, ghost cells are required
      ! because of call to donorcell in force...
      allocate (f1avg(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (f2avg(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (f3avg(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (f1(np,mp,l))
      allocate (f2(np,mp,l))
      allocate (f3(np,mp,l))
      allocate (rg_over_prrcp(np,mp,l)) !array for higrad
      allocate (cp_over_cv(np,mp,l)) !arrays for higrad
      allocate (c13(1-ih:np+ih,1-ih:mp+ih))
      allocate (c23(1-ih:np+ih,1-ih:mp+ih))
      allocate (gi(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (gmul(l))
      allocate (h(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (zs(1-ih:np+ih,1-ih:mp+ih))
      allocate (zsio(n,m))
      allocate (x(1-ih:n+ih))
      allocate (y(1-ih:m+ih))
      allocate (z(l))
      if(ibctopbot.eq.1) allocate (zedge(l+1))
c     allocate (tau(l,1-ih:np+ih,1-ih:mp+ih))
      allocate (relaxxv(1-ih:np+ih,1-ih:mp+ih,l,4))
c     allocate (relx(1-ih:np+ih)) 
c     allocate (rely(1-ih:mp+ih))
      !allocate(zcoords(1-ih:np+ih,1-ih:mp+ih,L));     zcoords=0.
      !allocate(zcoordsf(1-ih:np+ih,1-ih:mp+ih,L+1));     zcoordsf=0.
      !if(irad.GE.1.and.irod.eq.1) then             !KOO eq->GE
      !allocate(Ef(1-ih:np+ih,1-ih:mp+ih,LL));     Ef=0.
      !allocate(Es(1-ih:np+ih,1-ih:mp+ih,LL));     Es=0.
      !allocate(xflux(1-ih:np+ih,1-ih:mp+ih));     xflux=0.
      !allocate(yflux(1-ih:np+ih,1-ih:mp+ih));     yflux=0.
      !allocate(zflux(1-ih:np+ih,1-ih:mp+ih));     zflux=0.
      !allocate(lix(1-ih:n+ih));     lix=0.
      !allocate(liy(1-ih:m+ih));     liy=0.
      !allocate(liz(LL));     liz=0.
      !allocate(rnetsol(1-ih:np+ih,1-ih:mp+ih,l)); rnetsol=0.
      !allocate(rnetgas(1-ih:np+ih,1-ih:mp+ih,l)); rnetgas=0.
      !endif
      allocate(tauw(1-ih:np+ih,1-ih:mp+ih))
      allocate(uzs(1-ih:np+ih,1-ih:mp+ih))
      allocate(vzs(1-ih:np+ih,1-ih:mp+ih))
      allocate(actualfueldepth(nfuel,1-ih:np+ih,1-ih:mp+ih,l))
      allocate(rhof(nfuel,1-ih:np+ih,1-ih:mp+ih,l))
      allocate(sizescale(nfuel,1-ih:np+ih,1-ih:mp+ih,l))   !FP :ss as array
      allocate(rhofinitial(nfuel,1-ih:np+ih,1-ih:mp+ih,l))
      allocate(rhomicro(nfuel,1-ih:np+ih,1-ih:mp+ih,l))
      allocate(fib(1-ih:np+ih,1-ih:mp+ih,l))
      if(iturb.ge.1) then
        allocate(fka(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(fkb(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(rkc(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(sqrtk(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(cd(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(sa(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(saxy(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(saz(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(sb(1-ih:np+ih,1-ih:mp+ih,l))
        !allocate(sc(1-ih:np+ih,1-ih:mp+ih,l))
      endif
      if(inonlocal.eq.0) then
        allocate(u(1-ih:np+ih,1-ih:mp+ih,0:l))
        allocate(v(1-ih:np+ih,1-ih:mp+ih,0:l))
        allocate(w(1-ih:np+ih,1-ih:mp+ih,0:l))
        allocate(theta(1-ih:np+ih,1-ih:mp+ih,0:l))
        allocate(tke_a(1-ih:np+ih,1-ih:mp+ih,0:l))
        allocate(tke_b(1-ih:np+ih,1-ih:mp+ih,0:l))
        allocate(rtke_abc(1-ih:np+ih,1-ih:mp+ih,0:l))
        allocate(ox(1-ih:np+ih,1-ih:mp+ih,0:l))
        if (irhovapor.eq.1) allocate(vap(1-ih:np+ih,1-ih:mp+ih,0:l))
        allocate(K_axy(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(K_az(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(K_b(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(sqrtG_Kxy(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(sqrtG_Kz(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(sqrtG_KG33(1-ih:np+ih,1-ih:mp+ih,l))
        allocate(AB11c(np+1,mp,l),AB12c(np,mp+1,l),AB13c(np,mp,l+1))
        allocate(AB21c(np+1,mp,l),AB22c(np,mp+1,l),AB23c(np,mp,l+1))
        allocate(AB31c(np+1,mp,l),AB32c(np,mp+1,l),AB33c(np,mp,l+1))
        allocate(Ci(np+1,mp,l),Cj(np,mp+1,l),Ck(np,mp,l+1))
      else ! nonlocal arrays 
        allocate(d13(1-ih:np+ih+1,1-ih:mp+ih,l+1))
        allocate(d23(1-ih:np+ih,1-ih:mp+ih+1,l+1))
        allocate(d13a(1-ih:np+ih+1,1-ih:mp+ih,l+1))
        allocate(d12a(1-ih:np+ih+1,1-ih:mp+ih+1,l))
        allocate(d23a(1-ih:np+ih,1-ih:mp+ih+1,l+1))
        allocate(d11a(1-ih:np+ih+1,1-ih:mp+ih,l))
        allocate(d22a(1-ih:np+ih,1-ih:mp+ih+1,l))
        allocate(d33a(1-ih:np+ih,1-ih:mp+ih,l+1))
        allocate(d13b(1-ih:np+ih+1,1-ih:mp+ih,l+1))
        allocate(d12b(1-ih:np+ih+1,1-ih:mp+ih+1,l))
        allocate(d23b(1-ih:np+ih,1-ih:mp+ih+1,l+1))
        allocate(d11b(1-ih:np+ih+1,1-ih:mp+ih,l))
        allocate(d22b(1-ih:np+ih,1-ih:mp+ih+1,l))
        allocate(d33b(1-ih:np+ih,1-ih:mp+ih,l+1))
        allocate(d13t(1-ih:np+ih+1,1-ih:mp+ih,l+1))
        allocate(d12t(1-ih:np+ih+1,1-ih:mp+ih+1,l))
        allocate(d23t(1-ih:np+ih,1-ih:mp+ih+1,l+1))
        allocate(d11t(1-ih:np+ih+1,1-ih:mp+ih,l))
        allocate(d22t(1-ih:np+ih,1-ih:mp+ih+1,l))
        allocate(d33t(1-ih:np+ih,1-ih:mp+ih,l+1))
      endif
      if(irod.eq.1) then
        allocate(ifirestart(1-ih:np+ih,1-ih:mp+ih,l)); ifirestart=0.
        allocate(convht(nfuel,1-ih:np+ih,1-ih:mp+ih,l));     convht=0. 
        allocate(frhosiesrad(nfuel,1-ih:np+ih,1-ih:mp+ih,l)); frhosiesrad=0.
        allocate(temps(nfuel,1-ih:np+ih,1-ih:mp+ih,l));     temps=0.
        allocate(tempg(1-ih:np+ih,1-ih:mp+ih,l));     tempg=0.
        allocate (tambientarray(1-ih:np+ih,1-ih:mp+ih,l));  tambientarray=0.
        allocate (firad(1-ih:np+ih,1-ih:mp+ih,l));     firad=0.
        if (iradeastflux.eq.1) then
          allocate (eastFlux(1-ih:np+ih,1-ih:mp+ih,l))
          eastFlux=0.
        endif
        allocate (foxb(1-ih:np+ih,1-ih:mp+ih,l));      foxb=0.
        allocate (rmoist(nfuel,1-ih:np+ih,1-ih:mp+ih,l));    rmoist=0.
        allocate (rhos(nfuel,1-ih:np+ih,1-ih:mp+ih,l));      rhos=0.
        allocate (rhowater(nfuel,1-ih:np+ih,1-ih:mp+ih,l));  rhowater=0.
        allocate (cpsolid(nfuel,1-ih:np+ih,1-ih:mp+ih,l));   cpsolid=0.
        allocate (sies(nfuel,1-ih:np+ih,1-ih:mp+ih,l));      sies=0.
        allocate (psiwmax(nfuel,1-ih:np+ih,1-ih:mp+ih,l));   psiwmax=0.
        if(irhovapor.eq.1) allocate(frhovaporb(1-ih:np+ih,1-ih:mp+ih,l))
      endif
      allocate (xv(1-ih:np+ih,1-ih:mp+ih,l,5))
      allocate (xvfuel(nfuel,1-ih:np+ih,1-ih:mp+ih,l,5))
      allocate (tkewght(nfuel))
      allocate (xvb(1-ih:np+ih,1-ih:mp+ih,l,nv))
      allocate (xe(1-ih:np+ih,1-ih:mp+ih,l,nv))
      if (uswitch.eq.2.or.vswitch.eq.2) 
     +   allocate (uprofile(1-ih:np+ih,1-ih:mp+ih,l)) !FP
      allocate (pr(1-ih:np+ih,1-ih:mp+ih,l))
      if(ilspgf.ge.1) then
        allocate (sintheta(1-ih:np+ih,1-ih:mp+ih,l))
        allocate (sinthetaf(1-ih:np+ih,1-ih:mp+ih,l))
        allocate (flspgf(np,mp))
      endif
      allocate (pre(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (u1(1-ih:np+ih+1,1-ih:mp+ih,l))
      allocate (u2(1-ih:np+ih,1-ih:mp+ih+1,l))
      allocate (u3(1-ih:np+ih,1-ih:mp+ih,l+1))
! FIXME KOO 
! FP: this fix was incorrect : no need for ghost cells in flux arrays
      allocate (fd1(np+1,mp, l))
      allocate (fd2(np,mp+1, l))
      allocate (fd3(np,mp, l+1))
!      allocate (fd1(1-ih:np+ih+1,1-ih:mp+ih, l))
!      allocate (fd2(1-ih:np+ih,1-ih:mp+ih+1, l))
!      allocate (fd3(1-ih:np+ih,1-ih:mp+ih, l+1))

      if (iord.eq.2) then
        allocate (v1(1-ih:np+ih+1,1-ih:mp+ih, l))
        allocate (v2(1-ih:np+ih,1-ih:mp+ih+1, l))
        allocate (v3(1-ih:np+ih,1-ih:mp+ih, l+1))
        allocate (pmx(np,mp, l))
        allocate (pmn(np,mp, l))
        allocate (cp(1-ih:np+ih,1-ih:mp+ih, l))
        allocate (cn(1-ih:np+ih,1-ih:mp+ih, l))
      endif
      allocate (wt(nts+1))
      allocate (wtf(nts+1))
      if (inonlocal.eq.1) then
         allocate (fg(1-ih:np+ih,1-ih:mp+ih,l));        fg=0.
         allocate (fhc(1-ih:np+ih,1-ih:mp+ih,l));       fhc=0.
         allocate (fhcb(1-ih:np+ih,1-ih:mp+ih,l));       fhcb=0.
         allocate (psig(1-ih:np+ih,1-ih:mp+ih,l));      psig=0.
         ! arrays below were arrays in local code earlier
         allocate (psif(nfuel,1-ih:np+ih,1-ih:mp+ih,l));      psif=0.
         allocate (psiw(nfuel,1-ih:np+ih,1-ih:mp+ih,l));      psiw=0.
         allocate (frhowater(nfuel,1-ih:np+ih,1-ih:mp+ih,l)); frhowater=0.
         allocate (frho(1-ih:np+ih,1-ih:mp+ih,l));      frho=0.
         if(irhovapor.eq.1)
     .    allocate(frhovapor(1-ih:np+ih,1-ih:mp+ih,l)) 
         allocate (fox(1-ih:np+ih,1-ih:mp+ih,l));       fox=0.
         allocate(frhof(nfuel,1-ih:np+ih,1-ih:mp+ih,l))
         allocate (frhosies(nfuel,1-ih:np+ih,1-ih:mp+ih,l));  frhosies=0.
         allocate (cf(nfuel,1-ih:np+ih,1-ih:mp+ih,l));         cf=0.
         allocate(thetasolid(nfuel,1-ih:np+ih,1-ih:mp+ih,l)); thetasolid=0.
         allocate (fi(nfuel,1-ih:np+ih,1-ih:mp+ih,l))
         allocate (ff(nfuel,1-ih:np+ih,1-ih:mp+ih,l));        ff=0.
         allocate (fw(nfuel,1-ih:np+ih,1-ih:mp+ih,l));        fw=0.
         allocate (qflux(1-ih:np+ih,1-ih:mp+ih,l));     qflux=0.
      endif
!234567 JMC these arrays are for atmospheric stability runs
      if(iwindfieldin.eq.1)then
        allocate (xvbdataold(1-ih:np+ih,1-ih:mp+ih,l,nv))
        allocate (xvbdatanew(1-ih:np+ih,1-ih:mp+ih,l,nv))
      endif
      if(iwallclock.EQ.1) then
        allocate(cputime(20,2))
        cputime(1:20,1)=0
        cputime(1:20,2)=0
      endif
      return
      end
