c234567********************************************************
      subroutine setup
      use gridsetup
      use radiation
      implicit none

      time=0         !time is the counter in seconds

      nn=n           !nn is number of x-cells for radiation domain
      mm=m           !number of y-cells for radiation domain
      nx=1           !nx is the first real x-cell in the radiation domain
      ny=1           !nx is the first real y-cell in the radiation domain
      nz=11          !nz is the first real z-cell in the radiation domain
      LL=L+(nz-1)        !LL is number of z-cells for radiation domain
      nnx=n+nx-1     !nx is the last real x-cell in the radiation domain
      mny=m+ny-1     !ny is the last real x-cell in the radiation domain
      Lnz=L+nz-1     !nz is the last real x-cell in the radiation domain

c parallelization constants
      np=n/nprocx    !number of cells per processor in the x direction
      mp=m/nprocy    !number of cells per processor in the y direction
      nproc=nprocx*nprocy     !number of processors
      nm=n*m         !number of real cells in horizontal plane
      nml=n*m*l      !number of total cells
      ml=m*l         !number of cells in plane perpendicular to x direction
      npmp=np*mp     !number of cells in horizontal plane per processor
      npmpl=np*mp*l           !total number of cells per processor
      mpl=mp*l     !number of cells in plane perp. to x per processor
      dt=dts*real(nts)        !dt is the large time step
      dtp=dt/real(ntp)        !dt is the physics time step for FIRETEC
      dxi=1./dx 
      dyi=1./dy 
      dzi=1./dz 
      dti=1./dt 
      gc1s=dts*dxi   !gc1s is dt/dx for the MOA subcycles      (s/m)
      gc2s=dts*dyi   !dt/dy for the MOA subcycles      (s/m) 
      gc3s=dts*dzi   !gc3s is dt/dz for the MOA subcycles      (s/m)
      gh1s=0.5*gc1s  !gh1s is .5*dt/dx for the MOA subcycles   (s/m)
      gh2s=0.5*gc2s  !gh2s is .5*dt/dx for the MOA subcycles   (s/m)
      gh3s=0.5*gc3s  !gh3s is .5*dt/dx for the MOA subcycles   (s/m)
      gc1=dt*dxi     !gc1 is dt/dx                             (s/m)
      gc2=dt*dyi     !gc2 is dt/dy
      gc3=dt*dzi     !gc3 is dt/dz                             (s/m)
      gh1=0.5*gc1    !gh1 is .5*dt/dx                          (s/m)
      gh2=0.5*gc2    !gh2 is .5*dt/dy                          (s/m)
      gh3=0.5*gc3    !gh3 is .5*dt/dz                          (s/m)
      ibxo=1-ibcx 
      ibyo=1-ibcy
      return
      end
