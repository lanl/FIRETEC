!-----------------------------------------------------------------------
! rinitrad initializes all arrays and variables for radiation 
! used throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitrad
  use gridlist_variables, only : l,dx,dy,nprocx,nprocy,wp
  use vegetation_gridlist_variables, only : nfuel,rhoMicro
  use gridsetup, only : np,mp
  use metric_variables, only : zedge
  use radiation_gridlist_variables, only : irad
  use radiation_general_variables, only : papvthresh
  use radiation_montecarlo_variables, only : volumeAZ, &
    zCartEdgeAZGlobal
  use fuel_variables, only : actualFuelDepth,kfindex
  use msga_variables, only : mpi_wp,mpi_comm_world,numprocs,ierror
  use fuel_variables, only : min_rhoFuel,sizeScale
  Implicit None

  ! Local Variables
  integer :: i,j,k,it
  integer :: ia,ja,iprocx,jprocy,iproc,nsize
  real(wp),external :: zcart
  real(wp),allocatable :: zCartEdgeAZLocal(:,:,:,:)
  real(wp),allocatable :: zCartEdgeAZTmp(:,:,:,:,:)

  ! Executable Code
  !Allocate Monte Carlo  method radiation arrays
  if(irad.eq.2)then
    papvthresh=2._wp*min_rhoFuel/(maxval(sizeScale)*rhoMicro)
    call defineRadArray

    allocate(zCartEdgeAZLocal(4,np,mp,0:l+nfuel))
    ! AZ  define zcartedge from top of the first cell (k=1) (real mesh) to the top 
    do k=1,l
      do j=1,mp
        do i=1,np
          zCartEdgeAZLocal(1,i,j,k+nfuel)=0.25_wp* &
            (zcart(zedge(k+1),i-1,j-1)+zcart(zedge(k+1),i-1,j)+ &
             zcart(zedge(k+1),i,j-1)+zcart(zedge(k+1),i,j))  !top of k cell, x-, y-
          zCartEdgeAZLocal(2,i,j,k+nfuel)=0.25_wp* &
            (zcart(zedge(k+1),i+1,j-1)+zcart(zedge(k+1),i+1,j)+ &
             zcart(zedge(k+1),i,j-1)+zcart(zedge(k+1),i,j))  !top of k cell, x+, y-
          zCartEdgeAZLocal(3,i,j,k+nfuel)=0.25_wp* &
            (zcart(zedge(k+1),i-1,j+1)+zcart(zedge(k+1),i-1,j)+ &
             zcart(zedge(k+1),i,j+1)+zcart(zedge(k+1),i,j))  !top of k cell, x-, y+
          zCartEdgeAZLocal(4,i,j,k+nfuel)=0.25_wp* &
            (zcart(zedge(k+1),i+1,j+1)+zcart(zedge(k+1),i+1,j)+ &
             zcart(zedge(k+1),i,j+1)+zcart(zedge(k+1),i,j))  !top of k cell, x+, y+
        enddo !i
      enddo !j
    enddo !k
    ! Ground Level AZ
    do j=1,mp
      do i=1,np
        zCartEdgeAZLocal(1,i,j,0)=0.25_wp*(zcart(zedge(1),i-1,j-1)+ &
          zcart(zedge(1),i-1,j)+zcart(zedge(1),i,j-1)+ &
          zcart(zedge(1),i,j))  !top of k cell, x-, y-
        zCartEdgeAZLocal(2,i,j,0)=0.25_wp*(zcart(zedge(1),i+1,j-1)+ &
          zcart(zedge(1),i+1,j)+zcart(zedge(1),i,j-1)+ &
          zcart(zedge(1),i,j))  !top of k cell, x+, y-
        zCartEdgeAZLocal(3,i,j,0)=0.25_wp*(zcart(zedge(1),i-1,j+1)+ &
          zcart(zedge(1),i-1,j)+zcart(zedge(1),i,j+1)+ &
          zcart(zedge(1),i,j))  !top of k cell, x-, y+
        zCartEdgeAZLocal(4,i,j,0)=0.25_wp*(zcart(zedge(1),i+1,j+1)+ & 
          zcart(zedge(1),i+1,j)+zcart(zedge(1),i,j+1)+ &
          zcart(zedge(1),i,j))  !top of k cell, x+, y+
      enddo !i
    enddo !j
          
    ! AZ for individual fuel zones in bottom cell
    do k=1,nfuel
      do j=1,mp
        do i=1,np
          do it=1,4
            zCartEdgeAZLocal(it,i,j,k)=max( &
              zCartEdgeAZLocal(it,i,j,k-1)+ &
              0.0001_wp,min(zCartEdgeAZLocal(it,i,j,nfuel+1)- &
              0.0001_wp*(1+nfuel-k),zCartEdgeAZLocal(it,i,j,0)+ &
              actualFuelDepth(kfindex(i,j,k),i,j)))
          enddo !it 
        enddo !i
      enddo !j
    enddo !k

    ! AZ  volume
    do k=1,l+nfuel
      do j=1,mp
        do i=1,np
          volumeAZ(i,j,k)=dx*dy*0.25_wp*( &
            sum(zCartEdgeAZLocal(:,i,j,k))- &
            sum(zCartEdgeAZLocal(:,i,j,k-1)))
        enddo !i
      enddo !j
    enddo !k

    ! this section to build zCartEdgeAZGlobal (n,m,l+nfuel) can probably be improved   
    allocate(zCartEdgeAZTmp(4,np,mp,0:l+nfuel,numprocs))
    nsize=4*np*mp*(l+nfuel+1)
    call mpi_allgather(zCartEdgeAZLocal,nsize,mpi_wp, &
      zCartEdgeAZTmp,nsize,mpi_wp,mpi_comm_world,ierror)
    do iprocx=1,nprocx
      do jprocy=1,nprocy
        iproc=1+(iprocx-1)+(jprocy-1)*nprocx
        do k=0,l+nfuel
          do j=1,mp
            do i=1,np
              ia=(iprocx-1)*np+i
              ja=(jprocy-1)*mp+j
              zCartEdgeAZGlobal(:,ia,ja,k)=zCartEdgeAZTmp(:,i,j,k,iproc)
            enddo !i
          enddo !j
        enddo !k
      enddo !jprocy
    enddo !iprocx
    deallocate(zCartEdgeAZLocal,zCartEdgeAZTmp)
  endif    

end subroutine rinitrad

