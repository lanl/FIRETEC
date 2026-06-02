!----------------------------------------------------------------
! setup_MPI will set up the MPI environment
!----------------------------------------------------------------
subroutine setup_MPI
  use gridlist_variables, only : nprocx,nprocy,ibcx,ibcy,n,m
  use gridsetup, only : np,mp
  use msga_variables
  Implicit None

  ! Local Variables

  ! Executable Code 
  ! Set up MPI communication
  call MPI_Init(ierror)
  call MPI_Comm_rank(mpi_comm_world,mpi_rank,ierror)
  call MPI_Comm_size(mpi_comm_world,numprocs,ierror)

  ! Parallelization constants
  np=n/nprocx    !number of cells per processor in the x direction
  mp=m/nprocy    !number of cells per processor in the y direction
  npos = mod((mpi_rank+nprocx),nprocx)+1
  mpos = mpi_rank/nprocx+1
  
  ! set up edge information (0=>not an edge, 1=>edge)
  rightedge=0
  if(mod((mpi_rank+1),nprocx)==0) rightedge=1
  leftedge=0
  if(mod((mpi_rank+1),nprocx)==1.or.nprocx==1) leftedge=1
  botedge=0
  if ((mpi_rank+1).le.nprocx) botedge=1
  topedge=0
  if (((numprocs)-(mpi_rank+1))<nprocx) topedge=1

  ! set up neighbor information
  peleft=mpi_rank-1
  if (peleft<((mpos-1)*nprocx)) peleft=mpi_rank+(nprocx-1)
  peright=mpi_rank+1
  if (peright>(mpos*nprocx-1)) peright=mpi_rank-(nprocx-1)
  peabove=mpi_rank+nprocx
  if (peabove>(numprocs-1)) peabove=npos-1
  pebelow=mpi_rank-nprocx
  if (pebelow<0) pebelow=mpi_rank+((nprocy-1)*nprocx)
  if (npos<nprocx) then
     perightabove=peabove+1
     perightbelow=pebelow+1
  else
     perightabove=peabove-(nprocx-1)
     perightbelow=pebelow-(nprocx-1)
  end if
  if (npos/=1) then
     peleftabove=peabove-1
     peleftbelow=pebelow-1
  else
     peleftabove=peabove+(nprocx-1)
     peleftbelow=pebelow+(nprocx-1)
  end if
  if (mpi_rank==0) &
    print *,' Architecture setup ',mpi_rank,peleft,peright,pebelow, &
    peabove,perightabove,perightbelow,peleftabove,peleftbelow
  if(ibcx>0)then
    rightdedge=rightedge
    leftdedge=leftedge
  endif
  if(ibcy>0)then
    botdedge=botedge
    topdedge=topedge
  endif

end subroutine setup_MPI
