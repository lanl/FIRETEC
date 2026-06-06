!----------------------------------------------------------------
! tinit initializes all arrays and variables in relation to the 
! to and sides BC absorbers
!----------------------------------------------------------------
subroutine tinit
  use gridlist_variables, only : n,m,l,nr,ibcx,ibcy,tow,zab,zabt, &
    iab,ibctopopen
  use gridsetup, only : np,mp
  use metric_variables, only : zb,z
  use constants, only : pi
  use msga_variables, only : npos,mpos,mpi_rank
  use xvall, only : relaxxv
  Implicit None
  
  ! Local Variables
  integer :: i,j,k,ia,ja,it
  real :: zl,t1
  real,external :: zcart
  real,allocatable :: relb(:)

  ! Executable Code
  allocate(relb(nr+1))
  ! Define relaxing boundary (relb)
  if(ibcx.ge.1.or.ibcy.ge.1)then
    if(ibcx.ne.ibcy) then
      if(mpi_rank.eq.0) print*,'ibcx and ibcy are not the same reverts relaxation scheme ',max(ibcx,ibcy)
    endif
    select case(ibcx)
      case(3)
        do i=1,nr+1
          relb(i)=exp(-8.*real(i-1)/real(nr-1))
        enddo
      case(2)
        do i=1,nr+1
          relb(i)=1.5*0.25/(tow*(1.+(real(i)/(0.5*real(nr)-1.))**2.))
        enddo
      case(1)
        do i=1,nr+1
          relb(i)=0.05*exp(-real(i)*(real(i)-1./real(i))**2./(2.*real(nr-1)**2.))
        enddo
    end select
    relb(1:2)=1. ! Mass conservation and Neumann boundary conditions
  endif
 
  do k=1,l
    do j=1,mp
      do i=1,np
        zl=zcart(z(k),i,j)
        ! Upper dampers
        if(zl.ge.zab)then
          if(ibctopopen.eq.0)then
            t1=max(0.,zl-zab)
            select case(iab)
              case(1)
                relaxxv(i,j,k,:)=t1/((zb-zab)+1.e-10)/tow**2.
              case(2)
                relaxxv(i,j,k,:)=0.05*exp(-t1/(-(zb-zab)/log(tow**(-2.))))
              case(3)
                relaxxv(i,j,k,1:3)=0.66*0.25*(1.-cos(pi*(zl-zab)/(zb-zab)))/tow
            end select
          else
            relaxxv(i,j,k,3)=0.66*0.25*(1.-cos(pi*(zl-zab)/(zb-zab)))/tow
          endif
        endif
        if(iab.eq.4.and.zl.ge.zabt) &
          relaxxv(i,j,k,4)=0.66*0.25*(1.-cos(pi*(zl-zab)/(zb-zabt)))/tow
        ! Side dampers
        if(ibcx.ge.1)then
          ia=(npos-1)*np+i
          if(ia.le.nr)then
            do it=1,4
              relaxxv(i,j,k,it)=max(relaxxv(i,j,k,it),relb(ia))
            enddo
          elseif(ia.gt.n-nr)then
            do it=1,4
              relaxxv(i,j,k,it)=max(relaxxv(i,j,k,it),relb(n-ia+1))
            enddo
          endif
        endif
        if(ibcy.ge.1)then
          ja=(mpos-1)*mp+j
          if(ja.le.nr)then
            do it=1,4
              relaxxv(i,j,k,it)=max(relaxxv(i,j,k,it),relb(ja))
            enddo
          elseif(ja.gt.m-nr)then
            do it=1,4
              relaxxv(i,j,k,it)=max(relaxxv(i,j,k,it),relb(m-ja+1))
            enddo
          endif
        endif
      enddo
    enddo
  enddo
  deallocate(relb)

end subroutine tinit
