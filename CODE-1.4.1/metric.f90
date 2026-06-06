subroutine metric()
  use gridlist_variables, only : aa1,l,ih,dx,dy
  use metric_variables, only : zb,z,zs,gmul,gi,c13,c23,aa2,aa3
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank,botdedge,topdedge, &
    leftdedge,rightdedge
  use gridlist_variables, only : ibcx,ibcy
  use metric_variables_old
  Implicit None

  ! Local variables
  integer :: i,j,k
  real :: f
  real :: gdeform,gdeform1
  real,allocatable :: go(:,:)

  ! Executable code
  allocate (go(1-ih:np+ih,1-ih:mp+ih))
  f=0. ! 0<=f<=1, f=0=pure cubic fit, f=1=pure quadratic fit
  aa2=f*(1.-aa1)/zb
  aa3=(1.-aa2*zb-aa1)/zb**2.
  if (mpi_rank.eq.0) write(6,*) 'aa1 is ', aa1 
 
  do k=1,l
    gdeform=aa3*z(k)**3+aa2*z(k)**2+aa1*z(k)
    gdeform1=3*aa3*z(k)**2+2*aa2*z(k)+aa1
    gmul(k)=(zb-gdeform)/gdeform1
    do j=1,mp
      do i=1,np
        gi(i,j,k)=zb/(zb-zs(i,j))/gdeform1
        h(i,j,k)=1./gi(i,j,k)
        go(i,j)=zb/(zb-zs(i,j))
      enddo
    enddo
  enddo

  call update(gi,gi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
  call update(h,h,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
  call update(go,go,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0)
  
  do j=1,mp
    do i=1+leftdedge,np-rightdedge
      c13(i,j)=.5*(1./go(i+1,j)-1./go(i-1,j))*go(i,j)/dx
    enddo
  enddo
  if(leftdedge.eq.1)then
    do j=1,mp
      c13(1,j)=0.
    enddo
  endif
  if(rightdedge.eq.1)then
    do j=1,mp
      c13(np,j)=0.
    enddo
  endif
  
  do i=1,np
    do j=1+botdedge,mp-topdedge
      c23(i,j)=.5*(1./go(i,j+1)-1./go(i,j-1))*go(i,j)/dy
    enddo
  enddo
  if(botdedge.eq.1)then
    do i=1,np
      c23(i,1)=0.
    enddo
  endif
  if(topdedge.eq.1)then
    do i=1,np
      c23(i,mp)=0. 
    enddo
  endif

  call update(c13,c13,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0)
  call update(c23,c23,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0)
  deallocate (go)
end subroutine metric

!----------------------------------------------------------------
! zcart computes the cartesian vertical coordinate
!----------------------------------------------------------------
real function zcart(sigmax,i,j)
  use gridlist_variables, only : aa1
  use metric_variables, only : aa2,aa3,zb,zs
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: i,j
  real,intent(in) :: sigmax
  
  real :: gdeform

  ! Executable Code
  gdeform=aa3*sigmax**3+aa2*sigmax**2+aa1*sigmax
  zcart=gdeform*(zb-zs(i,j))/zb+zs(i,j)

end function zcart
