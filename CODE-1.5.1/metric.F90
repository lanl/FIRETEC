subroutine metric()
  use gridlist_variables, only : aa1,l,ih,dx,dy,wp
  use metric_variables, only : zb,z,zs,gmul,h,gi,c13,c23,aa2,aa3
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank,botedge,topedge, &
    leftedge,rightedge
  Implicit None

  ! Local variables
  integer :: i,j,k
  real(wp) :: f
  real(wp) :: gdeform,gdeform1
  real(wp),allocatable :: go(:,:)

  ! Executable code
  allocate (go(1-ih:np+ih,1-ih:mp+ih)); go = 0._wp
  f=0._wp ! 0<=f<=1, f=0=pure cubic fit, f=1=pure quadratic fit
  aa2=f*(1._wp-aa1)/zb
  aa3=(1._wp-aa2*zb-aa1)/zb**2
  if (mpi_rank.eq.0) write(6,*) 'aa1 is ', aa1 
 
  do k=1,l
    gdeform=aa3*z(k)**3+aa2*z(k)**2+aa1*z(k)
    gdeform1=3*aa3*z(k)**2+2*aa2*z(k)+aa1
    gmul(k)=(zb-gdeform)/gdeform1
    do j=1,mp
      do i=1,np
        gi(i,j,k)=zb/(zb-zs(i,j))/gdeform1
        h(i,j,k)=1._wp/gi(i,j,k)
        go(i,j)=zb/(zb-zs(i,j))
      enddo
    enddo
  enddo

  call update(gi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0,0)
  call update(h,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0,0)
  call update(go,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0,0)
  
  do j=1,mp
    do i=1+leftedge,np-rightedge
      c13(i,j)=0.5_wp*(1._wp/go(i+1,j)-1._wp/go(i-1,j))*go(i,j)/dx
    enddo
  enddo
  if(leftedge.eq.1)then
    do j=1,mp
      c13(1,j)=0._wp
    enddo
  endif
  if(rightedge.eq.1)then
    do j=1,mp
      c13(np,j)=0.0_wp
    enddo
  endif
  
  do i=1,np
    do j=1+botedge,mp-topedge
      c23(i,j)=0.5_wp*(1._wp/go(i,j+1)-1._wp/go(i,j-1))*go(i,j)/dy
    enddo
  enddo
  if(botedge.eq.1)then
    do i=1,np
      c23(i,1)=0.0_wp
    enddo
  endif
  if(topedge.eq.1)then
    do i=1,np
      c23(i,mp)=0.0_wp 
    enddo
  endif

  call update(c13,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0,0)
  call update(c23,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1,0,0)
  deallocate (go)
end subroutine metric

!----------------------------------------------------------------
! zcart computes the cartesian vertical coordinate
!----------------------------------------------------------------
function zcart(sigmax,i,j) result(z)
  use gridlist_variables, only : wp,aa1
  use metric_variables, only : aa2,aa3,zb,zs
  Implicit None

  ! Local Variables
  integer,intent(in) :: i,j
  real(wp),intent(in) :: sigmax

  real(wp) :: gdeform
  real(wp) :: z

  ! Executable Code
  gdeform=aa3*sigmax**3+aa2*sigmax**2+aa1*sigmax
  z=gdeform*(zb-zs(i,j))/zb+zs(i,j)

end function zcart
