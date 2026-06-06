! Routine computes the generation rate of soot moments in a cell
subroutine soot(ift,i,j,k,ff)
  
  end subroutine soot

  !*****************************************************************
  ! Quenched zone mechanisms
  !*****************************************************************
!  if(xv(i,j,k,iM0).gt.0.and.xv(i,j,k,iM1).gt.0)then
!    d_p = (6.0*xv(i,j,k,iM1)/(xv(i,j,k,iM0)*
! +    rhost*pi))**(r13) ! Particle diameter (m)
!    Kn = kB*(xv(i,j,k,4)/xv(i,j,k,nv))/
! +    (SQRT(2.0)*pi*d_p**3.0*pr(i,j,k)) ! Knudsen Number
!    beta_c = 8.0*kB*(xv(i,j,k,4)/xv(i,j,k,nv))/
! +    (3.0*Visc)*(1.0+1.257*Kn) ! Frequency of collision continuum regime
!    Cg = beta_c*(xv(i,j,k,iM0))**2.0 ! Coagulation rate in the air zone
!    xv(i,j,k,iM0) = xv(i,j,k,iM0)
! +    *EXP(-Cg/(xv(i,j,k,iM0))*dtp)
!  endif
!  
!  end subroutine soot
      
!***********************************************************************      
real function double_lin_interp(table,x,y,it0,it1)
  ! Computes a linearly interpolated value from a 2D table
  ! table is 2D data array
  ! x is the coordinate (0-1) for the first value
  ! y is the coordinate (0-1) for the second value
  ! it0 is the index value on the table for the first value
  ! it1 is the index value on the table for the second value

  implicit none
  real,intent(in) :: table(11,11)
  real   :: x,y
  integer:: it0,it1
  real   :: v0,v1,v2,v3
  real   :: x0,x1

  v0 = table(it0,it1)
  v1 = table(it0,it1+1)
  v2 = table(it0+1,it1)
  v3 = table(it0+1,it1+1)

  x0 = v0+x*(v2-v0)
  x1 = v1+x*(v3-v1)
  double_lin_interp= x0+y*(x1-x0)

  return
end function
