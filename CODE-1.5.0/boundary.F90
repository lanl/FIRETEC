!-----------------------------------------------------------------------
! boundary updates and buffers simulation boundary cells with
! environmental or wind-run variables
!-----------------------------------------------------------------------
subroutine boundary(force,xv,xe,relaxxv,dt,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : prec
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(prec),intent(in) :: dt
  real(prec),intent(in) :: relaxxv(il:iu,jl:ju,lls,4)
  real(prec),intent(in) :: xe(il:iu,jl:ju,lls,nvp)
  real(prec),intent(in) :: xv(il:iu,jl:ju,lls,nvp)
  real(prec),intent(inout) :: force(il:iu,jl:ju,lls,nvp)
  integer :: i,j,k,kv

  ! Executable Code
  ! Relaxation application
  do kv=1,nvp
    do k=1,lls
      do j=jl,ju
        do i=il,iu
          force(i,j,k,kv)=force(i,j,k,kv)+relaxxv(i,j,k,min(kv,4))* &
            (xe(i,j,k,kv)-xv(i,j,k,kv))/dt
        enddo
      enddo
    enddo
  enddo

end subroutine boundary
