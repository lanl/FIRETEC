!-----------------------------------------------------------------------
! boundary updates and buffers simulation boundary cells with
! environmental or wind-run variables
!-----------------------------------------------------------------------
subroutine boundary(force,xv,xe,relaxxv,dt,ilh,iuh,jlh,juh, &
    il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : wp
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,il,iu,jl,ju,lls,nvp
  real(wp),intent(in) :: dt
  real(wp),intent(in) :: relaxxv(il:iu,jl:ju,lls,4)
  real(wp),intent(in) :: xe(il:iu,jl:ju,lls,nvp)
  real(wp),intent(in) :: xv(ilh:iuh,jlh:juh,lls,nvp)
  real(wp),intent(inout) :: force(ilh:iuh,jlh:juh,lls,nvp)
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
