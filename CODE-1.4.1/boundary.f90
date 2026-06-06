!----------------------------------------------------------------
! boundary updates and buffers simulation boundary cells with
! environmental or wind-run variables
!----------------------------------------------------------------
subroutine boundary(xv,xe,relaxxv,il,iu,jl,ju,lls,nvp)
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real,intent(in) :: relaxxv(il:iu,jl:ju,lls,4)
  real,intent(in) :: xe(il:iu,jl:ju,lls,nvp)
  real,intent(inout) :: xv(il:iu,jl:ju,lls,nvp)
  integer :: i,j,k,kv

  ! Executable Code
  ! Relaxation application
  do kv=1,nvp
    do k=1,lls
      do j=jl,ju
        do i=il,iu
          xv(i,j,k,kv)=xv(i,j,k,kv)*(1.-relaxxv(i,j,k,min(kv,4))) &
            +xe(i,j,k,kv)*relaxxv(i,j,k,min(kv,4))
        enddo
      enddo
    enddo
  enddo

end subroutine boundary
