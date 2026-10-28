!-----------------------------------------------------------------------
! efSource calculates forcing terms for source emissions with
! fixed emission factors as specified by the user
!-----------------------------------------------------------------------
subroutine efSource(forcing,fuelForcing,ilh,iuh,jlh,juh,ku, &
    il,iu,jl,ju,lls,nfuel,nEmit)
  use gridlist_variables, only : wp
  use emission_factor_variables, only : efEmit
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,ku,il,iu,jl,ju,lls,nfuel,nEmit
  real(wp),intent(inout) :: forcing(ilh:iuh,jlh:juh,ku,nEmit)
  real(wp),intent(in)    :: fuelForcing(nfuel,il:iu,jl:ju,lls)

  integer :: ift,i,j,k,isp

  ! Executable Code
  do ift=1,nfuel
    do i=il,iu
      do j=jl,ju
        do k=1,lls
          do isp=1,nEmit
            forcing(i,j,k,isp)=forcing(i,j,k,isp) &
              -fuelForcing(ift,i,j,k)*efEmit(isp) 
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine efSource
