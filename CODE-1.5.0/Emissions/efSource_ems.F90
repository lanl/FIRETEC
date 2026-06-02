!-----------------------------------------------------------------------
! efSource calculates forcing terms for source emissions with
! fixed emission factors as specified by the user
!-----------------------------------------------------------------------
subroutine efSource(forcing,fuelConsumed,waterEvap,il,iu,jl,ju,lls, &
  nfuel,nEmit,iH2O)
  use gridlist_variables, only : prec
  use emission_general_variables, only : spEmit
  use emission_factor_variables, only : efEmit
  !use xvall, only : iH2O
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nfuel,nEmit,iH2O
  real(prec),intent(inout) :: forcing(il:iu,jl:ju,lls,nEmit)
  real(prec),intent(in)    :: fuelConsumed(nfuel,il:iu,jl:ju,lls)
  real(prec),intent(in)    :: waterEvap(nfuel,il:iu,jl:ju,lls)

  integer :: ift,i,j,k,isp

  ! Executable Code
  do ift=1,nfuel
    do i=il,iu
      do j=jl,ju
        do k=1,lls
          do isp=1,nEmit
            forcing(i,j,k,isp)=forcing(i,j,k,isp) &
              +fuelConsumed(ift,i,j,k)*efEmit(isp) 
          enddo
        enddo
      enddo
    enddo
  enddo

  if(ANY(spEmit.eq."H2O"))then
    do ift=1,nfuel
      do i=il,iu
        do j=jl,ju
          do k=1,lls
            forcing(i,j,k,iH2O)= &
              forcing(i,j,k,iH2O)+waterEvap(ift,i,j,k)
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine efSource
