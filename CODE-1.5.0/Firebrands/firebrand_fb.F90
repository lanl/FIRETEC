!-----------------------------------------------------------------------
! firebrand is the master driver routine for the firebrand physics 
! module executed during the time iterations
! Original Author : Eunmo Koo (2011)
! Last Updated : Alex Josephson (July 2025)
!-----------------------------------------------------------------------
subroutine firebrand(u,v,w,ilh,iuh,jlh,juh,klh,lls,rhog,tempg, &
    rhof,temps,il,iu,jl,ju,lfuel,nfuel)
  use gridlist_variables, only : prec
  use firebrand_general_variables, only : xvbrand,xvbrand_list, &
    numBrands,nvbrand
  use msga_variables, only : mpi_integer,mpi_sum,mpi_comm_world, &
    mpi_rank,ierror
  Implicit None

  ! Local Variables 
  integer,intent(in) :: ilh,iuh,jlh,juh,klh
  integer,intent(in) :: il,iu,jl,ju
  integer,intent(in) :: lls,lfuel,nfuel
  real(prec),intent(in) :: u(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: v(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: w(ilh:iuh,jlh:juh,klh:lls)
  real(prec),intent(in) :: rhog(il:iu,jl:ju,lls)
  real(prec),intent(in) :: tempg(il:iu,jl:ju,lls)
  real(prec),intent(in) :: rhof(nfuel,il:iu,jl:ju,lfuel)
  real(prec),intent(in) :: temps(nfuel,il:iu,jl:ju,lfuel)

  integer :: totNB
 
  ! Executable Code
  call launchFB(u,v,w,rhog,ilh,iuh,jlh,juh,klh,lls, &
    rhof,temps,il,iu,jl,ju,lfuel,nfuel)

  call evolveFB(u,v,w,rhog,tempg,ilh,iuh,jlh,juh,klh,il,iu,jl,ju,lls)

  call igniteFB()

  ! Only print to standard output when there are firebrands to print
  call mpi_allreduce(numBrands,totNB,1,mpi_integer,mpi_sum, &
    mpi_comm_world,ierror)
  if(totNB.ne.0) then
    if(mpi_rank.eq.0) print*,'Total Firebrands',totNB
    call rmaxmin_1D(xvbrand,xvbrand_list,1,numBrands,nvbrand)
  endif

end subroutine firebrand
