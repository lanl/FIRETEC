!-----------------------------------------------------------------------
! pointSource calculates forcing terms for point source emissions
! as specified by the user
!-----------------------------------------------------------------------
subroutine pointSource(emitForcing,il,iu,jl,ju,lls, &
    nEmit,nAero,nMAero)
  use gridlist_variables, only : prec
  use emission_general_variables, only : spEmit,spAero
  use emission_pointSource_variables, only : nEmitPoints,nAeroPoints, &
    emitPointLocation,aeroPointLocation,emitPointRate,aeroPointRate, &
    emitPointSpecies,aeroPointSpecies
  use gridlist_variables, only : dx,dy
  use metric_variables, only : zedge
  use msga_variables, only : npos,mpos
  use gridsetup, only : np,mp
  use zcart_function
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nEmit,nAero,nMAero
  real(prec),intent(inout) ::  &
    emitForcing(il:iu,jl:ju,lls,nEmit+nAero*nMAero)

  integer :: iemit
  integer :: i,j,k,isp
  integer :: str,stp

  ! Executable Code
  if(nEmitPoints.gt.0)then
    do iemit=1,nEmitPoints
      if(npos.eq.int(emitPointLocation(iemit,1)/dx/np)+1.and. &
         mpos.eq.int(emitPointLocation(iemit,2)/dy/mp)+1)then
        i=mod(int(emitPointLocation(iemit,1)/dx),np)+1
        j=mod(int(emitPointLocation(iemit,2)/dy),mp)+1
        do k=1,lls
          if(zcart(zedge(k),i,j).gt.emitPointLocation(iemit,3)) exit
        enddo
        do isp=1,nEmit
          if(emitPointSpecies(iemit).eq.spEmit(isp)) exit
        enddo
        emitForcing(i,j,k,isp)=emitForcing(i,j,k,isp)+ &
          emitPointRate(iemit)/(dx*dy* &
          (zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)))
      endif
    enddo
  endif

  if(nAeroPoints.gt.0)then
    do iemit=1,nAeroPoints
      if(npos.eq.int(aeroPointLocation(iemit,1)/dx/np)+1.and. &
         mpos.eq.int(aeroPointLocation(iemit,2)/dy/mp)+1)then
        i=mod(int(aeroPointLocation(iemit,1)/dx),np)+1
        j=mod(int(aeroPointLocation(iemit,2)/dy),mp)+1
        do k=1,lls
          if(zcart(zedge(k),i,j).gt.aeroPointLocation(iemit,3)) exit
        enddo
        do isp=1,nAero
          if(aeroPointSpecies(iemit).eq.spAero(isp)) exit
        enddo
        str=nEmit+1+(isp-1)*nMAero
        stp=str+nMAero-1
        emitForcing(i,j,k,str:stp)=emitForcing(i,j,k,str:stp)+ &
          aeroPointRate(iemit,:)/(dx*dy* &
          (zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j)))
      endif
    enddo
  endif

end subroutine pointSource
