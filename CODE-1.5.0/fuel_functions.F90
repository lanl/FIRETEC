!-----------------------------------------------------------------------
! computePsi(temps) computes the pdf with several methods
!-----------------------------------------------------------------------

module computePsi_function
  use gridlist_variables, only : prec

  contains
  real(prec) function computePsi(temps)
    use gridlist_variables, only : ifuel
    use fuel_variables, only : tcrit,tfstep
    Implicit None
  !
    ! Local Variables
    real(prec),intent(in) :: temps
  !
    real(prec) :: rlcrit,rlramp,rljoint1,scalefactor,rpsi
    real(prec) :: rltop,rljoint2,rpsi2,rl
    real(prec) :: xerf,terf,er,sumf
    real(prec) :: tramp=500. ! the tail of the ramp straight line (if extrapolated to psi-0)
    real(prec) :: psijoint=0.30 ! height of the ramp joint
    real(prec) :: c1=0.5
    real(prec) :: c2=0.0079
    real(prec) :: c3=1.
    real(prec) :: p=0.47047
    real(prec) :: a1=0.3480242
    real(prec) :: a2=-0.0958798
    real(prec) :: a3=0.7478556
  !
    ! Executable Code  
    if(ifuel.eq.1)then  !Use RRL's method of calculating computePsi!
      rlcrit=(tcrit-tfstep)
      rlramp=(tramp-tfstep)
      rljoint1=psijoint*2.*(rlcrit-rlramp)+rlramp
      scalefactor=sqrt(rlramp**2-(rljoint1-rlramp)**2)/psijoint
      rpsi=2.0/scalefactor*(rlcrit-rlramp)*rljoint1+scalefactor*psijoint
      rltop=2.*rlcrit
      rljoint2=rltop-rljoint1
      rpsi2=scalefactor-rpsi

      ! temps is used here
      rl=temps-tfstep
      if (rl.le.0.0 ) then
        computePsi=0.0
      elseif (rl.le.rljoint1) then
        computePsi=1./scalefactor*(rpsi-sqrt(rpsi**2-rl**2))
      elseif (rl.le.rljoint2) then
        computePsi=.5/(rlcrit-rlramp)*(rl-rlramp)
      elseif (rl.le.rltop) then
        computePsi=1./scalefactor*(rpsi2+sqrt(rpsi**2-(rl-rltop)**2))
      elseif (rl.gt.rltop) then
        computePsi=1.
      endif
  !
    elseif(ifuel.eq.2)then  !Use Michael Clark's method of calculating computePsi!
      xerf = c2*(temps-tcrit)
      terf = 1./(1.+p*abs(xerf))
      sumf = a1*terf+a2*terf*terf+a3*terf*terf*terf
      er   = 1-sumf*exp(-xerf*xerf)
      computePsi = c1*(c3+er)
      if (xerf.lt.0.0) computePsi=-computePsi+c3

    elseif(ifuel.eq.3)then  !Use JAS's method of calculating computePsi!
      if(temps.le.tfstep)then 
        computePsi=0.0
      elseif(temps.le.(2*(tcrit-tfstep)+tfstep))then
        er = erf(c2*(temps-tcrit))                         
        computePsi = c1*(c3+er)
      else 
        computePsi = 1.0
      endif  
    endif    !end if(ifuel==1,2,or 3)        

  end function computePsi
end module computePsi_function
