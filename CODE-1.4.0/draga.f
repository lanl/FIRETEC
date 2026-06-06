      !subroutine draga(xv,il,iu,jl,ju,lls,nvp)
      subroutine draga()
      use turba
      use updatedFields
      use gridsetup
      use msga
      use constants
      use xvo
      use fireteca, only:nfuel
      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      !integer,intent(in) :: il,iu,jl,ju,lls,nvp

      !real xv(il:iu,jl:ju,lls,nvp)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,ift
      real :: sp,rpi
      real :: av,drag

      rpi=1./3.14159
        
      do k=1,lfuel
        do j=1,mp
          do i=1,np
            sp=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
            if (idrag.eq.1) drag=-sqrtk(i,j,k)*xvb(i,j,k,5)/sa(i,j,k)
            if (idrag.eq.2) drag=-0.93*sqrtk(i,j,k)*xvb(i,j,k,5)/sa(i,j,k)
            do ift=1,nfuel
              if (idrag.eq.1) av=2./sizescale(ift,i,j,k)*rpi*rhof(ift,i,j,k)/rhomicro(ift,i,j,k)
              if (idrag.eq.2) av=2./sizescale(ift,i,j,k)*cd(i,j,k)*rhof(ift,i,j,k)/rhomicro(ift,i,j,k)
              drag = drag - av*sp*xvb(i,j,k,5)
            enddo
          fka(i,j,k)=fka(i,j,k)+2.*drag*dt
          enddo
        enddo
      enddo
      do k=lfuel+1,l
        do j=1,mp
          do i=1,np
            if (idrag.eq.1) then 
              drag=-sqrtk(i,j,k)*xvb(i,j,k,5)/sa(i,j,k)
            else if (idrag.eq.2) then
              ! FP added 0.93 as in Deardorf in the decay term
              drag=-0.93*sqrtk(i,j,k)*xvb(i,j,k,5)/sa(i,j,k)
            endif
          fka(i,j,k)=fka(i,j,k)+2.*drag*dt
          enddo
        enddo
      enddo
      return
      end
