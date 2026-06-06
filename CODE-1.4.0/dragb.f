      !subroutine dragb(xv,il,iu,jl,ju,lls,nvp)
      subroutine dragb()
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
      real :: drag,sp,av,alphas,rpi
      rpi=1./3.14159
      do k=1,lfuel
        do j=1,mp
          do i=1,np
            sp=sqrt(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
            if (idrag.eq.1) drag = -sqrtk(i,j,k)*xvb(i,j,k,6)/sb(i,j,k)
            if (idrag.eq.2) drag = -0.93*sqrtk(i,j,k)*xvb(i,j,k,6)/sb(i,j,k) ! FP added 0.93 as in Deardorf
            do ift=1,nfuel
              av=2./sizescale(ift,i,j,k)*rhof(ift,i,j,k)/rhomicro(ift,i,j,k)  !correct for cylinder if ss is radius
              if (idrag.eq.1) then
                drag=drag+rpi*av*sp*(.25*xvb(i,j,k,5)-xvb(i,j,k,6))
              elseif (idrag.eq.2) then
                drag=drag+cd(i,j,k)*av*sp*(xvb(i,j,k,5)-xvb(i,j,k,6))
              endif
            enddo
            fkb(i,j,k)=fkb(i,j,k)+2.*drag*dt
          enddo
        enddo
      enddo
      do k=lfuel+1,l
        do j=1,mp
          do i=1,np
            if (idrag.eq.1) then
              drag=-sqrtk(i,j,k)*xvb(i,j,k,6)/sb(i,j,k)
            elseif (idrag.eq.2) then
              drag=-0.93*sqrtk(i,j,k)*xvb(i,j,k,6)/sb(i,j,k)
            endif
            fkb(i,j,k)=fkb(i,j,k)+2.*drag*dt
          enddo
        enddo
      enddo
      return
      end
