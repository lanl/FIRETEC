      subroutine dragm2(xv,il,iu,jl,ju,lls,nvp)
      use fireteca
      use gridsetup
      use turba
      use workavg
      use gridsetup
      use constants
      use msga
      use xvo
      use metryic
      use nonlocal
      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real xv(il:iu,jl:ju,lls,nvp)

      real,external :: zcart
      integer :: i,j,k,ift !,ia,ja
      real :: ztopcell,zbottomcell,fuelcorrection
      real :: dragxt,dragyt,dragzt,sp
      real :: u,v,w
      real :: av


      do k=1,lfuel
        do j=1,mp
          do i=1,np

            ztopcell=zcart(zedge(k+1),i,j)-zs(i,j)
            zbottomcell=zcart(zedge(k),i,j)-zs(i,j)

            dragxt=0.
            dragyt=0.
            dragzt=0.
            if(zbottomcell.le.fueldepth.and.sum(rhof(:,i,j,k)).gt.min_rhof) then
              u=xv(i,j,k,1)/xv(i,j,k,nv)
              v=xv(i,j,k,2)/xv(i,j,k,nv)
              w=xv(i,j,k,3)/xv(i,j,k,nv)
              sp=sqrt(u*u+v*v+w*w)
         
!FP integrated a drag that depends on k=1 (constant drag for ground and reduced drag on z)
! cd was previously 0.2 and is now 0.15 so that drag remain the same for a kermes oak garrigue of 0.75 m
! , increases for lighter canopy and decrease for more dense
       ! total area
              do ift=1,nfuel
                av=2./sizescale(ift,i,j,k)*rhof(ift,i,j,k)/rhomicro(ift,i,j,k)  !correct for cylinder if ss is radius
                if (k.eq.1) then
                  fuelcorrection=sqrt(actualfueldepth(ift,i,j,k)/ztopcell)
                  xvfuel(ift,i,j,k,1)=xv(i,j,k,1)*fuelcorrection
                  xvfuel(ift,i,j,k,2)=xv(i,j,k,2)*fuelcorrection
                  xvfuel(ift,i,j,k,3)=xv(i,j,k,3)*fuelcorrection
                  xvfuel(ift,i,j,k,4)=xv(i,j,k,5)*fuelcorrection**2.0
                  xvfuel(ift,i,j,k,5)=xv(i,j,k,6)*fuelcorrection**2.0
                  dragxt=dragxt+(.5*cd(i,j,k)*av+0.1/nfuel)*sp
                  dragyt=dragyt+(.5*cd(i,j,k)*av+0.1/nfuel)*sp
                  dragzt=dragzt+.5*cd(i,j,k)*sp*av
                else
                  xvfuel(ift,i,j,k,1)=xv(i,j,k,1)
                  xvfuel(ift,i,j,k,2)=xv(i,j,k,2)
                  xvfuel(ift,i,j,k,3)=xv(i,j,k,3)
                  xvfuel(ift,i,j,k,4)=xv(i,j,k,5)
                  xvfuel(ift,i,j,k,5)=xv(i,j,k,6)
                  dragxt=dragxt+.5*cd(i,j,k)*sp*av
                  dragyt=dragyt+.5*cd(i,j,k)*sp*av
                  dragzt=dragzt+.5*cd(i,j,k)*sp*av
                endif
              enddo
              !formerly implicit
              !xv(i,j,k,1)=xv(i,j,k,1)/(1+dragxt*dt)
              !xv(i,j,k,2)=xv(i,j,k,2)/(1+dragyt*dt)
              !xv(i,j,k,3)=xv(i,j,k,3)/(1+dragzt*dt)
              f1avg(i,j,k)=f1avg(i,j,k)-2.*dragxt*dt*xv(i,j,k,1)
              f2avg(i,j,k)=f2avg(i,j,k)-2.*dragyt*dt*xv(i,j,k,2)
              f3avg(i,j,k)=f3avg(i,j,k)-2.*dragzt*dt*xv(i,j,k,3)
            elseif (inonlocal.eq.1) then
              xvfuel(:,i,j,k,1)=xv(i,j,k,1)
              xvfuel(:,i,j,k,2)=xv(i,j,k,2)
              xvfuel(:,i,j,k,3)=xv(i,j,k,3)
              xvfuel(:,i,j,k,4)=xv(i,j,k,5)
              xvfuel(:,i,j,k,5)=xv(i,j,k,6)
            endif
          enddo
        enddo
      enddo
      return
      end

