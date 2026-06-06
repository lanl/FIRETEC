      subroutine dragm(xv,il,iu,jl,ju,lls,nvp)
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
      use fueldrag
      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real xv(il:iu,jl:ju,lls,nvp)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      real,external :: zcart
      integer :: i,j,k,ift,ia,ja
      real :: rroot2,groundcellinterp,rrhomicro
      real :: zla,ztopcell,zbottomcell
      real :: dragxt,dragyt,dragzt,sp,sph,fueltop,zla2,zla3,ztopcell2
      real :: u,v,w,u2,v2,w2,u3,v3,w3,u12,v12,w12,u23,v23,w23,uf,vf,wf
      real :: rinterp12,rinterp1,rinterp2,rinterp23
      real :: u1interp,v1interp,w1interp
      real :: u2interp,v2interp,w2interp,ufueltop,vfueltop,zair
      real :: grasslength,u0grassbend
      real :: bendh,bendx,bendy,bendz,shieldeffect
      real :: avcrosswind,avlengthwind  
      real :: rhofcell

      !(06/16) - variable declaration for drag coefficient   for urban
      !settings
      real :: fwallx,fwally,fwallz    ! friction factor on wall
      real :: Awall,Awallz    ! Area of walls in cell
      real :: Vc       ! volume of cell
      real :: Xx       
      real :: Redx, Redy, Redz      ! reynolds number for drag coefficient
      real :: walllength,floorlength       
      real :: visair   !viscosity of air
      !real, dimension(4) :: cddh, cddhz - declaired in variable subroutine

      rroot2=1./sqrt(2.)
      groundcellinterp=0.5

      do k=1,l
        do j=1,mp
          do i=1,np

            zla=zcart(z(k),i,j)-zs(i,j)
            ztopcell=zcart(zedge(k+1),i,j)-zs(i,j)
            zbottomcell=zcart(zedge(k),i,j)-zs(i,j)

            if((zbottomcell.le.fueldepth).and.
     &        (sum(rhof(:,i,j,k)).gt.min_rhof)) then
              u=xv(i,j,k,1)/xv(i,j,k,nv)
              v=xv(i,j,k,2)/xv(i,j,k,nv)
              w=xv(i,j,k,3)/xv(i,j,k,nv)
              sp=sqrt(u*u+v*v+w*w)
              rho=xv(i,j,k,nv) !density of air
              sph=sqrt(u*u+v*v)
              if(k.eq.1) then
                fueltop = 0.0
                do ift=1,nfuel    !JAS obtain a value of fueltop for this cell
                  if(actualfueldepth(ift,i,j,k).gt.fueltop)then
                    fueltop=actualfueldepth(ift,i,j,k)
                    if(actualfueldepth(ift,i,j,k).gt.ztopcell) then !rrl check  
                      fueltop=ztopcell    !rrl check
                    endif    !rrl check actualfueldepth(ift,i,j,k).gt.ztopcell
                  endif ! actualfueldepth(ift,i,j,k).gt.fueltop
                enddo !ift

                zla2=zcart(z(2),i,j)-zs(i,j)
                zla3=zcart(z(3),i,j)-zs(i,j)
                ztopcell2=zcart(zedge(2+1),i,j)-zs(i,j)

                u2=(xv(i,j,2,1)/xv(i,j,2,nv))
                v2=(xv(i,j,2,2)/xv(i,j,2,nv))
                w2=(xv(i,j,2,3)/xv(i,j,2,nv))
            
                u3=(xv(i,j,3,1)/xv(i,j,3,nv))
                v3=(xv(i,j,3,2)/xv(i,j,3,nv))
                w3=(xv(i,j,3,3)/xv(i,j,3,nv))

                rinterp12=(ztopcell-zla)/(zla2-zla)
                rinterp23=(ztopcell2-zla2)/(zla3-zla2)
                u12=rinterp12*(u2-u)+u
                v12=rinterp12*(v2-v)+v
                w12=rinterp12*(w2-w)+w
                u23=rinterp23*(u3-u2)+u2
                v23=rinterp23*(v3-v2)+v2
                w23=rinterp23*(w3-w2)+w2

                rinterp1=(zla-ztopcell2)/(ztopcell-ztopcell2)
                rinterp2=(zla2-ztopcell2)/(ztopcell-ztopcell2)

                u1interp=rinterp1*(u12-u23)+u23
                v1interp=rinterp1*(v12-v23)+v23
                w1interp=rinterp1*(w12-w23)+w23

                u2interp=rinterp2*(u12-u23)+u23
                v2interp=rinterp2*(v12-v23)+v23
                w2interp=rinterp2*(w12-w23)+w23

                zair=ztopcell-fueltop

                ufueltop=(u1interp*ztopcell-u2interp*.5*(zair)/(zla2-fueltop)*(zair))
     &            /((ztopcell-fueltop/2)-(.5*zair/(zla2-fueltop)*zair))
                vfueltop=(v1interp*ztopcell-v2interp*.5*(zair)/(zla2-fueltop)*(zair))
     &            /((ztopcell-fueltop/2)-(.5*zair/(zla2-fueltop)*zair))

                uf=u12/(abs(u12)+.0000001)*max(.0000001,abs(ufueltop))
                vf=v12/(abs(v12)+.0000001)*max(.0000001,abs(vfueltop))
                wf=w1interp*(fueltop)/ztopcell

                !This is where the fueldrag module routines should be put to use
! KOO
                call calcZoneHts(zbottomcell,ztopcell,rhof(:,i,j,k),actualfueldepth(:,i,j,k))
                call calcZoneRhos(zbottomcell,ztopcell,rhof(:,i,j,k),actualfueldepth(:,i,j,k))
                call calcZoneVels(zbottomcell,ztopcell,uf,vf)
                call calcXvfuels(xvfuel(:,i,j,k,1),xvfuel(:,i,j,k,2),
     &            tkewght(:),rhof(:,i,j,k),actualfueldepth(:,i,j,k))
                do ift=1,nfuel
                  xvfuel(ift,i,j,k,1)=xvfuel(ift,i,j,k,1)*xv(i,j,k,nv)
                  xvfuel(ift,i,j,k,2)=xvfuel(ift,i,j,k,2)*xv(i,j,k,nv)
                  xvfuel(ift,i,j,k,3)=wf*xv(i,j,k,nv)
                  xvfuel(ift,i,j,k,4)=xv(i,j,k,5)*tkewght(ift)
                  xvfuel(ift,i,j,k,5)=xv(i,j,k,6)*tkewght(ift)
                  !xvfuel(ift,i,j,k,4)=xv(i,j,k,5)
                  !xvfuel(ift,i,j,k,5)=xv(i,j,k,6)
                  !sp=sqrt(uf*uf+vf*vf+wf*wf) !JAS original when uf was uf/2
                  sp=sqrt((0.5*uf)**2+(0.5*vf)**2+wf*wf)
                enddo
            
              else    !k.ne.1
                do ift=1,nfuel
                  xvfuel(ift,i,j,k,1)=xv(i,j,k,1)
                  xvfuel(ift,i,j,k,2)=xv(i,j,k,2)
                  xvfuel(ift,i,j,k,3)=xv(i,j,k,3)
                  xvfuel(ift,i,j,k,4)=xv(i,j,k,5)
                  xvfuel(ift,i,j,k,5)=xv(i,j,k,6)
                enddo
              endif  !k.eq.1...else k.ne.1
          
              dragxt=0.   
              dragyt=0.
              dragzt=0.

              do ift=1,nfuel
                rrhomicro=1./rhomicro(ift,i,j,k)
               !write(*,*) "ift=", ift, nfuel  
!(06/16) if statement to specify vegetation or urban

                if(actualfueldepth(ift,i,j,k).ge.zbottomcell.and.k.eq.1) then
                  grasslength=max(actualfueldepth(ift,i,j,k),ztopcell)-zbottomcell
                  u0grassbend=3.                !units of velocity squared
  
                  bendh=cos(atan(((0.5*uf)**2+(0.5*vf)**2)/(u0grassbend+wf**2)))
                  bendx=1-((1-bendh)*((u**2)/(u**2+v**2+.000001)))
                  bendy=1-((1-bendh)*((v**2)/(u**2+v**2+.000001)))
                  bendz=sin(atan((u**2+v**2)/(u0grassbend+w**2)))

                  shieldeffect=1.
                  avcrosswind=2./sizescale(ift,i,j,k)*rhof(ift,i,j,k)*rrhomicro/3.14159   !correct for cylinder if sizescale is radius
                  avlengthwind=rhof(ift,i,j,k)*rrhomicro/grasslength !correct for cylinder if length grasslength)

                  dragxt=dragxt+.5*cd(i,j,k)*sp*(rroot2*bendx*avcrosswind*shieldeffect)
                  dragyt=dragyt+.5*cd(i,j,k)*sp*(rroot2*bendy*avcrosswind*shieldeffect)
                  dragzt=dragzt+.5*cd(i,j,k)*sp*((rroot2+(1.-rroot2)*bendz)*avcrosswind)

                else !actualfueldepth.ge.zbottom.and.k.eq.1 is not true
                  avcrosswind=2./sizescale(ift,i,j,k)*rhof(ift,i,j,k)*rrhomicro/3.14159   !correct for cylinder if sizescale is radius
                  dragxt=dragxt+.5*cd(i,j,k)*sp*(rroot2*(avcrosswind))
                  dragyt=dragyt+.5*cd(i,j,k)*sp*(rroot2*(avcrosswind))
                  dragzt=dragzt+.5*cd(i,j,k)*sp*(rroot2*(avcrosswind))
                endif !actualfueldepth.ge.zbottom.and.k.eq.1 or otherwise...
            !write(*,*) "fwallx=", fwallx, "cddh=", cddh(ift,i,j,k)
              enddo !ift
              xv(i,j,k,1)=xv(i,j,k,1)/(1+dragxt*dt)
              xv(i,j,k,2)=xv(i,j,k,2)/(1+dragyt*dt)
              xv(i,j,k,3)=xv(i,j,k,3)/(1+dragzt*dt)

            elseif(inonlocal.eq.1) then !This is necessary for nonlocal gas combustion above fueldepth
              do ift=1,nfuel
                xvfuel(ift,i,j,k,1)=xv(i,j,k,1)  
                xvfuel(ift,i,j,k,2)=xv(i,j,k,2)
                xvfuel(ift,i,j,k,3)=xv(i,j,k,3)
                xvfuel(ift,i,j,k,4)=xv(i,j,k,5)
                xvfuel(ift,i,j,k,5)=xv(i,j,k,6)
              enddo !ift
            endif !inonlocal.eq.1
          enddo !i
        enddo !j
      enddo !k
  
       !JAS debug block
        !if(mpi_rank.eq.0)then
        !   do ift=1,nfuel
        ! write(6,*) 'xvfuel(',ift,'1,20,40,1,1) = ', xvfuel(ift,20,40,1,1)
        ! write(6,*) 'xvfuel(',ift,'1,20,40,1,2) = ', xvfuel(ift,20,40,1,2)
        ! write(6,*) 'xvfuel(',ift,'1,20,40,1,3) = ', xvfuel(ift,20,40,1,3)
        !   enddo !ift
        !endif
        !if((time.le.0.21).and.(time.ge.0.18))then
       !open(unit=12,file='xvfuel.dat',status='unknown',form='unformatted')
       !call writeio4d(xvfuel(:,:,:,:,1),.false.,12,'xvfuel.dat',1-ih,np+ih,1-ih,mp+ih,l,nfuel)
       !call writeio4d(xvfuel(:,:,:,:,2),.false.,12,'xvfuel.dat',1-ih,np+ih,1-ih,mp+ih,l,nfuel)
       !call writeio4d(xvfuel(:,:,:,:,3),.false.,12,'xvfuel.dat',1-ih,np+ih,1-ih,mp+ih,l,nfuel)
       !call writeio4d(xvfuel(:,:,:,:,4),.false.,12,'xvfuel.dat',1-ih,np+ih,1-ih,mp+ih,l,nfuel)
       !call writeio4d(xvfuel(:,:,:,:,5),.false.,12,'xvfuel.dat',1-ih,np+ih,1-ih,mp+ih,l,nfuel)
           !write(12) xvfuel 
           !close(12)
        !endif
 
      return
      end

