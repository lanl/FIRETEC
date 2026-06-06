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
       Implicit None

       !JAS 3/7/06 added explicit declarations to comply with implicit none
       integer,intent(in) :: il,iu,jl,ju,lls,nvp

       real xv(il:iu,jl:ju,lls,nvp)

       !JAS 3/7/06 added explicit declarations to comply with implicit none
       real,external :: zcart
       integer :: i,j,k,ia,ja
       real :: rroot2,groundcellinterp,rrhomicro
       real :: zla,ztopcell,zbottomcell
       real :: dragxt,dragyt,dragzt,sp,sph,grassht,zla2,zla3,ztopcell2
       real :: u,v,w,u2,v2,w2,u3,v3,w3,u12,v12,w12,u23,v23,w23,uf,vf,wf
       real :: rinterp12,rinterp1,rinterp2,rinterp23
       real :: u1interp,v1interp,w1interp
       real :: u2interp,v2interp,w2interp,ufueltop,vfueltop,zair
       real :: speedrho2,speedfuelrho2,grasslength,u0grassbend
       real :: bendh,bendx,bendy,bendz,shieldeffect
       real :: avcrosswind,avlengthwind  

c           actualgrassheight=.7
c       actualgrassheight=.35
       rroot2=1./sqrt(2.)
       groundcellinterp=0.5

       do k=1,lfuel
       do j=1,mp
       do i=1,np
       rrhomicro=1./rhomicro(i,j,k)

       zla=zcart(z(k),i,j)-zs(i,j)
       ztopcell=zcart(zedge(k+1),i,j)-zs(i,j)
       zbottomcell=zcart(zedge(k),i,j)-zs(i,j)

       dragxt=0.
       dragyt=0.
       dragzt=0.
       if(zbottomcell.le.fueldepth
     &         .and.rhof(i,j,k).gt.min_rhof) then
          u=xv(i,j,k,1)/xv(i,j,k,nv)
          v=xv(i,j,k,2)/xv(i,j,k,nv)
          w=xv(i,j,k,3)/xv(i,j,k,nv)
          sp=sqrt(u*u+v*v+w*w)
          rho=xv(i,j,k,nv)
          sph=sqrt(u*u+v*v)
          if (k.eq.1) then
c             grassht=actualgrassheight
             grassht=actualfueldepth(i,j,k)
             zla2=zcart(z(2),i,j)-zs(i,j)
             zla3=zcart(z(3),i,j)-zs(i,j)
             ztopcell2=zcart(zedge(2+1),i,j)-zs(i,j)

             u2=(xv(i,j,2,1)/xv(i,j,2,nv))
             v2=(xv(i,j,2,2)/xv(i,j,2,nv))
             w2=(xv(i,j,2,3)/xv(i,j,2,nv))
C**JLWams from Rod
             u3=(xv(i,j,3,1)/xv(i,j,3,nv))
             v3=(xv(i,j,3,2)/xv(i,j,3,nv))
             w3=(xv(i,j,3,3)/xv(i,j,3,nv))

cJLWams      u3=(xv(i,j,2,1)/xv(i,j,2,nv))
cJLWams      v3=(xv(i,j,2,2)/xv(i,j,2,nv))
cJLWams      w3=(xv(i,j,2,3)/xv(i,j,2,nv))

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


             zair=ztopcell-grassht

             ufueltop=(u1interp*ztopcell
     &                 -u2interp*.5*(zair)
     &                     /(zla2-grassht)*(zair))
     &               /((ztopcell-grassht/2)-
     &                (.5*zair/(zla2-grassht)*zair))
             vfueltop=(v1interp*ztopcell
     &                 -v2interp*.5*(zair)
     &                     /(zla2-grassht)*(zair))
     &               /((ztopcell-grassht/2)-
     &                (.5*zair/(zla2-grassht)*zair))


c             ufueltop=(u*ztopcell
c     &                 -u2*.5*(zair)
c     &                     /(zla2-grassht)*(zair))
c     &               /((ztopcell-grassht/2)-
c     &                (.5*zair/(zla2-grassht)*zair))
cc             vfueltop=(v*ztopcell
c     &                 -v2*.5*(zair)
c     &                     /(zla2-grassht)*(zair))
c     &               /((ztopcell-grassht/2)-
c     &                (.5*zair/(zla2-grassht)*zair))
c

             uf=u12/(abs(u12)+.0000001)
     &                *max(.0000001,abs(ufueltop)/2.)
             vf=v12/(abs(v12)+.0000001)
     &                *max(.0000001,abs(vfueltop)/2.)
             wf=w1interp*(grassht)/ztopcell

             if (uf.eq..0000001.or.temps(i,j,k).gt.320.) then

              ia=(npos-1)*np+i
              ja=(mpos-1)*mp+j
c             if (ia.lt.30.or.ja.lt.40) 
c    &         write (*,*) 'lower boundary',ia,ja,k,u
c    &       ,u2,uf,temps(i,j,k),xv(i,j,k,4)/xv(i,j,k,nv)
c    &        ,v
c    &       ,v2,vf,convht(i,j,k),frhosiesrad(i,j,k),
c    &        thetasolid(i,j,k)
        endif

             xvfuel(i,j,k,1)=uf*xv(i,j,k,nv)
             xvfuel(i,j,k,2)=vf*xv(i,j,k,nv)
             xvfuel(i,j,k,3)=wf*xv(i,j,k,nv)
             sp=sqrt(uf*uf+vf*vf+wf*wf)

             speedrho2=(xv(i,j,k,1)**2+
     &                  xv(i,j,k,2)**2+
     &                  xv(i,j,k,3)**2)
             speedfuelrho2=(xvfuel(i,j,k,1)**2+
     &                      xvfuel(i,j,k,2)**2+
     &                      xvfuel(i,j,k,3)**2)

c             xvfuel(i,j,k,4)=xv(i,j,k,5)*speedfuelrho2/speedrho2
c             xvfuel(i,j,k,4)=xv(i,j,k,5)*speedfuelrho2/speedrho2
             xvfuel(i,j,k,4)=xv(i,j,k,5)
             xvfuel(i,j,k,5)=xv(i,j,k,6)

          else
             xvfuel(i,j,k,1)=xv(i,j,k,1)
             xvfuel(i,j,k,2)=xv(i,j,k,2)
             xvfuel(i,j,k,3)=xv(i,j,k,3)
             xvfuel(i,j,k,4)=xv(i,j,k,5)
             xvfuel(i,j,k,5)=xv(i,j,k,6)
             sp=sqrt(u*u+v*v+w*w)
          endif
c          if (actualgrassheight.ge.zbottomcell.and.k.eq.1) then
          if (actualfueldepth(i,j,k).ge.zbottomcell.and.k.eq.1) then
         
             grasslength=max(actualfueldepth(i,j,k),ztopcell)
     &                -zbottomcell
c             write (*,*) actualfueldepth(i,j,k), grasslength

             u0grassbend=3.                !units of velocity squared
cJLWams      u0grassbend=300000.           !units of velocity squared

             bendh=cos(atan((uf**2+vf**2)/(u0grassbend+wf**2)))

             bendx=1-((1-bendh)*((u**2)/(u**2+v**2+.000001)))
             bendy=1-((1-bendh)*((v**2)/(u**2+v**2+.000001)))

             bendz=sin(atan((u**2+v**2)/(u0grassbend+w**2)))

c            bending=.5+.5*exp(-.125*sp**1.5)
c             shieldeffect=.25+.75*exp(-.125*sph**1.5)
c            bending=1.
             shieldeffect=1.
c            av=2.*(rhof(i,j,k))*0.004/ss

             !avcrosswind=2./ss*rhof(i,j,k)*rrhomicro/3.14159   !correct for cylinder if ss is radius
             avcrosswind=2./sizescale(i,j,k)*rhof(i,j,k)*rrhomicro/3.14159   !correct for cylinder if ss is radius
		!FP
             avlengthwind=rhof(i,j,k)*rrhomicro/grasslength  !correct for cylinder if length grasslength)

c         drag=.25*.375/sb(i,j,k)*cd(i,j,k)*sp*av
             dragxt=.5*cd(i,j,k)*sp
     &           *(rroot2*bendx*avcrosswind*shieldeffect)
c     &           +(rroot2+(1.-rroot2)*sqrt(1-bendx**2))*avlengthwind) 
             dragyt=.5*cd(i,j,k)*sp
     &           *(rroot2*bendy*avcrosswind*shieldeffect)
c     &           +(rroot2+(1.-rroot2)*sqrt(1-bendy**2))*avlengthwind) 
             dragzt=.5*cd(i,j,k)*sp
     &           *((rroot2+(1.-rroot2)*bendz)*avcrosswind)
c     &          sqrt(rroot2+(1.-rroot2)*(1.-bendz**2))*avlengthwind) 
          else
			!FP
             !avcrosswind=2./ss*rhof(i,j,k)*rrhomicro/3.14159   !correct for cylinder if ss is radius
             avcrosswind=2./sizescale(i,j,k)*rhof(i,j,k)*rrhomicro/3.14159   !correct for cylinder if ss is radius
             dragxt=.5*cd(i,j,k)*sp
     &           *(rroot2*(avcrosswind))
             dragyt=.5*cd(i,j,k)*sp
     &           *(rroot2*(avcrosswind))
             dragzt=.5*cd(i,j,k)*sp
     &           *(rroot2*(avcrosswind))
          endif

c         dragx=-rho*u*dragxt
c         dragy=-rho*v*dragyt
c         dragz=-rho*w*dragzt
c        f1avg(i,j,k)=f1avg(i,j,k)+2.*dragx*dt
c        f2avg(i,j,k)=f2avg(i,j,k)+2.*dragy*dt
c        f3avg(i,j,k)=f3avg(i,j,k)+2.*dragz*dt
         xv(i,j,k,1)=xv(i,j,k,1)/(1+dragxt*dt)
         xv(i,j,k,2)=xv(i,j,k,2)/(1+dragyt*dt)
         xv(i,j,k,3)=xv(i,j,k,3)/(1+dragzt*dt)
c        if(k.eq.1.and.mpi_rank.eq.0) 
c    .   print*,i,2.*dragx*dt
       elseif (inonlocal.eq.1) then
         xvfuel(i,j,k,1)=xv(i,j,k,1)
         xvfuel(i,j,k,2)=xv(i,j,k,2)
         xvfuel(i,j,k,3)=xv(i,j,k,3)
         xvfuel(i,j,k,4)=xv(i,j,k,5)
         xvfuel(i,j,k,5)=xv(i,j,k,6)
       endif
       enddo
       enddo
       enddo
   
       return
       end

