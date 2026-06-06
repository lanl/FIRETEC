C-grid R_ij calculation
      subroutine rij(xv,sax,say,saz,sb,il,iu,jl,ju,lls,nvp)
      use metryic
      use gridsetup
      use turbb
      use bc
      use msga
      use io

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp  
      real xv(il:iu,jl:ju,lls,nvp)
      real sax(il:iu,jl:ju,lls)
      real say(il:iu,jl:ju,lls)
      real saz(il:iu,jl:ju,lls)
      real sb(il:iu,jl:ju,lls)
      real,allocatable::  u(:,:,:),
     .                    v(:,:,:), 
     .                    w(:,:,:),
     .                  rho(:,:,:),
     .                 cdrg(:,:,:),
     .                  uka(:,:,:),
     .                  ukb(:,:,:) 
      real,allocatable:: rka(:,:,:),
     .                   rkb(:,:,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1
      integer :: ip,im,jp,jm,islipf
      real :: hdxi,hdyi,hdzi,hx,hy,bx,by,dxil,dyil
      real :: vtang,btt,dlt,gmm,gia,gmut,gii,g13,g23,cc,dva,dvb
      real :: uxa,uya,uza,vxa,vya,vza,wxa,wya,wza
      real :: uxt,uzt,vyt,vzt,wzt,vtax,vtaz,vtay,vtb,rhoa
c sax, say and saz (not isotropic...)
      
c vt=0.09*s(length scale)*sqrt(rk which is k_a or k_b))
cdrg is drag due to bottom boundary
components of r_ij
      real,allocatable:: tmp(:,:,:)
      if(iturb.ge.1) then
      allocate (u(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (rho(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (cdrg(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (uka(1-ih:np+ih,1-ih:mp+ih,l+1))
      allocate (ukb(1-ih:np+ih,1-ih:mp+ih,l+1))
      allocate (rka(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (rkb(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (tmp(1-ih:np+ih,1-ih:mp+ih,l))
      endif
 
compute some local constants
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi
      islipf=1

      do j=1,mp
      do i=1,np
      uzs(i,j)=0.
      vzs(i,j)=0.
      enddo
      enddo


      d13=0.; d23=0.; d12a=0.; d12b=0.  ! wss change 12/18/00
      do k=1,l
      do j=1,mp
      do i=1,np
      u(i,j,k)=xv(i,j,k,1)/xv(i,j,k,nv)
      v(i,j,k)=xv(i,j,k,2)/xv(i,j,k,nv)
      w(i,j,k)=xv(i,j,k,3)/xv(i,j,k,nv)
      rka(i,j,k)=xv(i,j,k,5)/xv(i,j,k,nv)
      rkb(i,j,k)=xv(i,j,k,6)/xv(i,j,k,nv)
      rho(i,j,k)=xv(i,j,k,nv)
      enddo
      enddo
      enddo
      call updated(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(w,w,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rka,rka,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rkb,rkb,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rho,rho,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      
c surface fluxes --------> for rod's stuff cdrg=0.
      do j=1,mp
      do i=1,np
      cdrg(i,j,1)=0.
      hx=-c13(i,j)/gi(i,j,1)*zb
      hy=-c23(i,j)/gi(i,j,1)*zb
      vtang=sqrt( (u(i,j,1)+hx*w(i,j,1))**2/(1.+hx**2)
     .           +(v(i,j,1)+hy*w(i,j,1))**2/(1.+hy**2) )
      tauw(i,j)=cdrg(i,j,1)*vtang*rho(i,j,1)
      enddo
      enddo
compute velocity z-derivatives.  note du/dz is in d13, dv/dz is in d23
      do k=2,l
        do j=1,mp
          do i=1,np
            d13(i,j,k)=dzi*(u(i,j,k)-u(i,j,k-1))
            d23(i,j,k)=dzi*(v(i,j,k)-v(i,j,k-1))
          end do
        end do
      end do
Constraints for "free-slip" boundaries:
c        in 2D: [ d13*(1.-hx**2)=(d11-d33)*hx ]
c in 3D:  (1-hx**2)*d13-    hx*hy*d23=(d11-d33)*hx+d12*hy
c in 3D:     -hx*hy*d13+(1-hy**2)*d23=(d22-d33)*hy+d12*hx
c where di3 denotes stress' elements; these conditions lead
c to auxiliary conditions on dui/dz denoted by di3 temporarily
c note that top boundary is assumed to be flat
      if(islipf.eq.1) then
      jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
      julim = mp + (ibcy-j3)*topdedge                       !add d rrl
      do j=jllim,julim
         if (topdedge.eq.1 .and. j.eq.mp) then                          !add d rrl
            if (ibcy.eq.1) then
               jp1 = mp + 2
            else
               jp1 = mp
            end if
         else
            jp1 = j + j3
         end if
         if (botdedge.eq.1 .and. j.eq.1) then                         !add d rrl
            if (ibcy.eq.1) then
               jm1 = -1
            else
               jm1 = 1
            end if
         else
            jm1 = j - j3
         end if

c      jp1=(j+j3-j/m*(m-1))*ibcy+(1-ibcy)*min0(j+j3,m)
c      jm1=(j-j3+(m-j)/(m-j3)*(m-j3))*ibcy+(1-ibcy)*max0(j-j3,1)
      dyil=hdyi
      if ((ibcy.eq.0) .and. (j.eq.1 .and. botdedge.eq.1)) dyil=dyi                         !add d rrl
      if ((ibcy.eq.0) .and. (j.eq.mp .and. topdedge.eq.1)) dyil=dyi                         !add d rrl
c      if((ibcy.eq.0) .and. (j.eq.1.or.j.eq.m) ) dyil=dyi
        do i=1,np
           if (rightdedge.eq.1 .and. i.eq.np) then                         !add d rrl
              if (ibcx.eq.1) then
                 ip1 = np + 2
              else
                 ip1 = np
              end if
           else
              ip1 = i + 1
           end if
           if (leftdedge.eq.1 .and. i.eq.1) then                         !add d rrl
              if (ibcx.eq.1) then
                 im1 = -1
              else
                 im1 = 1
              end if
           else
              im1 = i - 1
           end if
c        ip1=(i+1-i/n*(n-1))*ibcx+(1-ibcx)*min0(i+1,n)
c        im1=(i-1+(n-i)/(n-1)*(n-1))*ibcx+(1-ibcx)*max0(i-1,1)
        dxil=hdxi
        if ((ibcx.eq.0) .and. (i.eq.1 .and. leftdedge.eq.1)) dxil=dxi                         !add d rrl
        if ((ibcx.eq.0) .and. (i.eq.np .and. rightdedge.eq.1)) dxil=dxi                         !add d rrl
c        if((ibcx.eq.0) .and. (i.eq.1.or.i.eq.n) ) dxil=dxi
          hx=-c13(i,j)/gi(i,j,1)*zb
          hy=-c23(i,j)/gi(i,j,1)*zb
          btt=1.+hx**2+hy**2
          gmm=1-hx**2
          dlt=1-hy**2
          cc=1./gi(i,j,1)/btt
      d13(i,j,1)=cc*(-gmm*(w(ip1,j,1)-w(im1,j,1))*dxil
     .          -hx/cc*(3.*w(i,j,2)-2.*w(i,j,1)-w(i,j,3))*dzi
     . +( 2.*hx*(u(ip1,j,1)-u(im1,j,1))*dxil
     .   +   hy*( (u(i,jp1,1)-u(i,jm1,1))*dyil
     .           +(v(ip1,j,1)-v(im1,j,1))*dxil )
     .   +hx*hy*(w(i,jp1,1)-w(i,jm1,1))*dyil      ))
     .+tauw(i,j)*(u(i,j,1)+hx*w(i,j,1))/sqrt(btt)
      d23(i,j,1)=cc*(-dlt*(w(i,jp1,1)-w(i,jm1,1))*dyil
     .          -hy/cc*(3.*w(i,j,2)-2.*w(i,j,1)-w(i,j,3))*dzi
     . +( 2.*hy*(v(i,jp1,1)-v(i,jm1,1))*dyil
     .   +   hx*( (u(i,jp1,1)-u(i,jm1,1))*dyil
     .           +(v(ip1,j,1)-v(im1,j,1))*dxil )
     .   +hx*hy*(w(ip1,j,1)-w(im1,j,1))*dxil      ))
     .+tauw(i,j)*(v(i,j,1)+hy*w(i,j,1))/sqrt(btt)
c uzs and vzs are in the center
      uzs(i,j)=d13(i,j,1)
      vzs(i,j)=d23(i,j,1)
      d13(i,j,1)=2.*d13(i,j,1)-d13(i,j,2)
      d23(i,j,1)=2.*d23(i,j,1)-d23(i,j,2)
      d13(i,j,l+1)=-d13(i,j,l)
      d23(i,j,l+1)=-d23(i,j,l)
c     d13(i,j,1)=cc*(-gmm*(w(ip1,j,1)-w(im1,j,1))*dxil
c    .          -hx/cc*(3.*w(i,j,2)-2.*w(i,j,1)-w(i,j,3))*dzi
c    . +( 2.*hx*(u(ip1,j,1)-u(im1,j,1))*dxil
c    .   +   hy*( (u(i,jp1,1)-u(i,jm1,1))*dyil
c    .           +(v(ip1,j,1)-v(im1,j,1))*dxil )
c    .   +hx*hy*(w(i,jp1,1)-w(i,jm1,1))*dyil      ))
c    .+tauw(i,j)*(u(i,j,1)+hx*w(i,j,1))/sqrt(btt)
c     if(mpi_rank.eq.0) then
c    .   +   hy*( (u(i,jp1,1)-u(i,jm1,1))*dyil
c    .           +(v(ip1,j,1)-v(im1,j,1))*dxil )
c    .   +hx*hy*(w(i,jp1,1)-w(i,jm1,1))*dyil)
c     endif
      enddo
      enddo
close free-slip constraint
      endif


compute divergence of the wind vector at (i,j,k+-1/2)
      if(j3.eq.1) then
      do k=2,l
      do j=1,mp
      jp=j+1
      jm=j-1
      by=1.
      if(botdedge.eq.1.and.j.eq.1) jm=-ibcy+(1-ibcy)                         !add d rrl
      if(topdedge.eq.1.and.j.eq.mp) jp= ibcy*(mp+2)+mp*(1-ibcy)                         !add d rrl
      if((botdedge.eq.1.and.j.eq.1).and.ibcy.eq.0) by=2.                         !add d rrl
      if((topdedge.eq.1.and.j.eq.mp).and.ibcy.eq.0) by=2.                         !add d rrl
      do i=1,np
      ip=i+1
      im=i-1
      bx=1.
      if(leftdedge.eq.1.and.i.eq.1) im=-ibcx+(1-ibcx)                         !add d rrl
      if(rightdedge.eq.1.and.i.eq.np) ip= ibcx*(np+2)+np*(1-ibcx)                         !add d rrl
      if((leftdedge.eq.1.and.i.eq.1).and.ibcx.eq.0) bx=2.                         !add d rrl
      if((rightdedge.eq.1.and.i.eq.np).and.ibcx.eq.0) bx=2.                         !add d rrl
      gia=0.5*(gi(i,j,k)+gi(i,j,k-1))

      gmut=0.5*(gmul(k)+gmul(k-1))
      uxt=0.25*(u(ip,j,k)+u(ip,j,k-1)-u(im,j,k)-u(im,j,k-1))*dxi*bx
      vyt=0.25*(v(i,jp,k)+v(i,jp,k-1)-v(i,jm,k)-v(i,jm,k-1))*dyi*by
      uzt=c13(i,j)*gmut*d13(i,j,k)
      vzt=c23(i,j)*gmut*d23(i,j,k)
      wzt=gia*(w(i,j,k)-w(i,j,k-1))*dzi
       vtax=0.09*0.5*(sax(i,j,k)+sax(i,j,k-1))
     . *sqrt(0.5*(rka(i,j,k)+rka(i,j,k-1)))
       vtay=0.09*0.5*(say(i,j,k)+say(i,j,k-1))
     . *sqrt(0.5*(rka(i,j,k)+rka(i,j,k-1)))
       vtaz=0.09*0.5*(saz(i,j,k)+saz(i,j,k-1))
     . *sqrt(0.5*(rka(i,j,k)+rka(i,j,k-1)))
       vtb=0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))
     . *sqrt(0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
       rhoa=0.5*(rho(i,j,k)+rho(i,j,k-1))/gia
      uka(i,j,k)=1./3.*rhoa*(vtax*uxt+vtay*vyt+vtaz*(uzt+vzt+wzt))
     .          +1./3.*rhoa*0.5*(rka(i,j,k)+rka(i,j,k-1)) 
      ukb(i,j,k)=1./3.*rhoa*vtb*(uxt+vyt+uzt+vzt+wzt)
     .          +1./3.*rhoa*0.5*(rkb(i,j,k)+rkb(i,j,k-1)) 
      enddo
      enddo
      enddo
      do j=1,mp
      do i=1,np
      uka(i,j,1)=2.*uka(i,j,2)-uka(i,j,3)
      uka(i,j,l+1)=2.*uka(i,j,l)-uka(i,j,l-1)
      ukb(i,j,1)=2.*ukb(i,j,2)-ukb(i,j,3)
      ukb(i,j,l+1)=2.*ukb(i,j,l)-ukb(i,j,l-1)
      enddo
      enddo
      else
      do k=2,l
      do i=1,np
      im=i-1
      ip=i+1
      if(leftdedge.eq.1.and.i.eq.1) im=-ibcx+(1-ibcx)                         !add d rrl
      if(rightdedge.eq.1.and.i.eq.np) ip=ibcx*(np+2)+np*(1-ibcx)                         !add d rrl
      bx=1.
      by=1.
      if((leftdedge.eq.1.and.i.eq.1).and.ibcx.eq.0) bx=2.                         !add d rrl
      if((rightdedge.eq.1.and.i.eq.np).and.ibcx.eq.0) bx=2.                         !add d rrl
      gia=0.5*(gi(i,1,k)+gi(i,1,k-1))
      gmut=0.5*(gmul(k)+gmul(k-1))
      uxt=0.25*(u(ip,1,k)+u(ip,1,k-1)-u(im,1,k)-u(im,1,k-1))*dxi*bx
      uzt=c13(i,1)*gmut*d13(i,1,k)
      wzt=gia*(w(i,1,k)-w(i,1,k-1))*dzi
      vtax=0.09*0.5*(sax(i,1,k)+sax(i,1,k-1))
     . *sqrt(0.5*(rka(i,1,k)+rka(i,1,k-1)))
      vtaz=0.09*0.5*(saz(i,1,k)+saz(i,1,k-1))
     . *sqrt(0.5*(rka(i,1,k)+rka(i,1,k-1)))
      vtb=0.09*0.5*(sb(i,1,k)+sb(i,1,k-1))
     . *sqrt(0.5*(rkb(i,1,k)+rkb(i,1,k-1)))
       rhoa=0.5*(rho(i,1,k)+rho(i,1,k-1))/gia
      uka(i,1,k)=1./3.*rhoa*(vtax*uxt+vtaz*(uzt+wzt))
     .         +1./3.*rhoa*0.5*(rka(i,1,k)+rka(i,1,k-1)) 
      ukb(i,1,k)=1./3.*rhoa*vtb*(uxt+uzt+wzt)
     .         +1./3.*rhoa*0.5*(rkb(i,1,k)+rkb(i,1,k-1))
      enddo
      enddo
      do i=1,np
      uka(i,1,1)=2.*uka(i,1,2)-uka(i,1,3)
      uka(i,1,l+1)=2.*uka(i,1,l)-uka(i,1,l-1)
      ukb(i,1,1)=2.*ukb(i,1,2)-ukb(i,1,3)
      ukb(i,1,l+1)=2.*ukb(i,1,l)-ukb(i,1,l-1)
      enddo
      endif

      call updated(uzs,uzs,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(vzs,vzs,np,mp,1,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(uka,uka,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(ukb,ukb,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1)

compute d11 at (i +- 1/2, j, k)

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d13,d13,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(d13,d13,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      end if


      do k=1,L
        do j=1,mp
           illim = 1+leftdedge                         !add d rrl
           iulim = np
           do i=illim,iulim
            g13=gmul(k)*0.5*(c13(i-1,j)+c13(i,j))
            uza=0.25*(d13(i-1,j,k)+d13(i,j,k)+
     1                d13(i-1,j,k+1)+d13(i,j,k+1))
      dva=0.25*(uka(i-1,j,k)+uka(i,j,k)+uka(i-1,j,k+1)+uka(i,j,k+1))
      dvb=0.25*(ukb(i-1,j,k)+ukb(i,j,k)+ukb(i-1,j,k+1)+ukb(i,j,k+1))
            uxa=dxi*(u(i,j,k)-u(i-1,j,k))
            vtax=0.09*0.5*(sax(i,j,k)+sax(i-1,j,k))
     1         *sqrt(0.5*(rka(i,j,k)+rka(i-1,j,k)))
            vtaz=0.09*0.5*(saz(i,j,k)+saz(i-1,j,k))
     1         *sqrt(0.5*(rka(i,j,k)+rka(i-1,j,k)))
            vtb=0.09*0.5*(sb(i,j,k)+sb(i-1,j,k))
     1         *sqrt(0.5*(rkb(i,j,k)+rkb(i-1,j,k)))
       rhoa=0.5*(rho(i,j,k)+rho(i-1,j,k))
     1/(0.5*(gi(i,j,k)+gi(i-1,j,k)))
!  TODO : all these dij should be checked because apparently different from arps!
       
            d11a(i,j,k) = 2.*(-rhoa*(vtax*uxa + vtaz*g13*uza) + dva)
            d11b(i,j,k) = 2.*(-vtb*rhoa*(uxa + g13*uza) + dvb)
          end do
        end do
      end do

      if (rightdedge.eq.0)                          !add d rrl
     .call updated(d11a,d11a,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      if (rightdedge.eq.1)                          !add d rrl
     .call updated(d11a,d11a,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      if (rightdedge.eq.0)                          !add d rrl
     .call updated(d11b,d11b,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      if (rightdedge.eq.1)                          !add d rrl
     .call updated(d11b,d11b,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)

      if(leftdedge.eq.1) then                         !add d rrl
      do k=1,l
      do j=1,mp
      d11a(1,j,k)=(1-ibcx)*d11a(2,j,k) + ibcx*d11a(-1,j,k)
      d11b(1,j,k)=(1-ibcx)*d11b(2,j,k) + ibcx*d11b(-1,j,k)
      enddo
      enddo
      endif

      if(rightdedge.eq.1) then                         !add d rrl
      do k=1,l
      do j=1,mp
      d11a(np+1,j,k)=(1-ibcx)*d11a(np,j,k) + ibcx*d11a(np+3,j,k)
      d11b(np+1,j,k)=(1-ibcx)*d11b(np,j,k) + ibcx*d11b(np+3,j,k)
      enddo
      enddo
      endif

      do k=1,l
      do j=1,mp
      do i=1,np+1
        d11t(i,j,k)=d11a(i,j,k) 
        if(iturb.eq.2) 
     .     d11t(i,j,k)=d11t(i,j,k)+d11b(i,j,k)+0.2*d11b(i,j,k)
      enddo
      enddo
      enddo

      if (rightdedge.eq.0)                          !add d rrl
     .call updated(d11a,d11a,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      if (rightdedge.eq.1)                          !add d rrl
     .call updated(d11a,d11a,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      if (rightdedge.eq.0)                          !add d rrl
     .call updated(d11b,d11b,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      if (rightdedge.eq.1)                          !add d rrl
     .call updated(d11b,d11b,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      if (rightdedge.eq.0)                          !add d rrl
     .call updated(d11t,d11t,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      if (rightdedge.eq.1)                          !add d rrl
     .call updated(d11t,d11t,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)

      
compute d33 at (i, j, k +- 1/2)
      do j=1,mp
        do i=1,np
          do k=2,L
          gia=0.5*(gi(i,j,k)+gi(i,j,k-1))
            wza = dzi*(w(i,j,k)-w(i,j,k-1))
              vtaz=0.09*0.5*(saz(i,j,k)+saz(i,j,k-1))
     .           *sqrt(0.5*(rka(i,j,k)+rka(i,j,k-1)))
              vtb=0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))
     .           *sqrt(0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
            rhoa=0.5*(rho(i,j,k)+rho(i,j,k-1))/gia
            d33a(i,j,k) = 2.*(-vtaz*gia*rhoa*wza + uka(i,j,k) )
            d33b(i,j,k) = 2.*(-vtb*gia*rhoa*wza + ukb(i,j,k) )
c         if(mpi_rank.eq.0) then
c50       format(e30.10,2x,e30.10)
c         endif
          end do
       end do
      end do
create boundary conditions at k=1 and k=L+1
      do j=1,mp
        do i=1,np
          wza = dzi*(3.*w(i,j,2)-2.*w(i,j,1)-w(i,j,3))
            vtaz=0.09*saz(i,j,1)*sqrt(rka(i,j,1))
            vtb=0.09*sb(i,j,1)*sqrt(rkb(i,j,1))
          rhoa=rho(i,j,1)/gi(i,j,1)
          d33a(i,j,1) = 2.*( -vtaz*gi(i,j,1)*rhoa*wza + uka(i,j,1) )
          d33b(i,j,1) = 2.*( -vtb*gi(i,j,1)*rhoa*wza + ukb(i,j,1) )
          wza = dzi*(-3.*w(i,j,l-1)+2.*w(i,j,l)+w(i,j,l-2))
            vtaz=0.09*saz(i,j,l)*sqrt(rka(i,j,l))
            vtb=0.09*sb(i,j,l)*sqrt(rkb(i,j,l))
          rhoa=rho(i,j,l)/gi(i,j,l)
          d33a(i,j,l+1) = 2.*( -vtaz*gi(i,j,l)*rhoa*wza + uka(i,j,l+1))
          d33b(i,j,l+1) = 2.*( -vtb*gi(i,j,l)*rhoa*wza + ukb(i,j,l+1))
        end do
      end do

      do k=1,l+1
      do j=1,mp
      do i=1,np
        d33t(i,j,k)=d33a(i,j,k)
        if(iturb.eq.2) 
     .     d33t(i,j,k)=d33t(i,j,k)+d33b(i,j,k)+0.2*d33b(i,j,k)
      enddo
      enddo
      enddo
      call updated(d33a,d33a,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(d33b,d33b,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(d33t,d33t,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1)

      if (topdedge.eq.0) then                         !add d rrl
         call updated(d23,d23,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      else
         call updated(d23,d23,np,mp+1,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      end if

compute d22 at (i, j +- 1/2, k)
         jllim = 1+botdedge                         !add d rrl
         julim = mp
        if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=jllim,julim
              g23=0.5*gmul(k)*(c23(i,j-1)+c23(i,j))
              vza=0.25*(d23(i,j-1,k  )+d23(i,j,k  )+
     .                  d23(i,j-1,k+1)+d23(i,j,k+1) )
              dva=0.25*(uka(i,j-1,k  )+uka(i,j,k  )+
     .                  uka(i,j-1,k+1)+uka(i,j,k+1) )
              dvb=0.25*(ukb(i,j-1,k  )+ukb(i,j,k  )+
     .                  ukb(i,j-1,k+1)+ukb(i,j,k+1) )
              vya=dyi*(v(i,j,k)-v(i,j-1,k))
              vtay=0.09*0.5*(say(i,j,k)+say(i,j-1,k))
     .          *sqrt(0.5*(rka(i,j,k)+rka(i,j-1,k)))
              vtaz=0.09*0.5*(saz(i,j,k)+saz(i,j-1,k))
     .          *sqrt(0.5*(rka(i,j,k)+rka(i,j-1,k)))
              vtb=0.09*0.5*(sb(i,j,k)+sb(i,j-1,k))
     .          *sqrt(0.5*(rkb(i,j,k)+rkb(i,j-1,k)))
              rhoa=0.5*(rho(i,j,k)+rho(i,j-1,k))/
     .             (0.5*(gi(i,j,k)+gi(i,j-1,k)))
              d22a(i,j,k) = 2.*( -rhoa*(vya*vtay + g23*vza*vtaz) + dva )
              d22b(i,j,k) = 2.*( -vtb*rhoa*(vya + g23*vza) + dvb )
            end do
          end do
        end do

        if (topdedge.eq.0) then                         !add d rrl
         call updated(d22a,d22a,np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,1)
         call updated(d22b,d22b,np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,1)
        else
         call updated(d22a,d22a,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1)
         call updated(d22b,d22b,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1)
        end if


        if(botdedge.eq.1) then                         !add d rrl
        do k=1,l
        do i=1,np
        d22a(i,1,k)=(1-ibcy)*d22a(i,2,k) + ibcy*d22a(i,-1,k)
        d22b(i,1,k)=(1-ibcy)*d22b(i,2,k) + ibcy*d22b(i,-1,k)
        enddo
        enddo
        endif

        if(topdedge.eq.1) then                         !add d rrl
        do k=1,l
        do i=1,np
        d22a(i,mp+1,k)=(1-ibcy)*d22a(i,mp,k) + ibcy*d22a(i,mp+3,k)
        d22b(i,mp+1,k)=(1-ibcy)*d22b(i,mp,k) + ibcy*d22b(i,mp+3,k)
        enddo
        enddo
        endif

        else

        do k=1,L
          do i=1,np
              d22a(i,1,k) = 2.*uka(i,1,k) 
              d22a(i,2,k)=d22a(i,1,k)
              d22b(i,1,k) = 2.*ukb(i,1,k) 
              d22b(i,2,k)=d22b(i,1,k)
            end do
        end do

        endif

        do k=1,l
        do j=1,mp+1
        do i=1,np
        d22t(i,j,k)=d22a(i,j,k)
        if(iturb.eq.2) 
     .     d22t(i,j,k)=d22t(i,j,k)+d22b(i,j,k)+0.2*d22b(i,j,k)
        enddo
        enddo
        enddo

        if (topdedge.eq.0) then                         !add d rrl
         call updated(d22a,d22a,np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,1)
        else
         call updated(d22a,d22a,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1)
        end if


        if (topdedge.eq.0) then                         !add d rrl
         call updated(d22b,d22b,np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,1)
        else
         call updated(d22b,d22b,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1)
        end if


        if (topdedge.eq.0) then                         !add d rrl
         call updated(d22t,d22t,np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,1)
        else
         call updated(d22t,d22t,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1)
        end if

compute d12 at (i +- 1/2, j+- 1/2, k)
        if(j3.eq.1) then
        jllim = 1+botdedge                         !add d rrl
        julim = mp
        illim = 1+leftdedge                         !add d rrl
        iulim = np
        do k=1,L
          do j=jllim,julim
            do i=illim,iulim
      uya=0.
              g13=gmul(k)*0.25*(  c13(i-1,j-j3) + c13(i,j-j3)
     .                          + c13(i-1,j   ) + c13(i,j   )  )
              g23=gmul(k)*0.25*(  c23(i-1,j-j3) + c23(i,j-j3)
     .                          + c23(i-1,j   ) + c23(i,j   )  )
              uya=hdyi*((u(i-1,j,k)-u(i-1,j-1,k))+(u(i,j,k)-u(i,j-1,k)))
              vxa=hdxi*((v(i,j-1,k)-v(i-1,j-1,k))+(v(i,j,k)-v(i-1,j,k)))
              uza=0.125*(  d13(i-1,j-j3,k  )+d13(i,j-j3,k  )
     .                   + d13(i-1,j   ,k  )+d13(i,j   ,k  )
     .                   + d13(i-1,j-j3,k+1)+d13(i,j-j3,k+1)
     .                   + d13(i-1,j   ,k+1)+d13(i,j   ,k+1) )
              vza=0.125*(  d23(i-1,j-j3,k  )+d23(i,j-j3,k  )
     .                   + d23(i-1,j   ,k  )+d23(i,j   ,k  )
     .                   + d23(i-1,j-j3,k+1)+d23(i,j-j3,k+1)
     .                   + d23(i-1,j   ,k+1)+d23(i,j   ,k+1) )
       vtax=0.09*0.25*(sax(i,j,k)+sax(i,j-1,k)+sax(i-1,j,k)+sax(i-1,j-1,k))
     .   *sqrt(0.25*(rka(i,j,k)+rka(i,j-1,k)
     .              +rka(i-1,j,k)+rka(i-1,j-1,k)))
       vtay=0.09*0.25*(say(i,j,k)+say(i,j-1,k)+say(i-1,j,k)+say(i-1,j-1,k))
     .   *sqrt(0.25*(rka(i,j,k)+rka(i,j-1,k)
     .              +rka(i-1,j,k)+rka(i-1,j-1,k)))
       vtaz=0.09*0.25*(saz(i,j,k)+saz(i,j-1,k)+saz(i-1,j,k)+saz(i-1,j-1,k))
     .   *sqrt(0.25*(rka(i,j,k)+rka(i,j-1,k)
     .              +rka(i-1,j,k)+rka(i-1,j-1,k)))
       vtb=0.09*0.25*(sb(i,j,k)+sb(i,j-1,k)+sb(i-1,j,k)+sb(i-1,j-1,k))
     .   *sqrt(0.25*(rkb(i,j,k)+rkb(i,j-1,k)
     .              +rkb(i-1,j,k)+rkb(i-1,j-1,k)))
       rhoa=0.25*(rho(i,j,k)+rho(i-1,j,k)+rho(i-1,j-1,k)+rho(i,j-1,k))/
     .     (0.25*(gi(i,j,k)+gi(i-1,j,k)+gi(i-1,j-1,k)+gi(i,j-1,k)))
              d12a(i,j,k) = -rhoa*(uya*vtay + g23*uza*vtaz + vxa*vtax + g13*vza*vtaz)
              d12b(i,j,k) = -vtb*rhoa*(uya + g23*uza + vxa + g13*vza)
            end do
          end do
        end do

      if (rightdedge.eq.0 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12a,d12a,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12a,d12a,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                         !add d rrl
         call updated(d12a,d12a,np,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(d12a,d12a,np+1,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

      if (rightdedge.eq.0 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12b,d12b,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12b,d12b,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                         !add d rrl
         call updated(d12b,d12b,np,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(d12b,d12b,np+1,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

        if(leftdedge.eq.1) then                         !add d rrl
        do k=1,l
        do j=1,mp
        d12a(1,j,k)=(1-ibcx)*d12a(2,j,k) + ibcx*d12a(-1,j,k)
        d12b(1,j,k)=(1-ibcx)*d12b(2,j,k) + ibcx*d12b(-1,j,k)
        enddo
        enddo
        endif
        if(rightdedge.eq.1) then                         !add d rrl
        do k=1,l
        do j=1,mp
        d12a(np+1,j,k)=(1-ibcx)*d12a(np,j,k) + ibcx*d12a(np+3,j,k)
        d12b(np+1,j,k)=(1-ibcx)*d12b(np,j,k) + ibcx*d12b(np+3,j,k)
        enddo
        enddo
        endif
        if(botdedge.eq.1) then                         !add d rrl
        do k=1,l
        do i=1,np
        d12a(i,1,k)=(1-ibcy)*d12a(i,2,k) + ibcy*d12a(i,-1,k)
        d12b(i,1,k)=(1-ibcy)*d12b(i,2,k) + ibcy*d12b(i,-1,k)
        enddo
        enddo
        endif
        if(topdedge.eq.1) then                         !add d rrl
        do k=1,l
        do i=1,np
        d12a(i,mp+1,k)=(1-ibcy)*d12a(i,mp,k) + ibcy*d12a(i,mp+3,k)
        d12b(i,mp+1,k)=(1-ibcy)*d12b(i,mp,k) + ibcy*d12b(i,mp+3,k)
        enddo
        enddo
        endif
        else
        do k=1,l
        do i=1+leftdedge,np                         !add d rrl
        g13=gmul(k)*0.5*(c13(i-1,1)+c13(i,1))
        vxa=dxi*(v(i,1,k)-v(i-1,1,k))
        vza=0.25*(d23(i,1,k)+d23(i-1,1,k)+d23(i,1,k+1)+d23(i-1,1,k+1))
        vtax=0.09*0.5*(sax(i,1,k)+sax(i-1,1,k))*
     .      sqrt(0.5*(rka(i,1,k)+rka(i-1,1,k)))
        vtaz=0.09*0.5*(saz(i,1,k)+saz(i-1,1,k))*
     .      sqrt(0.5*(rka(i,1,k)+rka(i-1,1,k)))
        vtb=0.09*0.5*(sb(i,1,k)+sb(i-1,1,k))*
     .      sqrt(0.5*(rkb(i,1,k)+rkb(i-1,1,k)))
        rhoa=0.5*(rho(i,1,k)+rho(i-1,1,k))
        d12a(i,1,k)=-rhoa*(vxa*vtax+g13*vza*vtaz)
        d12a(i,2,k)=d12a(i,1,k)
        d12b(i,1,k)=-vtb*rhoa*(vxa+g13*vza)
        d12b(i,2,k)=d12b(i,1,k)
        enddo
        enddo

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d12a,d12a,np,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(d12a,d12a,np+1,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d12b,d12b,np,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(d12b,d12b,np+1,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

        if(leftdedge.eq.1) then                         !add d rrl
        do k=1,l
        d12a(1,1,k)=(1-ibcx)*d12a(2,1,k) + ibcx*d12a(-1,1,k)
        d12a(1,2,k)=d12a(1,1,k)
        d12b(1,1,k)=(1-ibcx)*d12b(2,1,k) + ibcx*d12b(-1,1,k)
        d12b(1,2,k)=d12b(1,1,k)
        enddo
        endif

        if(rightdedge.eq.1) then                         !add d rrl
        do k=1,l
        d12a(np+1,1,k)=(1-ibcx)*d12a(np,1,k) + ibcx*d12a(np+3,1,k)
        d12a(np+1,2,k)=d12a(np+1,1,k)
        d12b(np+1,1,k)=(1-ibcx)*d12b(np,1,k) + ibcx*d12b(np+3,1,k)
        d12b(np+1,2,k)=d12b(np+1,1,k)
        enddo
        endif

        endif

        do k=1,l
        do j=1,mp+1
        do i=1,np+1
        d12t(i,j,k)=d12a(i,j,k)
        if(iturb.eq.2) 
     .     d12t(i,j,k)=d12t(i,j,k)+d12b(i,j,k)+0.2*d12b(i,j,k)
        enddo
        enddo
        enddo

      if (rightdedge.eq.0 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12a,d12a,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12a,d12a,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                         !add d rrl
         call updated(d12a,d12a,np,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(d12a,d12a,np+1,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

      if (rightdedge.eq.0 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12b,d12b,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12b,d12b,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                         !add d rrl
         call updated(d12b,d12b,np,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(d12b,d12b,np+1,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

      if (rightdedge.eq.0 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12t,d12t,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                         !add d rrl
         call updated(d12t,d12t,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                         !add d rrl
         call updated(d12t,d12t,np,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(d12t,d12t,np+1,mp+1,l,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

        if(j3.eq.1) then
        jllim = 1+botdedge                         !add d rrl
        julim = mp
        do k=2,L
          do j=julim,jllim,-1
            do i=1,np
              gii=0.25*(gi(i,j-j3,k)+gi(i,j,k)
     1                 +gi(i,j,k-1)+gi(i,j-j3,k-1))
              g23=0.25*(gmul(k-1)+gmul(k))*(c23(i,j-j3)+c23(i,j))
              wya=hdyi*((w(i,j,k-1)-w(i,j-1,k-1))+(w(i,j,k)-w(i,j-1,k)))
              vza=0.5*(d23(i,j-j3,k)+d23(i,j,k))
              wza=hdzi*((w(i,j-j3,k)-w(i,j-j3,k-1))+
     1                  (w(i,j   ,k)-w(i,j   ,k-1)))
        vtay=0.09*0.25*(say(i,j,k)+say(i,j-1,k)+say(i,j,k-1)+say(i,j-1,k-1))*
     .     sqrt(0.25*(rka(i,j,k)+rka(i,j-1,k)
     .               +rka(i,j,k-1)+rka(i,j-1,k-1)))
        vtaz=0.09*0.25*(saz(i,j,k)+saz(i,j-1,k)+saz(i,j,k-1)+saz(i,j-1,k-1))*
     .     sqrt(0.25*(rka(i,j,k)+rka(i,j-1,k)
     .               +rka(i,j,k-1)+rka(i,j-1,k-1)))
        vtb=0.09*0.25*(sb(i,j,k)+sb(i,j-1,k)+sb(i,j,k-1)+sb(i,j-1,k-1))*
     .     sqrt(0.25*(rkb(i,j,k)+rkb(i,j-1,k)
     .               +rkb(i,j,k-1)+rkb(i,j-1,k-1)))
        rhoa=0.25*(rho(i,j,k)+rho(i,j-1,k)+rho(i,j,k-1)+rho(i,j-1,k-1))/
     .       gii
              d23a(i,j,k) = -rhoa*( gii*vza*vtaz + wya*vtay + g23*wza*vtaz )
              d23b(i,j,k) = -rhoa*vtb*( gii*vza + wya + g23*wza )
            end do
          end do
        end do

      if (topdedge.eq.0) then                         !add d rrl
         call updated(d23,d23,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      else
         call updated(d23,d23,np,mp+1,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      end if


        do j=julim,jllim,-1
          do i=1,np
            gii=0.5*(gi(i,j-j3,1)+gi(i,j,1))
            g23=.25*(3.*gmul(1)-gmul(2))*(c23(i,j-j3)+c23(i,j))
            wya=hdyi*(3.*(w(i,j,1)-w(i,j-j3,1))-(w(i,j,2)-w(i,j-j3,2)))
            vza=0.5*(d23(i,j-j3,1)+d23(i,j,1))
            wza=(3.*w(i,j   ,2)-2.*w(i,j   ,1)-w(i,j   ,3))*hdzi
     1         +(3.*w(i,j-j3,2)-2.*w(i,j-j3,1)-w(i,j-j3,3))*hdzi
            vtay=0.09*0.5*(say(i,j,1)+say(i,j-1,1))*
     1         sqrt(0.5*(rka(i,j,1)+rka(i,j-1,1)))
            vtaz=0.09*0.5*(saz(i,j,1)+saz(i,j-1,1))*
     1         sqrt(0.5*(rka(i,j,1)+rka(i,j-1,1)))
            vtb=0.09*0.5*(sb(i,j,1)+sb(i,j-1,1))*
     1         sqrt(0.5*(rkb(i,j,1)+rkb(i,j-1,1)))
            rhoa=0.5*(rho(i,j,1)+rho(i,j-1,1))/
     1           (0.5*(gi(i,j,1)+gi(i,j-1,1)))
            d23a(i,j,1) = -rhoa*( gii*vza*vtaz+ wya*vtay+ g23*wza*vtaz)
            d23a(i,j,L+1) = -0.5*(d23(i,j-j3,L+1)+d23(i,j,L+1))
            d23b(i,j,1) = -vtb*rhoa*( gii*vza+ wya+ g23*wza)
            d23b(i,j,L+1) = -0.5*(d23(i,j-j3,L+1)+d23(i,j,L+1))
          end do
        end do

        if (topdedge.eq.0) then                         !add d rrl
         call updated(d23a,d23a,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
         call updated(d23b,d23b,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
        else
         call updated(d23a,d23a,np,mp+1,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
         call updated(d23b,d23b,np,mp+1,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
        end if

        if(botdedge.eq.1) then                         !add d rrl
        do k=1,l+1
        do i=1,np
        d23a(i,1,k)=(1-ibcy)*d23a(i,2,k)+ibcy*d23a(i,-1,k)
        d23b(i,1,k)=(1-ibcy)*d23b(i,2,k)+ibcy*d23b(i,-1,k)
        enddo
        enddo
        endif

        if(topdedge.eq.1) then                         !add d rrl
        do k=1,l+1
        do i=1,np
        d23a(i,mp+1,k)=(1-ibcy)*d23a(i,mp,k)+ibcy*d23a(i,mp+3,k)
        d23b(i,mp+1,k)=(1-ibcy)*d23b(i,mp,k)+ibcy*d23b(i,mp+3,k)
        enddo
        enddo
        endif

        else
        do k=2,l
        do i=1,np
        vtaz=0.09*0.5*(saz(i,1,k)+saz(i,1,k-1))*
     .  sqrt(0.5*(rka(i,1,k)+rka(i,1,k-1)))
        vtb=0.09*0.5*(sb(i,1,k)+sb(i,1,k-1))*
     .sqrt(0.5*(rkb(i,1,k)+rkb(i,1,k-1)))
        gia=0.5*(gi(i,1,k)+gi(i,1,k-1))
        rhoa=0.5*(rho(i,1,k)+rho(i,1,k-1))/gia
        d23a(i,1,k)=-vtaz*rhoa*gia*d23(i,1,k)
        d23a(i,2,k)=d23a(i,1,k)
        d23b(i,1,k)=-vtb*rhoa*gia*d23(i,1,k)
        d23b(i,2,k)=d23b(i,1,k)
        enddo
        enddo
        do i=1,np
        vtaz=0.09*saz(i,1,1)*sqrt(rka(i,1,1))
        vtb=0.09*sb(i,1,1)*sqrt(rkb(i,1,1))
        d23a(i,1,1)=-vtaz*gi(i,1,1)*d23(i,1,1)
        d23a(i,2,1)=d23a(i,1,1)
        d23a(i,1,l+1)=-d23(i,1,l+1)
        d23a(i,2,l+1)=-d23(i,1,l+1)
        d23b(i,1,1)=-vtb*gi(i,1,1)*d23(i,1,1)
        d23b(i,2,1)=d23b(i,1,1)
        d23b(i,1,l+1)=-d23(i,1,l+1)
        d23b(i,2,l+1)=-d23(i,1,l+1)
        enddo

        endif

      do k=1,l+1
      do j=1,mp+1
      do i=1,np
      d23t(i,j,k)=d23a(i,j,k)
      if(iturb.eq.2) 
     .     d23t(i,j,k)=d23t(i,j,k)+d23b(i,j,k)+0.2*d23b(i,j,k)
      enddo
      enddo
      enddo

      if (topdedge.eq.0) then                         !add d rrl
         call updated(d23a,d23a,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      else
         call updated(d23a,d23a,np,mp+1,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      end if
      if (topdedge.eq.0) then                         !add d rrl
         call updated(d23b,d23b,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      else
         call updated(d23b,d23b,np,mp+1,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      end if
      if (topdedge.eq.0) then                         !add d rrl
         call updated(d23t,d23t,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      else
         call updated(d23t,d23t,np,mp+1,l+1,1-ih,np+ih,1-ih,mp+ih+1,1)
      end if

      illim = 1+leftdedge                         !add d rrl
      iulim = np
      do k=2,L
        do j=1,mp
          do i=iulim,illim,-1
            gii=0.25*(gi(i-1,j,k)+gi(i,j,k)
     .               +gi(i-1,j,k-1)+gi(i,j,k-1))
            g13=0.25*(gmul(k-1)+gmul(k))*(c13(i-1,j)+c13(i,j))
            wxa=hdxi*((w(i,j,k-1)-w(i-1,j,k-1))+(w(i,j,k)-w(i-1,j,k)))
            uza=0.5*(d13(i-1,j,k)+d13(i,j,k))
            wza=hdzi*((w(i-1,j,k)-w(i-1,j,k-1))+(w(i,j,k)-w(i,j,k-1)))
       vtax=0.09*0.25*(sax(i,j,k)+sax(i-1,j,k)+sax(i-1,j,k-1)+sax(i,j,k-1))
     .   *sqrt(0.25*(rka(i,j,k)+rka(i-1,j,k)
     .              +rka(i-1,j,k-1)+rka(i,j,k-1)))
       vtaz=0.09*0.25*(saz(i,j,k)+saz(i-1,j,k)+saz(i-1,j,k-1)+saz(i,j,k-1))
     .   *sqrt(0.25*(rka(i,j,k)+rka(i-1,j,k)
     .              +rka(i-1,j,k-1)+rka(i,j,k-1)))
       vtb=0.09*0.25*(sb(i,j,k)+sb(i-1,j,k)+sb(i-1,j,k-1)+sb(i,j,k-1))
     .   *sqrt(0.25*(rkb(i,j,k)+rkb(i-1,j,k)
     .              +rkb(i-1,j,k-1)+rkb(i,j,k-1)))
       rhoa=0.25*(rho(i,j,k)+rho(i-1,j,k)+rho(i-1,j,k-1)+rho(i,j,k-1))/
     .      gii
            d13a(i,j,k) = -rhoa*(gii*uza*vtaz + wxa*vtax + g13*wza*vtaz)
            d13b(i,j,k) = -vtb*rhoa*(gii*uza + wxa + g13*wza)
   
150   format(e30.16)
151   format(e30.16,2x,e30.16)
c     endif
          end do
        end do
      end do

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d13a,d13a,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(d13a,d13a,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d13b,d13b,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(d13b,d13b,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      do j=1,mp
        do i=iulim,illim,-1
            gii=0.5*(gi(i-1,j,1)+gi(i,j,1))
            g13=.25*(3.*gmul(1)-gmul(2))*(c13(i-1,j)+c13(i,j))
            wxa=hdxi*(3.*(w(i,j,1)-w(i-1,j,1))-(w(i,j,2)-w(i-1,j,2)))
            uza=0.5*(d13(i-1,j,1)+d13(i,j,1))
            wza=(3.*w(i  ,j,2)-2.*w(i  ,j,1)-w(i  ,j,3))*hdzi
     1         +(3.*w(i-1,j,2)-2.*w(i-1,j,1)-w(i-1,j,3))*hdzi
            vtax=0.09*0.5*(sax(i,j,1)+sax(i-1,j,1))
     1             *sqrt(0.5*(rka(i,j,1)+rka(i-1,j,1)))
            vtaz=0.09*0.5*(saz(i,j,1)+saz(i-1,j,1))
     1             *sqrt(0.5*(rka(i,j,1)+rka(i-1,j,1)))
            vtb=0.09*0.5*(sb(i,j,1)+sb(i-1,j,1))
     1             *sqrt(0.5*(rkb(i,j,1)+rkb(i-1,j,1)))
            rhoa=0.5*(rho(i,j,1)+rho(i-1,j,1))
            d13a(i,j,1) = -rhoa*( gii*uza*vtaz+ wxa*vtax+ g13*wza*vtaz)
            d13b(i,j,1) = -vtb*rhoa*( gii*uza+ wxa+ g13*wza)
            d13a(i,j,L+1) =-0.5*(d13(i-1,j,L+1)-d13(i,j,L+1))
            d13b(i,j,L+1) =-0.5*(d13(i-1,j,L+1)-d13(i,j,L+1))
c           if(i.eq.43.and.mpi_rank.eq.0) print*,d13b(i,j,k),'d13b'
        end do
      end do

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d13a,d13a,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(d13a,d13a,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d13b,d13b,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(d13b,d13b,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      if(leftdedge.eq.1) then                         !add d rrl
      do k=1,l+1
      do j=1,mp
      d13a(1,j,k)=(1-ibcx)*d13a(2,j,k) + ibcx*d13a(-1,j,k)
      d13b(1,j,k)=(1-ibcx)*d13b(2,j,k) + ibcx*d13b(-1,j,k)
      enddo
      enddo
      endif

      if(rightdedge.eq.1) then                         !add d rrl
      do k=1,l+1
      do j=1,mp
      d13a(np+1,j,k)=(1-ibcx)*d13a(np,j,k) + ibcx*d13a(np+3,j,k)
      d13b(np+1,j,k)=(1-ibcx)*d13b(np,j,k) + ibcx*d13b(np+3,j,k)
      enddo
      enddo
      endif

      do k=1,l+1
      do j=1,mp
      do i=1,np+1
      d13t(i,j,k)=d13a(i,j,k)
      if(iturb.eq.2) 
     .     d13t(i,j,k)=d13t(i,j,k)+d13b(i,j,k)+0.2*d13b(i,j,k)
      enddo
      enddo
      enddo

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d13a,d13a,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(d13a,d13a,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d13b,d13b,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(d13b,d13b,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      if (rightdedge.eq.0) then                         !add d rrl
         call updated(d13t,d13t,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(d13t,d13t,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih,1)
      end if
c     if(mpi_rank.eq.0) print*,'d13t',d13a(51,1,2)

      deallocate (u)
      deallocate (v)
      deallocate (w)
      deallocate (rho)
      deallocate (cdrg)
      deallocate (uka)
      deallocate (ukb)
      deallocate (rka)
      deallocate (rkb)
      deallocate (tmp)

      return
      end
