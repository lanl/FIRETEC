      subroutine lapdi(xv,fi,s,d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use metryic
      use turba
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real      xv(il:iu, jl:ju,lls,nvp),
     .          s(il:iu, jl:ju,lls),
     .          fi(il:iu, jl:ju,lls)
      real            d13(il:iu+1, jl:ju,lls+1),
     .                d12(il:iu+1, jl:ju+1,lls),
     .                d23(il:iu, jl:ju+1,lls+1),
     .                d11(il:iu+1, jl:ju,lls),
     .                d22(il:iu, jl:ju+1,lls),
     .                d33(il:iu, jl:ju,lls+1)
            real,allocatable::
     .           r(:, :,:),
     .          ri(:, :,:),
     .         rkb(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1,kmz
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: rkba
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2,hxp1,hxp2,hyp1,hyp2,hzp1,hzp2
      real :: hxn,hxnp1,hymp1,hzLp1
      real :: coef,d12a,d13a,d21a,d23a,d31a,d32a

      if(iturb.ge.1) then
      allocate (ri(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (rkb(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
      endif


 
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi

      do k=1,l 
      do j=1,mp 
      do i=1,np 
      ri(i,j,k)=xv(i,j,k,4)/xv(i,j,k,nv)
      rkb(i,j,k)=xv(i,j,k,5)/xv(i,j,k,nv)
      enddo
      enddo
      enddo
      call updated(ri,ri,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rkb,rkb,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(ri(i,j,k)-ri(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
          pz(i,j,L+1)=-pz(i,j,L)
        end do
      end do

      do j=1,mp
      do i=1,np
        srff(i,j)=0.
      enddo
      enddo

      do j=1,mp
      do i=1,np
      g33=(c13(i,j)*gmul(1))**2+(c23(i,j)*gmul(1))**2+gi(i,j,1)**2
      srff(i,j) = sqrt(g33)*srff(i,j)
      enddo
      enddo
      jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
      julim = mp + (ibcy-j3)*topdedge                       !add d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then                       !add d rrl
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then                       !add d rrl
               jm1 = -1
            else
               jm1 = j - 1
            end if
         else
            jp1=1
            jm1=1
         end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
         illim = 1  + (1-ibcx)*leftdedge                       !add d rrl
         iulim = np + (ibcx-1)*rightdedge                       !add d rrl
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then                       !add d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then                       !add d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(ri(ip1,j,1)-ri(im1,j,1))
      py=hdyi*(ri(i,jp1,1)-ri(i,jm1,1))*j3
   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)                       !add d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)                       !add d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
            julim = mp + (ibcy-j3)*topdedge                       !add d rrl
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then                       !add d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then                       !add d rrl
                     jm1 = -1
                  else
                     jm1 = j - 1
                  end if
               else
                  jp1=1
                  jm1=1
               end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(ri(i,jp1,1)-ri(i,jm1,1))*j3
 2111      pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)                       !add d rrl
              julim = mp*topdedge +  1*(1-topdedge)                       !add d rrl
              do 2112 j=jllim,julim,mp-j3
 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
     .                  pz(i,j,2)
           endif
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)                       !add d rrl
         julim = mp*topdedge +  1*(1-topdedge)                       !add d rrl
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                       !add d rrl

            iulim = np + (ibcx-1)*rightdedge                       !add d rrl
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then                       !add d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then                       !add d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(ri(ip1,j,1)-ri(im1,j,1))
 2121      pz(i,j,1)=-2.*(srff(i,j)+g13*px)/(g33-g23*g23)-pz(i,j,2)
 212    continue
      endif

      if (rightdedge.eq.0 .and. topdedge.eq.0) then                       !add d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                       !add d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                       !add d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

       
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np                       !add d rrl
          do k=1,L
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(ri(i,j,k)-ri(i-1,j,k))
            hx(i,j,k)=( pxa + g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k)

        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c
compute y-flux at (i,j+-1/2,k)
      if (j3.eq.1) then
        do i=1,np
          do j=1+botdedge,mp                       !add d rrl
            do k=1,L
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(ri(i,j,k)-ri(i,j-j3,k))
              hy(i,j,k)=( pya + g23*pza )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                       !add d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif

compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= gii*hza + g13*hxa
          end if
        end do
      end do
! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+

     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                       !add d rrl
        do k=2,L
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.25*(hy(i,mp,k-1)+hy1+hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        end do
      end if
      endif
create boundary conditions at k=1; for k=L see divergence below
c surface fluxes:
      do j=1,mp
      do i=1,np
      hz(i,j,1)=-hz(i,j,2)
      end do
      end do
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
  
compute u
      do k=1,l
      do j=1,mp
      do i=1+leftdedge,np                       !add d rrl
      rkba=0.5*(rkb(i,j,k)+rkb(i-1,j,k))
c     coef=0.09*0.5*(sc(i,j,k)+sc(i-1,j,k))/sqrt(0.2*rkba)
      coef=0.09*0.5*(s(i,j,k)+s(i-1,j,k))/sqrt(rkba)
      d12a=0.
      if(j3.eq.1) d12a=0.5*(d12(i,j,k)+d12(i,j+1,k))
      d13a=0.5*(d13(i,j,k)+d13(i,j,k+1))
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i-1,j+1,k)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)                        !add d rrl
     .                        hyp1=hy(i,j,k)*(ibcy-1)
     .                            +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)                        !add d rrl
     .                        hyp2=hy(i-1,j,k)*(ibcy-1)
     .                            +hy(i-1,mp+2,k)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i-1,j,k)+hyp2)
      if(k.ne.l) hzp1=hz(i,j,k+1)
      if(k.ne.l) hzp2=hz(i-1,j,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i-1,j,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i-1,j,k)+hzp2)
      u(i,j,k)=0.2*coef*(d11(i,j,k)*hx(i,j,k)+d12a*hya+
     .               d13a*hza)
      enddo
      enddo
      enddo
      call updated(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          u(1,j,k) = (ibcx-1)*u(2,j,k) + ibcx*u(0,j,k)
        end do
      end do
      endif

compute v
      if(j3.eq.1) then
      do k=1,l
      do j=1+botdedge,mp                       !add d rrl
      do i=1,np
      rkba=0.5*(rkb(i,j,k)+rkb(i,j-1,k))
c     coef=0.09*0.5*(sc(i,j,k)+sc(i,j-1,k))/sqrt(0.2*rkba)
      coef=0.09*0.5*(s(i,j,k)+s(i,j-1,k))/sqrt(rkba)
      d21a=0.5*(d12(i,j,k)+d12(i+1,j,k))
      d23a=0.5*(d23(i,j,k)+d23(i,j,k+1))
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j-1,k)
      if(rightdedge.eq.1.and.i.eq.np) hxp1=hx(i,j,k)*(ibcx-1)                       !add d rrl
     .                                   +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np) hxp2=hx(i,j-1,k)*(ibcx-1)                       !add d rrl
     .                                   +hx(np+2,j-1,k)*ibcx
      hxa=0.25*(hx(i,j,k)+hxp1+hx(i,j-1,k)+hxp2)
      if(k.ne.l) hzp1=hz(i,j,k+1)
      if(k.ne.l) hzp2=hz(i,j-1,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i,j-1,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i,j-1,k)+hzp2)
      v(i,j,k)=0.2*coef*(d21a*hxa+d22(i,j,k)*hy(i,j,k)+
     .                   d23a*hza)
      enddo
      enddo
      enddo
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(botdedge.eq.1) then                       !add d rrl
      do k=1,l
      do i=1,np
      v(i,1,k)= (ibcy-1)*v(i,2,k) + ibcy*v(i,0,k)
      end do
      end do
      endif
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif
compute w
      do k=2,l
      do j=1,mp
      do i=1,np
      kmz=k-1
      rkba=0.5*(rkb(i,j,k)+rkb(i,j,k-1))
c     coef=0.09*0.5*(sc(i,j,k)+sc(i,j,kmz))/sqrt(0.2*rkba)
      coef=0.09*0.5*(s(i,j,k)+s(i,j,kmz))/sqrt(rkba)
      d31a=0.5*(d13(i,j,k)+d13(i+1,j,k))
      d32a=0.
      if(j3.eq.1) d32a=0.5*(d23(i,j,k)+d23(i,j+1,k))
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j,k-1)
      if(rightdedge.eq.1.and.i.eq.np) hxp1=hx(i,j,k)*(ibcx-1)                       !add d rrl
     .                                   +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np) hxp2=hx(i,j,k-1)*(ibcx-1)                       !add d rrl
     .                                   +hx(np+2,j,k-1)*ibcx
      hxa=0.25*(hx(i,j,k)+hx(i,j,k-1)+hxp1+hxp2)
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i,j+1,k-1)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)                        !add d rrl
     .                             hyp1=hy(i,j,k)*(ibcy-1)
     .                                 +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)                        !add d rrl
     .                             hyp2=hy(i,j,k-1)*(ibcy-1)
     .                                 +hy(i,mp+2,k-1)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i,j,k-1)+hyp2)
      w(i,j,k)=0.2*coef*(d31a*hxa+d32a*hya+
     .             d33(i,j,k)*hz(i,j,k))
      enddo
      enddo
      enddo

      do j=1,mp
      do i=1,np
      w(i,j,1)=-w(i,j,2)
      enddo
      enddo
       
      do k=1,l
      do j=1,mp
      do i=1,np
      hx(i,j,k)=u(i,j,k)
      hy(i,j,k)=v(i,j,k)
      hz(i,j,k)=w(i,j,k)
      enddo
      enddo
      enddo
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute Laplacian term by term
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.
      end do
      end do
      end do
compute d/dx(dh/dx)
      do k=1,L
        do j=1,mp
          do i=1,np-rightdedge                       !add d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                       !add d rrl
         do k=1,L
            do j=1,mp
               hxn=hx(np,j,k)
               hxnp1 = (ibcx-1)*hxn + ibcx*hx(np+2,j,k)
               r(np,j,k) = dxi*(hxnp1-hxn)
            end do
         end do
      end if
c for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                       !add d rrl
          do k=1,L
             do i=1,np
create boundary conditions at j=m+1
                hymp1= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
       end if
      endif
compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,L-1
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
          hzLp1 =-hz(i,j,L)
          r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do

      call updated(fi,fi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      do k=1,l
      do j=1,mp
      do i=1,np
       r(i,j,k)=r(i,j,k)*gi(i,j,k)
      fi(i,j,k)=fi(i,j,k)+2.*r(i,j,k)*dt
      enddo
      enddo
      enddo
      call updated(fi,fi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
150   format(5x,e30.15)
151   format(e30.15,2x,e30.15)

      deallocate (ri)
      deallocate (r)
      deallocate (rkb)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
 
      return
      end subroutine lapdi
!**************************end subroutine lapdi*****************************!
      !subroutine lapdis(xv,fi,s,
      !.d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls,nvp)
      !JAS 3/7/06 changed routine signature here and in calling routine
      !(turb) as "s" and "d-arrays" are not used!
      subroutine lapdis(xv,fi,
     .il,iu,jl,ju,lls,nvp)
      use gridsetup
      use metryic
      use xve, only:xe
      use turba
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real      xv(il:iu, jl:ju,lls,nvp),
      !.          s(il:iu, jl:ju,lls),
     .          fi(il:iu, jl:ju,lls)
      !real            d13(il:iu+1, jl:ju,lls+1),
      !.                d12(il:iu+1, jl:ju+1,lls),
      !.                d23(il:iu, jl:ju+1,lls+1),
      !.                d11(il:iu+1, jl:ju,lls),
      !.                d22(il:iu, jl:ju+1,lls),
      !.                d33(il:iu, jl:ju,lls+1)
            real,allocatable::
     .           r(:, :,:),
     .          ri(:, :,:),
     .         rkb(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2
      real :: hxn,hxnp1,hymp1,hzLp1
      real :: coef,dcr1

      if(iturb.ge.1) then
      allocate (ri(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (rkb(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l+1))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
      endif


 
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi

      do k=1,l 
      do j=1,mp 
      do i=1,np 
      ri(i,j,k)=xv(i,j,k,4)/xv(i,j,k,nv)
      rkb(i,j,k)=xv(i,j,k,5)/xv(i,j,k,nv)  !FP
      enddo
      enddo
      enddo
      call updated(ri,ri,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rkb,rkb,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(ri(i,j,k)-ri(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
c          pz(i,j,L+1)=-pz(i,j,L)
          pz(i,j,L+1)=dzi*(xe(i,j,l,4)/xe(i,j,l,nv)-ri(i,j,l))   !rrl
        end do
      end do

      do j=1,mp
      do i=1,np
        srff(i,j)=0.
      enddo
      enddo

      do j=1,mp
      do i=1,np
      g33=(c13(i,j)*gmul(1))**2+(c23(i,j)*gmul(1))**2+gi(i,j,1)**2
      srff(i,j) = sqrt(g33)*srff(i,j)
      enddo
      enddo
      jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
      julim = mp + (ibcy-j3)*topdedge                       !add d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then                       !add d rrl
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then                       !add d rrl
               jm1 = -1
            else
               jm1 = j - 1
            end if
         else
            jp1=1
            jm1=1
         end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
         illim = 1  + (1-ibcx)*leftdedge                       !add d rrl
         iulim = np + (ibcx-1)*rightdedge                       !add d rrl
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then                       !add d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then                       !add d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(ri(ip1,j,1)-ri(im1,j,1))
      py=hdyi*(ri(i,jp1,1)-ri(i,jm1,1))*j3
c   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
   21 pz(i,j,1)=-1.*(srff(i,j)+g13*px+g23*py)/g33   !rrl (insulated ground)
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)                       !add d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)                       !add d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
            julim = mp + (ibcy-j3)*topdedge                       !add d rrl
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then                       !add d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then                       !add d rrl
                     jm1 = -1
                  else
                     jm1 = j - 1
                  end if
               else
                  jp1=1
                  jm1=1
               end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(ri(i,jp1,1)-ri(i,jm1,1))*j3
c 2111      pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
 2111      pz(i,j,1)=-1.*(srff(i,j)+g23*py)/(g33-g13*g13)   !rrl (insulated ground)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)                       !add d rrl
              julim = mp*topdedge +  1*(1-topdedge)                       !add d rrl
              do 2112 j=jllim,julim,mp-j3
c 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
c     .                  pz(i,j,2)
 2112         pz(i,j,1)=-1.*srff(i,j)/(gi(i,j,1)**2)   !rrl (insulated ground)
           endif
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)                       !add d rrl
         julim = mp*topdedge +  1*(1-topdedge)                       !add d rrl
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                       !add d rrl

            iulim = np + (ibcx-1)*rightdedge                       !add d rrl
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then                       !add d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then                       !add d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(ri(ip1,j,1)-ri(im1,j,1))
 2121      pz(i,j,1)=-1.*(srff(i,j)+g13*px)/(g33-g23*g23)  !rrl (insulated ground)
 212    continue
      endif

      if (rightdedge.eq.0 .and. topdedge.eq.0) then                       !add d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                       !add d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                       !add d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

       
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np                       !add d rrl
          do k=1,L
c            coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i-1,j,k))*   
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i-1,j,k)))
c     .      *0.5*(xv(i-1,j,k,nv)+xv(i,j,k,nv))
            if (ist.eq.0) then
              coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i-1,j,k))*   
     .           sqrt(0.5*(rkb(i,j,k)+rkb(i-1,j,k)))
     .        *0.5*(xv(i-1,j,k,nv)+xv(i,j,k,nv))
            else !ist.eq.1
              coef=0.09*rturbprandtl*0.5*(sa(i,j,k)+sa(i-1,j,k))*   
     .           0.5*(sqrtk(i,j,k)+sqrtk(i-1,j,k))
     .         *0.5*(xv(i-1,j,k,nv)+xv(i,j,k,nv))
            endif
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(ri(i,j,k)-ri(i-1,j,k))
            hx(i,j,k)=coef*( pxa + g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k)

        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c
compute y-flux at (i,j+-1/2,k)
      if (j3.eq.1) then
        do i=1,np
          do j=1+botdedge,mp                       !add d rrl
            do k=1,L
c            coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j-1,k))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j-1,k)))
c     .      *0.5*(xv(i,j-1,k,nv)+xv(i,j,k,nv))
            if (ist.eq.0) then
                coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i,j-1,k))*
     .           sqrt(0.5*(rkb(i,j,k)+rkb(i,j-1,k)))
     .         *0.5*(xv(i,j-1,k,nv)+xv(i,j,k,nv))
            else !ist.eq.1
                coef=0.09*rturbprandtl*0.5*(sa(i,j,k)+sa(i,j-1,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j-1,k))
     .         *0.5*(xv(i,j-1,k,nv)+xv(i,j,k,nv))
            endif
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(ri(i,j,k)-ri(i,j-j3,k))
              hy(i,j,k)=coef*( pya + g23*pza )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                       !add d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif

compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
c     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
             if (ist.eq.0) then
               coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i,j,k-1))*
     .           sqrt(0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
     .          *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
             else !ist.eq.1
               coef=0.09*rturbprandtl*0.5*(sa(i,j,k)+sa(i,j,k-1))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j,k-1))
     .          *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
             endif
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
c     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
             if (ist.eq.0) then
                coef=0.09*rturbprandtl*0.5*(sb(np,j,k)+sb(np,j,k-1))*
     .           sqrt(0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
     .         *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
            else !ist.eq.1
                coef=0.09*rturbprandtl*0.5*(sa(np,j,k)+sa(np,j,k-1))*
     .           0.5*(sqrtk(np,j,k)+sqrtk(np,j,k-1))
     .         *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
            endif
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
     .      write(6,*) 'diffusion too strong in theta',i,j,k
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do
      end do
      k=1
      do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
c     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
c             coef=0.66*0.09*sb(i,j,k)*
c     .           sqrt(1.2*rkb(i,j,k))
c     .          *xv(i,j,k,nv)
            if (ist.eq.0) then
             coef=0.09*rturbprandtl*sb(i,j,k)*
     .           sqrt(rkb(i,j,k))
     .          *xv(i,j,k,nv)
            else !ist.eq.1
             coef=0.09*rturbprandtl*sa(i,j,k)*
     .           sqrtk(i,j,k)
     .          *xv(i,j,k,nv)
            endif
c234567
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
!     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
     .      write(6,*) 'diffusion too strong in theta',i,j,k
            g13=gmul(k)*c13(i,j)
            gii=gi(i,j,k)
            hxa=0.5*(hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
c     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
cc             coef=0.66*0.09*sb(i,j,k)*
c     .           sqrt(1.2*rkb(i,j,k))
c     .      *xv(i,j,k,nv)
            if (ist.eq.0) then
             coef=0.09*rturbprandtl*sb(i,j,k)*
     .           sqrt(rkb(i,j,k))
     .      *xv(i,j,k,nv)
            else !ist.eq1
             coef=0.09*rturbprandtl*sa(i,j,k)*
     .           sqrtk(i,j,k)
     .      *xv(i,j,k,nv)
            endif

            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
     .      write(6,*) 'diffusion too strong in theta',i,j,k
             g13=gmul(k)*c13(np,j)
             gii=gi(np,j,k)
c             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.5*(hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do
        k=l+1
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
c     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
c             coef=0.66*0.09*sb(i,j,k-1)*
c     .           sqrt(1.2*rkb(i,j,k-1))
c     .      *xv(i,j,k-1,nv)
            if (ist.eq.0) then
             coef=0.09*rturbprandtl*sb(i,j,k-1)*
     .           sqrt(rkb(i,j,k-1))
     .      *xv(i,j,k-1,nv)
            else !ist.eq.1
             coef=0.09*rturbprandtl*sa(i,j,k-1)*
     .           sqrtk(i,j,k-1)
     .      *xv(i,j,k-1,nv)
            endif

            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=gmul(k-1)*c13(i,j)
            gii=gi(i,j,k-1)
            hxa=0.5*(hx(i,j,k-1)+hx(i+1,j,k-1))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
c     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
c             coef=0.66*0.09*sb(i,j,k-1)*
c     .           sqrt(1.2*rkb(i,j,k-1))
c     .      *xv(i,j,k-1,nv)
             if (ist.eq.0) then
             coef=0.09*rturbprandtl*sb(i,j,k-1)*
     .           sqrt(rkb(i,j,k-1))
     .      *xv(i,j,k-1,nv)
             else !ist.eq.1
             coef=0.09*rturbprandtl*sa(i,j,k-1)*
     .           sqrtk(i,j,k-1)
     .      *xv(i,j,k-1,nv)
             endif

            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))

             g13=gmul(k-1)*c13(np,j)
             gii=gi(np,j,k-1)
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
c             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=.5*(hx(np,j,k-1)+hx1)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do

! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+

     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                       !add d rrl
        do k=2,L
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.25*(hy(i,mp,k-1)+hy1+hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        end do
      end if
      k=1
         do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=gmul(k)*c23(i,j)
              hya=0.5*(
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
        end do
      if (topdedge.eq.1) then                       !add d rrl
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=gmul(k)*c23(i,mp)
c            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.5*(hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
      end if
      k=l+1
         do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=gmul(k-1)*c23(i,j)
              hya=0.5*(hy(i,j,k-1)+hy(i,j+j3,k-1))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
      if (topdedge.eq.1) then                       !add d rrl
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=gmul(k-1)*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
c            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.5*(hy(i,mp,k-1)+hy1)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
      end if

      endif
create boundary conditions at k=1; for k=L see divergence below
c surface fluxes:
c      do j=1,mp
c      do i=1,np
c      hz(i,j,1)=-hz(i,j,2)
c      end do
c      end do
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
  

compute Laplacian term by term
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.
      end do
      end do
      end do
compute d/dx(dh/dx)
      do k=1,L
        do j=1,mp
          do i=1,np-rightdedge                       !add d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                       !add d rrl
         do k=1,L
            do j=1,mp
               hxn=hx(np,j,k)
               hxnp1 = (ibcx-1)*hxn + ibcx*hx(np+2,j,k)
               r(np,j,k) = dxi*(hxnp1-hxn)
            end do
         end do
      end if
c for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                       !add d rrl
          do k=1,L
             do i=1,np
create boundary conditions at j=m+1
                hymp1= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
       end if
      endif
compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,L
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
          hzLp1 =-hz(i,j,L)
          r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do

      call updated(fi,fi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      do k=1,l
      do j=1,mp
      do i=1,np
       r(i,j,k)=r(i,j,k)*gi(i,j,k)
      fi(i,j,k)=fi(i,j,k)+2.*r(i,j,k)*dt
      enddo
      enddo
      enddo
      call updated(fi,fi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
150   format(5x,e30.15)
151   format(e30.15,2x,e30.15)

      deallocate (ri)
      deallocate (r)
      deallocate (rkb)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
 
      return
      end
