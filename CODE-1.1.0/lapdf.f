      subroutine lapdf(rk,fk,s,d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls)
      use gridsetup
      use turba
      use metryic
      use msga
      
      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls

      real      rk(il:iu, jl:ju,lls),
     .           s(il:iu, jl:ju,lls),
     .          fk(il:iu, jl:ju,lls)
      real            d13(il:iu+1, jl:ju,lls+1), 
     .                d12(il:iu+1, jl:ju+1,lls),
     .                d23(il:iu, jl:ju+1,lls+1),
     .                d11(il:iu+1, jl:ju,lls),
     .                d22(il:iu, jl:ju+1,lls),
     .                d33(il:iu, jl:ju,lls+1)
      real,allocatable::
     .           r(:, :,:),
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
      real :: saa,sqrtka 
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2,hxp1,hxp2,hyp1,hyp2,hzp1,hzp2
      real :: hxn,hxnp1,hymp1,hzLp1
      real :: coef,coefa,d12a,d13a,d21a,d23a,d31a,d32a


      if(iturb.ge.1) then
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
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
      call updated(rk,rk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(rk(i,j,k)-rk(i,j,k-1))
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
c This section needs to be addressed with new ibcx=0 bc rrl       
      jllim = 1  + (j3-ibcy)*botdedge                 !added d rrl
      julim = mp + (ibcy-j3)*topdedge                 !added d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then                 !added d rrl
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then                 !added d rrl
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
         illim = 1  + (1-ibcx)*leftdedge                 !added d rrl
         iulim = np + (ibcx-1)*rightdedge                 !added d rrl
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then                 !added d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then                 !added d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
      py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)                 !added d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)                 !added d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                 !added d rrl
            julim = mp + (ibcy-j3)*topdedge                 !added d rrl
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then                 !added d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then                 !added d rrl
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
               py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
 2111      pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)                 !added d rrl
              julim = mp*topdedge +  1*(1-topdedge)                 !added d rrl
              do 2112 j=jllim,julim,mp-j3
 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
     .                  pz(i,j,2)
           endif
c what is the value of px?  
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)                 !added d rrl
         julim = mp*topdedge +  1*(1-topdedge)                 !added d rrl
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                 !added d rrl
            iulim = np + (ibcx-1)*rightdedge                 !added d rrl
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then                 !added d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then                 !added d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
 2121      pz(i,j,1)=-2.*(srff(i,j)+g13*px)/(g33-g23*g23)-pz(i,j,2)
c  in this case what is the value of py
 212    continue
      endif
c the following lines use rightedge and topedge and need examination.  They might be ok.  rrl
      if (rightdedge.eq.0 .and. topdedge.eq.0) then                 !added d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                 !added d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                 !added d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if
 
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np                 !added d rrl
          do k=1,L 
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(rk(i,j,k)-rk(i-1,j,k))
            hx(i,j,k)=( pxa + g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then                 !added d rrl
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
          do j=1+botdedge,mp                 !added d rrl
            do k=1,L
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(rk(i,j,k)-rk(i,j-j3,k))
              hy(i,j,k)=( pya + g23*pza )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                 !added d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
      endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
 
compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                 !added d rrl
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 !added d rrl
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
            do j=1,mp-topdedge                 !added d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                 !added d rrl
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
      do i=1+leftdedge,np                 !added d rrl
      saa=0.5*(s(i-1,j,k)+s(i,j,k))
      sqrtka=0.5*(sqrtk(i-1,j,k)+sqrtk(i,j,k))
      coef=0.
      if(sqrtka.ne.0.) coef=saa/sqrtka
      coefa=coef*dt*(dxi**2+j3*dyi**2+dzi**2)
      d12a=0.
      if(j3.eq.1) d12a=0.5*(d12(i,j,k)+d12(i,j+1,k))
      d13a=0.5*(d13(i,j,k)+d13(i,j,k+1))
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i-1,j+1,k)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)                  !added d rrl
     .               hyp1=hy(i,j,k)*(ibcy-1)
     .                   +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)                  !added d rrl
     .               hyp2=hy(i-1,j,k)*(ibcy-1)
     .                   +hy(i-1,mp+2,k)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i-1,j,k)+hyp2)
      if(k.lt.l) hzp1=hz(i,j,k+1)
      if(k.lt.l) hzp2=hz(i-1,j,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i-1,j,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i-1,j,k)+hzp2)
c     if(abs(d11(i,j,k)*hx(i,j,k))*coefa.gt.0.5) print*,'error'
c     if(abs(d12a*hy(i,j,k))*coefa.gt.0.5) print*,'error'
c     if(abs(d13a*hz(i,j,k))*coefa.gt.0.5) print*,'error'
      u(i,j,k)=coef*(d11(i,j,k)*hx(i,j,k)+d12a*hya+d13a*hza)
      enddo
      enddo
      enddo
      call updated(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      if(leftdedge.eq.1) then                 !added d rrl
      do j=1,mp
        do k=1,L
          u(1,j,k) = (ibcx-1)*u(2,j,k) + ibcx*u(0,j,k)
        end do
      end do
      endif

compute v
      if(j3.eq.1) then
      do k=1,l
      do j=1+botdedge,mp                 !added d rrl
      do i=1,np
      saa=0.5*(s(i,j-1,k)+s(i,j,k))
      sqrtka=0.5*(sqrtk(i,j-1,k)+sqrtk(i,j,k))
      coef=0.
      if(sqrtka.ne.0.) coef=saa/sqrtka
      d21a=0.5*(d12(i,j,k)+d12(i+1,j,k))
      d23a=0.5*(d23(i,j,k)+d23(i,j,k+1))
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j-1,k)
      if(rightdedge.eq.1.and.i.eq.np) hxp1=hx(i,j,k)*(ibcx-1)                 !added d rrl
     .                                   +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np) hxp2=hx(i,j-1,k)*(ibcx-1)                 !added d rrl
     .                                   +hx(np+2,j-1,k)*ibcx
      hxa=0.25*(hx(i,j,k)+hxp1+hx(i,j-1,k)+hxp2)
      if(k.ne.l) hzp1=hz(i,j,k+1)
      if(k.ne.l) hzp2=hz(i,j-1,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i,j-1,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i,j-1,k)+hzp2)
      v(i,j,k)=coef*(d21a*hxa+d22(i,j,k)*hy(i,j,k)+d23a*hza)
      enddo
      enddo
      enddo
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(botdedge.eq.1) then                 !added d rrl
      do k=1,l
      do i=1,np
      v(i,1,k)= (ibcy-1)*v(i,2,k) + ibcy*v(i,0,k)
      end do
      end do
      endif
      endif
compute w
      do k=2,l
      do j=1,mp
      do i=1,np
      saa=0.5*(s(i,j,k)+s(i,j,k-1))
      sqrtka=0.5*(sqrtk(i,j,k)+sqrtk(i,j,k-1))
      if(sqrtka.eq.0.) coef=0.
      if(sqrtka.ne.0.) coef=saa/sqrtka
      d31a=0.5*(d13(i,j,k)+d13(i+1,j,k))
      d32a=0. 
      if(j3.eq.1) d32a=0.5*(d23(i,j,k)+d23(i,j+1,k))
      coefa=coef*dt*(dxi**2+j3*dyi**2+dzi**2)
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j,k-1)
      if(rightdedge.eq.1.and.i.eq.np)                 !added d rrl 
     .                    hxp1=hx(i,j,k)*(ibcx-1)
     .                        +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np)                  !added d rrl
     .                    hxp2=hx(i,j,k-1)*(ibcx-1)
     .                        +hx(np+2,j,k-1)*ibcx
      hxa=0.25*(hx(i,j,k)+hx(i,j,k-1)+hxp1+hxp2)
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i,j+1,k-1)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)                  !added d rrl
     .                    hyp1=hy(i,j,k)*(ibcy-1)
     .                        +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)                  !added d rrl
     .                    hyp2=hy(i,j,k-1)*(ibcy-1)
     .                        +hy(i,mp+2,k-1)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i,j,k-1)+hyp2)
c     if(abs(d31a*hx(i,j,k))*coefa.gt.0.5) print*,'error'
c     if(abs(d32a*hy(i,j,k))*coefa.gt.0.5) print*,'error'
c     if(abs(d33(i,j,k)*hz(i,j,k))*coefa.gt.0.5) print*,'error'
      w(i,j,k)=coef*(d31a*hxa+d32a*hya+d33(i,j,k)*hz(i,j,k))
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
c     w(i,j,k)=0.
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
          do i=1,np-rightdedge                 !added d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                 !added d rrl
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
            do j=1,mp-topdedge                 !added d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                 !added d rrl
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
 
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.09*2.*r(i,j,k)*gi(i,j,k)/3.
      fk(i,j,k)=fk(i,j,k)+2.*r(i,j,k)*dt
      end do
      end do
      end do
      call updated(fk,fk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      deallocate (r)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
 
      return
      end subroutine lapdf
!*************************end subroutine lapdf*****************************!
      !subroutine lapdfs(rk,fk,s,
      !.d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls,iflg)
      !JAS 3/7/06 changed routine signature here 
      ! as "s" and "d-arrays" are not used!
      subroutine lapdfs(rk,fk,
     .il,iu,jl,ju,lls,iflg)
      use gridsetup
      use turba
      use xvo
      use xve, only:xe
      use metryic
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls

      real      rk(il:iu, jl:ju,lls),
      !.           s(il:iu, jl:ju,lls),
     .          fk(il:iu, jl:ju,lls)
      !real            d13(il:iu+1, jl:ju,lls+1), 
      !.                d12(il:iu+1, jl:ju+1,lls),
      !.                d23(il:iu, jl:ju+1,lls+1),
      !.                d11(il:iu+1, jl:ju,lls),
      !.                d22(il:iu, jl:ju+1,lls),
      !.                d33(il:iu, jl:ju,lls+1)
      real,allocatable::
     .           r(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)!,
!     .         rkb(:,:,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1,iflg
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2
      real :: hxn,hxnp1,hymp1
      real :: coefx,coefy,coefz,aswitch ! ,dcr1

      if(iturb.ge.1) then
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l+1))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
!      if (iturb.eq.2) then
!        allocate (rkb(1-ih:np+ih, 1-ih:mp+ih,l))
!        rkb(:,:,:)=xvb(:,:,:,6)/xvb(:,:,:,nv)
!      endif
      endif
 
c      aswitch=-1.*real(iflg-6)
      aswitch=0.
c      write (*,*) 'iflg=',iflg,' aswitch= ',aswitch
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi
      call updated(rk,rk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c      write (*,*) 'inside lapdfs',xe(5,5,l,iflg),xe(5,5,l,nv)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(rk(i,j,k)-rk(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
          pz(i,j,L+1)=dzi*(xe(i,j,l,iflg)/xe(i,j,l,nv)-rk(i,j,l))
c      write (*,*) 'after pz top',i,j,dzi,xe(i,j,l,iflg),
c     &                xe(i,j,l,nv),
c     &                rk(i,j,l)
c23456
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
       
      jllim = 1  + (j3-ibcy)*botdedge                 !added d rrl
      julim = mp + (ibcy-j3)*topdedge                 !added d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then                 !added d rrl
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then                 !added d rrl
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
         illim = 1  + (1-ibcx)*leftdedge                 !added d rrl
         iulim = np + (ibcx-1)*rightdedge                 !added d rrl
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then                 !added d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then                 !added d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
      py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
c   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
   21 pz(i,j,1)=-1.*(srff(i,j)+g13*px+g23*py)/g33
     &          +aswitch*2.*dzi*rk(i,j,1)
c   21 pz(i,j,1)=-1.*(srff(i,j)+g13*px+g23*py)/g33
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)                 !added d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)                 !added d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                 !added d rrl
            julim = mp + (ibcy-j3)*topdedge                 !added d rrl
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then                 !added d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then                 !added d rrl
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
               py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
c 2111  pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
c234567
 2111  pz(i,j,1)=-1.*(srff(i,j)+g23*py)/(g33-g13*g13)
     &           +aswitch*2.*dzi*rk(i,j,1)
c 2111  pz(i,j,1)=-1.*(srff(i,j)+g23*py)/(g33-g13*g13)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)                 !added d rrl
              julim = mp*topdedge +  1*(1-topdedge)                 !added d rrl
              do 2112 j=jllim,julim,mp-j3
c 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
c     .                  pz(i,j,2)
 2112         pz(i,j,1)=-1.*srff(i,j)/(gi(i,j,1)**2)
     &                  +aswitch*2.*dzi*rk(i,j,1)
c 2112         pz(i,j,1)=-1.*srff(i,j)/(gi(i,j,1)**2)
           endif
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)                 !added d rrl
         julim = mp*topdedge +  1*(1-topdedge)                 !added d rrl
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                 !added d rrl
            iulim = np + (ibcx-1)*rightdedge                 !added d rrl
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then                 !added d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then                 !added d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
c 2121      pz(i,j,1)=-2.*(srff(i,j)+g13*px)/(g33-g23*g23)-pz(i,j,2)
 2121      pz(i,j,1)=-1.*(srff(i,j)+g13*px)/(g33-g23*g23)
     &                  +aswitch*2.*dzi*rk(i,j,1)
c 2121      pz(i,j,1)=-1.*(srff(i,j)+g13*px)/(g33-g23*g23)
 212    continue
      endif

      if (rightdedge.eq.0 .and. topdedge.eq.0) then                 !added d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                 !added d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                 !added d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if
 
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np
          do k=1,L 
c            coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i-1,j,k))*   !rrl
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i-1,j,k)))  !rrl 
            coefx=0.09*0.5*(saxy(i,j,k)+saxy(i-1,j,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i-1,j,k))
     .      *0.5*(xvb(i-1,j,k,nv)+xvb(i,j,k,nv))
            coefz=0.09*0.5*(saz(i,j,k)+saz(i-1,j,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i-1,j,k))
     .      *0.5*(xvb(i-1,j,k,nv)+xvb(i,j,k,nv))
            !dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(rk(i,j,k)-rk(i-1,j,k))
            hx(i,j,k)=(coefx* pxa + coefz*g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then
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
          do j=1+botdedge,mp                 !added d rrl
            do k=1,L
c           coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j-1,k))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j-1,k)))
            coefy=0.09*0.5*(saxy(i,j,k)+saxy(i,j-1,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j-1,k))
     .     *0.5*(xvb(i,j-1,k,nv)+xvb(i,j,k,nv))
            coefz=0.09*0.5*(saz(i,j,k)+saz(i,j-1,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j-1,k))
     .     *0.5*(xvb(i,j-1,k,nv)+xvb(i,j,k,nv))
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(rk(i,j,k)-rk(i,j-j3,k))
              hy(i,j,k)=(coefy*pya + g23*pza*coefz )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                 !added d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
      endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
 
compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                 !added d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
            coefz=0.09*0.5*(saz(i,j,k)+saz(i,j,k-1))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j,k-1))
     .      *0.5*(xvb(i,j,k-1,nv)+xvb(i,j,k,nv))
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coefz*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 !added d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
            coefz=0.09*0.5*(saz(np,j,k)+saz(np,j,k-1))*
     .           0.5*(sqrtk(np,j,k)+sqrtk(np,j,k-1))
     .      *0.5*(xvb(np,j,k-1,nv)+xvb(np,j,k,nv))
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coefz*gii*hza + g13*hxa
          end if
        end do
      end do
      k=1
        do j=1,mp
          do i=1,np-1*rightdedge                 !added d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
            coefz=0.09*saz(i,j,k)*
     .           sqrtk(i,j,k)
     .          *xvb(i,j,k,nv)
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=gmul(k)*c13(i,j)
            gii=gi(i,j,k)
            hxa=0.5*(hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coefz*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 !added d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
            coefz=0.09*saz(np,j,k)*
     .          sqrtk(np,j,k)
     .         *xvb(np,j,k,nv)
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=gmul(k)*c13(np,j)
             gii=gi(np,j,k)
c             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.5*(hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coefz*gii*hza + g13*hxa
          end if
        end do
        k=l+1
        do j=1,mp
          do i=1,np-1*rightdedge                 !added d rrl
c             coef=0.66*0.09*sb(i,j,k-1)*
c     .           sqrt(1.2*rkb(i,j,k-1))
            coefz=0.09*saz(i,j,k-1)*
     .           sqrtk(i,j,k-1)
     .      *xvb(i,j,k-1,nv)
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=gmul(k-1)*c13(i,j)
            gii=gi(i,j,k-1)
            hxa=0.5*(hx(i,j,k-1)+hx(i+1,j,k-1))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coefz*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 !added d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
            coefz=0.09*saz(np,j,k-1)*
     .           sqrtk(np,j,k-1)
     .      *xvb(np,j,k-1,nv)
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=gmul(k-1)*c13(np,j)
             gii=gi(np,j,k-1)
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
c             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.5*(hx(np,j,k-1)+hx1)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coefz*gii*hza + g13*hxa
          end if
        end do

! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                 !added d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                 !added d rrl
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
            do j=1,mp-topdedge                 !added d rrl
              g23=gmul(k)*c23(i,j)
              hya=0.5*(
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
        end do
      if (topdedge.eq.1) then                 !added d rrl
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
            do j=1,mp-topdedge                 !added d rrl
              g23=gmul(k-1)*c23(i,j)
              hya=0.5*(hy(i,j,k-1)+hy(i,j+j3,k-1))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
      if (topdedge.eq.1) then                 !added d rrl
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
c     do j=1,mp
c     do i=1,np
c     hz(i,j,1)=-hz(i,j,2)
c     end do
c     end do
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
          do i=1,np-rightdedge                 !added d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                 !added d rrl
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
            do j=1,mp-topdedge                 !added d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                 !added d rrl
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
c         hzLp1 =-hz(i,j,L)
c         r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do
 
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=r(i,j,k)*gi(i,j,k)
      fk(i,j,k)=fk(i,j,k)+2.*r(i,j,k)*dt
      end do
      end do
      end do
      call updated(fk,fk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      deallocate (r)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
      !deallocate (rkb)
 
      return
      end
