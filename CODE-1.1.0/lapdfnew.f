      subroutine lapdfsnew(fk,il,iu,jl,ju,lls,iv)
      use gridsetup
      use turba
      use xvo
      use xve, only:xe
      use metryic
      use msga

      Implicit None

      integer,intent(in) :: il,iu,jl,ju,lls
      integer,intent(in) :: iv ! to know with field is diffused :
                    !4 for theta
                    !5 for tka
                    !6 for tkb
                    !7 for O2

      real      fk(il:iu, jl:ju,lls) ! rhs of variable to diffuse
      real,allocatable::
     .           rk(:, :,:),  ! =xv(:,:,:iv)/xv(:,:,:nv), field to diffuve
     .           r(:, :,:),
     .          hx(:, :,:),  ! xflux
     .          hy(:, :,:),  ! yflux
     .          hz(:, :,:),  ! zflux
     .          pz(:, :,:)  ! z derivatives
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: hxa,hya,hza
      real :: hxnp1,hymp1
      real :: coef
      real :: cdiff  ! diffusion constant
      real ::aswitch ! TODO : additional source of diffusion for dz/dx in k=1 (remove?)

!      aswitch=-1.*real(iflg-6)
      aswitch=0.
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi

      allocate (rk(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l+1))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
! computation of diffusion constant
	if (iv.eq.4.or.iv.eq.7.or.(iv.eq.8.and.irhovapor.eq.1)) then ! theta, O2, rhovapor
	  cdiff=diffcst * rturbprandtl ! =0.09*2
	else ! KA,KB
	  cdiff=diffcst ! =0.09	  
	endif
      
! definition of quantity to diffuse      
      do k=1,l
          do j=1,mp
            do i=1,np
              rk(i,j,k)=xvb(i,j,k,iv)/xvb(i,j,k,nv)
            enddo
         enddo
      enddo
      call updated(rk,rk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

! compute z-derivatives at (i,j,k+-1/2) and its bc
      do k=2,l
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(rk(i,j,k)-rk(i,j,k-1))
          end do
        end do
      end do
      ! top bc:
      do j=1,mp
        do i=1,np
          pz(i,j,l+1)=dzi*(xe(i,j,l,iv)/xe(i,j,l,nv)-rk(i,j,l))
        end do
      end do

      !  bottom bc  : general case (including cyclic because rk is updated)  
      ! NB TODO : this is wrong because the BC is not at the cell center but
      ! at the boundary
      ! However when no topo, pz(i,j,1)=0 which is correct
      jllim = 1  + botdedge                 
      julim = mp + topdedge                 
      do j=jllim,julim
         if (j3.eq.1) then
            jp1 = j + 1
            jm1 = j - 1
          else	
            jp1=1
            jm1=1
         end if
         illim = 1  + leftdedge                 
         iulim = np + rightdedge                 
         do i=illim,iulim
            ip1 = i + 1
            im1 = i - 1
          
      		g13=c13(i,j)*gmul(1)
      		g23=c23(i,j)*gmul(1)
      		g33=g13**2+g23**2+gi(i,j,1)**2
      		px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
      		py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
            pz(i,j,1)=-(g13*px+g23*py)/g33
!     &          +aswitch*2.*dzi*rk(i,j,1)
          enddo
      enddo
      ! bottom bc : here we treat left and/or right edge
      if(ibcx.eq.0) then  ! this if test is not required
         illim = 1*leftdedge   + np*(1-leftdedge)                 
         iulim = np*rightdedge +  1*(1-rightdedge)                
         ! here we have between 0 step (not an edge), 1 (left or right) and 2 steps
         ! (left and right) 
         do i=illim,iulim,np-1
            jllim = 1  + botdedge                 
            julim = mp - topdedge                 
            do j=jllim,julim
               if (j3.eq.1) then
                     jp1 = j + 1
                     jm1 = j - 1
               else
                  jp1=1
                  jm1=1
               end if
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
               pz(i,j,1)=-g23*py/(g33-g13*g13)
!     &           +aswitch*2.*dzi*rk(i,j,1)
            enddo
          ! bottom bc : here we treat corners when ibcx=ibcy=0       
           if(ibcy.eq.0.and.j3.eq.1) then  
             jllim = 1*botdedge  + mp*(1-botdedge)                 
             julim = mp*topdedge +  1*(1-topdedge)                 
           ! here we have between 0 step (not an edge), 1 (bottom or top) and 2 steps
           ! (bottom and top) 
           do j=jllim,julim,mp-j3
              pz(i,j,1)=0.0
!     &                  +aswitch*2.*dzi*rk(i,j,1)
           enddo
          endif !(ibcy=0)
       enddo
      endif  !ibcx.eq.0
      ! bottom bc : here we treat bottom and/or top edge
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)                 
         julim = mp*topdedge +  1*(1-topdedge)                 
         ! here we have between 0 step (not an edge), 1 (bottom or top) and 2 steps
         ! (bottom and top) 

         do j=jllim,julim,mp-j3
            illim = 1  + leftdedge                 
            iulim = np - rightdedge                
            do i=illim,iulim
               ip1 = i + 1
               im1 = i - 1
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
               pz(i,j,1)=-(g13*px)/(g33-g23*g23)
!     &                  +aswitch*2.*dzi*rk(i,j,1)
            enddo
          enddo
      endif    
        
      if (rightdedge.eq.0 .and. topdedge.eq.0) then                
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then                 
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then                 
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if
! end computation of pz
 
! compute x-flux at (i-1/2,j,k) hx
      do j=1,mp
        do i=1+leftdedge,np
          do k=1,l 
            gii=0.5*(gi(i,j,k)+gi(i-1,j,k))
            !coef=coefxy
            coef=cdiff*0.5*(saxy(i,j,k)+saxy(i-1,j,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i-1,j,k))
     .      *0.5*(xvb(i-1,j,k,nv)+xvb(i,j,k,nv))/gii
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            pxa=dxi*(rk(i,j,k)-rk(i-1,j,k))
            hx(i,j,k)=coef * (pxa + g13 * pza)
          end do
        end do
      end do
	 ! create boundary conditions at i=1 when required (leftdege.eq.1)
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then
      do j=1,mp
        do k=1,l
          !hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k) !ibcx=0 here (leftDegde)
          hx(1,j,k) = - hx(2,j,k)
        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

! compute y-flux at (i,j-1/2,k) hy
      if (j3.eq.1) then
        do i=1,np
         do j=1+botdedge,mp                 
           do k=1,l
            gii=0.5*(gi(i,j,k)+gi(i,j-1,k))
            !coef=coefxy
            coef=cdiff*0.5*(saxy(i,j,k)+saxy(i,j-1,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j-1,k))
     .     *0.5*(xvb(i,j-1,k,nv)+xvb(i,j,k,nv))/gii
            
            g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
            pza=0.25*(pz(i,j-1,k)+pz(i,j,k)+pz(i,j-1,k+1)+pz(i,j,k+1) )
            pya=dyi*(rk(i,j,k)-rk(i,j-1,k))
            hy(i,j,k)=coef * (pya + g23*pza)
           end do
         end do
        end do
        !create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                 
        do k=1,l
          do i=1,np
            !hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
            hy(i,1,k)= -hy(i,2,k)
          end do
        end do
        endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif ! end j3.eq.1
      
 
! compute z-flux at (i,j,k-1/2) hz 
!i)(include the dh/dx and dh/dz terms)
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                 
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            !coef=coefz
            coef=cdiff*0.5*(saz(i,j,k)+saz(i,j,k-1))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j,k-1))
     .      *0.5*(xvb(i,j,k-1,nv)+xvb(i,j,k,nv))/gii
       
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii * pz(i,j,k)
            hz(i,j,k)= coef * gii * hza + g13 * hxa
          end do
          ! bc on right edge when ibcx=0
          if (rightdedge.eq.1) then                 
            gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
            coef=cdiff*0.5*(saz(np,j,k)+saz(np,j,k-1))*
     .           0.5*(sqrtk(np,j,k)+sqrtk(np,j,k-1))
     .      *0.5*(xvb(np,j,k-1,nv)+xvb(np,j,k,nv))/gii
     
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             hxa=0.0 !(ibcx=0)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef * gii * hza + g13 * hxa
          end if
        end do
      end do
      ! compute zflux at k=1-1/2
      ! FP simplified a lot this section => hz(i,j,1)=0 by definition of the BC!
      k=1
        do j=1,mp
          do i=1,np-1*rightdedge
            hz(i,j,k)=0.0                 
!            gii=gi(i,j,k)
!             coef, g13 are approximated by coef(i,j,k) and g13(i,j,k)
!            coef=cdiff*saz(i,j,k)*
!     .           sqrtk(i,j,k)
!     .          *xvb(i,j,k,nv)/gii
!            g13=gmul(k)*c13(i,j)
!            hxa=0.5*(hx(i,j,k)+hx(i+1,j,k))
!            hza=gii*pz(i,j,k)
!            hz(i,j,k)= coef * gii * hza + g13 * hxa
          end do
          if (rightdedge.eq.1) then                 
           !  b.c. on right  i=n+1 when ibcx=0
!            gii=gi(np,j,k)
!            coef=cdiff*saz(np,j,k)*
!     .          sqrtk(np,j,k) *xvb(np,j,k,nv)/gii
            
!            g13=gmul(k)*c13(np,j)
!            hxa=0.0 !(ibcx=0)
!            hza=gii * pz(np,j,k)
!            hz(np,j,k)= coef*gii*hza + g13*hxa
            hz(np,j,k)=0.0
          end if
        end do
        ! compute zflux at k=l+1/2
        k=l+1
        do j=1,mp
          do i=1,np-1*rightdedge                 
            gii=gi(i,j,k-1)
            coef=cdiff*saz(i,j,k-1)*
     .           sqrtk(i,j,k-1)
     .      *xvb(i,j,k-1,nv)/gii
            g13=gmul(k-1)*c13(i,j)
            hxa=0.5*(hx(i,j,k-1)+hx(i+1,j,k-1))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 
            !  b.c. on right  i=n+1 when ibcx=0
            gii=gi(np,j,k-1)
            coef=cdiff*saz(np,j,k-1)*
     .           sqrtk(np,j,k-1)
     .      *xvb(np,j,k-1,nv)/gii
             g13=gmul(k-1)*c13(np,j)
             hxa=0.0 !(ibcx=0)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do

! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                 
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
! FP commented the following lines that was doing nothing (bcy=0)        
!      if (topdedge.eq.1) then !ibcy=0 => hya=0...                 
!        do k=2,L
!          do i=1,np
!            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
!            hya=0.0
!            hz(i,mp,k)=hz(i,mp,k) + g23*hya
!          end do
!        end do
!      end if
        k=1
          do i=1,np
            do j=1,mp-topdedge                 
              g23=gmul(k)*c23(i,j)
              hya=0.5*(hy(i,j,k)+hy(i,j+1,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
        end do
! FP : next lines not required cause ibcy=0        
!      if (topdedge.eq.1) then ! ibcy=0                
!          do i=1,np
!corporate b.c. for hy on j=m+1
!            g23=gmul(k)*c23(i,mp)
!            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
!            hya=0.5*(hy(i,mp,k)+hy2)
!            hz(i,mp,k)=hz(i,mp,k) + g23*hya
!          end do
!      end if
       k=l+1
          do i=1,np
            do j=1,mp-topdedge                 
              g23=gmul(k-1)*c23(i,j)
              hya=0.5*(hy(i,j,k-1)+hy(i,j+j3,k-1))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
!      if (topdedge.eq.1) then                 
!          do i=1,np
!corporate b.c. for hy on j=m+1
!            g23=gmul(k-1)*c23(i,mp)
!            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
!c            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
!            hya=0.5*(hy(i,mp,k-1)+hy1)
!            hz(i,mp,k)=hz(i,mp,k) + g23*hya
!          end do
!      end if

      endif  ! if j3=1

      !create boundary conditions at k=1; for k=l see divergence below
      ! surface fluxes:
!     do j=1,mp
!     do i=1,np
!     hz(i,j,1)=-hz(i,j,2)
!     end do
!     end do
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

! compute Laplacian term by term
      do k=1,l
        do j=1,mp
          do i=1,np
             r(i,j,k)=0.
          end do
         end do
      end do
! compute d/dx(dh/dx)
      do k=1,l
        do j=1,mp
          do i=1,np-rightdedge                 
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
      ! create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                 
         do k=1,l
            do j=1,mp
               hxnp1 = -hx(np,j,k)
               r(np,j,k) = dxi*(hxnp1-hx(np,j,k))
            end do
         end do
      end if

! for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,l
          do i=1,np
            do j=1,mp-topdedge                 
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
        ! create boundary conditions at j=m+1
        if (topdedge.eq.1) then                 
          do k=1,l
             do i=1,np

                hymp1= - hy(i,mp,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
        end if
      endif  !j3.eq.1

! compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,l
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
! corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
! TODO : why is this bc commented?
!         hzLp1 =-hz(i,j,L)
!         r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
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

      deallocate (rk)
      deallocate (r)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      
      return
      end
