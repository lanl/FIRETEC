!**********************************************************************!
!     this routine computes the advection on the large time steps and
!     the nv variables. It can use either donorcell or fct
!**********************************************************************!
module advection_functions
  use gridlist_variables, only : wp

  contains

  real(wp) function pp(y)
    Implicit None
    real(wp),intent(in) :: y

    pp = max(0._wp,y)

  end function pp
  !*************************************
  real(wp) function pn(y)
    use gridlist_variables, only : wp
    Implicit None
    real(wp),intent(in) :: y

    pn = -min(0._wp,y)

  end function pn
  !*************************************
  real(wp) function donor(y1,y2,a10)
    Implicit None
    real(wp),intent(in) :: y1,y2,a10

    donor = pp(a10)*y1-pn(a10)*y2

  end function donor
  !*************************************
  real(wp) function vdiv1(a1,a2,a3,r)
    use gridlist_variables, only : wp
    Implicit None
    real(wp),intent(in) :: a1,a2,a3,r

    vdiv1 = 0.25_wp*a2*(a3-a1)/r

  end function vdiv1
  !*************************************
  real(wp) function vdiv2(aa,b1,b2,b3,b4,r)
    Implicit None
    real(wp),intent(in) :: aa,b1,b2,b3,b4,r

    vdiv2 = 0.25_wp*aa*(b1+b2-b3-b4)/r

  end function vdiv2
  !*************************************
  real(wp) function rat2(z1,z2)
    use gridlist_variables, only : wp
    Implicit None
    real(wp),intent(in) :: z1,z2

    rat2 = 0.5_wp*(z2-z1)

  end function rat2
  !*************************************
  real(wp) function rat4(z0,z1,z2,z3)
    Implicit None
    real(wp),intent(in) :: z0,z1,z2,z3

    rat4 = (z3+z2-z1-z0)*0.25_wp

  end function rat4
  !*************************************
  real(wp) function vdyf(x1,x2,aa,r)
    Implicit None
    real(wp),intent(in) :: x1,x2,aa,r
    !real(wp),external :: rat2

    vdyf = (abs(aa)-aa**2/r)*rat2(x1,x2)

  end function vdyf
  !*************************************
  real(wp) function vcorr(aa,b,y0,y1,y2,y3,r)
    Implicit None
    real(wp),intent(in) :: aa,b,y0,y1,y2,y3,r
    !real(wp),external :: rat4

    vcorr = -0.125_wp*aa*b/r*rat4(y0,y1,y2,y3)

  end function vcorr
  !*************************************

end module advection_functions

!-----------------------------------------------------------------------
! from u,v,o divided by gi and cell centered, this routine set 
! u1,u2,u3, the contravariant face coordinates*dt/dx*1/gi required 
!-----------------------------------------------------------------------

subroutine setAdvectiveVelocities(u,v,o,il,iu,jl,ju,lls,dta)
  use gridlist_variables, only : l,ih,wp,ibctopopen
  use forcings, only : u1,u2,u3
  use gridsetup, only : np,mp,dxi,dyi,dzi
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real(wp),intent(in) :: dta
  real(wp),intent(inout) :: u(il:iu,jl:ju,lls),v(il:iu,jl:ju,lls)
  real(wp),intent(inout) :: o(il:iu,jl:ju,lls)

  integer :: i,j,k
  real(wp) :: gcx,gcy,gcz

  ! Executable Code
  gcx=dta*dxi
  gcy=dta*dyi
  gcz=dta*dzi

  ! compute contravariant cell face in the x direction
  call update(u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,2,0,1)
  do k=1,l
    do j=1,mp
      do i=1,np+1
        u1(i,j,k)=0.5_wp*(u(i,j,k)+u(i-1,j,k))*gcx           
      enddo
    enddo
  enddo

  call update(v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,3,0,2)
  do k=1,l
    do j=1,mp+1
      do i=1,np
        u2(i,j,k)=0.5_wp*(v(i,j,k)+v(i,j-1,k))*gcy
      enddo
    enddo
  enddo

  ! compute contravariant cell face in the z direction
  do k=2,l
    do j=1,mp
      do i=1,np
        u3(i,j,k)=0.5_wp*(o(i,j,k)+o(i,j,k-1))*gcz
      enddo
    enddo
  enddo
  u3(:,:,1)=0._wp
  if(ibctopopen.eq.0)then
    u3(:,:,l+1)=0._wp
  else
    do j=1,mp
      do i=1,np
        u3(i,j,l+1)=o(i,j,l)*gcz
      enddo
    enddo
  endif

end subroutine setAdvectiveVelocities

!-----------------------------------------------------------------------
! mpdata is the 2nd order advection scheme of Smolarkiewicz (1984)
!-----------------------------------------------------------------------
subroutine mpdata(x,xrho,h,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : l,ih,iord,nonos,idiv,nfct,nonosold, &
    wp,ibctopopen
  use forcings, only : u1,u2,u3
  use gridsetup, only : np,mp
  use msga_variables, only : leftedge,rightedge,botedge,topedge, &
    npos,mpos
  use xvall, only : xe
  use constants, only : tolerance
  use advection_functions
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(wp),intent(in) :: h(il:iu,jl:ju,lls)
  real(wp),intent(inout) :: x(il:iu,jl:ju,lls,nvp)
  real(wp),intent(inout) :: xrho(il:iu,jl:ju,0:lls,nvp-1)
  
  integer :: kv,i,j,k,ip,im,jp,jm,kp,km,itrfct,iuv
  integer :: itr,illim,iulim,jllim,julim,ia,ja
  real(wp) :: rhoin,rhoout,ain,aout,rmxuse,rmnuse,a1p,a2p,a1n,a2n
  real(wp) :: tmpp,tmpn,v1d,v2d,v3d,c1,c2
  real(wp),allocatable:: v1mp(:,:,:,:),v2mp(:,:,:,:),v3mp(:,:,:,:)
  real(wp),allocatable:: f1(:,:,:,:),f2(:,:,:,:),f3(:,:,:,:)
  real(wp),allocatable:: f1o(:,:,:,:),f2o(:,:,:,:),f3o(:,:,:,:)
  real(wp),allocatable:: cpmp(:,:,:,:),cnmp(:,:,:,:)
  real(wp),allocatable:: mxo(:,:,:,:),mno(:,:,:,:)
  real(wp),allocatable:: mx(:,:,:,:),mn(:,:,:,:)

  ! Executable Code
  allocate (v1mp(1:np+1,1:mp  ,l  ,nvp))
  allocate (v2mp(1:np  ,1:mp+1,l  ,nvp))
  allocate (v3mp(1:np  ,1:mp  ,l+1,nvp))
  allocate (f1(il:iu+1,jl:ju  ,lls  ,nvp)  )
  allocate (f2(il:iu  ,jl:ju+1,lls  ,nvp)  )
  allocate (f3(il:iu  ,jl:ju  ,lls+1,nvp)  )
  allocate (f1o(il:iu+1,jl:ju  ,lls  ,nvp) )
  allocate (f2o(il:iu  ,jl:ju+1,lls  ,nvp) )
  allocate (f3o(il:iu  ,jl:ju  ,lls+1,nvp) )
  allocate (cpmp(il:iu,jl:ju,lls,nvp)      )
  allocate (cnmp(il:iu,jl:ju,lls,nvp)      )
  allocate (mx(il:iu,jl:ju,lls,nvp-1)      )
  allocate (mn(il:iu,jl:ju,lls,nvp-1)      )
  allocate (mxo(il:iu,jl:ju,lls,nvp)       )
  allocate (mno(il:iu,jl:ju,lls,nvp)       )

  do kv=1,nvp
    v1mp(1:np+1,1:mp  ,1:l  ,kv)=u1
    v2mp(1:np  ,1:mp+1,1:l  ,kv)=u2
    v3mp(1:np  ,1:mp  ,1:l+1,kv)=u3
  enddo

  if(nonos.eq.1) then
    do kv=1,nvp-1
      do k=1,l
        km=max(k-1,1)
        kp=min(k+1,l)
        do j=1,mp
          jm=j-1
          jp=j+1
          do i=1,np
            im=i-1
            ip=i+1
            mx(i,j,k,kv)=max(xrho(i,j,k,kv),xrho(im,j,k,kv), &
              xrho(ip,j,k,kv),xrho(i,jm,k,kv),xrho(i,jp,k,kv), &
              xrho(i,j,kp,kv),xrho(i,j,km,kv))
            mn(i,j,k,kv)=min(xrho(i,j,k,kv),xrho(im,j,k,kv), &
              xrho(ip,j,k,kv),xrho(i,jm,k,kv),xrho(i,jp,k,kv), &
              xrho(i,j,kp,kv),xrho(i,j,km,kv))
          enddo ! i=1:np
        enddo ! j=1:mp
      enddo ! k=1:l(l)
    enddo ! kv=1:nv-1 
  endif ! nonos==1

  if(nonosold.eq.1) then
    do kv=1,nvp
      do k=1,l
        km=max(k-1,1)
        kp=min(k+1,l)
        do j=1,mp
          jm=j-1
          jp=j+1
          do i=1,np
            im=i-1
            ip=i+1
            mxo(i,j,k,kv)=max(x(im,j,k,kv),x(i,j,k,kv),x(ip,j,k,kv), &
              x(i,jm,k,kv),x(i,jp,k,kv),x(i,j,kp,kv),x(i,j,km,kv))
            mno(i,j,k,kv)=min(x(im,j,k,kv),x(i,j,k,kv),x(ip,j,k,kv), &
              x(i,jm,k,kv),x(i,jp,k,kv),x(i,j,kp,kv),x(i,j,km,kv))
          enddo ! i=1:np
        enddo ! j=1:mp
      enddo ! k=1:nm3(l)
    enddo ! kv=1:nv
  endif ! nonosold==1

  c1=1._wp
  c2=0._wp

  do itr=1,iord ! loop to LINE 1402  
    do kv=1,nvp
      illim=1+leftedge                       
      do k=1,l
        do j=1,mp
          do i=illim,np
            f1(i,j,k,kv)=donor(c1*x(i-1,j,k,kv)+c2, &
              c1*x(i,j,k,kv)+c2,v1mp(i,j,k,kv))
          enddo
        enddo
      enddo
      if (leftedge.eq.1) then                       
        do k=1,l
          do j=1,mp
            f1(1,j,k,kv)=donor(c1*xe(1,j,k,kv)+c2, &
               c1*x(1,j,k,kv)+c2, &
               v1mp(1,j,k,kv))
          enddo
        enddo
      endif
      if (rightedge.eq.1) then                       
        do k=1,l
          do j=1,mp
            f1(np+1,j,k,kv)= donor(c1*x(np,j,k,kv)+c2, &
                            c1*xe(np,j,k,kv)+c2, &
                            v1mp(np+1,j,k,kv))
          enddo
        enddo
      endif

      jllim=1+botedge                       
      do k=1,l
        do j=jllim,mp
          do i=1,np
            f2(i,j,k,kv)=donor(c1*x(i,j-1,k,kv)+c2, &
                         c1*x(i,j,k,kv)+c2, &
                         v2mp(i,j,k,kv))
          enddo
        enddo
      enddo
      if (botedge.eq.1) then                       
        do k=1,l
          do i=1,np
            f2(i,1,k,kv)=donor(c1*xe(i,1,k,kv)+c2, &
                         c1*x(i,1,k,kv)+c2, &
                         v2mp(i,1,k,kv))
          enddo
        enddo
      endif
      if (topedge.eq.1) then                       
        do k=1,l
          do i=1,np
             f2(i,mp+1,k,kv)= donor(c1*x(i,mp,k,kv)+c2, &
                             c1*xe(i,mp,k,kv)+c2, &
                             v2mp(i,mp+1,k,kv))
          enddo
        enddo
      endif

      do k=2,l
        do j=1,mp
          do i=1,np
            f3(i,j,k,kv)=donor(c1*x(i,j,k-1,kv)+c2,c1*x(i,j,k,kv)+c2, &
                              v3mp(i,j,k,kv))
          enddo
        enddo
      enddo
      f3(1:np,1:mp,  1,kv)=0._wp
      if(ibctopopen.eq.0)then
        f3(1:np,1:mp,l+1,kv)=0._wp
      else
        do j=1,mp
          do i=1,np
            f3(i,j,l+1,kv)=donor(c1*x(i,j,l,kv)+c2, &
              c1*x(i,j,l,kv)+c2,max(0.0_wp,v3mp(i,j,l+1,kv)))
          enddo
        enddo
      endif

      if (rightedge.eq.0) then                       
        call update(f1(:,:,:,kv),np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,1,1) 
      else
        call update(f1(:,:,:,kv),np+1,mp,l, &
                    1-ih,np+ih+1,1-ih,mp+ih,2,0,1) 
      endif
      if(topedge.eq.0) then                       
        call update(f2(:,:,:,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,2,2) 
      else
        call update(f2(:,:,:,kv),np,mp+1,l, &
                    1-ih,np+ih,1-ih,mp+ih+1,3,0,2) 
      endif

      do k=1,l
        do j=1,mp
          do i=1,np
            x(i,j,k,kv)=x(i,j,k,kv)-(f1(i+1,j,k,kv)-f1(i,j,k,kv) &
                                    +f2(i,j+1,k,kv)-f2(i,j,k,kv) &
                                    +f3(i,j,k+1,kv)-f3(i,j,k,kv)) &
                      /h(i,j,k)
          enddo
        enddo
      enddo
      call loop_update_flags(kv,iuv)
      call update(x(:,:,:,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0,iuv) 
    enddo ! kv=1:nv

    if(itr.ne.iord)then  !! Skip if iord.eq.1: start of mpdata when iord=2: escaping to LINE 1403
      c1=0._wp
      c2=1._wp
      do kv=1,nvp
        f1(1:np+1,1:mp,1:l,kv)=v1mp(:,:,:,kv)
        f2(1:np,1:mp+1,1:l,kv)=v2mp(:,:,:,kv)
        f3(1:np,1:mp,1:l+1,kv)=v3mp(:,:,:,kv)
        v1mp(:,:,:,kv)=0._wp
        v2mp(:,:,:,kv)=0._wp
        v3mp(:,:,:,kv)=0._wp
        
        if (rightedge.eq.0) then                       
          call update(f1(:,:,:,kv),np,mp,l, &
           1-ih,np+ih+1,1-ih,mp+ih,1,1,1) 
        else
          call update(f1(:,:,:,kv),np+1,mp,l, &
           1-ih,np+ih+1,1-ih,mp+ih,1,0,1) 
        endif
        if (topedge.eq.0) then 
          call update(f2(:,:,:,kv),np,mp,l, &
                      1-ih,np+ih,1-ih,mp+ih+1,1,2,2) 
        else
          call update(f2(:,:,:,kv),np,mp+1,l, &
                      1-ih,np+ih,1-ih,mp+ih+1,1,0,2)  
        endif
        call update(f3(:,:,:,kv),np,mp,l+1, &
                      1-ih,np+ih,1-ih,mp+ih,1,0,0) 

        !!! compute antidiffusive velocities in x direction
        illim = 1 + leftedge                       
        iulim = np
        jllim = 1  + botedge                       
        julim = mp - topedge                       
        do k=2,l-1
          do j=jllim,julim
            do i=illim,iulim
              v1mp(i,j,k,kv)=vdyf(  &
                x(i-1,j,k,kv),x(i,j,k,kv),f1(i,j,k,kv), &
               0.5_wp*(h(i-1,j,k)+h(i,j,k))) &
               +vcorr(f1(i,j,k,kv),f2(i-1,j,k,kv)+f2(i-1,j+1,k,kv) &
               +f2(i,j+1,k,kv)+f2(i,j,k,kv), &
               x(i-1,j-1,k,kv),x(i,j-1,k,kv), &
               x(i-1,j+1,k,kv),x(i,j+1,k,kv), &
               0.5_wp*(h(i-1,j,k)+h(i,j,k))) &
               +vcorr(f1(i,j,k,kv), &
               f3(i-1,j,k,kv)+f3(i-1,j,k+1,kv) &
               +f3(i,j,k+1,kv)+f3(i,j,k,kv), &
               x(i-1,j,k-1,kv),x(i,j,k-1,kv), &
               x(i-1,j,k+1,kv),x(i,j,k+1,kv), &
               0.5_wp*(h(i-1,j,k)+h(i,j,k)))
            enddo
          enddo
        enddo
        if(idiv.eq.1) then
          do k=2,l-1
            do j=jllim,julim
              do i=illim,iulim
                v1d=-vdiv1(f1(i-1,j,k,kv),f1(i,j,k,kv),f1(i+1,j,k,kv), &
                 0.5_wp*(h(i-1,j,k)+h(i,j,k))) &
                 -vdiv2(f1(i,j,k,kv),f2(i-1,j+1,k,kv), &
                 f2(i,j+1,k,kv),f2(i-1,j,k,kv),&
                 f2(i,j,k,kv),0.5_wp*(h(i-1,j,k)+h(i,j,k))) &
                 -vdiv2(f1(i,j,k,kv),f3(i-1,j,k+1,kv), &
                 f3(i,j,k+1,kv),f3(i-1,j,k,kv), &
                 f3(i,j,k,kv),0.5_wp*(h(i-1,j,k)+h(i,j,k)))
                v1mp(i,j,k,kv)=v1mp(i,j,k,kv)+ &
                  (pp(v1d)*x(i-1,j,k,kv)-pn(v1d)*x(i,j,k,kv))
              enddo
            enddo
          enddo
        endif ! idiv.eq.1

        !!! compute antidiffusive velocities in y direction
        illim = 1  + leftedge                       
        iulim = np - rightedge                       
        jllim = 1  + botedge                       
        julim = mp
        do k=2,l-1
          do j=jllim,julim
            do i=illim,iulim
              v2mp(i,j,k,kv)=vdyf( &
                x(i,j-1,k,kv),x(i,j,k,kv),f2(i,j,k,kv), &
               0.5_wp*(h(i,j-1,k)+h(i,j,k))) &
               +vcorr(f2(i,j,k,kv), &
               f1(i,j-1,k,kv)+f1(i,j,k,kv) &
               +f1(i+1,j,k,kv)+f1(i+1,j-1,k,kv), &
               x(i-1,j-1,k,kv),x(i-1,j,k,kv), &
               x(i+1,j-1,k,kv),x(i+1,j,k,kv), &
               0.5_wp*(h(i,j-1,k)+h(i,j,k))) &
               +vcorr(f2(i,j,k,kv), &
               f3(i,j-1,k,kv)+f3(i,j,k,kv) &
               +f3(i,j,k+1,kv)+f3(i,j-1,k+1,kv), &
               x(i,j-1,k-1,kv),x(i,j,k-1,kv), &
               x(i,j-1,k+1,kv),x(i,j,k+1,kv), &
               0.5_wp*(h(i,j-1,k)+h(i,j,k)))
            enddo
          enddo
        enddo
        if(idiv.eq.1) then
          do k=2,l-1
            do j=jllim,julim
              do i=illim,iulim
                v2d=-vdiv1(f2(i,j-1,k,kv),       &
                 f2(i,j,k,kv),f2(i,j+1,k,kv),   &
                 0.5_wp*(h(i,j-1,k)+h(i,j,k)))      &
                 -vdiv2(f2(i,j,k,kv),           &
                 f1(i+1,j-1,k,kv),f1(i+1,j,k,kv), &
                 f1(i,j-1,k,kv),f1(i,j,k,kv),   &
                 0.5_wp*(h(i,j-1,k)+h(i,j,k)))      &
                 -vdiv2(f2(i,j,k,kv),f3(i,j-1,k+1,kv), &
                 f3(i,j,k+1,kv),f3(i,j-1,k,kv), &
                 f3(i,j,k,kv),    0.5_wp*(h(i,j-1,k)+h(i,j,k)))
               v2mp(i,j,k,kv)=v2mp(i,j,k,kv)+(pp(v2d)*x(i,j-1,k,kv)- &
                 pn(v2d)*x(i,j,k,kv))
              enddo
            enddo
          enddo
        endif ! idiv.eq.1 

        !!! compute antidiffusive velocities in z direction
        illim = 1  + leftedge                       
        iulim = np - rightedge                       
        jllim = 1  + botedge                       
        julim = mp - topedge                       
        do k=2,l
          do j=jllim,julim
            do i=illim,iulim
              v3mp(i,j,k,kv)=vdyf( &
                x(i,j,k-1,kv),x(i,j,k,kv),f3(i,j,k,kv), &
                0.5_wp*(h(i,j,k-1)+h(i,j,k)))               &
                +vcorr(f3(i,j,k,kv),                    &
                f1(i,j,k-1,kv)+f1(i,j,k,kv)             &
                +f1(i+1,j,k,kv)+f1(i+1,j,k-1,kv),       &
                x(i-1,j,k-1,kv),x(i-1,j,k,kv),          &
                x(i+1,j,k-1,kv),x(i+1,j,k,kv),          &
                0.5_wp*(h(i,j,k-1)+h(i,j,k)))               &
                +vcorr(f3(i,j,k,kv),                    &
                f2(i,j,k-1,kv)+f2(i,j+1,k-1,kv)         &
                +f2(i,j+1,k,kv)+f2(i,j,k,kv),           &
                x(i,j-1,k-1,kv),x(i,j-1,k,kv),          &
                x(i,j+1,k-1,kv),x(i,j+1,k,kv),          &
                0.5_wp*(h(i,j,k-1)+h(i,j,k)))
            enddo
          enddo
        enddo
        if(idiv.eq.1) then
          do k=2,l
            do j=jllim,julim
              do i=illim,iulim
                ia=(npos-1)*np + i
                ja=(mpos-1)*mp + j

                v3d=-vdiv1(f3(i,j,k-1,kv), &
                 f3(i,j,k,kv),f3(i,j,k+1,kv), &
                 0.5_wp*(h(i,j,k-1)+h(i,j,k))) &
                 -vdiv2(f3(i,j,k,kv), &
                 f1(i+1,j,k-1,kv),f1(i+1,j,k,kv), &
                 f1(i,j,k-1,kv), &
                 f1(i,j,k,kv),  0.5_wp*(h(i,j,k-1)+h(i,j,k))) &
                 -vdiv2(f3(i,j,k,kv), &
                 f2(i,j+1,k-1,kv),f2(i,j+1,k,kv), &
                 f2(i,j,k-1,kv), &
                 f2(i,j,k,kv),    0.5_wp*(h(i,j,k-1)+h(i,j,k)))
                v3mp(i,j,k,kv)=v3mp(i,j,k,kv) &
                 +(pp(v3d)*x(i,j,k-1,kv)-pn(v3d)*x(i,j,k,kv))
              enddo
            enddo
          enddo
        endif ! idiv==1
      enddo ! kv=1:nv

      if(nonosold.eq.1) then  ! non-osscilatory option
        do kv=1,nvp
          do k=1,l
            km=max0(k-1,1)
            kp=min0(k+1,l)
            do j=1,mp
              jm=j-1
              jp=j+1
              do i=1,np
                im=i-1
                ip=i+1
                mxo(i,j,k,kv)=max(x(im,j,k,kv),x(i,j,k,kv), &
                  x(ip,j,k,kv),mxo(i,j,k,kv),               &
                  x(i,jm,k,kv),x(i,jp,k,kv),                &
                  x(i,j,kp,kv),x(i,j,km,kv))
                mno(i,j,k,kv)=min(x(im,j,k,kv),x(i,j,k,kv), &
                  x(ip,j,k,kv),mno(i,j,k,kv), &
                  x(i,jm,k,kv),x(i,jp,k,kv), &
                  x(i,j,kp,kv),x(i,j,km,kv))
              enddo
            enddo
          enddo

          illim = 1
          iulim = np+rightedge                       
          jllim = 1
          julim = mp
          do k=1,l
            do j=jllim,julim
              do i=illim,iulim
                f1(i,j,k,kv)=donor(c2,c2,v1mp(i,j,k,kv))
              enddo
            enddo
          enddo

          illim = 1
          iulim = np
          jllim = 1
          julim = mp+topedge                       
          do k=1,l
            do j=jllim,julim
              do i=illim,iulim
                f2(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))
              enddo
            enddo
          enddo

          do k=1,l+1
            do j=1,mp
              do i=1,np
                f3(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))
              enddo
            enddo
          enddo

          if (rightedge.eq.0) then                       
            call update(f1(:,:,:,kv), &
              np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,1,1) 
          else
            call update(f1(:,:,:,kv), &
              np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0,1) 
          endif

          if (topedge.eq.0) then                       
            call update(f2(:,:,:,kv), &
              np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,2,2) 
          else
            call update(f2(:,:,:,kv), &
              np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,3,0,2) 
          endif

          do k=1,l
            do j=1,mp
              do i=1,np
                cpmp(i,j,k,kv)=(mxo(i,j,k,kv)-x(i,j,k,kv))*h(i,j,k)/ &
                  (pn(f1(i+1,j,k,kv))+pp(f1(i,j,k,kv))              &
                  +pn(f2(i,j+1,k,kv))+pp(f2(i,j,k,kv))               &
                  +pn(f3(i,j,k+1,kv))+pp(f3(i,j,k,kv))+tolerance)

                cnmp(i,j,k,kv)=(x(i,j,k,kv)-mno(i,j,k,kv))*h(i,j,k)/ &
                  (pp(f1(i+1,j,k,kv))+pn(f1(i,j,k,kv)) &
                  +pp(f2(i,j+1,k,kv))+pn(f2(i,j,k,kv)) &
                  +pp(f3(i,j,k+1,kv))+pn(f3(i,j,k,kv))+tolerance)
              enddo
            enddo
          enddo

          call loop_update_flags(kv,iuv)
          call update(cpmp(:,:,:,kv), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0,iuv) 
          call update(cnmp(:,:,:,kv), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0,iuv) 

          illim = 1+leftedge                       
          do k=1,l
            do j=1,mp
              do i=illim,np
                v1mp(i,j,k,kv)= pp(v1mp(i,j,k,kv))&
                 *min(1._wp,cpmp(i,j,k,kv),cnmp(i-1,j,k,kv))&
                 -pn(v1mp(i,j,k,kv))&
                 *min(1._wp,cpmp(i-1,j,k,kv),cnmp(i,j,k,kv))
              enddo
            enddo
          enddo

          jllim = 1+botedge                       
          do k=1,l
            do j=jllim,mp
              do i=1,np
                v2mp(i,j,k,kv)= pp(v2mp(i,j,k,kv))&
                 *min(1._wp,cpmp(i,j,k,kv),cnmp(i,j-1,k,kv))&
                 -pn(v2mp(i,j,k,kv))&
                 *min(1._wp,cpmp(i,j-1,k,kv),cnmp(i,j,k,kv))
              enddo
            enddo
          enddo

          do k=2,l
            do j=1,mp
              do i=1,np
                v3mp(i,j,k,kv)= pp(v3mp(i,j,k,kv)) &
                 *min(1._wp,cpmp(i,j,k,kv),cnmp(i,j,k-1,kv)) &
                 -pn(v3mp(i,j,k,kv)) &
                 *min(1._wp,cpmp(i,j,k-1,kv),cnmp(i,j,k,kv))
              enddo
            enddo
          enddo
        enddo ! kv=1:nv
      endif  !! nonosold.eq.1 

      if(nonos.eq.1) then
        do itrfct=1,nfct !!nfct=1 is default
          if(itrfct.eq.1) then
            do kv=1,nvp
              illim = 1
              iulim = np+rightedge                       
              jllim = 1
              julim = mp
              do k=1,l
                do j=jllim,julim
                  do i=illim,iulim
                    f1(i,j,k,kv)=donor(c2,c2,v1mp(i,j,k,kv))
                  enddo
                enddo
              enddo

              illim = 1
              iulim = np
              jllim = 1
              julim = mp + topedge                       
              do k=1,l
                do j=jllim,julim
                  do i=illim,iulim
                    f2(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))
                  enddo
                enddo
              enddo

              do k=1,l+1
                do j=1,mp
                  do i=1,np
                    f3(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))
                  enddo
                enddo
              enddo

              if (rightedge.eq.0) then                       
                call update(f1(:,:,:,kv), &
                  np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,1,1) 
              else
                call update(f1(:,:,:,kv), &
                  np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0,1) 
              endif
              if (topedge.eq.0) then                       
                call update(f2(:,:,:,kv), &
                  np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,2,2) 
              else
                call update(f2(:,:,:,kv), &
                  np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,3,0,2) 
              endif
            enddo !kv=1:nv
          endif ! itrfct.eq.1

          !!! correction coefficients for variables nv-1
          f3(1:np,1:mp,1,:)=-f3(1:np,1:mp,2,:)
          if(ibctopopen.eq.0) f3(1:np,1:mp,l+1,:)=-f3(1:np,1:mp,l,:)
          do kv=1,nvp-1
            do k=1,l
              do j=1,mp
                do i=1,np
                  if(abs(mx(i,j,k,kv)).lt.tolerance) mx(i,j,k,kv)=0._wp
                  if(abs(mn(i,j,k,kv)).lt.tolerance) mn(i,j,k,kv)=0._wp
                  rhoin=(pn(f1(i+1,j,k,nvp))+pp(f1(i,j,k,nvp)) &
                    +pn(f2(i,j+1,k,nvp))+pp(f2(i,j,k,nvp)) &
                    +pn(f3(i,j,k+1,nvp))+pp(f3(i,j,k,nvp)))
                  rhoout=-(pp(f1(i+1,j,k,nvp))+pn(f1(i,j,k,nvp)) &
                    +pp(f2(i,j+1,k,nvp))+pn(f2(i,j,k,nvp)) &
                    +pp(f3(i,j,k+1,nvp))+pn(f3(i,j,k,nvp)))
                  ain=(pn(f1(i+1,j,k,kv))+pp(f1(i,j,k,kv)) &
                    +pn(f2(i,j+1,k,kv))+pp(f2(i,j,k,kv)) &
                    +pn(f3(i,j,k+1,kv))+pp(f3(i,j,k,kv)))
                  aout=-(pp(f1(i+1,j,k,kv))+pn(f1(i,j,k,kv)) &
                    +pp(f2(i,j+1,k,kv))+pn(f2(i,j,k,kv)) &
                    +pp(f3(i,j,k+1,kv))+pn(f3(i,j,k,kv)))
                  cpmp(i,j,k,kv)= &
                    pp(mx(i,j,k,kv)*x(i,j,k,nvp)-x(i,j,k,kv))*h(i,j,k) &
                    /(ain-pp(mx(i,j,k,kv))*rhoout &
                    +pn(mx(i,j,k,kv))*rhoin+tolerance)
                  cnmp(i,j,k,kv)= &
                    pp(x(i,j,k,kv)-mn(i,j,k,kv)*x(i,j,k,nvp))*h(i,j,k) &
                    /(-aout+pp(mn(i,j,k,kv))*rhoin &
                    -pn(mn(i,j,k,kv))*rhoout+tolerance)
                enddo ! i=1:np
              enddo ! j=1:mp
            enddo ! k=1:l

            call loop_update_flags(kv,iuv)
            call update(cpmp(:,:,:,kv), &
              np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0,iuv) 
            call update(cnmp(:,:,:,kv), &
              np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0,iuv) 
          enddo ! kv=1:nv-1

        !!! correction coeffiecients for rho
          do kv=1,nvp-1
            do k=1,l
              do j=1,mp
                do i=1,np
                  rmxuse=-1._wp*mx(i,j,k,kv)
                  if(mx(i,j,k,kv).lt.tolerance) rmxuse=0.0_wp
                  rmnuse=-1._wp*mn(i,j,k,kv)
                  if(mn(i,j,k,kv).lt.tolerance) rmnuse=0.0_wp
                  a1p=cpmp(i,j,k,kv)+max(0._wp,sign(1._wp,mx(i,j,k,kv)))
                  a2p=cnmp(i,j,k,kv)+max(0._wp,sign(1._wp,rmnuse))
                  a1n=cnmp(i,j,k,kv)+max(0._wp,sign(1._wp,mn(i,j,k,kv)))
                  a2n=cpmp(i,j,k,kv)+max(0._wp,sign(1._wp,rmxuse))
                  tmpp=min(a1p,a2p)
                  tmpn=min(a1n,a2n)
                  if(kv.eq.1) cpmp(i,j,k,nvp)=tmpp
                  if(kv.eq.1) cnmp(i,j,k,nvp)=tmpn
                  cpmp(i,j,k,nvp)=min(tmpp,cpmp(i,j,k,nvp))
                  cnmp(i,j,k,nvp)=min(tmpn,cnmp(i,j,k,nvp))
                enddo !i=1:np
              enddo !j=1:mp
            enddo !k=1:l
          enddo !kv=1:nv-1 

          call update(cpmp(:,:,:,nvp), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0,0) 
          call update(cnmp(:,:,:,nvp), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0,0) 

          do kv=1,nvp
            do k=1,l
              do j=1,mp
                do i=1+leftedge,np                        
                  v1mp(i,j,k,kv)= pp(v1mp(i,j,k,kv)) &
                   *min(1._wp,cpmp(i,j,k,kv),cnmp(i-1,j,k,kv)) &
                   -pn(v1mp(i,j,k,kv)) &
                   *min(1._wp,cpmp(i-1,j,k,kv),cnmp(i,j,k,kv))
                enddo
              enddo
            enddo

            do k=1,l
              do j=1+botedge,mp                        
                do i=1,np
                  v2mp(i,j,k,kv)= pp(v2mp(i,j,k,kv)) &
                    *min(1._wp,cpmp(i,j,k,kv),cnmp(i,j-1,k,kv)) &
                    -pn(v2mp(i,j,k,kv)) &
                    *min(1._wp,cpmp(i,j-1,k,kv),cnmp(i,j,k,kv))
                enddo
              enddo
            enddo

            do k=2,l
              do j=1,mp
                do i=1,np
                  v3mp(i,j,k,kv)= pp(v3mp(i,j,k,kv)) &
                    *min(1._wp,cpmp(i,j,k,kv),cnmp(i,j,k-1,kv)) &
                    -pn(v3mp(i,j,k,kv)) &
                    *min(1._wp,cpmp(i,j,k-1,kv),cnmp(i,j,k,kv))
                enddo
              enddo
            enddo
          enddo
          if(itrfct.lt.nfct) then
            do kv=1,nvp
              do k=1,l
                do j=1,mp
                  do i=1,np+1
                    f1o(i,j,k,kv)=donor(c2,c2,v1mp(i,j,k,kv))
                    f1(i,j,k,kv)=f1(i,j,k,kv)-f1o(i,j,k,kv)
                  enddo
                enddo
              enddo

              do k=1,l
                do j=1,mp+1
                  do i=1,np
                    f2o(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))
                    f2(i,j,k,kv)=f2(i,j,k,kv)-f2o(i,j,k,kv)
                  enddo
                enddo
              enddo

              do k=1,l+1
                do j=1,mp
                  do i=1,np
                    f3o(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))
                    f3(i,j,k,kv)=f3(i,j,k,kv)-f3o(i,j,k,kv)
                  enddo
                enddo
              enddo

              if(rightedge.eq.0) then                       
                call update(f1o(:,:,:,kv), &
                  np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,1,1) 
              else
                call update(f1o(:,:,:,kv), &
                  np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0,1) 
              endif
 
              if(topedge.eq.0) then                       
                call update(f2o(:,:,:,kv), &
                  np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,2,2) 
              else
                call update(f2o(:,:,:,kv), &
                  np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,3,0,2) 
              endif

              do k=1,l
                do j=1,mp
                  do i=1,np
                    x(i,j,k,kv)=x(i,j,k,kv)-( &
                      f1o(i+1,j,k,kv)-f1o(i,j,k,kv) &
                     +f2o(i,j+1,k,kv)-f2o(i,j,k,kv) &
                     +f3o(i,j,k+1,kv)-f3o(i,j,k,kv) ) &
                     /h(i,j,k)
                  enddo !i
                enddo !j
              enddo !k 
              call loop_update_flags(kv,iuv)
              call update(x(:,:,:,kv),np,mp,l, &
                1-ih,np+ih,1-ih,mp+ih,1,0,iuv) 
            enddo !kv=1:nv
          endif ! itrfct < nfct
        enddo ! itrfct=1:nfct
      endif ! nonos.eq.1 
    endif ! itr.ne.iord  
  enddo ! itr=1,iord  

  deallocate (v1mp,v2mp,v3mp)
  deallocate (f1,f2,f3)
  deallocate (f1o,f2o,f3o)
  deallocate (cpmp,cnmp)
  deallocate (mx,mn,mxo,mno)

end subroutine mpdata
