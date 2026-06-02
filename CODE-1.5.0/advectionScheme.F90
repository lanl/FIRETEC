!**********************************************************************!
!     this routine computes the advection on the large time steps and
!     the nv variables. It can use either donorcell or fct
!
!**********************************************************************!

module advection_functions
  use gridlist_variables, only : prec

  contains

  real(prec) function pp(y)
    Implicit None
    real(prec),intent(in) :: y

    pp = max(0.0,y)

  end function pp
  !*************************************
  real(prec) function pn(y)
    Implicit None
    real(prec),intent(in) :: y

    pn = -min(0.0,y)

  end function pn
  !*************************************
  real(prec) function donor(y1,y2,a10)
    Implicit None
    real(prec),intent(in) :: y1,y2,a10
    !real(prec),internal :: pp,pn
    ! donor=y1 for positive velocity

    donor = pp(a10)*y1-pn(a10)*y2

  end function donor
  !*************************************
  real(prec) function vdiv1(a1,a2,a3,r)
    Implicit None
    real(prec),intent(in) :: a1,a2,a3,r

    vdiv1 = 0.25*a2*(a3-a1)/r

  end function vdiv1
  !*************************************
  real(prec) function vdiv2(aa,b1,b2,b3,b4,r)
    Implicit None
    real(prec),intent(in) :: aa,b1,b2,b3,b4,r

    vdiv2 = 0.25*aa*(b1+b2-b3-b4)/r

  end function vdiv2
  !*************************************
  real(prec) function rat2(z1,z2)
    Implicit None
    real(prec),intent(in) :: z1,z2

    rat2 = (z2-z1)*0.5

  end function rat2
  !*************************************
  real(prec) function rat4(z0,z1,z2,z3)
    Implicit None
    real(prec),intent(in) :: z0,z1,z2,z3

    rat4 = (z3+z2-z1-z0)*0.25

  end function rat4
  !*************************************
  real(prec) function vdyf(x1,x2,aa,r)
    Implicit None
    real(prec),intent(in) :: x1,x2,aa,r
    !real(prec),external :: rat2

    vdyf = (abs(aa)-aa**2/r)*rat2(x1,x2)

  end function vdyf
  !*************************************
  real(prec) function vcorr(aa,b,y0,y1,y2,y3,r)
    Implicit None
    real(prec),intent(in) :: aa,b,y0,y1,y2,y3,r
    !real(prec),external :: rat4

    vcorr = -0.125*aa*b/r*rat4(y0,y1,y2,y3)

  end function vcorr
  !*************************************

end module advection_functions

!-----------------------------------------------------------------------
! from u,v,o divided by gi and cell centered, this routine set 
! u1,u2,u3, the contravariant face coordinates*dt/dx*1/gi required 
! for advective schemes when iorder=2, u1,u2,u3 update 
! (for antidiffusive vel in mpdata...)
!-----------------------------------------------------------------------

subroutine setAdvectiveVelocities(u,v,o,il,iu,jl,ju,lls,dta,order)
  use gridlist_variables, only : l,ih,ibclatopen,prec
  use forcings, only : u1,u2,u3
  use xvall, only : iuvel,ivvel,xe,nv
  use metric_variables, only : gi
  use gridsetup, only : np,mp,dxi,dyi,dzi
  use msga_variables, only : leftdedge,rightdedge,botdedge,topdedge
  use metric_variables_old
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,order
  real(prec),intent(in) :: dta
  real(prec),intent(inout) :: u(il:iu,jl:ju,lls),v(il:iu,jl:ju,lls)
  real(prec),intent(inout) :: o(il:iu,jl:ju,lls)

  integer :: i,j,k
  real(prec) :: gcx,gcy,gcz

  ! Executable Code
  gcx=dta*dxi
  gcy=dta*dyi
  gcz=dta*dzi

  ! compute contravariant cell face in the x direction
  call update(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,2,0)
  do k=1,l
    do j=1,mp
      do i=1,np+1
        u1(i,j,k)=0.5*(u(i,j,k)+u(i-1,j,k))*gcx           
      enddo
    enddo
  enddo
  if(leftdedge.eq.1) then
    do k=1,l
      do j=1,mp
        ! 0 gradient between xe, cause xe(0 is not defined)
        u1(1,j,k)=xe(1,j,k,iuvel)/xe(1,j,k,nv)*gcx/gi(i,j,k)
      enddo
    enddo
  endif
  if(rightdedge.eq.1) then
    if (ibclatopen.eq.0) then
      do k=1,l
        do j=1,mp
          u1(np+1,j,k)=xe(np,j,k,1)/xe(np,j,k,nv)*gcx/gi(i,j,k)
          ! here we assume 0 gradient between xe
        enddo
      enddo
    elseif(ibclatopen.eq.1) then
      do k=1,l
        do j=1,mp
          u1(np+1,j,k)=max(u1(np,j,k),0.0)
        enddo
      enddo
    elseif(ibclatopen.eq.2) then
      do k=1,l
        do j=1,mp
          u1(np+1,j,k)=u1(np,j,k)
        enddo
      enddo
    endif 
  endif

  ! compute contravariant cell face in the y direction
  call update(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,3,0)
  do k=1,l
    do j=1,mp+1
      do i=1,np
        u2(i,j,k)=0.5*(v(i,j,k)+v(i,j-1,k))*gcy
      enddo
    enddo
  enddo
  if(botdedge.eq.1)then
    if(ibclatopen.eq.0) then
      do k=1,l
        do i=1,np
          ! here we assume 0 gradient between xe
          u2(i,1,k)=xe(i,1,k,ivvel)/xe(i,1,k,nv)*gcy/gi(i,j,k)
        enddo
      enddo
    elseif(ibclatopen.eq.1) then
      do k=1,l
        do i=1,np
          u2(i,1,k)=min(0.0,u2(i,2,k))
        enddo
      enddo
    elseif(ibclatopen.eq.2) then
      do k=1,l
        do i=1,np
          u2(i,1,k)=u2(i,2,k)
        enddo
      enddo
    endif
  endif
  if(topdedge.eq.1)then
    if(ibclatopen.eq.0) then
      do k=1,l
        do i=1,np
          ! here we assume 0 gradient between xe
          u2(i,mp+1,k)=xe(i,mp,k,ivvel)/xe(i,mp,k,nv)*gcy/gi(i,j,k)
        enddo
      enddo
    elseif(ibclatopen.eq.1) then
      do k=1,l
        do i=1,np
          u2(i,mp+1,k)=max(0.0,u2(i,mp,k))
        enddo
      enddo
    elseif(ibclatopen.eq.2) then
      do k=1,l
        do i=1,np
          u2(i,mp+1,k)=u2(i,mp,k)
        enddo
      enddo
    endif
  endif

  ! compute contravariant cell face in the z direction
  do k=2,l
    do j=1,mp
      do i=1,np
        u3(i,j,k)=0.5*(o(i,j,k)+o(i,j,k-1))*gcz
      enddo
    enddo
  enddo
  do j=1,mp
    do i=1,np
      u3(i,j,1)=0.0
      u3(i,j,l+1)=0.0
    enddo
  enddo

  ! the following updated are required for antidiffusive vel in
  ! mpdata
  if (order.eq.2) then
    call update(u1,u1,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1,1)
    call update(u2,u2,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1,2)
    call update(u3,u3,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1,0)
  endif

end subroutine setAdvectiveVelocities

!-----------------------------------------------------------------------
! mpdata is the 2nd order advection scheme of Smolarkiewicz (1984)
!-----------------------------------------------------------------------
subroutine mpdata(x,h,xe,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : l,ih,iord,nonos,idiv,nfct,nonosold,prec
  use forcings, only : u1,u2,u3
  use gridsetup, only : np,mp
  use msga_variables, only : leftdedge,rightdedge,botdedge,topdedge, &
    npos,mpos
  use advection_functions
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(prec),intent(in) :: h(il:iu,jl:ju,lls)
  real(prec),intent(in) :: xe(il:iu,jl:ju,lls,nvp)
  real(prec),intent(inout) :: x(il:iu,jl:ju,lls,nvp)
  
  !real(prec),external :: donor,vcorr,vdiv1,vdiv2,vdyf,pn,pp

  integer :: n3,kv,i,j,k,ip,im,jp,jm,kp,km,itrfct
  integer :: itr,illim,iulim,jllim,julim,ia,ja
  real(prec) :: ep,rhoin,rhoout,ain,aout,rmxuse,rmnuse,a1p,a2p,a1n,a2n
  real(prec) :: tmpp,tmpn,tmp,v1d,v2d,c1,c2
  real(prec) :: one=1.0
  real(prec),allocatable:: v1mp(:,:,:,:),v2mp(:,:,:,:),v3mp(:,:,:,:)
  real(prec),allocatable:: f1(:,:,:,:),f2(:,:,:,:),f3(:,:,:,:)
  real(prec),allocatable:: f1o(:,:,:,:),f2o(:,:,:,:),f3o(:,:,:,:)
  real(prec),allocatable:: cpmp(:,:,:,:),cnmp(:,:,:,:)
  real(prec),allocatable:: mxo(:,:,:,:),mno(:,:,:,:)
  real(prec),allocatable:: mx(:,:,:,:),mn(:,:,:,:)
  real(prec),allocatable::  a(:,:,:,:)

  ! Executable Code
  n3=l+1
#ifdef DBL
  ep=1.d-10
#else
  ep=1.0e-10
#endif
  allocate (v1mp(il:iu+1,jl:ju  ,lls  ,nvp)); v1mp=0.0
  allocate (v2mp(il:iu  ,jl:ju+1,lls  ,nvp)); v2mp=0.0
  allocate (v3mp(il:iu  ,jl:ju  ,lls+1,nvp)); v3mp=0.0
  allocate (f1(il:iu+1,jl:ju  ,lls  ,nvp)  ); f1=  0.0
  allocate (f2(il:iu  ,jl:ju+1,lls  ,nvp)  ); f2=  0.0
  allocate (f3(il:iu  ,jl:ju  ,lls+1,nvp)  ); f3=  0.0
  allocate (f1o(il:iu+1,jl:ju  ,lls  ,nvp) ); f1o= 0.0
  allocate (f2o(il:iu  ,jl:ju+1,lls  ,nvp) ); f2o= 0.0
  allocate (f3o(il:iu  ,jl:ju  ,lls+1,nvp) ); f3o= 0.0
  allocate (cpmp(il:iu,jl:ju,lls,nvp)      ); cpmp=0.0
  allocate (cnmp(il:iu,jl:ju,lls,nvp)      ); cnmp=0.0
  allocate (mx(il:iu,jl:ju,lls,nvp-1)      ); mx=  0.0
  allocate (mn(il:iu,jl:ju,lls,nvp-1)      ); mn=  0.0
  allocate (mxo(il:iu,jl:ju,lls,nvp)       ); mxo= 0.0
  allocate (mno(il:iu,jl:ju,lls,nvp)       ); mno= 0.0
  allocate (a(il:iu,jl:ju,lls,nvp-1)       ); a=   0.0

  do kv=1,nvp
    call update(x(:,:,:,kv),x(:,:,:,kv),np,mp,l, &
      1-ih,np+ih,1-ih,mp+ih,1,0)  

    do k=1,l+1
      do j=1,mp
        do i=1,np
          v3mp(i,j,k,kv) = u3(i,j,k)
        enddo
      enddo
    enddo

    do k=1,l
      do j=1,mp
        do i=1,np+1
          v1mp(i,j,k,kv) = u1(i,j,k)
        enddo
      enddo

      do i=1,np
        do j=1,mp+1
          v2mp(i,j,k,kv) = u2(i,j,k)
        enddo
      enddo
    enddo
  enddo

  if(nonos.eq.1) then
    do kv=1,nvp-1
      do k=1,l
        do j=1,mp
          do i=1,np
            a(i,j,k,kv)=x(i,j,k,kv)/x(i,j,k,nvp)
          enddo ! i=1:np
        enddo ! j=1:mp
      enddo ! k=1:l
      call update(a(1-ih,1-ih,1,kv),a(1-ih,1-ih,1,kv),np,mp,l, &
        1-ih,np+ih,1-ih,mp+ih,1,0) 

      do k=1,l
        km=max(k-1,1)
        kp=min(k+1,l)
        do j=1,mp
          if (botdedge.eq.1 .and. j.eq.1) then 
            jm = 1
          else
            jm = j - 1
          endif
          if (topdedge.eq.1 .and. j.eq.mp) then 
            jp = mp
          else
            jp = j + 1
          endif
          do i=1,np
            if (leftdedge.eq.1 .and. i.eq.1) then 
              im = 1
            else
              im = i - 1
            endif
            if (rightdedge.eq.1 .and. i.eq.np) then 
              ip = np
            else
              ip = i + 1
            endif
            mx(i,j,k,kv)=max(a(im,j,k,kv),a(i,j,k,kv),a(ip,j,k,kv), &
              a(i,jm,k,kv),a(i,jp,k,kv),a(i,j,kp,kv),a(i,j,km,kv))
            mn(i,j,k,kv)=min(a(im,j,k,kv),a(i,j,k,kv),a(ip,j,k,kv), &
              a(i,jm,k,kv),a(i,jp,k,kv),a(i,j,kp,kv),a(i,j,km,kv))
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
          if (botdedge.eq.1 .and. j.eq.1) then 
            jm = 1
          else
            jm = j - 1
          end if
          if (topdedge.eq.1 .and. j.eq.mp) then 
            jp = mp
          else
            jp = j + 1
          end if
          do i=1,np
            if (leftdedge.eq.1 .and. i.eq.1) then 
              im = 1
            else
              im = i - 1
            endif
            if (rightdedge.eq.1 .and. i.eq.np) then 
              ip = np
            else
              ip = i + 1
            endif
            mxo(i,j,k,kv)=max(x(im,j,k,kv),x(i,j,k,kv),x(ip,j,k,kv), &
              x(i,jm,k,kv),x(i,jp,k,kv),x(i,j,kp,kv),x(i,j,km,kv))
            mno(i,j,k,kv)=min(x(im,j,k,kv),x(i,j,k,kv),x(ip,j,k,kv), &
              x(i,jm,k,kv),x(i,jp,k,kv),x(i,j,kp,kv),x(i,j,km,kv))
          enddo ! i=1:np
        enddo ! j=1:mp
      enddo ! k=1:nm3(l)
    enddo ! kv=1:nv
  endif ! nonosold==1

  c1=1.
  c2=0.

  do itr=1,iord ! loop to LINE 1402  
    do kv=1,nvp
      call update(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l, &
                   1-ih,np+ih,1-ih,mp+ih,1,0) 

      illim=1+leftdedge*1                       
      do k=1,l
        do j=1,mp
          do i=illim,np
            f1(i,j,k,kv)=donor(c1*x(i-1,j,k,kv)+c2, &
              c1*x(i,j,k,kv)+c2,v1mp(i,j,k,kv))
          enddo
        enddo
      enddo

      if (rightdedge.eq.0) then                       
        call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l, &
                     1-ih,np+ih+1,1-ih,mp+ih,2,0) 
      else
        call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l, &
                     1-ih,np+ih+1,1-ih,mp+ih,2,0) 
      endif

      if (leftdedge.eq.1) then                       
        do k=1,l
          do j=1,mp
            f1(1 ,j,k,kv)=donor(c1*xe(1,j,k,kv)+c2, &
               c1*x(1,j,k,kv)+c2, &
               v1mp(1,j,k,kv))
          enddo
        enddo
      endif

      if (rightdedge.eq.1) then                       
        do k=1,l
          do j=1,mp
            f1(np+1,j,k,kv)= donor(c1*x(np,j,k,kv)+c2, &
                            c1*xe(np,j,k,kv)+c2, &
                            v1mp(np+1,j,k,kv))
          enddo
        enddo
      endif

      if (topdedge.eq.0) then                       
        call update(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp,l, &
                    1-ih,np+ih,1-ih,mp+ih+1,3,0) 
      else
        call update(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv),np,mp+1, &
                    l,1-ih,np+ih,1-ih,mp+ih+1,3,0) 
      endif

      if (botdedge.eq.0) then                       
        do k=1,l
          do j=1,mp
            do i=1,np
              f2(i,j,k,kv)=donor(c1*x(i,j-1,k,kv)+c2, &
                                c1*x(i,j,k,kv)+c2, &
                                v2mp(i,j,k,kv))
            enddo
          enddo
        enddo
      endif
      if (botdedge.eq.1) then                       
        do k=1,l
          do j=2,mp
            do i=1,np
              f2(i,j,k,kv)=donor(c1*x(i,j-1,k,kv)+c2, &
                                c1*x(i,j,k,kv)+c2, &
                                v2mp(i,j,k,kv))
            enddo
          enddo
        enddo
      endif

      if (topdedge.eq.0) then                       
        call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l, &
                    1-ih,np+ih,1-ih,mp+ih+1,3,0) 
      else
        call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l, &
                    1-ih,np+ih,1-ih,mp+ih+1,3,0) 
      endif

      if (botdedge.eq.1) then                       
        do k=1,l
          do i=1,np
            f2(i,1 ,k,kv)= donor(c1*xe(i,1,k,kv)+c2, &
                          c1*x(i,1,k,kv)+c2, &
                          v2mp(i,1,k,kv))
          enddo
        enddo
      endif
      if (topdedge.eq.1) then                       
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

      do j=1,mp
        do i=1,np
          f3(i,j, 1,kv)=0.
          f3(i,j,n3,kv)=0.
        enddo
      enddo

      if (rightdedge.eq.0) then                       
        call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l, &
                    1-ih,np+ih+1,1-ih,mp+ih,2,0) 
      else
        call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l, &
                    1-ih,np+ih+1,1-ih,mp+ih,2,0) 
      endif
      if(topdedge.eq.0) then                       
        call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l, &
                    1-ih,np+ih,1-ih,mp+ih+1,3,0) 
      else
        call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l, &
                    1-ih,np+ih,1-ih,mp+ih+1,3,0) 
      endif

      do k=1,l
        do j=1,mp
          do i=1,np
            x(i,j,k,kv)=x(i,j,k,kv)-( f1(i+1,j,k,kv)-f1(i,j,k,kv) &
                                    +f2(i,j+1,k,kv)-f2(i,j,k,kv) &
                                    +f3(i,j,k+1,kv)-f3(i,j,k,kv) ) &
                      /h(i,j,k)
          enddo
        enddo
      enddo
    enddo ! kv=1:nv

    if(itr.ne.iord)then  !! Skip if iord.eq.1: start of mpdata when iord=2: escaping to LINE 1403
      c1=0.
      c2=1.
      do kv=1,nvp
        do k=1,l
          do j=1,mp
            do i=1,np+1
              f1(i,j,k,kv)=v1mp(i,j,k,kv)
              v1mp(i,j,k,kv)=0.
            enddo
          enddo
        enddo

        do k=1,l
          do j=1,mp+1
            do i=1,np
              f2(i,j,k,kv)=v2mp(i,j,k,kv)
              v2mp(i,j,k,kv)=0.
            enddo
          enddo
        enddo
 
        do k=1,n3
          do j=1,mp
            do i=1,np
              f3(i,j,k,kv)=v3mp(i,j,k,kv)
              v3mp(i,j,k,kv)=0.
            enddo
          enddo
        enddo

        !!! compute antidiffusive velocities in x direction
        call update(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l, &
                    1-ih,np+ih,1-ih,mp+ih,1,0) 

        if (topdedge.eq.0) then 
          call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l, &
                      1-ih,np+ih,1-ih,mp+ih+1,1,0) 
        else
          call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l, &
                      1-ih,np+ih,1-ih,mp+ih+1,1,0)  
        endif
        call update(f3(1-ih,1-ih,1,kv),f3(1-ih,1-ih,1,kv),np,mp,l+1, &
                      1-ih,np+ih,1-ih,mp+ih,1,0) 

        illim = 1 + leftdedge                       
        iulim = np
        jllim = 1  + botdedge                       
        julim = mp - topdedge                       
        do k=2,l-1
          do j=jllim,julim
            do i=illim,iulim
              v1mp(i,j,k,kv)=vdyf(  &
                x(i-1,j,k,kv),x(i,j,k,kv),f1(i,j,k,kv), &
               .5*(h(i-1,j,k)+h(i,j,k))) &
               +vcorr(f1(i,j,k,kv),f2(i-1,j,k,kv)+f2(i-1,j+1,k,kv) &
               +f2(i,j+1,k,kv)+f2(i,j,k,kv), &
               x(i-1,j-1,k,kv),x(i,j-1,k,kv), &
               x(i-1,j+1,k,kv),x(i,j+1,k,kv), &
               .5*(h(i-1,j,k)+h(i,j,k))) &
               +vcorr(f1(i,j,k,kv), &
               f3(i-1,j,k,kv)+f3(i-1,j,k+1,kv) &
               +f3(i,j,k+1,kv)+f3(i,j,k,kv), &
               x(i-1,j,k-1,kv),x(i,j,k-1,kv), &
               x(i-1,j,k+1,kv),x(i,j,k+1,kv), &
               .5*(h(i-1,j,k)+h(i,j,k)))
            enddo
          enddo
        enddo

        if (rightdedge.eq.0) then                       
          call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l, &
           1-ih,np+ih+1,1-ih,mp+ih,1,0) 
        else
          call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l, &
           1-ih,np+ih+1,1-ih,mp+ih,1,0) 
        endif

        if(idiv.eq.1) then
          illim = 1  + leftdedge                       
          iulim = np
          jllim = 1  + botdedge                       
          julim = mp - topdedge                       
          do k=2,l-1
            do j=jllim,julim
              do i=illim,iulim
                v1d=-vdiv1(f1(i-1,j,k,kv),f1(i,j,k,kv),f1(i+1,j,k,kv), &
                 .5*(h(i-1,j,k)+h(i,j,k))) &
                 -vdiv2(f1(i,j,k,kv),f2(i-1,j+1,k,kv), &
                 f2(i,j+1,k,kv),f2(i-1,j,k,kv),&
                 f2(i,j,k,kv),.5*(h(i-1,j,k)+h(i,j,k))) &
                 -vdiv2(f1(i,j,k,kv),f3(i-1,j,k+1,kv), &
                 f3(i,j,k+1,kv),f3(i-1,j,k,kv), &
                 f3(i,j,k,kv),.5*(h(i-1,j,k)+h(i,j,k)))
                v1mp(i,j,k,kv)=v1mp(i,j,k,kv)+ &
                  (pp(v1d)*x(i-1,j,k,kv)-pn(v1d)*x(i,j,k,kv))
              enddo
            enddo
          enddo
        endif ! idiv.eq.1

        !!! compute antidiffusive velocities in y direction
        illim = 1  + leftdedge                       
        iulim = np - rightdedge                       
        jllim = 1  + botdedge                       
        julim = mp
        do k=2,l-1
          do j=jllim,julim
            do i=illim,iulim
              v2mp(i,j,k,kv)=vdyf( &
                x(i,j-1,k,kv),x(i,j,k,kv),f2(i,j,k,kv), &
               .5*(h(i,j-1,k)+h(i,j,k))) &
               +vcorr(f2(i,j,k,kv), &
               f1(i,j-1,k,kv)+f1(i,j,k,kv) &
               +f1(i+1,j,k,kv)+f1(i+1,j-1,k,kv), &
               x(i-1,j-1,k,kv),x(i-1,j,k,kv), &
               x(i+1,j-1,k,kv),x(i+1,j,k,kv), &
               .5*(h(i,j-1,k)+h(i,j,k))) &
               +vcorr(f2(i,j,k,kv), &
               f3(i,j-1,k,kv)+f3(i,j,k,kv) &
               +f3(i,j,k+1,kv)+f3(i,j-1,k+1,kv), &
               x(i,j-1,k-1,kv),x(i,j,k-1,kv), &
               x(i,j-1,k+1,kv),x(i,j,k+1,kv), &
               .5*(h(i,j-1,k)+h(i,j,k)))
            enddo
          enddo
        enddo

        if(idiv.eq.1) then
          illim = 1  + leftdedge                       
          iulim = np - rightdedge                       
          jllim = 1  + botdedge                       
          julim = mp
          do k=2,l-1
            do j=jllim,julim
              do i=illim,iulim
                v2d=-vdiv1(f2(i,j-1,k,kv),       &
                 f2(i,j,k,kv),f2(i,j+1,k,kv),   &
                 .5*(h(i,j-1,k)+h(i,j,k)))      &
                 -vdiv2(f2(i,j,k,kv),           &
                 f1(i+1,j-1,k,kv),f1(i+1,j,k,kv), &
                 f1(i,j-1,k,kv),f1(i,j,k,kv),   &
                 .5*(h(i,j-1,k)+h(i,j,k)))      &
                 -vdiv2(f2(i,j,k,kv),f3(i,j-1,k+1,kv), &
                 f3(i,j,k+1,kv),f3(i,j-1,k,kv), &
                 f3(i,j,k,kv),    .5*(h(i,j-1,k)+h(i,j,k)))
               v2mp(i,j,k,kv)=v2mp(i,j,k,kv)+(pp(v2d)*x(i,j-1,k,kv)- &
                 pn(v2d)*x(i,j,k,kv))
              enddo
            enddo
          enddo
        endif ! idiv.eq.1 

        !!! compute antidiffusive velocities in z direction
        illim = 1  + leftdedge                       
        iulim = np - rightdedge                       
        jllim = 1  + botdedge                       
        julim = mp - topdedge                       

        do k=2,l
          do j=jllim,julim
            do i=illim,iulim
              v3mp(i,j,k,kv)=vdyf( &
                x(i,j,k-1,kv),x(i,j,k,kv),f3(i,j,k,kv), &
                .5*(h(i,j,k-1)+h(i,j,k)))               &
                +vcorr(f3(i,j,k,kv),                    &
                f1(i,j,k-1,kv)+f1(i,j,k,kv)             &
                +f1(i+1,j,k,kv)+f1(i+1,j,k-1,kv),       &
                x(i-1,j,k-1,kv),x(i-1,j,k,kv),          &
                x(i+1,j,k-1,kv),x(i+1,j,k,kv),          &
                .5*(h(i,j,k-1)+h(i,j,k)))               &
                +vcorr(f3(i,j,k,kv),                    &
                f2(i,j,k-1,kv)+f2(i,j+1,k-1,kv)         &
                +f2(i,j+1,k,kv)+f2(i,j,k,kv),           &
                x(i,j-1,k-1,kv),x(i,j-1,k,kv),          &
                x(i,j+1,k-1,kv),x(i,j+1,k,kv),          &
                .5*(h(i,j,k-1)+h(i,j,k)))
            enddo
          enddo
        enddo

        if(idiv.eq.1) then
          illim = 1  + leftdedge                      
          iulim = np - rightdedge                    
          jllim = 1  + botdedge                       
          julim = mp - topdedge                      
          do k=2,l
            do j=jllim,julim
              do i=illim,iulim
                ia=(npos-1)*np + i
                ja=(mpos-1)*mp + j

                v2d=-vdiv1(f3(i,j,k-1,kv), &
                 f3(i,j,k,kv),f3(i,j,k+1,kv), &
                 .5*(h(i,j,k-1)+h(i,j,k))) &
                 -vdiv2(f3(i,j,k,kv), &
                 f1(i+1,j,k-1,kv),f1(i+1,j,k,kv), &
                 f1(i,j,k-1,kv), &
                 f1(i,j,k,kv),  .5*(h(i,j,k-1)+h(i,j,k))) &
                 -vdiv2(f3(i,j,k,kv), &
                 f2(i,j+1,k-1,kv),f2(i,j+1,k,kv), &
                 f2(i,j,k-1,kv), &
                 f2(i,j,k,kv),    .5*(h(i,j,k-1)+h(i,j,k)))
                v3mp(i,j,k,kv)=v3mp(i,j,k,kv) &
                 +(pp(v2d)*x(i,j,k-1,kv)-pn(v2d)*x(i,j,k,kv))
              enddo
            enddo
          enddo

        endif ! idiv==1
        if (rightdedge.eq.0) then                       
          call update(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv), &
            np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1,0)  
        else
          call update(v1mp(1-ih,1-ih,1,kv),v1mp(1-ih,1-ih,1,kv), &
            np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1,0)  
        endif

        if (topdedge.eq.0) then                       
          call update(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,0)  
        else
          call update(v2mp(1-ih,1-ih,1,kv),v2mp(1-ih,1-ih,1,kv), &
            np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,3,0)   
        endif
      enddo ! kv=1:nv

      if(nonosold.eq.1) then  ! non-osscilatory option
        do kv=1,nvp
          do k=1,l
            km=max0(k-1,1  )
            kp=min0(k+1,l)
            do j=1,mp
              if (botdedge.eq.1 .and. j.eq.1) then   
                jm = 1
              else
                jm = j - 1
              endif
              if (topdedge.eq.1 .and. j.eq.mp) then  
                jp = mp
              else
                jp = j + 1
              endif
              do i=1,np
                if (leftdedge.eq.1 .and. i.eq.1) then 
                  im = 1
                else
                  im = i - 1
                endif
                if (rightdedge.eq.1 .and. i.eq.np) then 
                  ip = np
                else
                  ip = i + 1
                endif
                mxo(i,j,k,kv)=max(x(im,j,k,kv),x(i,j,k,kv), &
                  x(ip,j,k,kv),mxo(i,j,k,kv),                &
                  x(i,jm,k,kv),x(i,jp,k,kv),                 &
                  x(i,j,kp,kv),x(i,j,km,kv))
                mno(i,j,k,kv)=min(x(im,j,k,kv),x(i,j,k,kv), &
                  x(ip,j,k,kv),mno(i,j,k,kv), &
                  x(i,jm,k,kv),x(i,jp,k,kv), &
                  x(i,j,kp,kv),x(i,j,km,kv))
              enddo
            enddo
          enddo

          illim = 1
          iulim = np + rightdedge                       
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
          julim = mp + topdedge                       
          do k=1,l
            do j=jllim,julim
              do i=illim,iulim
                f2(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))
              enddo
            enddo
          enddo

          do k=1,n3
            do j=1,mp
              do i=1,np
                f3(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))
              enddo
            enddo
          enddo

          if (rightdedge.eq.0) then                       
            call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv), &
              np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0) 
          else
            call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv), &
              np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0) 
          endif

          if (topdedge.eq.0) then                       
            call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv), &
              np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,0) 
          else
            call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv), &
              np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,3,0) 
          endif

          do k=1,l
            do j=1,mp
              do i=1,np
                cpmp(i,j,k,kv)=(mxo(i,j,k,kv)-x(i,j,k,kv))*h(i,j,k)/ &
                  ( pn(f1(i+1,j,k,kv))+pp(f1(i,j,k,kv))              &
                  +pn(f2(i,j+1,k,kv))+pp(f2(i,j,k,kv))               &
                  +pn(f3(i,j,k+1,kv))+pp(f3(i,j,k,kv))+ep)

                cnmp(i,j,k,kv)=(x(i,j,k,kv)-mno(i,j,k,kv))*h(i,j,k)/ &
                  ( pp(f1(i+1,j,k,kv))+pn(f1(i,j,k,kv)) &
                  +pp(f2(i,j+1,k,kv))+pn(f2(i,j,k,kv)) &
                  +pp(f3(i,j,k+1,kv))+pn(f3(i,j,k,kv))+ep)
              enddo
            enddo
          enddo

          call update(cpmp(1-ih,1-ih,1,kv),cpmp(1-ih,1-ih,1,kv), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0) 
          call update(cnmp(1-ih,1-ih,1,kv),cnmp(1-ih,1-ih,1,kv), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0) 

          illim = 1 + leftdedge                       
          do k=1,l
            do j=1,mp
              do i=illim,np
                v1mp(i,j,k,kv)= pp(v1mp(i,j,k,kv))&
                 *min(1.,cpmp(i,j,k,kv),cnmp(i-1,j,k,kv))&
                 -pn(v1mp(i,j,k,kv))&
                 *min(1.,cpmp(i-1,j,k,kv),cnmp(i,j,k,kv))
              enddo
            enddo
          enddo

          jllim = 1 + botdedge                       
          do k=1,l
            do j=jllim,mp
              do i=1,np
                v2mp(i,j,k,kv)= pp(v2mp(i,j,k,kv))&
                 *min(1.,cpmp(i,j,k,kv),cnmp(i,j-1,k,kv))&
                 -pn(v2mp(i,j,k,kv))&
                 *min(1.,cpmp(i,j-1,k,kv),cnmp(i,j,k,kv))
              enddo
            enddo
          enddo

          do k=2,l
            do j=1,mp
              do i=1,np
                v3mp(i,j,k,kv)= pp(v3mp(i,j,k,kv)) &
                 *min(1.,cpmp(i,j,k,kv),cnmp(i,j,k-1,kv)) &
                 -pn(v3mp(i,j,k,kv)) &
                 *min(1.,cpmp(i,j,k-1,kv),cnmp(i,j,k,kv))
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
              iulim = np + rightdedge                       
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
              julim = mp + topdedge                       
              do k=1,l
                do j=jllim,julim
                  do i=illim,iulim
                    f2(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))
                  enddo
                enddo
              enddo

              do k=1,n3
                do j=1,mp
                  do i=1,np
                    f3(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))
                  enddo
                enddo
              enddo

              if (rightdedge.eq.0) then                       
                call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv), &
                  np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0) 
              else
                call update(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv), &
                  np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0) 
              endif

              if (topdedge.eq.0) then                       
                call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv), &
                  np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,0) 
              else
                call update(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv), &
                  np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,3,0) 
              endif
            enddo !kv=1:nv
          endif ! itrfct.eq.1

        !!! correction coefficients for variables nv-1
          do kv=1,nvp-1
            do k=1,l
              do j=1,mp
                do i=1,np
                  if(abs(mx(i,j,k,kv)).lt.ep) mx(i,j,k,kv)=0.
                  if(abs(mn(i,j,k,kv)).lt.ep) mn(i,j,k,kv)=0.
                  if(k.eq.1) f3(i,j,k,kv)=-f3(i,j,k+1,kv)
                  if(k.eq.1) f3(i,j,k,nvp)=-f3(i,j,k+1,nvp)
                  if(k.eq.l) f3(i,j,k+1,kv)=-f3(i,j,k,kv)
                  if(k.eq.l) f3(i,j,k+1,nvp)=-f3(i,j,k,nvp)
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
                   +pn(mx(i,j,k,kv))*rhoin+ep)
                 cnmp(i,j,k,kv)= &
                   pp(x(i,j,k,kv)-mn(i,j,k,kv)*x(i,j,k,nvp))*h(i,j,k) &
                   /(-aout+pp(mn(i,j,k,kv))*rhoin &
                   -pn(mn(i,j,k,kv))*rhoout+ep)
                enddo ! i=1:np
              enddo ! j=1:mp
            enddo ! k=1:l

            call update(cpmp(1-ih,1-ih,1,kv),cpmp(1-ih,1-ih,1,kv), &
              np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0) 
            call update(cnmp(1-ih,1-ih,1,kv),cnmp(1-ih,1-ih,1,kv), &
              np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0) 
          enddo ! kv=1:nv-1

        !!! correction coeffiecients for rho
          do kv=1,nvp-1
            do k=1,l
              do j=1,mp
                do i=1,np
                  rmxuse=-1.*mx(i,j,k,kv)
                  if(mx(i,j,k,kv).lt.1.e-8) rmxuse=0.
                  rmnuse=-1.*mn(i,j,k,kv)
                  if(mn(i,j,k,kv).lt.1.e-8) rmnuse=0.
                  a1p=cpmp(i,j,k,kv)+max(0.,sign(one, mx(i,j,k,kv)))
                  a2p=cnmp(i,j,k,kv)+max(0.,sign(one,rmnuse))
                  a1n=cnmp(i,j,k,kv)+max(0.,sign(one, mn(i,j,k,kv)))
                  a2n=cpmp(i,j,k,kv)+max(0.,sign(one,rmxuse))
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

          call update(cpmp(1-ih,1-ih,1,nvp),cpmp(1-ih,1-ih,1,nvp), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0) 
          call update(cnmp(1-ih,1-ih,1,nvp),cnmp(1-ih,1-ih,1,nvp), &
            np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0) 

          do kv=1,nvp
            do k=1,l
              do j=1,mp
                do i=1+leftdedge,np                        
                  v1mp(i,j,k,kv)= pp(v1mp(i,j,k,kv)) &
                   *min(1.,cpmp(i,j,k,kv),cnmp(i-1,j,k,kv)) &
                   -pn(v1mp(i,j,k,kv)) &
                   *min(1.,cpmp(i-1,j,k,kv),cnmp(i,j,k,kv))
                enddo
              enddo
            enddo

            do k=1,l
              do j=1+botdedge,mp                        
                do i=1,np
                  v2mp(i,j,k,kv)= pp(v2mp(i,j,k,kv)) &
                    *min(1.,cpmp(i,j,k,kv),cnmp(i,j-1,k,kv)) &
                    -pn(v2mp(i,j,k,kv)) &
                    *min(1.,cpmp(i,j-1,k,kv),cnmp(i,j,k,kv))
                enddo
              enddo
            enddo

            do k=2,l
              do j=1,mp
                do i=1,np
                  v3mp(i,j,k,kv)= pp(v3mp(i,j,k,kv)) &
                    *min(1.,cpmp(i,j,k,kv),cnmp(i,j,k-1,kv)) &
                    -pn(v3mp(i,j,k,kv)) &
                    *min(1.,cpmp(i,j,k-1,kv),cnmp(i,j,k,kv))
                enddo
              enddo
            enddo
          enddo
          if(itrfct.lt.nfct) then
            do kv=1,nvp
              do k=1,l
                do j=1,mp
                  do i=1,np+1
                    tmp=f1(i,j,k,kv)
                    f1o(i,j,k,kv)=donor(c2,c2,v1mp(i,j,k,kv))
                    f1(i,j,k,kv)=tmp-f1o(i,j,k,kv)
                  enddo
                enddo
              enddo

              do k=1,l
                do j=1,mp+1
                  do i=1,np
                    tmp=f2(i,j,k,kv)
                    f2o(i,j,k,kv)=donor(c2,c2,v2mp(i,j,k,kv))
                    f2(i,j,k,kv)=tmp-f2o(i,j,k,kv)
                  enddo
                enddo
              enddo

              do k=1,n3
                do j=1,mp
                  do i=1,np
                    tmp=f3(i,j,k,kv)
                    f3o(i,j,k,kv)=donor(c2,c2,v3mp(i,j,k,kv))
                    f3(i,j,k,kv)=tmp-f3o(i,j,k,kv)
                  enddo
                enddo
              enddo

              if(rightdedge.eq.0) then                       
                call update(f1o(1-ih,1-ih,1,kv),f1o(1-ih,1-ih,1,kv), &
                  np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0) 
              else
                call update(f1o(1-ih,1-ih,1,kv),f1o(1-ih,1-ih,1,kv), &
                  np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0) 
              endif
 
              if(botdedge.eq.0) then                       
                call update(f2o(1-ih,1-ih,1,kv),f2o(1-ih,1-ih,1,kv), &
                  np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,0) 
              else
                call update(f2o(1-ih,1-ih,1,kv),f2o(1-ih,1-ih,1,kv), &
                  np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,3,0) 
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
            enddo !kv=1:nv
          endif ! itrfct < nfct
        enddo ! itrfct=1:nfct
      endif ! nonos.eq.1 
    endif ! itr.ne.iord  
  enddo ! itr=1,iord  

  do kv=1,nvp
    call update(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l, &
                     1-ih,np+ih,1-ih,mp+ih,1,0) 
  enddo
  deallocate (v1mp,v2mp,v3mp)
  deallocate (f1,f2,f3)
  deallocate (f1o,f2o,f3o)
  deallocate (cpmp,cnmp)
  deallocate (mx,mn,mxo,mno)
  deallocate (a)

end subroutine mpdata
