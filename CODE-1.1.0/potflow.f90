  subroutine potflow(xv,il,iu,jl,ju,lls,nvp)
  use gridsetup
  use msga
  use xve
  use metryic
  implicit none
  integer::il,iu,jl,ju,lls,nvp
  real::xv(il:iu,jl:ju,lls,nvp)
  integer::i,j,k
  integer::iprocx,jprocy,ia,ja,nsize
  real::g13,g23
  real,allocatable::utilda(:,:,:), &
                    vtilda(:,:,:), &
                    wtilda(:,:,:), &
                    uconstl(:,:),uconstr(:,:), &
                    vconstb(:,:),vconstt(:,:)
  !fixed buffer error
   real,allocatable::tempu(:,:),tempv(:,:);

 
  allocate(utilda(il:iu,jl:ju,lls), &
           vtilda(il:iu,jl:ju,lls), &
           wtilda(il:iu,jl:ju,lls))
  utilda=0.
  vtilda=0.
  wtilda=0.

! ensure mass continuity is enforced at model boundaries
  call vbcadlite(xv,il,iu,jl,ju,lls,nvp)
  
  allocate(uconstl(m,lls),uconstr(m,lls)); uconstl=0.; uconstr=0.
  allocate(vconstb(n,lls),vconstt(n,lls)); vconstb=0.; vconstt=0.
  jprocy=mpi_rank/nprocx + 1
  iprocx=(mpi_rank+1)-(jprocy-1)*nprocx
  do k=1,L
    do j=1,mp
      ja=(jprocy-1)*mp + j
      if(leftedge==1)uconstl(ja,k)=xv(1,j,k,1)/gi(1,j,k)
      if(rightedge==1)uconstr(ja,k)=xv(np,j,k,1)/gi(np,j,k)
    enddo
    do i=1,np
      ia=(iprocx-1)*np + i
      if(botedge==1)vconstb(ia,k)=xv(i,1,k,2)/gi(i,1,k)
      if(topedge==1)vconstt(ia,k)=xv(i,mp,k,2)/gi(i,mp,k)
    enddo
  enddo

  !fixed buffer error
   allocate(tempu(m,lls))
   tempu=0.
   allocate(tempv(n,lls))
   tempv=0.
  nsize=m*lls
  call mpi_allreduce(uconstl,tempu,nsize,mpi_real,mpi_sum,mpi_comm_world,ierror)
  uconstl=tempu
  call mpi_allreduce(uconstr,tempu,nsize,mpi_real,mpi_sum,mpi_comm_world,ierror)
  uconstr=tempu
  nsize=n*lls
  call mpi_allreduce(vconstb,tempv,nsize,mpi_real,mpi_sum,mpi_comm_world,ierror)
  vconstb=tempv
  call mpi_allreduce(vconstt,tempv,nsize,mpi_real,mpi_sum,mpi_comm_world,ierror)
  vconstt=tempv

! Interior mass flows (not fluxes !) are linear combination of flows 
! at opposite ends of each boundary.
      jprocy=mpi_rank/nprocx + 1
      iprocx=(mpi_rank+1)-(jprocy-1)*nprocx
      do j=1,mp
        do i=1,np
          ia=(iprocx-1)*np + i
          ja=(jprocy-1)*mp + j
          utilda(i,j,:)=((n-ia)*uconstl(ja,:)  + (ia-1)*uconstr(ja,:))/real(n-1)
          vtilda(i,j,:)=((m-ja)*vconstb(ia,:)  + (ja-1)*vconstt(ia,:))/real(m-1)
        enddo
      enddo

! remove areal weight to get horizontal velocities  
  utilda=utilda*gi
  vtilda=vtilda*gi
  
! Compute cartesian velocity using constraint that velocity is tangent to 
! model topographical surface.
! Use wtilda to temporarily hold cartesian vertical velocity
  do k=1,L
    do j=1,mp
      do i=1,np
        g13=gmul(k)*c13(i,j)
        g23=gmul(k)*c23(i,j)
        wtilda(i,j,k)=-(g13*utilda(i,j,k)+g23*vtilda(i,j,k))/gi(i,j,k)
      enddo
    enddo
  enddo

  xv(:,:,:,1)=utilda
  xv(:,:,:,2)=vtilda
  xv(:,:,:,3)=wtilda
  
 
! call vbcad to ensure mass continuity at model boundaries
!  call vbcad(xv,1-ih,np+ih,1-ih,mp+ih,L,nv)
  deallocate(utilda,vtilda,wtilda,uconstl,uconstr,vconstb,vconstt)
  return
  end
  
  subroutine vbcadlite(xv,il,iu,jl,ju,lls,nvp)
! works like vbcad. Ensures that net mass flow through model domain is zero.
! There is no option for flow through the upper and ground boundaries. 
  use gridsetup
  use msga
  use xve
  use metryic
  implicit none
  integer::il,iu,jl,ju,lls,nvp
  real::xv(il:iu,jl:ju,lls,nvp)
  integer::kv,iflg
  real(8)::areain,areaout,massout,massperarea,area,massoutcheck,ratio
! compute net mass flow out of model domain (needs to be zero
! find total area of model boundaries, mass outflow area, inflow area, 
! net mass outflow.
    massout=0.
    areaout=0.
    areain=0.
    area=0
    area=leftedge*dy*dz*sum(1.d0/gi(1,1:mp,1:L)) &
        +rightedge*dy*dz*sum(1.d0/gi(np,1:mp,1:L)) &
        +botedge*dx*dz*sum(1.d0/gi(1:np,1,1:L)) &
        +topedge*dx*dz*sum(1.d0/gi(1:np,mp,1:L))
    call mpi_allreduce(area,area,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
! outflow area and inflow area around model domain
    areaout=leftedge*dy*dz*sum(1.d0/gi(1,1:mp,1:L),xv(1,1:mp,1:L,1)<0.) &
        +rightedge*dy*dz*sum(1.d0/gi(np,1:mp,1:L),xv(np,1:mp,1:L,1)>0.) &
        +botedge*dx*dz*sum(1.d0/gi(1:np,1,1:L),xv(1:np,1,1:L,2)<0.) &
        +topedge*dx*dz*sum(1.d0/gi(1:np,mp,1:L),xv(1:np,mp,1:L,2)>0.)
    call mpi_allreduce(areaout,areaout,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    areain=area-areaout
! compute net mass flow out of box. mass flow out is positive by convention.
    massout=massout &
          -leftedge*dy*dz*sum((1.d0*xv(1,1:mp,1:L,1))/gi(1,1:mp,1:L)) &
          +rightedge*dy*dz*sum((1.d0*xv(np,1:mp,1:L,1))/gi(np,1:mp,1:L)) &
          -botedge*dx*dz*sum((1.d0*xv(1:np,1,1:L,2))/gi(1:np,1,1:L)) &
          +topedge*dx*dz*sum((1.d0*xv(1:np,mp,1:L,2))/gi(1:np,mp,1:L))
    call mpi_allreduce(massout,massout,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    
    iflg=0
    if(iflg==0)then !--------------------------------------------
! adjust outflow velocities with mean mass outflow per area out (=rho*V.n)
! so massout is zero. Inflow velocities are not adjusted.
      massperarea=massout/areaout
      if(leftedge==1)then
        where(xv(1,:,:,1)<0.); xv(1,:,:,1)=xv(1,:,:,1)+massperarea; endwhere
      endif
      if(rightedge==1)then
        where(xv(np,:,:,1)>0.); xv(np,:,:,1)=xv(np,:,:,1)-massperarea; endwhere
      endif
      if(botedge==1)then
        where(xv(:,1,:,2)<0.); xv(:,1,:,2)=xv(:,1,:,2)+massperarea; endwhere
      endif
      if(topedge==1)then
        where(xv(:,mp,:,2)>0.); xv(:,mp,:,2)=xv(:,mp,:,2)-massperarea; endwhere
      endif
    else
! adjust all boundary velocities so massout is zero. This method is not 
! recommended if either u or v velocity adjustment at the appropriate 
! boundary is not desired. For instance, if u=6 on x boundaries, v=0 on 
! y boundaries, and you want v to remain zero after 
! adjustment, use the above method instead. 
      massperarea=massout/area
      if(leftedge==1)xv(1,:,:,1)=xv(1,:,:,1)+massperarea
      if(rightedge==1)xv(np,:,:,1)=xv(np,:,:,1)-massperarea
      if(botedge==1)xv(:,1,:,2)=xv(:,1,:,2)+massperarea
      if(topedge==1)xv(:,mp,:,2)=xv(:,mp,:,2)-massperarea
    endif  !-----------------------------------------------------------
! diagnostic - massout should be zero (or <= floating point uncertainty 
! for original massout).
    massoutcheck=0.
    massoutcheck=massoutcheck &
          -leftedge*dy*dz*sum(1.d0*xv(1,1:mp,1:L,1)/gi(1,1:mp,1:L)) &
          +rightedge*dy*dz*sum(1.d0*xv(np,1:mp,1:L,1)/gi(np,1:mp,1:L)) &
          -botedge*dx*dz*sum(1.d0*xv(1:np,1,1:L,2)/gi(1:np,1,1:L)) &
          +topedge*dx*dz*sum(1.d0*xv(1:np,mp,1:L,2)/gi(1:np,mp,1:L))
    call mpi_allreduce(massoutcheck,massoutcheck,1, &
                       mpi_real8,mpi_sum,mpi_comm_world,ierror)
! successful adjustment of net mass flow is indicated when
! ratio of corrected massout over original massout is small ( < 1.0e-4)
  if(massout/=0.)then
    ratio=massoutcheck/massout
  else
    ratio=0.
  endif
      
! update xv array
      do kv=1,2
        call update(xv(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      enddo
  end subroutine vbcadlite
        
