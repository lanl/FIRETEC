! xevariation takes into account any evolution of xe 

subroutine xevariation()
  use gridlist_variables, only : n,m,l,ih,u0,uramp,uramptime,uswitch,v0,vramp,vramptime, &
    vswitch,zu,iwindfieldin,windspeedupfactor,itinterp,ibcells,jbcells,xvdataname
  use xvall, only : xvdataold,xvdatanew,xe,iuvel,ivvel,irho,nv
  use fuel_variables, only : rhoFuel
  use metric_variables, only : z
  use gridsetup, only : ittot,np,mp,time
  use msga_variables, only : mpi_rank,npos,mpos
  use xvbin     
  Implicit None

  ! Local Variables
  real,external :: zcart
  integer :: i,j,k,ia,ja
  integer :: kv

  ! Executable Code
  if(iwindfieldin.eq.1) then
    if(mod(ittot,itinterp).eq.0) then
      xvbdataold = xvbdatanew
      itabsold = itabsnew
      itabsnew=itabsold+itinterp
      call namefile(itabsnew/windspeedupfactor,xvdataname,fxvbdataname)
      if (mpi_rank.eq.0) write(6,*) 'reading file ',fxvbdataname
      call windfld_read(fxvbdataname,xvbdatanew)
    endif

  ! here we update xe 
    do k=1,l
      do j=1,mp
        do i=1,np
          ia=(npos-1)*np+i
          ja=(mpos-1)*mp+j
          if(ia.le.ibcells.or.ia.gt.n-ibcells &
            .or.ja.le.jbcells.or.ja.gt.m-jbcells)then
            do kv=1,nv
              xe(i,j,k,kv)= &
                lin_interp(real(ittot),real(itabsnew),real(itabsold) &
                ,xvbdatanew(i,j,k,kv),xvbdataold(i,j,k,kv))
            enddo
          endif
        enddo
      enddo
    enddo
    if(mpi_rank.eq.0)then
      write(6,*)'told,tot,new',real(itabsold),real(ittot),real(itabsnew)
      write(6,*)'xvbdata:x...',xvbdataold(1,1,1,1),xe(1,1,1,1) &
                 ,xvbdatanew(1,1,1,1)
      write(6,*)'xvbdata:y...',xvbdataold(1,1,1,2),xe(1,1,1,2), &
                xvbdatanew(1,1,1,2)
    endif
    do kv=1,nv
      call update(xe(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
    enddo
  else ! iwindfield.ne.1
    if(uswitch==1)then
      do k=1,l
        do j=1,mp
          do i=1,np
            xe(i,j,k,1)=min(uramp*xe(i,j,k,nv), &
              (time*uramp/uramptime+u0)*xe(i,j,k,nv)) &
              *max(0.,(1.-1.5*sum(rhoFuel(:,i,j,k))))  &
              *(min(100.,zcart(z(k),1,1))/zu)**(1./7.)
          enddo
        enddo
      enddo
    endif
    if(vswitch==1)then
      do k=1,l
        do j=1,mp
          do i=1,np
            xe(i,j,k,2)=min(vramp*xe(i,j,k,nv), &
              (time*vramp/vramptime+v0)*xe(i,j,k,nv)) &
              *max(0.,(1.-1.5*sum(rhoFuel(:,i,j,k)))) &
              *(min(100.,zcart(z(k),1,1))/zu)**(1./7.)
          enddo
        enddo
      enddo
    endif
    do kv=1,nv
      call update(xe(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
    enddo
  endif
end subroutine xevariation
