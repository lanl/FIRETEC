subroutine forcing
  use gridlist_variables, only : l,ih
  use forcings, only : force,u1,u2,u3
  use xvall, only : xv,iuvel,iwvel
  use gridsetup, only : np,mp
  use workavg
  use higrad  
  Implicit None

  ! Local Variables 
  integer :: i,j,k,kv
 
  ! Executable Code
  u1(:,:,:)=u1(:,:,:)*0.5
  u2(:,:,:)=u2(:,:,:)*0.5
  u3(:,:,:)=u3(:,:,:)*0.5
  do kv=iuvel,iwvel
    call update(force(:,:,:,kv),force(:,:,:,kv),np,mp,l, &
      1-ih,np+ih,1-ih,mp+ih,0,0)
    call donorcell(force(:,:,:,kv),1-ih,np+ih,1-ih,mp+ih,l)
    xv(1:np,1:mp,1:l,kv)=xv(1:np,1:mp,1:l,kv)+force(1:np,1:mp,1:l,kv)
    call update(xv(:,:,:,kv),xv(:,:,:,kv),np,mp,l, &
      1-ih,np+ih,1-ih,mp+ih,0,0) 
  enddo

end subroutine forcing
