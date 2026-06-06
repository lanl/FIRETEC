      subroutine force
      use gridsetup
      use xvo
      use metryic
      use advo
      use workavg
      use msga
      use higrad  
      use xve       !this use allows access to the xe array

      Implicit None
   
      integer :: i,j,k,kv
      !FP09/2019 : KOO initially added these lines to make CODE3
      !similar to CODE, but it is incorrect, especially when cyclic BC are on
      ! because it uses 0 instead of cyclic values, so these lines were
      ! removed again and f1avg were used in updated below (instead of
      ! tmp)
      !real,allocatable::tmp(:,:,:,:)
      !allocate (tmp(1-ih:np+ih,1-ih:mp+ih,l,nv))
      !tmp=0.0
    
      !This array multiplied by a scalar can be written as a single line 
      u1(:,:,:)=u1(:,:,:)*0.5
      u2(:,:,:)=u2(:,:,:)*0.5
      u3(:,:,:)=u3(:,:,:)*0.5
      call updated(f1avg,f1avg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
      call donorcell(f1avg,1-ih,np+ih,1-ih,mp+ih,l)
      call updated(f2avg,f2avg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
      call donorcell(f2avg,1-ih,np+ih,1-ih,mp+ih,l)
      call updated(f3avg,f3avg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
      call donorcell(f3avg,1-ih,np+ih,1-ih,mp+ih,l)

      do k=1,l
       do j=1,mp
        do i=1,np
         xvb(i,j,k,1)=xvb(i,j,k,1)+f1avg(i,j,k)
         xvb(i,j,k,2)=xvb(i,j,k,2)+f2avg(i,j,k)
         xvb(i,j,k,3)=xvb(i,j,k,3)+f3avg(i,j,k)
         f1avg(i,j,k)=0.0
         f2avg(i,j,k)=0.0
         f3avg(i,j,k)=0.0
        enddo
       enddo
      enddo

      
      !JAS adding updates to keep halos in sync with real cells 7/7/06
      do kv=1,3
       call updated(xvb(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
      enddo

      !deallocate (tmp)
      return
      end
