      subroutine force
      use gridsetup
      use xvo
      use metryic
      use advo
      use workavg
      use msga
      use higrad    !JAS allows usage of the donorcell routine
      use xve       !this use allows access to the xe array

      Implicit None
   
      integer :: i,j,k,kv
      real,allocatable::tmp(:,:,:,:)
      allocate (tmp(1-ih:np+ih,1-ih:mp+ih,l,nv))
      tmp=0.0
    
      !This array multiplied by a scalar can be written as a single line 
      uab(:,:,:)=uab(:,:,:)*0.5
      vab(:,:,:)=vab(:,:,:)*0.5
      oab(:,:,:)=oab(:,:,:)*0.5

      call donorcell(uab,vab,oab,f1avg,tmp,1-ih,np+ih,1-ih,mp+ih,l)
      call donorcell(uab,vab,oab,f2avg,tmp,1-ih,np+ih,1-ih,mp+ih,l)
      call donorcell(uab,vab,oab,f3avg,tmp,1-ih,np+ih,1-ih,mp+ih,l)

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
       call updated(xvb(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      enddo

      deallocate (tmp)
      return
      end
