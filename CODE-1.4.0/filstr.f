      subroutine filstr()
      use gridsetup
      use msga
      use xvo
      use xve

      Implicit None
      real,allocatable:: atmpx(:), atmpy(:), atmpz(:)
      integer :: i,j,k,kv
      real ::f1,f2

      f2=1.0/frqfilstr
      !f2=10.0/frqfilstr when filtering is every ten time step to speed
      !up filtering
      f1=1.0-f2
      allocate (atmpx(np))
      allocate (atmpy(mp))
      allocate (atmpz(l))

      do kv=1,3
        call updated(xvb(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
        do k=1,l
          do j=1,mp
            do i=1,np
              atmpx(i)=0.25*(xvb(i+1,j,k,kv)+2.*xvb(i,j,k,kv)+xvb(i-1,j,k,kv))
            enddo
            do i=1,np
              xvb(i,j,k,kv)=f1 * xvb(i,j,k,kv) + f2 * atmpx(i)
            enddo
          enddo
        enddo
        call updated(xvb(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
        do k=1,l
          do i=1,np
            do j=1,mp
              atmpy(j)=0.25*(xvb(i,j+1,k,kv)+2.*xvb(i,j,k,kv)+xvb(i,j-1,k,kv))
            enddo
            do j=1,mp
              xvb(i,j,k,kv)=f1 * xvb(i,j,k,kv) + f2 * atmpy(j)
            enddo
          enddo
        enddo
        call updated(xvb(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0)
        do j=1,mp
          do i=1,np
            do k=2,l-1
              atmpz(k)=0.25*(xvb(i,j,k+1,kv)+2.*xvb(i,j,k,kv)+xvb(i,j,k-1,kv))
            enddo
            !if (kz<=2)  !du/dz=0
            !   atmpz(1)=0.25*xvb(i,j,2,kv)+0.75*xvb(i,j,1,kv);
            !else  !w=0
            !   atmpz(1)=0.25*xvb(i,j,2,kv)-0.25*xvb(i,j,1,kv);
            !endif
            atmpz(1)=xvb(i,j,1,kv)
            atmpz(l)=xvb(i,j,l,kv)
            do k=1,l
              xvb(i,j,k,kv)=f1 * xvb(i,j,k,kv) + f2 * atmpz(k)
            enddo
          enddo
        enddo
      enddo
      deallocate (atmpx)
      deallocate (atmpy)
      deallocate (atmpz)
      return
      end
