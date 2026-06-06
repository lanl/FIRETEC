      subroutine filstr(a,il,iu,jl,ju,lls)
      use gridsetup
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declaration to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls
      real,dimension(il:iu, jl:ju, lls) :: a
      real,allocatable:: sxa(:), sy(:), sz(:)
      integer :: i,j,k,illim,iulim,jllim,julim
      real :: temp

      allocate (sxa(np))
      allocate (sy(mp+1))
      allocate (sz(l))

      call updated(a,a,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      illim = 1  + 1*leftdedge                       !add d rrl
      iulim = np - 1*rightdedge                       !add d rrl

      do k=1,l
       do j=1,mp
        do i=illim,iulim
         sxa(i)=0.25*(a(i+1,j,k)+2.*a(i,j,k)+a(i-1,j,k))
        enddo
        if (leftdedge.eq.1) then                       !add d rrl
           sxa(1)=ibcx*0.25*(a(2,j,k)+2.*a(1,j,k)+a(-1,j,k))
     .          +(1-ibcx)*a(1,j,k)
        end if
        if (rightdedge.eq.1) then                       !add d rrl
c           !temp=ibcx*0.25*(a(np+2,j,k)+2.*a(np+1,j,k)+a(np-1,j,k))
c     .          +(1-ibcx)*a(np+1,j,k)
           temp=ibcx*0.25*(a(np+2,j,k)+2.*a(np,j,k)+a(np-1,j,k))
     .          +(1-ibcx)*a(np,j,k)       !FPRRL
           sxa(np)=ibcx*temp+(1-ibcx)*a(np,j,k)
        end if
        do i=1,np
         a(i,j,k)=sxa(i)
        enddo
       enddo
      enddo

      call updated(a,a,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      if(j3.eq.1) then
         jllim=1  + j3*botdedge                       !add d rrl
         julim=mp - j3*topdedge                       !add d rrl

         do k=1,l
            do i=1,np
               do j=jllim,julim
                  sy(j)=0.25*(a(i,j+j3,k)+2.*a(i,j,k)+a(i,j-j3,k))
               enddo
               if (botdedge.eq.1) then                       !add d rrl
                  sy(1)=ibcy*0.25*(a(i,1+j3,k)+2.*a(i,1,k)+a(i,-j3,k))
     .                 +(1-ibcy)*a(i,1,k)
               end if
               if (topdedge.eq.1) then                       !add d rrl
c                  temp=ibcy*0.25*(a(i,mp+1+j3,k)+2.*a(i,mp+1,k)+
c     .                 a(i,mp-j3,k)) +(1-ibcy)*a(i,mp+1,k)
                  temp=ibcy*0.25*(a(i,mp+1+j3,k)+2.*a(i,mp,k)+
     .                 a(i,mp-j3,k)) +(1-ibcy)*a(i,mp,k)   !FPRRL
                  sy(mp)=ibcy*temp+(1-ibcy)*a(i,mp,k)
               end if
               do j=1,mp
                  a(i,j,k)=sy(j)
               enddo
            enddo
         enddo
      endif

      call updated(a,a,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      do j=1,mp
         do i=1,np
            do k=2,l-1
               sz(k)=0.25*(a(i,j,k+1)+2.*a(i,j,k)+a(i,j,k-1))
            enddo
            sz(1)=a(i,j,1)
            sz(l)=a(i,j,l)
            do k=1,l
               a(i,j,k)=sz(k)
            enddo
         enddo
      enddo

      call updated(a,a,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      deallocate (sxa)
      deallocate (sy)
      deallocate (sz)

      return
      end
