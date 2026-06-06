      subroutine turb()
      use xvo
      use gridsetup
      use metryic
      use updatedFields
      use turba
      use turbb
      use workavg
      use msga
      use fireteca,only:foxb,frhovaporb
      Implicit None
     
      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer :: i,j,k
      real :: dtm

      f1avg=0.
      f2avg=0.
      f3avg=0.
      fib=0.
      fka=0.
      fkb=0.
      foxb=0.
      if (irhovapor.eq.1) frhovaporb=0.


      dtm=1.0*dt
! 1/ compute forces on momentum
      !stressrij, requires fieldUpdate call (in compress.f)      
      !if (mpi_rank.eq.0) write(6,*) 'stressrij'

      call stressrij(f1avg,f2avg,f3avg,1-ih,np+ih,1-ih,mp+ih,l)

      ! compute drag contribution
      !if (mpi_rank.eq.0) write(6,*) 'drag'
      if (idrag.eq.1) then
          call dragm(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)
      else if (idrag.eq.2) then
          call dragm2(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)
      end if
     
! 2/ compute forces for ka and kb
      
! 2.1/compute ka terms, requires fieldUpdate call (in compress.f)
      !if (mpi_rank.eq.0) write(6,*) 'ka terms'
      ! shear production
      call rijgradu(fka,K_axy,K_az,xvb(:,:,:,5),1-ih,np+ih,1-ih,mp+ih,l)
      ! turbulent diffusion
      call diffuse(fka,tke_a,1-ih,np+ih,1-ih,mp+ih,l,0)
      ! turbulence dissipation+vegetation effect
      call draga()
! 2.2/compute kb terms, requires fieldUpdate call (in compress.f)
      if (iturb.eq.2) then
      !if (mpi_rank.eq.0) write(6,*) 'kb terms'
        ! shear production
        call rijgradu(fkb,K_b,K_b,xvb(:,:,:,6),1-ih,np+ih,1-ih,mp+ih,l)
        ! turbulent diffusion
        call diffuse(fkb,tke_b,1-ih,np+ih,1-ih,mp+ih,l,0)
        ! turbulence dissipation+vegetation effect
        call dragb()
      endif
     
! 3/ compute diffusion of other fields
      call diffuse(fib,theta,1-ih,np+ih,1-ih,mp+ih,l,1)
! compute diffusion of o2 (moved by fp from firetec)
      if (ilapdo.eq.1)
     . call diffuse(foxb,ox,1-ih,np+ih,1-ih,mp+ih,l,1)
      ! compute diffusion of rhovapor (moved by fp from firetec)
      if (irhovapor.eq.1)
     . call diffuse(frhovaporb,vap,1-ih,np+ih,1-ih,mp+ih,l,1)
! prepare contrubution of these terms to rhs
      if(irod.eq.1) then
        do k=1,l
          do j=1,mp
            do i=1,np
              f1avg(i,j,k)=f1avg(i,j,k)*dti
              f2avg(i,j,k)=f2avg(i,j,k)*dti
              f3avg(i,j,k)=f3avg(i,j,k)*dti
              fka(i,j,k)=fka(i,j,k)*dti
              fkb(i,j,k)=fkb(i,j,k)*dti
              fib(i,j,k)=fib(i,j,k)*dti
            enddo
          enddo
         enddo
         if (ilapdo.eq.1)
     .     foxb=foxb * dti
         if (irhovapor.eq.1)
     .     frhovaporb=frhovaporb * dti
      else  ! irod.eq.0
        do k=1,l
          do j=1,mp
            do i=1,np
              xvb(i,j,k,1)=xvb(i,j,k,1)+0.5*f1avg(i,j,k)
              xvb(i,j,k,2)=xvb(i,j,k,2)+0.5*f2avg(i,j,k)
              xvb(i,j,k,3)=xvb(i,j,k,3)+0.5*f3avg(i,j,k)
c              xvb(i,j,k,4)=xvb(i,j,k,4)+0.5*fib(i,j,k)
              xvb(i,j,k,5)=xvb(i,j,k,5)+0.5*fka(i,j,k)
              xvb(i,j,k,6)=xvb(i,j,k,6)+0.5*fkb(i,j,k)
             enddo
           enddo
         enddo
      endif 

150   format(5x,e30.15)

      return
      end 
