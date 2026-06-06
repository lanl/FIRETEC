!23456789112345678921234567893123456789412345678951234567896123456789712
c  JMC:5/17/5
c  This subroutine reads all information for a slice of the grid and
c  writes it to a file for later use.

      subroutine icfmewind(xvball,ilow,ihigh,jlow,jhigh,lhigh,
     &                     nvicfme)

      use gridsetup
      use msga
     
      Implicit None

      
      integer,intent(in) :: ilow,ihigh,jlow,jhigh,lhigh,nvicfme
      real :: xvball(ilow:ihigh,jlow:jhigh,lhigh,nvicfme)
      real,allocatable :: tmparray(:,:,:,:,:),xvbwind(:,:,:)
      integer :: nvpw,ipw,jpw,kpw

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: ia,ja,iproc,iprocx,jprocy,nsize

      if (mpi_rank.eq.0)
     & allocate (tmparray(ilow:ihigh,jlow:jhigh,lhigh,nvicfme,nproc),
     &           xvbwind(m,lhigh,nvicfme))
      nsize=(ihigh-ilow+1)*(jhigh-jlow+1)*lhigh*nvicfme
      call mpi_gather(xvball,nsize,mpi_real,tmparray,nsize,mpi_real,
     &     0,mpi_comm_world,ierror)

      if (mpi_rank.eq.0) then
        do iprocx=1,nprocx
          do jprocy=1,nprocy
            iproc=1+(iprocx-1)+(jprocy-1)*nprocx
            do nvpw=1,nvicfme
              do kpw=1,lhigh
                do jpw=1,mp
                  do ipw=1,np
                    ia=(iprocx-1)*np + ipw
                    ja=(jprocy-1)*mp + jpw

                    if(ia.eq.100)then
c                     write(*,*)'JMC i= ',ipw
c                     write(*,*)'JMC ja= ',ja
                      xvbwind(ja,kpw,nvpw)=
     &                         tmparray(ipw,jpw,kpw,nvpw,iproc)
                    endif

                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
 
        open (701, file='icfmeallprof.dat',form='unformatted',
     &        status='unknown')
        write(701) xvbwind
        close (701)
 
        deallocate (tmparray)
        deallocate (xvbwind)
      endif
      
      end
