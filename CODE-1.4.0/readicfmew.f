!23456789112345678921234567893123456789412345678951234567896123456789712
c  JMC:5/18/5
c  This subroutine reads all information from a slice of the grid and
c  places it in array xe.

      subroutine readicfmew(array,ilow,ihigh,jlow,jhigh,lhigh,nvicfme)

      use gridsetup
      use msga

      Implicit None

      
      integer,intent(in) :: ilow,ihigh,jlow,jhigh,lhigh,nvicfme
      real :: array(ilow:ihigh,jlow:jhigh,lhigh,nvicfme)
      real,allocatable :: tmparray(:,:,:,:,:),xvbwind(:,:,:)
      integer :: iprocx,jprocy,iproc,ia,ja,nsize,nvpw,ipw,jpw,kpw

      if (mpi_rank.eq.0) then
        allocate (tmparray(ilow:ihigh,jlow:jhigh,lhigh,nvicfme,nproc),
     &           xvbwind(m,lhigh,nvicfme))

        open (701, file='icfmeallprof.dat',form='unformatted',
     &        status='old')
        read(701) xvbwind
        close (701)
 

        do iprocx=1,nprocx
          do jprocy=1,nprocy
            iproc=1+(iprocx-1)+(jprocy-1)*nprocx
            do nvpw=1,nvicfme
              do kpw=1,lhigh
                do jpw=1,mp
                  do ipw=1,np
                    ia=(iprocx-1)*np + ipw
                    ja=(jprocy-1)*mp + jpw

                    tmparray(ipw,jpw,kpw,nvpw,iproc)=
     &                          xvbwind(ja,kpw,nvpw)

                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
 
      nsize=(ihigh-ilow+1)*(jhigh-jlow+1)*lhigh*nvicfme
      call MPI_Scatter(tmparray,nsize,mpi_real,array,nsize,mpi_real,
     &     0,mpi_comm_world,ierror)
 
      if (mpi_rank.eq.0)
     &  deallocate (tmparray,xvbwind)

      return 
      end
