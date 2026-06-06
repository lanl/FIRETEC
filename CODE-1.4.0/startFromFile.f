
      subroutine startFromFile()
      use xvo
      use xvbin
      use xve
      use gridsetup
      use metryic
      use msga
      use restarta
      use pres
      use filesub
      use io
      use windfield
      use fireteca
      use turba
      use constants
      Implicit None
      integer :: i,j,k,kv!,ia,ja
      integer :: idot,iend
      if (irst.eq.1) then ! normal restart, read a comp.out
        idot=index(restartfile,'out.')+4
        if(idot.eq.4)then
          write(6,*) 'restart: incorrect restart name'
          call mpi_finalize()
          stop
        endif
        iend=len_trim(restartfile)
        read(restartfile(idot:iend),*)itrestart
        if(mpi_rank.eq.0)write(6,*)'restarting with ',itrestart
        call irstreadio(restartfile)
      else ! restart from a windfieldstart (irst.eq.0 and iwindfieldin.eq.1)
        itrestart=0
        if(mpi_rank.eq.0)write(6,*)'starting with windfieldstart :',windfieldstartfile
        call irstreadio2(windfieldstartfile)
      endif 

! below this line is restart in case of iwindfieldin (old comes from restartfile and new is read in xvbdata)
      if(iwindfieldin.eq.1)then
        itabsold=itrestart
!        call namefile(itabsold,xvbdataname,fxvbdataname)
!        call windfld_read(fxvbdataname,xvbdataold)
        do kv=1,nv
          do k=1,l
            do j=1,mp
              do i=1,np
                xvbdataold(i,j,k,kv)=xvb(i,j,k,kv)
                if(kv.eq.7) xvbdataold(i,j,k,kv)=0.233*xvb(i,j,k,nv)
                if(irhovapor.eq.1.and.kv.eq.8)
     &            xvbdataold(i,j,k,kv)=specifichumidity*xvb(i,j,k,nv)
              enddo
            enddo
          enddo
        enddo

          
        itabsnew=itabsold+itinterp
        ! here fp added a windspeedupfactor in order to entail to do wind runs
        ! with larger timestep that fire runs
        call namefile(itabsnew/windspeedupfactor,xvbdataname,fxvbdataname)
        if (mpi_rank.eq.0) write(6,*) 'reading file ',fxvbdataname
        call windfld_read(fxvbdataname,xvbdatanew)

!        if(mpi_rank.eq.0)then
!          write(6,*)'xvbdata:x...',xvbdataold(5,5,5,1),xe(5,5,5,1)
!     &,xvbdatanew(5,5,5,1)
!          write(6,*)'xvbdata:y...',xvbdataold(5,5,5,2),xe(5,5,5,2)
!     &,xvbdatanew(4,2,5,2)
!        endif
      endif
      return
      end 
