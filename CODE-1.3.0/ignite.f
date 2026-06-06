c2345678***************************************************
      module ignite
      use msga
      Implicit None

        integer::igntype !flag for ignition type 0=none,1=rinitfire,2=terratorch,
                         !  3=multiple rinitfire,4=aerial,5=atv/driptorch
        integer::numignpts !number of points in ignition file
        integer,allocatable::ignpts(:,:)
        real::targettemp,startigntime,endigntime,ramprate
              !targettemp = target temperature for ignition
               !startigntime = igntime (s) at which to start ignition
                !ramprate = rate (k/s) at which to increase temperature
        real:: xfirelinelow,xfirelinehigh,yfirelinelow,yfirelinehigh
        real::flamedistance
       !JLD twolinesignition !
       real::startigntime1,endigntime1,startigntime2,endigntime2
       real::xlineor1,xlineend1,ylineor1,ylineend1
       real::xlineor2,xlineend2,ylineor2,ylineend2
        character(len=60)::ignfile='ignite.dat'    !name of ignition parameters file
        character(len=60) :: buffer

c  JLW 02/20/13 aerial ignition
        integer,allocatable::ignloc(:,:)
        real,allocatable ::igntime(:)
        integer:: naerial

c  JLW 02/20/13 atv/driptorch ignition
        real,allocatable ::atvstart(:,:),atvstop(:,:),atvtime(:,:)
        integer:: natv
       save

       contains
       !***************************************************!
       subroutine ign_setup()
       Implicit None

        integer::i,ierror
        namelist/rinitlist/ numignpts,targettemp,startigntime,ramprate
        namelist/torchlist/ targettemp,startigntime,endigntime, 
     &      xfirelinelow,xfirelinehigh,
     &      yfirelinelow,yfirelinehigh,flamedistance
        namelist/aeriallist/ naerial,targettemp,ramprate
        namelist/atvlist/ natv,targettemp,flamedistance
        namelist/twolineslist/ targettemp,flamedistance,
     &       startigntime1,endigntime1,
     &      xlineor1,xlineend1,ylineor1,ylineend1,
     &      startigntime2,endigntime2,
     &      xlineor2,xlineend2,ylineor2,ylineend2

        if(mpi_rank==0)then
         open(unit=1000,file=ignfile,form='formatted',status='old')
         read(1000,'(A15,I5)') buffer,igntype

         if(igntype==1)then  !read in a rinitfire style ignition file !
           read (1000,nml=rinitlist)
           allocate(ignpts(numignpts,2))
           do i=1,numignpts
            read(1000,*) ignpts(i,1),ignpts(i,2)
           enddo
           !write(6,nml=rinitlist)
           !write(6,*) ignpts
         elseif(igntype==2.or.igntype==7)then !read in a terratorch style ignition file !
           read (1000,nml=torchlist)
           !write(6,nml=torchlist)
         elseif(igntype==4)then !read in an aerial ignition style ignition file !
           read (1000,nml=aeriallist)
           allocate(ignloc(naerial,2))
           allocate(igntime(naerial))
           do i=1,naerial
             read(1000,*) ignloc(i,1),ignloc(i,2),igntime(i)
           enddo
         elseif(igntype==5)then !read in a atv/driptorch ignition file !
           read (1000,nml=atvlist)
           allocate(atvstart(natv,2))
           allocate(atvstop(natv,2))
           allocate(atvtime(natv,2))
           do i=1,natv
             read(1000,*) atvstart(i,1),atvstart(i,2),
     &                    atvstop(i,1),atvstop(i,2),
     &                    atvtime(i,1),atvtime(i,2)
           enddo
         elseif(igntype==6)then !JLD read in a twolines  style ignition file
            read (1000,nml=twolineslist)
            startigntime=min(startigntime1,startigntime2)  
         endif
         close(1000)
        endif     !end if(mpi_rank == 0) 
        
       call MPI_Bcast(igntype,1,mpi_integer,0,mpi_comm_world,ierror)
       if(igntype==1)then
        call MPI_Bcast(numignpts,1,mpi_integer,0,mpi_comm_world,ierror)
        if(mpi_rank.ne.0)then
         allocate(ignpts(numignpts,2))
        endif     
        call MPI_Bcast(targettemp,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(startigntime,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(ramprate,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(ignpts,2*numignpts,mpi_integer,0,
     &   mpi_comm_world,ierror)
        
       elseif(igntype==2.or.igntype==7)then
        call MPI_Bcast(targettemp,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(startigntime,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(endigntime,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(xfirelinelow,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(xfirelinehigh,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(yfirelinelow,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(yfirelinehigh,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(flamedistance,1,mpi_real,0,mpi_comm_world,ierror)

       elseif(igntype==4)then     ! JLW 02/20/13 aerial ignition
        call MPI_Bcast(targettemp,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(ramprate,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(naerial,1,mpi_integer,0,mpi_comm_world,ierror)
        if(mpi_rank.ne.0)then
          allocate(ignloc(naerial,2))
          allocate(igntime(naerial))
        endif
        call MPI_Bcast(ignloc,2*naerial,mpi_integer,0,
     &                 mpi_comm_world,ierror)
        call MPI_Bcast(igntime,naerial,mpi_real,0,
     &                 mpi_comm_world,ierror)
        if (mpi_rank.eq.0) print *,'ignite.f ',ignloc(naerial,1),
     &                          ignloc(naerial,2),igntime(naerial)

       elseif(igntype==5)then     ! JLW 02/20/13 atv/driptorch ignition
        call MPI_Bcast(targettemp,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(flamedistance,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(natv,1,mpi_integer,0,mpi_comm_world,ierror)
        if(mpi_rank.ne.0)then
          allocate(atvstart(natv,2))
          allocate(atvstop(natv,2))
          allocate(atvtime(natv,2))
        endif
        call MPI_Bcast(atvstart,2*natv,mpi_real,0,
     &                 mpi_comm_world,ierror)
        call MPI_Bcast(atvstop,2*natv,mpi_real,0,
     &                 mpi_comm_world,ierror)
        call MPI_Bcast(atvtime,2*natv,mpi_real,0,
     &                 mpi_comm_world,ierror)
        if (mpi_rank.eq.0) print *,'ignite.f ',atvstart(natv,1),
     &                        atvstart(natv,2),atvtime(natv,1)
        elseif(igntype==6)then
        call MPI_Bcast(targettemp,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(flamedistance,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(startigntime,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(startigntime1,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(endigntime1,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(xlineor1,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(xlineend1,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(ylineor1,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(ylineend1,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(startigntime2,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(endigntime2,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(xlineor2,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(xlineend2,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(ylineor2,1,mpi_real,0,mpi_comm_world,ierror)
        call MPI_Bcast(ylineend2,1,mpi_real,0,mpi_comm_world,ierror)
       endif

       end subroutine ign_setup

       !***************************************************!
       subroutine ign_cleanup()
        Implicit None

       if(igntype==1)then
        deallocate(ignpts)
       else if (igntype==4) then
! DO NOT DEALLOCATE - ARRAYS MAY BE NEEDED THOUGHOUT RUN
       else if (igntype==5) then
! DO NOT DEALLOCATE - ARRAYS MAY BE NEEDED THOUGHOUT RUN
       endif

       end subroutine ign_cleanup
       !***************************************************!

      end module ignite
c2345678***************************************************
