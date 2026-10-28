!-----------------------------------------------------------------------
! rinitEmission initializes all arrays and variables related to 
! smoke, emissions, and plume evolution throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitHeatSource()
  use heatsource_gridlist_variables, only : MLRHRRinputs_hs,TME,MLR,HRR, &
    MLRHRR_eof_hs,inputName,nsource_hs
  use gridlist_variables, only : irst, wp
  use gridsetup, only : time
  use msga_variables, only : mpi_rank
  Implicit None

  ! Local Variables
  real(wp) :: told,mlrold,hrrold
  real(wp) :: tnew,mlrnew,hrrnew
  integer :: iostat,itsource
  integer :: unit_hs
  character(len=257) :: filename

  ! Executable Code
  if(any(MLRHRRinputs_hs.eq.1)) then
    allocate(TME(nsource_hs,2))
    allocate(MLR(nsource_hs,2))
    allocate(HRR(nsource_hs,2))
    allocate(MLRHRR_eof_hs(nsource_hs))

    TME = 0.0_wp
    MLR = 0.0_wp
    HRR = 0.0_wp
    MLRHRR_eof_hs = .false.

    do itsource = 1,nsource_hs

      if(MLRHRRinputs_hs(itsource).ne.1) cycle

      unit_hs = 7825 + itsource
      call namefile(itsource,inputName,filename)
      open(unit_hs,file=filename,form='formatted',status='old', &
        iostat=iostat)
      if(iostat.ne.0) then
        if(mpi_rank.eq.0) print*, &
          'ERROR opening MLR/HRR input file for heat source ',itsource
        stop
      endif

      read(unit_hs,*,iostat=iostat) ! headers
      if(iostat.ne.0) then
        if(mpi_rank.eq.0) print*, &
          'ERROR reading MLR/HRR header for heat source ',itsource
        stop
      endif

      MLRHRR_eof_hs(itsource) = .false.

      if(irst.eq.1) then

        ! Restart case:
        ! Scan forward until told <= time <= tnew, if possible.
        ! At file bounds, clamp by duplicating the first or last record.

        read(unit_hs,*,iostat=iostat) told, mlrold, hrrold
        if(iostat.ne.0) then
          if(mpi_rank.eq.0) print*, &
            'ERROR: no MLR/HRR data for heat source ',itsource
          stop
        endif

        if(time.le.told) then

          ! Restart time is before or exactly at the first data record.
          TME(itsource,1) = told
          MLR(itsource,1) = mlrold
          HRR(itsource,1) = hrrold

          TME(itsource,2) = told
          MLR(itsource,2) = mlrold
          HRR(itsource,2) = hrrold

        else

          read(unit_hs,*,iostat=iostat) tnew, mlrnew, hrrnew

          if(iostat.lt.0) then

            ! Only one data record exists. Clamp to it permanently.
            TME(itsource,1) = told
            MLR(itsource,1) = mlrold
            HRR(itsource,1) = hrrold

            TME(itsource,2) = told
            MLR(itsource,2) = mlrold
            HRR(itsource,2) = hrrold

            MLRHRR_eof_hs(itsource) = .true.

          elseif(iostat.gt.0) then

            if(mpi_rank.eq.0) print*, &
              'ERROR reading MLR/HRR data for heat source ',itsource
            stop

          else

            do while(tnew.lt.time)
              told   = tnew
              mlrold = mlrnew
              hrrold = hrrnew

              read(unit_hs,*,iostat=iostat) tnew, mlrnew, hrrnew

              if(iostat.lt.0) then

                ! Restart time is after the last record.
                ! Clamp to the last known record permanently.
                tnew   = told
                mlrnew = mlrold
                hrrnew = hrrold

                MLRHRR_eof_hs(itsource) = .true.
                exit

              elseif(iostat.gt.0) then

                if(mpi_rank.eq.0) print*, &
                  'ERROR reading MLR/HRR data for heat source ',itsource
                stop

              endif
            enddo

            TME(itsource,1) = told
            MLR(itsource,1) = mlrold
            HRR(itsource,1) = hrrold

            TME(itsource,2) = tnew
            MLR(itsource,2) = mlrnew
            HRR(itsource,2) = hrrnew

          endif

        endif

      else

        ! Non-restart case: start from the first two records.
        read(unit_hs,*,iostat=iostat) TME(itsource,1), &
          MLR(itsource,1),HRR(itsource,1)
        if(iostat.ne.0) then
          if(mpi_rank.eq.0) print*, &
            'ERROR: no MLR/HRR data for heat source ',itsource
          stop
        endif

        read(unit_hs,*,iostat=iostat) TME(itsource,2), &
          MLR(itsource,2),HRR(itsource,2)

        if(iostat.lt.0) then

          ! Only one data record exists. Clamp to it permanently.
          TME(itsource,2) = TME(itsource,1)
          MLR(itsource,2) = MLR(itsource,1)
          HRR(itsource,2) = HRR(itsource,1)

          MLRHRR_eof_hs(itsource) = .true.

        elseif(iostat.gt.0) then

          if(mpi_rank.eq.0) print*, &
            'ERROR reading MLR/HRR data for heat source ',itsource
          stop

        endif

      endif

      if(mpi_rank.eq.0) print*,'Reading MLR/HRR data for heat source ', &
        itsource

    enddo
  endif
end subroutine rinitHeatSource
