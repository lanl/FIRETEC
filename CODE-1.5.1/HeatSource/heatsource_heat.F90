!-----------------------------------------------------------------------
! updates position of a circular heat source
!-----------------------------------------------------------------------
subroutine updateHeatSourcePosition(iheatsource,il,iu,jl,ju)
  use gridlist_variables, only : wp,dx,dy
  use heatsource_gridlist_variables, only : nsource_hs,starttime_hs, &
  endtime_hs,radius_hs,xcen_hs,ycen_hs,xl_hs,xu_hs,yl_hs,yu_hs
  use gridsetup, only : np,mp,time
  use msga_variables, only : npos,mpos
  Implicit none
 
  ! Local Variables
  integer,intent(in) :: iheatsource
  integer,intent(in) :: il,iu,jl,ju

  integer :: i,j,ia,ja,itsource
  real(wp) :: xlo,xhi,ylo,yhi,overlap_area=0._wp
  real(wp),external :: hs_rect_area,hs_circle_rect_area

  ! Executable Code
  ! here we update value of rhof in the heat source 
  ! FIXME :: No temporal translation of position at this time AJJ (15 July 24)
  if(iheatsource.eq.1) then ! Circular
    do j=jl,ju
      ja=(mpos-1)*mp+j
      do i=il,iu
        ia=(npos-1)*np+i
        do itsource=1,nsource_hs
          if(time.lt.starttime_hs(itsource) .or. &
            time.ge.endtime_hs(itsource)) cycle
          xlo = (ia-1)*dx
          xhi = ia*dx
          ylo = (ja-1)*dy
          yhi = ja*dy
          overlap_area = hs_circle_rect_area(xlo,xhi,ylo,yhi, &
            xcen_hs(itsource),ycen_hs(itsource),radius_hs(itsource))
        enddo
      enddo
    enddo 
  elseif(iheatsource.eq.2) then ! Rectangular
    do j=jl,ju
      ja=(mpos-1)*mp+j
      do i=il,iu
        ia=(npos-1)*np+i
        do itsource=1,nsource_hs
          if(time.lt.starttime_hs(itsource) .or. &
            time.ge.endtime_hs(itsource)) cycle
          xlo = (ia-1)*dx
          xhi = ia*dx
          ylo = (ja-1)*dy
          yhi = ja*dy
          overlap_area = hs_rect_area(xlo,xhi,ylo,yhi, &
            xl_hs(itsource),xu_hs(itsource), &
            yl_hs(itsource),yu_hs(itsource))
        enddo
      enddo
    enddo
  endif 

end subroutine updateHeatSourcePosition

!-----------------------------------------------------------------------
! heatSource calculates forcing terms for artificial heat as specified 
! by the user
!-----------------------------------------------------------------------
subroutine heatSource(iheatsource,force,il,iu,jl,ju,lls,nv)
  use heatsource_gridlist_variables
  use gridlist_variables, only : dx,dy,wp,ih
  use gridsetup, only : time
  use constants, only : pi,pref
  use thermo_variables, only : pr,cp_gas,cv_gas
  use gridsetup, only : np,mp
  use metric_variables, only : zedge,zs
  use msga_variables, only : npos,mpos,mpi_rank,mpi_comm_world,ierror,mpi_sum
  use msga_variables, only : mpi_wp
  use xvall, only : irho,itheta
  use fuel_variables, only : rnfuel,hf
  use constants, only : tolerance
  Implicit None

  ! Local Variables
  integer,intent(in) :: iheatsource
  integer,intent(in) :: il,iu,jl,ju,lls,nv
  real(wp),intent(inout) :: force(1-ih:np+ih,1-ih:mp+ih,lls,nv)

  integer :: i,j,k,itsource,iostat
  integer :: ia,ja
  real(wp) :: source_function=0.0_wp,source_time
  real(wp) :: energy,fi,frho
  real(wp) :: volume,source_volume,hrr_now=0._wp
  real(wp) :: xlo,xhi,ylo,yhi,zlo,zhi,zlo_agl,zhi_agl 
  real(wp) :: overlap_area=0._wp
  real(wp) :: overlap_height,overlap_volume,coverfrac
  real(wp) :: source_volume_local
  real(wp),external :: zcart,hs_rect_area,hs_circle_rect_area
  logical, save, allocatable :: heatsource_logged(:)
  logical, save, allocatable :: heatsource_off_logged(:)
  logical :: heatsource_active

  ! Executable Code
  if(.not.allocated(heatsource_logged)) then
    allocate(heatsource_logged(nsource_hs))
    heatsource_logged = .false.
  endif
  if(.not.allocated(heatsource_off_logged)) then
    allocate(heatsource_off_logged(nsource_hs))
    heatsource_off_logged = .false.
  endif

  do itsource=1,nsource_hs

    heatsource_active = time.ge.starttime_hs(itsource) .and. &
      time.lt.endtime_hs(itsource)

    if(heatsource_active) then

      if((.not.heatsource_logged(itsource)).and. &
        (iheatsource.eq.1.or.iheatsource.eq.2)) then
        if(mpi_rank.eq.0) then
          print*,'Heat source activated: itsource = ',itsource, &
            ' time = ',time,' iheatsource = ',iheatsource
        endif
        heatsource_logged(itsource) = .true.
      endif

    else

      if(time.ge.endtime_hs(itsource) .and. &
        heatsource_logged(itsource) .and. &
        (.not.heatsource_off_logged(itsource)) .and. &
        (iheatsource.eq.1.or.iheatsource.eq.2)) then
        if(mpi_rank.eq.0) then
          print*,'Heat source deactivated: itsource = ',itsource, &
            ' time = ',time,' iheatsource = ',iheatsource
        endif
        heatsource_off_logged(itsource) = .true.
      endif

      cycle

    endif

    source_function=0.0_wp
    source_time=time-starttime_hs(itsource)
    select case(isource_function_hs(itsource))
      case(0)
        source_function=1.0_wp
      case(1)
        if(sin(2.0_wp*pi*freq_hs(itsource)*source_time).ge.0.0_wp)then
          source_function = 1.0_wp
        endif
      case(2)
        source_function = 0.5_wp* &
          (1.0_wp+sin(2.0_wp*pi*freq_hs(itsource)*source_time))
      case(3)
        if((freq_hs(itsource).gt.0.0_wp).and. &
          (freq_hs(itsource)*source_time.le.1.0_wp))then
          source_function = 1.0_wp
        endif
      case default
        if(mpi_rank.eq.0) print*, &
          'ERROR: unsupported isource_function_hs for heat source ', &
          itsource,isource_function_hs(itsource)
        STOP
    end select

    if(MLRHRRinputs_hs(itsource).eq.1) then
      do while((.not.MLRHRR_eof_hs(itsource)).and. &
        TME(itsource,2).lt.time)

        TME(itsource,1) = TME(itsource,2)
        MLR(itsource,1) = MLR(itsource,2)
        HRR(itsource,1) = HRR(itsource,2)

        read(7825+itsource,*,iostat=iostat) &
          TME(itsource,2),MLR(itsource,2),HRR(itsource,2)

        if(iostat.lt.0) then ! iostat<0 means end of file

          TME(itsource,2) = TME(itsource,1)
          MLR(itsource,2) = MLR(itsource,1)
          HRR(itsource,2) = HRR(itsource,1)

          MLRHRR_eof_hs(itsource) = .true.

        elseif(iostat.gt.0) then

          if(mpi_rank.eq.0) print*, &
            'ERROR reading MLR/HRR data for heat source ',itsource
          stop

        else

          if(mpi_rank.eq.0) print*,'Updating MLR/HRR data'

        endif

      enddo

      if(time.le.TME(itsource,1)) then
        fm_hs(itsource) = MLR(itsource,1)
      elseif(abs(TME(itsource,2)-TME(itsource,1)).lt.tolerance) then
        fm_hs(itsource) = MLR(itsource,2)
      else
        fm_hs(itsource) = MLR(itsource,1) + &
          (time-TME(itsource,1))/(TME(itsource,2)-TME(itsource,1))* &
          (MLR(itsource,2)-MLR(itsource,1))
      endif
    endif

    if(MLRHRRinputs_hs(itsource).eq.1) then
      if(time.le.TME(itsource,1)) then
        hrr_now = HRR(itsource,1)
      elseif(abs(TME(itsource,2)-TME(itsource,1)).lt.tolerance) then
        hrr_now = HRR(itsource,2)
      else
        hrr_now = HRR(itsource,1) + &
          (time-TME(itsource,1))/(TME(itsource,2)-TME(itsource,1))* &
          (HRR(itsource,2)-HRR(itsource,1))
      endif
      source_volume_local = 0.0_wp
      do k=1,lls
        do j=jl,ju
          do i=il,iu
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            zlo = zcart(zedge(k),i,j)
            zhi = zcart(zedge(k+1),i,j)
            zlo_agl = zlo - zs(i,j)
            zhi_agl = zhi - zs(i,j)
            overlap_height = max(0.0_wp, &
              min(zhi_agl,depth_hs(itsource))-zlo_agl)
            if(overlap_height.le.0.0_wp) cycle
            if(iheatsource.eq.1) then
              xlo = (ia-1)*dx
              xhi = ia*dx
              ylo = (ja-1)*dy
              yhi = ja*dy
              overlap_area = hs_circle_rect_area(xlo,xhi,ylo,yhi, &
                xcen_hs(itsource),ycen_hs(itsource), &
                radius_hs(itsource))
            elseif(iheatsource.eq.2) then
              xlo = (ia-1)*dx
              xhi = ia*dx
              ylo = (ja-1)*dy
              yhi = ja*dy
              overlap_area = hs_rect_area(xlo,xhi,ylo,yhi, &
                xl_hs(itsource),xu_hs(itsource), &
                yl_hs(itsource),yu_hs(itsource))
            endif
            source_volume_local = source_volume_local + &
              overlap_area*overlap_height
          enddo
        enddo
      enddo
      call mpi_allreduce(source_volume_local,source_volume,1,mpi_wp, &
        mpi_sum,mpi_comm_world,ierror)
    endif

    if(iheatsource.eq.1)then !Circular source
      do k=1,lls
        do j=jl,ju
          do i=il,iu
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            xlo = (ia-1)*dx
            xhi = ia*dx
            ylo = (ja-1)*dy
            yhi = ja*dy
            zlo = zcart(zedge(k),i,j)
            zhi = zcart(zedge(k+1),i,j)
            zlo_agl = zlo - zs(i,j)
            zhi_agl = zhi - zs(i,j)
            overlap_area = hs_circle_rect_area(xlo,xhi,ylo,yhi, &
              xcen_hs(itsource),ycen_hs(itsource),radius_hs(itsource))
            overlap_height = max(0.0_wp, &
              min(zhi_agl,depth_hs(itsource))-zlo_agl)
            volume=dx*dy*(zhi-zlo)
            overlap_volume = overlap_area*overlap_height
            if(overlap_volume.gt.0.0_wp.and.volume.gt.0.0_wp) then
                coverfrac = overlap_volume/volume
                if(MLRHRRinputs_hs(itsource).eq.1) then
                  if(source_volume.le.0.0_wp) cycle
                  energy = hrr_now/source_volume
                else
                  energy= &
                    av_hs(itsource)*flux_hs(itsource)*source_function
                endif
                fi=coverfrac*energy/cp_gas(i,j,k)* &
                  (pref/pr(i,j,k))**(1.0_wp-cv_gas(i,j,k)/cp_gas(i,j,k))
                ! computation of energy really realised in gas phase
                force(i,j,k,itheta)=force(i,j,k,itheta)+fi

                !mass source from energy:
                if(imass_source_hs(itsource).eq.1)then
                  frho=coverfrac*energy/(hf/rnfuel)
                  force(i,j,k,irho)=force(i,j,k,irho)+frho
                elseif(imass_source_hs(itsource).eq.2)then
                  if(MLRHRRinputs_hs(itsource).eq.1) then
                    force(i,j,k,irho)=force(i,j,k,irho)+ &
                      coverfrac*fm_hs(itsource)/source_volume
                  else
                    force(i,j,k,irho)=force(i,j,k,irho)+ &
                      coverfrac*fm_hs(itsource)
                  endif
                endif !(imass_source_hs option)
            endif
          enddo
        enddo
      enddo
    elseif(iheatsource.eq.2)then !Rectangular source
      do k=1,lls
        do j=jl,ju
          do i=il,iu
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            xlo = (ia-1)*dx
            xhi = ia*dx
            ylo = (ja-1)*dy
            yhi = ja*dy
            zlo = zcart(zedge(k),i,j)
            zhi = zcart(zedge(k+1),i,j)
            zlo_agl = zlo - zs(i,j)
            zhi_agl = zhi - zs(i,j)
            overlap_area = hs_rect_area(xlo,xhi,ylo,yhi, &
              xl_hs(itsource),xu_hs(itsource), &
              yl_hs(itsource),yu_hs(itsource))
            overlap_height = max(0.0_wp, &
              min(zhi_agl,depth_hs(itsource))-zlo_agl)
            volume = dx*dy*(zhi-zlo)
            overlap_volume = overlap_area*overlap_height
            if(overlap_volume.gt.0.0_wp.and.volume.gt.0.0_wp)then
              coverfrac = overlap_volume/volume

              if(MLRHRRinputs_hs(itsource).eq.1) then
                if(source_volume.le.0.0_wp) cycle
                energy = hrr_now/source_volume
              else
                energy=av_hs(itsource)*flux_hs(itsource)*source_function
              endif
              fi=coverfrac*energy/cp_gas(i,j,k)* &
                (pref/pr(i,j,k))**(1.0_wp-cv_gas(i,j,k)/cp_gas(i,j,k))
              ! computation of energy really realised in gas phase
              force(i,j,k,itheta)=force(i,j,k,itheta)+fi

              !mass source from energy:
              if(imass_source_hs(itsource).eq.1)then
                frho=coverfrac*energy/(hf/rnfuel)
                force(i,j,k,irho)=force(i,j,k,irho)+frho
              elseif(imass_source_hs(itsource).eq.2)then
                if(MLRHRRinputs_hs(itsource).eq.1) then
                  force(i,j,k,irho)=force(i,j,k,irho)+ &
                    coverfrac*fm_hs(itsource)/source_volume
                else
                  force(i,j,k,irho)=force(i,j,k,irho)+ &
                    coverfrac*fm_hs(itsource)
                endif
              endif !(imass_source_hs option)
            endif
          enddo
        enddo
      enddo
    endif
  enddo

end subroutine heatSource

!-----------------------------------------------------------------------
! Area of overlap between two axis-aligned rectangles.
!-----------------------------------------------------------------------
real(wp) function hs_rect_area(x1,x2,y1,y2,x3,x4,y3,y4)
  use gridlist_variables, only : wp
  Implicit None

  real(wp),intent(in) :: x1,x2,y1,y2,x3,x4,y3,y4

  hs_rect_area = max(0.0_wp,min(x2,x4)-max(x1,x3))* &
    max(0.0_wp,min(y2,y4)-max(y1,y3))

end function hs_rect_area

!-----------------------------------------------------------------------
! Exact area of overlap between a circle and an axis-aligned rectangle.
!-----------------------------------------------------------------------
real(wp) function hs_circle_rect_area(x1,x2,y1,y2,xc,yc,r)
  use gridlist_variables, only : wp
  Implicit None

  real(wp),intent(in) :: x1,x2,y1,y2,xc,yc,r
  integer,parameter :: maxpts=10
  integer :: n,ii,jj
  real(wp) :: xs(maxpts),tmp,a,b,c,d,xl,xr,ym,ss,ytop,ybot
  real(wp) :: total,root

  hs_circle_rect_area = 0.0_wp
  if(r.le.0.0_wp) return

  a = max(min(x1,x2)-xc,-r)
  b = min(max(x1,x2)-xc, r)
  c = min(y1,y2)-yc
  d = max(y1,y2)-yc
  if(a.ge.b.or.c.ge.r.or.d.le.-r) return

  n = 2
  xs(1) = a
  xs(2) = b
  call add_y_break(c)
  call add_y_break(d)

  do ii=1,n-1
    do jj=ii+1,n
      if(xs(jj).lt.xs(ii)) then
        tmp = xs(ii)
        xs(ii) = xs(jj)
        xs(jj) = tmp
      endif
    enddo
  enddo

  total = 0.0_wp
  do ii=1,n-1
    xl = xs(ii)
    xr = xs(ii+1)
    if(xr.le.xl) cycle
    ym = 0.5_wp*(xl+xr)
    ss = sqrt(max(0.0_wp,r*r-ym*ym))
    ytop = min(d,ss)
    ybot = max(c,-ss)
    if(ytop.le.ybot) cycle
    total = total + int_y(ytop,xl,xr) - int_y(ybot,xl,xr)
  enddo

  hs_circle_rect_area = max(0.0_wp,total)

contains

  subroutine add_y_break(yy)
    real(wp),intent(in) :: yy
    if(abs(yy).lt.r) then
      root = sqrt(max(0.0_wp,r*r-yy*yy))
      call add_x(max(a,min(b,-root)))
      call add_x(max(a,min(b, root)))
    endif
  end subroutine add_y_break

  subroutine add_x(xx)
    real(wp),intent(in) :: xx
    integer :: kk
    do kk=1,n
      if(abs(xs(kk)-xx).le.100.0_wp*epsilon(1.0_wp)*max(1.0_wp,r)) &
        return
    enddo
    if(n.lt.maxpts) then
      n = n + 1
      xs(n) = xx
    endif
  end subroutine add_x

  real(wp) function int_y(yy,xx1,xx2)
    real(wp),intent(in) :: yy,xx1,xx2
    real(wp) :: smid
    smid = sqrt(max(0.0_wp,r*r-(0.5_wp*(xx1+xx2))**2))
    if(abs(yy-smid).le.100.0_wp*epsilon(1.0_wp)*max(1.0_wp,r)) then
      int_y = circle_int(xx2)-circle_int(xx1)
    elseif(abs(yy+smid).le.100.0_wp*epsilon(1.0_wp)*max(1.0_wp,r)) then
      int_y = -(circle_int(xx2)-circle_int(xx1))
    else
      int_y = yy*(xx2-xx1)
    endif
  end function int_y

  real(wp) function circle_int(xx)
    real(wp),intent(in) :: xx
    real(wp) :: xclamp
    xclamp = max(-r,min(r,xx))
    circle_int = 0.5_wp*(xclamp*sqrt(max(0.0_wp,r*r-xclamp*xclamp)) &
      + r*r*asin(xclamp/r))
  end function circle_int

end function hs_circle_rect_area
