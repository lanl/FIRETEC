!-----------------------------------------------------------------------
! updates position of a circular heat source
!-----------------------------------------------------------------------
subroutine updateHeatSourcePosition(iheatsource,rhof,nfuel,il,iu,jl,ju &
    ,lls)
  use gridlist_variables, only : prec,dx,dy
  use heatsource_gridlist_variables, only : nsource_hs,radius_hs, &
    xcen_hs,ycen_hs,xl_hs,xu_hs,yl_hs,yu_hs
  use gridsetup, only : np,mp
  use msga_variables, only : npos,mpos
  Implicit none
 
  ! Local Variables
  integer,intent(in) :: iheatsource
  integer,intent(in) :: nfuel,il,iu,jl,ju,lls
  real(prec),intent(inout) :: rhof(nfuel,il:iu,jl:ju,lls)

  integer :: i,j,ia,ja,itsource

  ! Executable Code
  ! here we update value of rhof in the heat source 
  ! FIXME :: No temporal translation of position at this time AJJ (15 July 24)
  if(iheatsource.eq.1) then ! Circular
    do j=jl,ju
      ja=(mpos-1)*mp+j
      do i=il,iu
        ia=(npos-1)*np+i
        do itsource=1,nsource_hs
          if((sqrt(((ia*dx) - xcen_hs(itsource))**2 + &
            ((ja*dy) - ycen_hs(itsource))**2).lt.radius_hs(itsource))) &
            rhof(:,i,j,1) = 0.01 ! burn fuel (drag is very small)
        enddo
      enddo
    enddo 
  elseif(iheatsource.eq.2) then ! Rectangular
    do j=jl,ju
      ja=(mpos-1)*mp+j
      do i=il,iu
        ia=(npos-1)*np+i
        do itsource=1,nsource_hs
          if((ia*dx).ge.xl_hs(itsource).and.(ia*dx).le.xu_hs(itsource) &
            .and.(ja*dy).ge.yl_hs(itsource) &
            .and.(ja*dy).le.yu_hs(itsource))&
            rhof(:,i,j,1)=0.01 ! burn fuel (drag is very small)
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
  use gridlist_variables, only : dx,dy,prec
  use gridsetup, only : time
  use constants, only : pi,pref
  use thermo_variables, only : pr,cp_gas,cv_gas
  use gridsetup, only : np,mp
  use metric_variables, only : zedge
  use msga_variables, only : npos,mpos,mpi_rank
  use xvall, only : irho,itemp
  use fuel_variables, only : rnfuel,hf
  use zcart_function, only : zcart
  use constants, only : tolerance
  Implicit None

  ! Local Variables
  integer,intent(in) :: iheatsource
  integer,intent(in) :: il,iu,jl,ju,lls,nv
  real(prec),intent(inout) :: force(il:iu,jl:ju,lls,nv)

  integer :: i,j,k,itsource,iostat
  integer :: ia,ja
  real(prec) :: source_function
  real(prec) :: energy,fi,frho
  real(prec) :: volume

  

  ! Executable Code
  do itsource=1,nsource_hs
    select case(isource_function)
      case(0)
        source_function=1.0
      case(1)
        if(sin(2.*pi*freq_hs(itsource)*time).ge.0.0)then
          source_function = 1.0
        else
          source_function = 0.0
        endif
      case(2)
        source_function = 0.5*(1.+sin(2.*pi*freq_hs(itsource)*time))
      case(3)
        if((freq_hs(itsource).gt.0).and. &
          (freq_hs(itsource)*time.le.1.0))then
          source_function = 1.0
        else
          source_function = 0.0
        endif
    end select

    if(MLRHRRinputs.eq.1) then
      if(TME(itsource,2).lt.time) then
        TME(itsource,1) = TME(itsource,2)
        MLR(itsource,1) = MLR(itsource,2)
        HRR(itsource,1) = HRR(itsource,2)
        read(7825+itsource,*,iostat=iostat) &
          TME(itsource,2),MLR(itsource,2),HRR(itsource,2)
        if(mpi_rank.eq.0) print*,'Updating MLR/HRR data'
        if(iostat.lt.0) then ! iostat<0 means end of file
          TME(itsource,2) = TME(itsource,1)
          MLR(itsource,2) = MLR(itsource,1)
          HRR(itsource,2) = HRR(itsource,1)
        endif
      endif

      if(abs((TME(itsource,1)-TME(itsource,2))/TME(itsource,1)).lt. &
        tolerance) then
        fm_hs(itsource) = MLR(itsource,2)
      else
        fm_hs(itsource) = MLR(itsource,1) + &
          (time-TME(itsource,1))/(TME(itsource,2)-TME(itsource,1))* &
          (MLR(itsource,2)-MLR(itsource,1))
      endif
    endif

    if(iheatsource.eq.1)then !Cicular source
      do k=1,lls
        do j=jl,ju
          do i=il,iu
            ia=(npos-1)*np+i
            ja=(mpos-1)*mp+j
            if((sqrt((ia*dx+dx/2. - xcen_hs(itsource))**2+ &
              (ja*dy+dy/2. - ycen_hs(itsource))**2).lt. &
              radius_hs(itsource)).and. &
              (zcart(zedge(k),i,j).lt.depth_hs(itsource))) then
                volume=dx*dy*(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))
                if(MLRHRRinputs.eq.1) then
                  if(abs((TME(itsource,1)-TME(itsource,2)) &
                    /TME(itsource,1)).lt.tolerance) then
                    energy = HRR(itsource,2)/volume
                  else
                    energy = (HRR(itsource,1) + &
                      (time-TME(itsource,1))/ &
                      (TME(itsource,2)-TME(itsource,1))* &
                      (HRR(itsource,2)-HRR(itsource,1)))/volume
                  endif
                else
                  energy= &
                    av_hs(itsource)*flux_hs(itsource)*source_function
                endif
                fi=energy/cp_gas(i,j,k)* &
                  (pref/pr(i,j,k))**(1.-cv_gas(i,j,k)/cp_gas(i,j,k))
                ! computation of energy really realised in gas phase
                force(i,j,k,itemp)=force(i,j,k,itemp)+fi

                !mass source from energy:
                if(imass_source.eq.1)then
                  frho=energy/(hf/rnfuel)
                  force(i,j,k,irho)=force(i,j,k,irho)+frho
                elseif(imass_source.eq.2)then
                  force(i,j,k,irho)=force(i,j,k,irho)+fm_hs(itsource)
                endif !(imass_source.ne.0)then
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

            if((ia-1)*dx.ge.xl_hs(itsource) &
              .and.(ja-1)*dy.ge.yl_hs(itsource) &
              .and.ia*(dx).le.xu_hs(itsource) &
              .and.ja*(dy).le.yu_hs(itsource) &
              .and.zcart(zedge(k),i,j).lt.depth_hs(itsource))then
              volume = dx*dy*(zcart(zedge(k+1),i,j)-zcart(zedge(k),i,j))

              if(MLRHRRinputs.eq.1) then
                if(abs((TME(itsource,1)-TME(itsource,2)) &
                 /TME(itsource,1)).lt.tolerance) then
                  energy = HRR(itsource,2)/volume
                else
                  energy = (HRR(itsource,1) + &
                    (time-TME(itsource,1))/ &
                    (TME(itsource,2)-TME(itsource,1)) * &
                    (HRR(itsource,2)-HRR(itsource,1)))/volume
                endif
              else
                energy=av_hs(itsource)*flux_hs(itsource)*source_function
              endif
              fi=energy/cp_gas(i,j,k)* &
                (pref/pr(i,j,k))**(1.-cv_gas(i,j,k)/cp_gas(i,j,k))
              ! computation of energy really realised in gas phase
              force(i,j,k,itemp)=force(i,j,k,itemp)+fi

              !mass source from energy:
              if(imass_source.eq.1)then
                frho=energy/(hf/rnfuel)
                force(i,j,k,irho)=force(i,j,k,irho)+frho
              elseif(imass_source.eq.2)then
                force(i,j,k,irho)=force(i,j,k,irho)+fm_hs(itsource)
              endif !(imass_source.ne.0)then
            endif
          enddo
        enddo
      enddo
    endif
  enddo

end subroutine heatSource

