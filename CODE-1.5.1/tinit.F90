!----------------------------------------------------------------
! tinit initializes all arrays and variables in relation to the 
! to and sides BC absorbers
!----------------------------------------------------------------
subroutine tinit
  use gridlist_variables, only : n,m,l,nr,tow,zab, &
    iab,ibctopopen,wp,ibc_top,ibc_bot,ibc_left,ibc_right
  use gridsetup, only : np,mp
  use metric_variables, only : zb,z
  use constants, only : pi,tolerance
  use msga_variables, only : npos,mpos,mpi_rank
  use xvall, only : relaxxv
  use gridlist_variables, only : wp
  Implicit None
  
  ! Local Variables
  integer :: i,j,k,ia,ja,it
  integer, dimension(4) :: bcs
  real(wp) :: zl,t1
  real(wp),external :: zcart
  real(wp),allocatable :: relb(:,:)

  ! Executable Code
  if(tow.le.0._wp)then
    print*,'ERROR: tow must be positive'
    STOP
  endif
  allocate(relb(4,nr+1)); relb = 0._wp
  ! Define relaxing boundary (relb)
  bcs(1) = ibc_top
  bcs(2) = ibc_bot
  bcs(3) = ibc_left
  bcs(4) = ibc_right

  do j=1,4
    if(bcs(j).ge.1.and.bcs(j).lt.5) then
      select case(bcs(j))
        case(4)
          do i=1,nr+1
            relb(j,i)=exp(-real(i,kind=wp)*(real(i,kind=wp)-1._wp/ &
              real(i,kind=wp))**2/(2._wp*real(nr-1,kind=wp)**2))
          enddo
        case(3)
          do i=1,nr+1
            relb(j,i)=exp(-8._wp*real(i-1,kind=wp)/ &
              real(nr-1,kind=wp))
          enddo
        case(2)
          do i=1,nr+1
            relb(j,i)=1.5_wp*0.25_wp/(tow*(1._wp+(real(i,kind=wp)/ &
              (0.5_wp*real(nr,kind=wp)-1._wp))**2))
          enddo
        case(1)
          do i=1,nr+1
            relb(j,i)=0.05_wp*exp(-real(i,kind=wp)*(real(i,kind=wp)- &
              1._wp/real(i,kind=wp))**2/(2._wp*real(nr-1,kind=wp)**2))
          enddo
      end select
    endif
  enddo
  relb(:,1:2)=1._wp ! Mass conservation and Neumann boundary conditions
 
  if(zab.ge.zb)then
    if(mpi_rank.eq.0)then
      write(*,*) 'WARNING: top absorption layer is inactive.'
      write(*,*) '         zab is at or above the domain top.'
      write(*,*) '         zab = ', zab, ' zb = ', zb
      write(*,*) '         No vertical/top damping will be applied.'
      write(*,*) '         To enable it, set zab < zb in the gridlist.'
    endif
  endif

  do k=1,l
    do j=1,mp
      do i=1,np
        zl=zcart(z(k),i,j)
        ! Upper dampers
        if(zl.ge.zab)then
          if(ibctopopen.eq.0)then
            t1=max(0._wp,zl-zab)
            select case(iab)
              case(1)
                relaxxv(i,j,k,:)=t1/((zb-zab)+tolerance)/tow**2._wp
              case(2)
                relaxxv(i,j,k,:)=0.05_wp* &
                  exp(-t1/(-(zb-zab)/log(tow**(-2._wp))))
              case(3)
                relaxxv(i,j,k,1:3)= &
                  0.66_wp*0.25_wp*(1._wp-cos(pi*(zl-zab)/(zb-zab)))/tow
              case default
                print*,'ERROR: unsupported top absorber option iab=',iab
                STOP
            end select
          else
            relaxxv(i,j,k,3)=0.66_wp*0.25_wp* &
              (1._wp-cos(pi*(zl-zab)/(zb-zab)))/tow
          endif
        endif
        
        ! Side dampers
        if(bcs(1).ge.1.and.bcs(1).lt.5)then !top boundary
          ja=(mpos-1)*mp+j
          if(ja.gt.m-nr)then
            do it=1,4
              relaxxv(i,j,k,it)=max(relaxxv(i,j,k,it),relb(1,m-ja+1))
            enddo
          endif
        endif
        
        if(bcs(2).ge.1.and.bcs(2).lt.5)then !bottom boundary
          ja=(mpos-1)*mp+j
          if(ja.le.nr)then
            do it=1,4
              relaxxv(i,j,k,it)=max(relaxxv(i,j,k,it),relb(2,ja))
            enddo
          endif
        endif

        if(bcs(3).ge.1.and.bcs(3).lt.5) then !left boundary
          ia=(npos-1)*np+i
          if(ia.le.nr)then
            do it=1,4
              relaxxv(i,j,k,it)=max(relaxxv(i,j,k,it),relb(3,ia))
            enddo
          endif
        endif

        if(bcs(4).ge.1.and.bcs(4).lt.5) then !right boundary
          ia=(npos-1)*np+i
          if(ia.gt.n-nr)then
            do it=1,4
              relaxxv(i,j,k,it)=max(relaxxv(i,j,k,it),relb(4,n-ia+1))
            enddo
          endif
        endif

      enddo
    enddo
  enddo
  deallocate(relb)

end subroutine tinit
