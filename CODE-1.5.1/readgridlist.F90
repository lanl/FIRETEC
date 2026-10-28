module QueryGridlist_mod
  use gridlist_variables, only : wp
  use iso_fortran_env, only : int64
  Implicit None

  contains
    
    !-------------------------------------------------------------------
    ! Queries gridlist for input variable values regardless of type
    !-------------------------------------------------------------------
    subroutine QueryGridlist_variable(variableName,variable,fileunit)
      Implicit None
    
      ! Local Variables
      character(len=*),intent(in) :: variableName
      class(*),intent(inout) :: variable
      integer,intent(in)  :: fileunit
    
      integer :: ierror,eqpos,nfound
      character(len=20)   :: word
      character(len=1000) :: queryName
      character(len=1000) :: text,valueText
    
      ! Executable Code
      queryName = variableName
      call lowercase_inplace(queryName)
      nfound = 0
      do ! Iterate through all lines of file
        read(fileunit,"(a)",iostat=ierror) text ! Read line into text
        if (ierror/=0) exit
        if(text/='')then
          eqpos=index(text,'=')
          word =trim(adjustl(text(:eqpos-1)))
          call lowercase_inplace(word)
          if(trim(word).eq.trim(queryName))then
            nfound = nfound + 1
            if(nfound.gt.1)then
              print*,'ERROR: duplicate gridlist key: ', &
                trim(variableName)
              STOP
            endif
            valueText=trim(text(eqpos+1:))
            select type (x=>variable)
              type is (integer)
                read(valueText,*) x
              type is (integer(int64))
                read(valueText,*) x
              type is (real(wp))
                read(valueText,*) x
              type is (character(len=*))
                read(valueText,*) x
              class default
                print*,"Unexpected type in reading gridlist"
            end select
          endif
        endif
      enddo
      rewind(fileunit)
    end subroutine QueryGridlist_variable
    
    !-------------------------------------------------------------------
    ! Queries gridlist for input array values regardless of type
    !-------------------------------------------------------------------
    subroutine QueryGridlist_array(variableName,array,size,fileunit)
      Implicit None

      ! Local Variables
      character(len=*),intent(in) :: variableName
      integer,intent(in)  :: size
      integer,intent(in)  :: fileunit
      class(*),intent(inout) :: array(size)

      integer :: ierror,eqpos,nfound
      character(len=20)   :: word
      character(len=1000) :: queryName
      character(len=1000) :: text,valueText

      ! Executable Code
      queryName = variableName
      call lowercase_inplace(queryName)
      nfound = 0
      do ! Iterate through all lines of file
        read(fileunit,"(a)",iostat=ierror) text ! Read line into text
        if (ierror/=0) exit
        if(text/='')then
          eqpos=index(text,'=')
          word =trim(adjustl(text(:eqpos-1)))
          call lowercase_inplace(word)
          if(trim(word).eq.trim(queryName))then
            nfound = nfound + 1
            if(nfound.gt.1)then
              print*,'ERROR: duplicate gridlist key: ', &
                trim(variableName)
              STOP
            endif
            valueText=trim(text(eqpos+1:))
            select type (x=>array)
              type is (real(wp))
                read(valueText,*) x
              type is (character(len=*))
                read(valueText,*) x
              class default
                print*,"Unexpected type in reading gridlist"
            end select
          endif
        endif
      enddo
      rewind(fileunit)
    end subroutine QueryGridlist_array

    subroutine lowercase_inplace(text)
      Implicit None

      ! Local Variables
      character(len=*),intent(inout) :: text
      integer :: i,c

      ! Executable Code
      do i=1,len_trim(text)
        c = iachar(text(i:i))
        if(c.ge.iachar('A').and.c.le.iachar('Z'))then
          text(i:i) = achar(c + iachar('a') - iachar('A'))
        endif
      enddo
    end subroutine lowercase_inplace

end module QueryGridlist_mod
