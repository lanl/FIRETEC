!----------------------------------------------------------------
! Queries gridlist for input integer values
!----------------------------------------------------------------
subroutine QueryGridlist_integer(variableName,variable,fileunit)
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(out) :: variable
  integer,intent(in)  :: fileunit

  integer :: ierror
  character(len=1)    :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    if(text/='')then
      read(text,*) word
      if(word.eq.variableName)then
        read(text,*) word,equal,variable
        exit
      endif
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryGridlist_integer

!----------------------------------------------------------------
! Queries gridlist for input real values
!----------------------------------------------------------------
subroutine QueryGridlist_real(variableName,variable,fileunit)
  use gridlist_variables, only : prec
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: fileunit
  real(prec),intent(out) :: variable

  integer :: ierror
  character(len=1)    :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    if(text/='')then
      read(text,*) word
      if(word.eq.variableName)then
        read(text,*) word,equal,variable
        exit
      endif
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryGridlist_real

!----------------------------------------------------------------
! Queries gridlist for input real values
!----------------------------------------------------------------
subroutine QueryGridlist_string(variableName,variable,fileunit)
  Implicit None

  ! Local Variables
  character(len=*),intent(in)  :: variableName
  character(len=*),intent(out) :: variable
  integer,intent(in)  :: fileunit

  integer :: ierror
  character(len=1)    :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    if(text/='')then
      read(text,*) word
      if(word.eq.variableName)then
        read(text,*) word,equal,variable
        exit
      endif
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryGridlist_string

!----------------------------------------------------------------
! Queries gridlist for input string arrays
!----------------------------------------------------------------
subroutine QueryGridlist_string_array(variableName,array, &
    arraysize,fileunit)
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: arraysize
  integer,intent(in)  :: fileunit
  character(len=*),intent(out) :: array(arraysize)

  integer :: ierror
  character(len=20)   :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    if(text/='')then
      read(text,*) word
      if(word.eq.variableName)then
        read(text,*) word,equal,array
        exit
      endif
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryGridlist_string_array

!----------------------------------------------------------------
! Queries gridlist for input real arrays
!----------------------------------------------------------------
subroutine QueryGridlist_real_array(variableName,array, &
    arraysize,fileunit)
    use gridlist_variables, only : prec
  Implicit None

  ! Local Variables
  character(len=*),intent(in) :: variableName
  integer,intent(in)  :: arraysize
  integer,intent(in)  :: fileunit
  real(prec),intent(out) :: array(arraysize)

  integer :: ierror
  character(len=1)    :: equal
  character(len=20)   :: word
  character(len=1000) :: text

  ! Executable Code
  do ! Iterate through all lines of file
    read(fileunit,"(a)",iostat=ierror) text ! Read line into text
    if (ierror/=0) exit
    if(text/='')then
      read(text,*) word
      if(word.eq.variableName)then
        read(text,*) word,equal,array
        exit
      endif
    endif
  enddo
  rewind(fileunit)
    
end subroutine QueryGridlist_real_array
