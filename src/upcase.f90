
subroutine upcase(char_string)

! convert character string to all upper case
! modeled after routine written by r. j. cochran
! lower case alphabet characters lie in the ascii
! character set range 97-122

! input:
! ncol_max    = maximum number of columns
! char_string = character string of maximum length 80

! output:
! char_string = character string converted to upper case

  implicit none
  character(len=*), intent(inout) :: char_string
  integer :: i, ia, num_char
  character(len=1) :: a

  interface
     subroutine strip_blank(record)
       implicit none
       character(len=80) :: record
     end subroutine strip_blank
  end interface

! strip leading blanks from record

  call strip_blank(char_string)

  num_char = len(char_string)

  do i=1,num_char
     a = char_string(i:i)
     ia = ichar(a)
     if (ia >= 97 .and. ia <= 122) then   ! convert lower case to upper
        ia = ia - 32
        a = char(ia)
        char_string(i:i) = a
     end if
  end do

  return

end subroutine upcase
