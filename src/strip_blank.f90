
subroutine strip_blank(record)

! strip leading blanks from data record
! ascii character 32 is blank or space

! input:
! record = character string of maximum length 80

! output:
! record = character string with leading blanks removed

  implicit none
  character(len=80) :: record
! local variables
  integer :: i, ia, num_char
  character(len=1) :: a

  num_char = len(record)

  do i=1,num_char
     a = record(i:i)
     ia = ichar(a)
     if (ia .ne. 32) exit
  end do

! there are i-1 blank characters

  record = record(i:num_char)

  return

end subroutine strip_blank
