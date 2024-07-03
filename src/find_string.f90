
subroutine find_string(record,substring)

  ! reads successive input records till you find a record containing substring
  ! use cautiously in order to not inadvertently reach an 'end of file' on read

  ! input:
  !   record    = input record
  !   substring = substring your are lookin for

  ! output:
  !  record = output record that contains substring

  use mod_files
  use mod_to_lower
  implicit none

  interface
     function compare_string(string,substring)
       implicit none
       logical :: compare_string
       character (len=*), intent (in) :: string, substring
     end function compare_string
  end interface

  integer :: i
  integer, parameter :: imax = 100
  character (len=*), intent (in) :: substring
  character (len=*), intent (out) :: record

  do i=1,imax  ! till you find substring, e.g. 'begin material'
     read(kinp,"(a)") record  ! read a new record
     write(kout,*) record
     record = to_lower(record)      ! convert to upper case and remove leading blanks
     if(compare_string(record,substring)) then
        exit
     end if
  end do

  if(i == imax) then
     write(6,"('****error in find_string, i > imax')")
     stop
  end if

  return

end subroutine find_string
