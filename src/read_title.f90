
subroutine read_title

  use mod_parsing
  use mod_files
  use mod_to_lower
  implicit none
  integer :: i, num_title
  character (len=ncol_max), dimension(10):: title

  ncol = 80
  do i=1,11                          ! maximum of 10 records, 80 col/record
     read( kinp,"(a)",end=99 ) record
     record = to_lower(record)       ! convert to lower case
     if (record(1:3) == 'end' .or. record(2:4) == 'end') then
        exit
     else
        backspace(kinp)
        read(kinp,"(80a)") title(i)
        write(kout,"(80a)") title(i)
     end if
  end do

  num_title = i - 1
  write(kout,"(i3,' title cards read')") num_title

   return

99 continue
   write(6,"('***execution terminated, no title information')")
   write(kout,"('***execution terminated, no title information')")
   stop

end subroutine read_title
