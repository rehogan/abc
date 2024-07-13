

function compare_string(string,substring)

! search 'string' for the presence of 'substring'
! modeled after f77 routine developed by r. j. cochran

  implicit none
  logical :: compare_string
  character (len=*), intent (in) :: string, substring

  if(index(string,substring) /= 0) then
     compare_string = .true.
  else
     compare_string = .false.
  end if

end function compare_string
