

module mod_to_lower
  implicit none
contains

  function to_lower(strIn) result(strOut)
    ! Adapted from original code from Clive Page
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j

    interface
       subroutine strip_blank(record)
         implicit none
         character(len=80) :: record
       end subroutine strip_blank
    end interface

    call strip_blank(strIn)        ! remove leading blanks first
    do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j>= iachar("A") .and. j<=iachar("Z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))+32)
       else
          strOut(i:i) = strIn(i:i)
       end if
    end do

  end function to_lower

end module mod_to_lower
