

module mod_interface_definitions
  interface
     function compare_string(string,substring)
       implicit none
       logical :: compare_string
       character (len=*), intent (in) :: string, substring
     end function compare_string

     subroutine record_type(icol,icol_max,record,rec_typ)
       implicit none
       integer, intent(in) ::  icol, icol_max
       character (len=icol_max) :: record
       character (len=1) :: rec_typ
     end subroutine record_type

  end interface

end module mod_interface_definitions
