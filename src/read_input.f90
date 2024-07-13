
subroutine read_input
  ! read all the input records based on detection of "begin xxxx" on
  ! previously read record
  ! possible key words include "TITLE", "GLOBAL PARAMETERS,
  ! "TIME INTEGRATION", "MATERIAL PROPERTIES"
  use mod_bc
  use mod_files
  use mod_parsing
  use mod_to_lower
  implicit none
! local variables
  integer :: j, mat

  interface
     subroutine read_title
       use mod_parsing
       use mod_files
       implicit none
       integer :: i, num_title
       character (len=ncol_max), dimension(10):: title
     end subroutine read_title

     function compare_string(string,substring)
       implicit none
       logical :: compare_string
       character (len=*), intent (in) :: string, substring
     end function compare_string

     subroutine read_conv_bc()
       use mod_constants
       use mod_bc
       use mod_files
       use mod_global_parameters
       use mod_parsing
       use mod_to_lower
       implicit none
     end subroutine read_conv_bc

     subroutine read_glob_param
       use mod_constants
       use mod_files
       use mod_global_parameters
       use mod_parsing
       use mod_to_lower
       implicit none
     end subroutine read_glob_param

     subroutine read_therm_prop
       use mod_constants
       use mod_nodes_elements
       use mod_files
       use mod_mat_prop
       use mod_parsing
       implicit none
     end subroutine read_therm_prop

     subroutine read_spec_t_bc()
       use mod_constants
       use mod_bc
       use mod_files
       use mod_global_parameters
       use mod_parsing
       use mod_to_lower
       implicit none
     end subroutine read_spec_t_bc

     subroutine read_encl_rad_bc
       use mod_constants
       use mod_files
       use mod_global_parameters
       use mod_parsing
       use mod_to_lower
       implicit none
       interface
          function compare_string(string,substring)
            implicit none
            logical :: compare_string
            character (len=*), intent (in) :: string, substring
          end function compare_string
       end interface
     end subroutine read_encl_rad_bc

  end interface

  ! initialize some counters in case this type bc is not read
  no_conv_bc = 0
  no_spec_t_bc = 0
10 continue
  read( kinp,"(a)",end=99 ) record
  ! when ! or # appears in columns 1 or 2, treat as comment
  if (record(1:1) == '!' .or. record(1:1) == '#' .or. &
       record(2:2) == '!' .or. record(2:2) == '#') then
     go to 10
  end if
  write(kout,*) record
  record = to_lower(record)               ! convert to lower case
  if(compare_string(record,'title')) then ! read title records, 10 max
     write(6,"('begin title block')")
     call read_title
     write(6,"('end title block'/)")
  else if(compare_string(record,'begin global par')) then ! read global parameters
     write(6,"('begin global parameters block')")
     call read_glob_param
     write(6,"('end global parameters block'/)")
  else if(compare_string(record,'begin thermal prop')) then ! read thermal props
     write(6,"('begin thermal property block')")
     call read_therm_prop
     write(6,"('end thermal property block'/)")
  else if(compare_string(record,'begin convection bc')) then ! read thermal props
     write(6,"('begin convection bc block')")
     call read_conv_bc
     write(6,"('end convection bc block'/)")
  else if(compare_string(record,'begin specified t bc')) then ! read spec T bc
     write(6,"('begin specified t bc block')")
     call read_spec_t_bc()
     write(6,"('end specified t bc block'/)")
  else if(compare_string(record,'begin encl rad bc')) then ! read encl rad bc
     write(6,"('begin encl rad bc block')")
     call read_encl_rad_bc
     write(6,"('end encl rad bc block'/)")  
  end if

  go to 10

99 continue

end subroutine read_input
