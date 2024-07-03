
subroutine read_glob_param

! read the global parameters block

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

  character(len=12) :: char1, char2, char3, char4   ! dummy character variables
  character(len=2) :: delim        ! delimiter

  do     ! till you find the record 'end global constants block'
     read(kinp,"(a)") record ! read input record
     write(kout,*) record
     record = to_lower(record)         ! convert to lower case
     if(compare_string(record,'end global')) then
        exit
     else if(compare_string(record,'verification')) then    ! read verif
        read(record,*) char1,delim,verif
     else if(compare_string(record,'debug')) then           ! read debug
        read(record,*) char1,delim,debug
     else if(compare_string(record,'analytical')) then      ! read anal
        read(record,*) char1,char2,char3,char4,delim,anal_ini_prof
     else if(compare_string(record,'sigma')) then           ! read sigma
        read(record,*) char1,delim,sigma
     else if(compare_string(record,'length conv')) then     ! read len_conv
        read(record,*) char1,char2,delim,len_conv
     else if(compare_string(record,'maximum non')) then     ! read max_nl_iter
        read(record,*) char1,char2,char3,delim,max_nl_iter
     else if(compare_string(record,'maximum glob')) then    ! read max_glob_iter
        read(record,*) char1,char2,char3,delim,max_glob_iter
     else if(compare_string(record,'nonlinear iter')) then  ! read nl_iter_tol
        read(record,*) char1,char2,char3,delim,nl_iter_tol
     else if(compare_string(record,'initial temp')) then    ! read unif_ini_temp
        read(record,*) char1,char2,delim,unif_ini_temp
     else if(compare_string(record,'room width')) then      ! read room width
        read(record,*) char1,char2,delim,rm_w
     else if(compare_string(record,'room height')) then     ! read room height
        read(record,*) char1,char2,delim,rm_h
     else if(compare_string(record,'room depth')) then      ! read room depth
        read(record,*) char1,char2,delim,rm_d
     else if(compare_string(record,'z coord')) then         ! read z_coord_shft
        read(record,*) char1,char2,char3,delim,z_coord_shft
     else if (compare_string(record,'problem units')) then
        read(record,*) char1,char2,delim,units
        if(compare_string(units,'eng')) units ='engineering'
        if (compare_string(units,'si') .or.                             &
             compare_string(units,'eng')) then
           continue
        else
           write(*,*) "***execution terminated: problem units"
           write(*,*) "must be 'si' or 'engineering'"
           stop
        end if
     end if

  end do

  if(verif == 'yes') then
     open(unit=kver, file=ver_file, status="replace", form='formatted', &
          iostat=ierror) ! verification file
  end if
  return

end subroutine read_glob_param
