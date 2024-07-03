subroutine read_spec_t_bc()

  ! read the convection boundary condition block
  ! convection coefficient is evaluated in subroutine natl_conv_corr
  ! upon entering this routine, the record 
  ! "begin convection bc block" has already been found

  use mod_constants
  use mod_bc
  use mod_files
  use mod_global_parameters
  use mod_parsing
  use mod_to_lower
  implicit none
  interface
     subroutine read_record(kifile,kofile,ncol,ncol_max,record,rec_typ)
       integer,                 intent(in)  :: kifile, kofile, ncol, ncol_max
       character(len=1),        intent(out) :: rec_typ
       character(len=ncol_max), intent(out) :: record
     end subroutine read_record

     function compare_string(string,substring)
       implicit none
       logical :: compare_string
       character (len=*), intent (in) :: string, substring
     end function compare_string
  end interface
  integer :: j

  no_spec_t_bc = 0
  j = 0
  do     ! till you find the record 'end specified t bc block' 
     read(kinp,"(a)") record ! read input record
!     write(kout,*) record
     record = to_lower(record)         ! convert to lower case
     call read_record(kinp,kout,ncol,ncol_max,record,rec_typ)
     if(rec_typ == 'e') then
        exit
     else if(rec_typ == 'd') then
        j = j + 1
        read(record,*) ni_bc_spec_t(j),t_st_bc(j)
     end if
  end do
  no_spec_t_bc = j
end subroutine read_spec_t_bc
