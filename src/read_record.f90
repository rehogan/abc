

subroutine read_record(kifile,kofile,ncol,ncol_max,record,rec_typ)

  ! read a data record. if it is a comment record, then skip it
  ! and read the next data or terminator record
  ! when data record is read, input unit is backspaced

  ! input:
  ! kifile   = input file number
  ! kofile   = output file number
  ! ncol     = number of columns of record to be examined
  !            ex. if ncol = 60, then only the first 60 columns are
  !            examined and there can be comments past column 60
  ! ncol_max = maxinum number of columns allowed

  ! output:
  ! record   = character string containing data record and is returned
  !           from this routine
  ! rec_typ  = record type as determined in routine rectype,
  !                       rec_typ = c, comment record
  !                               = d, data record
  !                               = e, terminator record such as "end"
  !                               = b, begin record such as "begin"

  ! written by:
  ! ben blackwell
  ! corrales, nm 87048
  ! 505-400-0205
  ! bblackwell13@comcast.net

  ! declarations

  use mod_to_lower
  implicit none
  interface
     subroutine record_type(icol,icol_max,record,rec_typ)
       implicit none
       integer, intent(in) ::  icol, icol_max
       character (len=icol_max) :: record
       character (len=1) :: rec_typ
     end subroutine record_type

  end interface
  integer,                 intent(in)  :: kifile, kofile, ncol, ncol_max
  character(len=1),        intent(out) :: rec_typ
  character(len=ncol_max), intent(out) :: record

  ! local variables
  integer :: ierror

  do      ! till we find data (d), begin (b) or end (e) record

     read(kifile,"(a)",iostat=ierror) record ! read input record
     if (ierror /= 0) then  ! we have an end of file
        rec_typ = 'e'
        exit
     end if
     write(kofile,*) record
     record = to_lower(record)   ! convert to lower case

     ! examine the first ncol columns of this record to determine if
     ! it is comment, data, or terminator

     call record_type(ncol,ncol_max,record,rec_typ)

     ! if record is comment, then ignore it and read next record

     if (rec_typ == 'c') then
        cycle

        ! if record is data, backspace input file and return to calling routine

     else if (rec_typ == 'd') then
        backspace (kifile)
        exit

        ! if record is end, return without backspacing input file

     else if (rec_typ == 'e') then
        exit

        ! if record is begin, return without backspacing input file

     else if (rec_typ == 'b') then
        exit
     endif

  end do

  return

end subroutine read_record
