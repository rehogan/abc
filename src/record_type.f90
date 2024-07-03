

subroutine record_type(icol,icol_max,record,rec_typ)

! this subroutine takes a record of length icol and decides record type
! the following characters define begin, comment, and end
! 'begin' == b, begin block
! '!'     == c, comment
! '#'     == c, comment
! '*'     == c, comment
! 20b     == c, comment  (1st 20 columns are blank)
! 'end'   == e, end block
! '#end'  == e, end block
! '# end' == e, end block
! if none of the above exist, then data record is data (d)

! this procedure works best if leading blanks have been removed using
! strip_blank.f90

  implicit none
  integer, intent(in) ::  icol, icol_max
  character (len=icol_max) :: record
  character (len=1) :: rec_typ

  if (record(1:3) == 'end' .or. record(1:4) == '#end' .or. &
       record(1:5) == '# end') then
     rec_typ = 'e'
  else if (record(1:1) == '!' .or. record(1:1) == '#' &
       .or. record(1:1) == '*') then
     rec_typ = 'c'
  else if (record(2:2) == '!' .or. record(2:2) == '#' ) then
     rec_typ = 'c'
  else if (record(1:20) == ' ') then
     rec_typ = 'c'
  else if (record(1:5) == 'begin' .or. record(1:6) == '#begin' &
       .or. record(1:7) == '# begin') then
     rec_typ = 'b'
  else
     rec_typ = 'd'
  endif

  return

end subroutine record_type
