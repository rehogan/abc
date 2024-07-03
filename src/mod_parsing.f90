

module mod_parsing
! things related to parsing of input records
  integer :: ncol
  integer, parameter :: ncol_max = 80
  character(len=ncol_max) :: record
  character(len=1) :: rec_typ
end module mod_parsing
