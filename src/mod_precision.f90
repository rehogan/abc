
module mod_precision

  ! select single (sp) and double precision (dp)

  implicit none
  save
  integer, parameter :: sp = kind(0.0  )    ! single precision
  integer, parameter :: dp = kind(0.0d0)    ! double precision

end module mod_precision
