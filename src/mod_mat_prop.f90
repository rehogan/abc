
module mod_mat_prop
! things related to temperature dependent material properties
! property tables are arranged as follows:
! the tables can be arranged in any order and each table can have a
! different number of entries (minimum of 2 for linear interpolation)
! tt_cp(1)              t_cp(1)     ! this is cp table
! tt_cp(2)              t_cp(2)
! ...                  ...
! tt_cp(num_cp_ent)     t_cp(num_cp_ent)

! tt_cond(1)            t_cond(1)   ! this is cond table
! tt_cond(2)            t_cond(2)
! ...                  ...
! tt_cond(num_cond_ent) t_cond(num_cond_ent)

! tt_emit(1)            t_emit(1)   ! this is emit table
! tt_emit(2)            t_emit(2)
! ...                  ...
! tt_emit(num_emit_ent) t_emit(num_emit_ent)

! dimensioned as (# materials,# table entries)

  use mod_precision
  use mod_constants
  implicit none
  integer :: num_mat
  integer, parameter                 :: max_no_mat = 10, max_no_ent = 6
  integer, allocatable, dimension(:) :: no_cp_ent,   no_cond_ent, &
                                        no_emit_ent, no_abs_ent
  real(dp),          allocatable, dimension(:) :: rho      ! density (kg/m^3)
  character(len=40), allocatable, dimension(:) :: mat_name ! array mat names
  integer,           allocatable, dimension(:) :: mat_num  ! material #

  ! specific heat, conductivity, emittance, absorbtance, internal energy
  ! dimensionsed as xx(i,j): i = entry #, j = material #
  real(dp), allocatable, dimension(:,:) :: TT_cp,   T_cp
  real(dp), allocatable, dimension(:,:) :: TT_cond, T_cond
  real(dp), allocatable, dimension(:,:) :: TT_emit, T_emit
  real(dp), allocatable, dimension(:,:) :: TT_abs,  T_abs
  real(dp), allocatable, dimension(:,:) :: T_int_en
end module mod_mat_prop
