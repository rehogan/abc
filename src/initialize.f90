
subroutine initialize
  ! initialize various arrays and some preliminary calculations
  ! that do not change once the geometry is read in
  use mod_precision
  use mod_constants
  use mod_global_parameters
  use mod_mat_prop
  use mod_nodes_elements
  use mod_time_params
  implicit none
  real(dp),  dimension(4) :: x_g, y_g ! global element coord
  integer                 :: i, ne

  interface

     subroutine get_en_cnt
       use mod_precision
       use mod_constants
       use mod_mat_prop
       use mod_nodes_elements
       use mod_time_params
       implicit none
     end subroutine get_en_cnt

  end interface

  ! for any specified temperature nodes, initialize the last iteration
  ! Temp t_li, presently hard wired
  ! temperature array, 2nd index is for
  ! times nm1, n, np1 (n-1, n, n+1)
  t(1:no_nodes,1:3) = unif_ini_temp   
  t_li(1:no_nodes) = unif_ini_temp    ! temp at last iteration (li)
  dt(1:no_nodes) = zero        ! temp correction for iteration
  t_li(1:12) = 300._dp         ! hard wired
  t_li(13:24) = 1300._dp
  t_li(11:12) = 400_dp
  t_li(13:14) = 800_dp
  t_li(11:12) = 373.371_dp
  t_li(13:14) = 880.736_dp
  t_li(3:22) = 1000._dp

  first_index = 1          ! first index of time(i): (n-1), (n), (n+1)

  ! initialize the energy content arrays
  call get_en_cnt

end subroutine initialize
