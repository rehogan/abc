

module mod_nodes_elements

  ! max_no_nodes = maximum number of nodes
  ! max_no_elem = maximum number of elements
  ! max_no_el_blk = maximum number of element blocks
  ! no_el_in_blk(i) = no elements in element block i

  use mod_precision
  implicit none
  save

  integer, parameter :: max_elem_in_blk = 100, max_no_el_blk = 10,  &
                        max_no_elem   = max_elem_in_blk*max_no_el_blk


  integer            :: max_elem_blk, no_el_blk, nfb, nhb, nhbp1, num_nodes_in_vec
  integer            :: no_nodes, no_elem   ! read in Abaqus grid file
  integer            :: ne_eb(max_no_elem,max_no_el_blk),&
                                                    no_el_in_blk(max_no_el_blk)

!  integer,          dimension(4,max_no_elem)     :: nn ! global node to local node map
  integer,  allocatable,  dimension(:,:)        :: nn ! global node to local node map
  integer,          dimension(max_no_el_blk)    :: mat_no  ! elem blk
  integer,          dimension(max_no_elem)      :: seq_el_no
  real(dp), allocatable,  dimension(:)      :: x, y, xi, yi, dt, t_li , b, z
  real(dp), allocatable,  dimension(:,:)    :: c
  ! T(j,i), Tdot(j,i), j = node index, i = time index (n-1, n, n+1)
  real(dp), allocatable,  dimension(:,:)    :: T , Tdot
  ! energy content en_cnt(i,j,k), i=el no, j=time index, k=local node no
  real(dp), allocatable,  dimension(:,:,:)  :: en_cnt
  ! volume of sub-control volume, v_scv(i,j), i = local node no, j = elem no
  real(dp), allocatable,  dimension(:,:)    :: v_scv
  real(dp)                                  :: zi_dum, dz

  character(len=3)   :: itinit
  character(len=72)  :: grid_title

end module mod_nodes_elements
