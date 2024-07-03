subroutine get_en_cnt
  ! get energy content of all sub-control volumes at the initial
  ! element temperatures
  use mod_precision
  use mod_constants
  use mod_mat_prop
  use mod_nodes_elements
  use mod_time_params
  implicit none

  interface
     subroutine get_prop(tbar, mat_no_eb, cond, dk_dt, cp, dcp_dt, int_en)
       use mod_precision
       use mod_constants
       use mod_mat_prop
       implicit none
       real(dp), intent(inout)  :: tbar
       integer,       intent(in)     :: mat_no_eb
       real(dp), intent(out)    :: cond, dk_dt, cp, dcp_dt, int_en
     end subroutine get_prop

     subroutine vol_of_scv_quad(ne,x_g,y_g)
       use mod_precision
       use mod_constants
       use mod_nodes_elements
       implicit none
       save
       integer,  intent(in)                  :: ne       ! elem no
       real(dp), intent(in),  dimension(4)   :: x_g, y_g ! global element coord
     end subroutine vol_of_scv_quad
  end interface

  integer                      :: i, j, jj, k, lnn, mat_no_eb, ne, nscv
  real(dp)                     :: rho_int_en, coef_1, coef_2, dum, den
  real(dp),       dimension(4) :: aa, vol
  real(kind=dp),  dimension(4) :: t_gn, x_g, y_g ! corner (x,y) ne
  real(kind=dp)                :: cond, dk_dt, cp, dcp_dt, tbar, int_en

  ! calculate and store volume of each sub-control volume
  allocate(v_scv(4,no_elem))
  dz = one                 ! length of element in z-direction
  do ne=1,no_elem
     x_g(:) = x(nn(:,ne))  ! vector (1:4) of x-coord of elem ne
     y_g(:) = y(nn(:,ne))  ! vector (1:4) of y-coord of elem ne
     call vol_of_scv_quad(ne,x_g,y_g)
  end do

  jj = 0
  do k=1,no_el_blk           ! loop on element blocks
     mat_no_eb = mat_no(k)   ! mat no for this element block
     den = rho(mat_no_eb)

     write(6,"('element block = ',i3,' material no = ',i3)") k, mat_no_eb

     do j=1,no_el_in_blk(k)  ! loop on # elements in element block k
        jj = jj + 1          ! sequential element no
        ne = seq_el_no(jj)
        ! array of 4 nodal T's for this elem no ne
        T_gn(:) = T(nn(:,ne),1)
        Tbar = pt25*(sum(t_gn(:)))  ! avg of 4 nodal T's for this elem no
        ! get therm prop at avg elem T for elem mat_no_eb
        call get_prop(tbar, mat_no_eb, dum, dum, dum, dum, int_en) 
        rho_int_en = den*int_en
        ! loop on # of sub-control volumes
        vol(:) = v_scv(:,ne)       
        en_cnt(ne,1,:) = rho_int_en*vol(:) 
        write(6,"(6i4)") jj,ne,(nn(lnn,ne),lnn=1,4)
     end do
  end do

  ! fill in energy content at other time levels
  en_cnt(:,2,:) = en_cnt(:,1,:)
  en_cnt(:,3,:) = en_cnt(:,1,:)

end subroutine get_en_cnt
