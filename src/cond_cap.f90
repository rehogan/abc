
subroutine cond_cap
  ! conduction + capacitance w/o any sources but with 
  ! temp dependent properties
  use mod_precision
  use mod_constants
  use mod_mat_prop
  use mod_nodes_elements
  use mod_time_params
  implicit none
  integer            :: i, j, jj, k, ne, mat_no_eb
  integer, parameter :: nnpe = 4
  real(kind=dp),  dimension(4)   :: t_gn, x_g, y_g ! corner (x,y) ne
  real(kind=dp),  dimension(4,4) :: el_cond  ! el cond matrix
  real(kind=dp),  dimension(4,4) :: el_jac   ! Jac matrix for elem e
  real(kind=dp),  dimension(4)   :: res      ! residual vector elem i
  real(kind=dp)                  :: cond, dk_dt, cp, dcp_dt, tbar, int_en

  interface
     subroutine get_prop(tbar, mat_no_eb, cond, dk_dt, cp, dcp_dt, int_en)
       use mod_precision
       use mod_constants
       use mod_mat_prop
       implicit none
       real(kind=dp), intent(inout)  :: tbar
       integer,       intent(in)     :: mat_no_eb
       real(kind=dp), intent(out)    :: cond, dk_dt, cp, dcp_dt, int_en
     end subroutine get_prop

     subroutine elem_cond_quad_iso(ne, x_g, y_g, t_gn, cond, dk_dt, el_jac, res)
       use mod_precision
       use mod_constants
       implicit none
       save

       integer,       intent(in)                  :: ne
       real(kind=dp), intent(in),  dimension(4)   :: x_g, y_g, t_gn ! global coord elem # ne
       real(kind=dp), intent(in)                  :: cond, dk_dt     ! k of elem & dk/dT
       real(kind=dp), intent(out), dimension(4,4) :: el_jac          ! Jac matrix for elem e
       real(kind=dp), intent(out), dimension(4)   :: res            ! residual vector elem i
     end subroutine elem_cond_quad_iso

     subroutine el_matrix_to_global_matrix(ne, nnpe, el_jac, res)

       ! ne    = element no
       ! nnpe  = no nodes per element, 3 for triangle, 4 for quad

       use mod_precision
       use mod_constants
       use mod_nodes_elements
       implicit none
       integer,        intent(in)                       :: ne, nnpe
       real(kind=dp),  intent(in), dimension(nnpe,nnpe) :: el_jac
       real(kind=dp),  intent(in), dimension(nnpe)      :: res
     end subroutine el_matrix_to_global_matrix

     subroutine elem_cap_quad(ne,rho,cp,int_en,el_cap,res)
       use mod_precision
       use mod_constants
       use mod_nodes_elements
       use mod_time_params
       implicit none
       save
       integer,  intent(in)                  :: ne
       real(dp), intent(in)                  :: rho    ! element cond
       real(dp), intent(in)                  :: cp     ! elem spec heat
       real(dp), intent(in)                  :: int_en ! elem internal energy
       real(dp), intent(out), dimension(4,4) :: el_cap
       real(dp), intent(out), dimension(4)   :: res
     end subroutine elem_cap_quad
  end interface

  jj = 0                                    ! initialize seq elem counter
  elem_blk: do k = 1,no_el_blk              ! loop on element blocks
     mat_no_eb = mat_no(k)                  ! mat no for this element block
     elem_in_blk: do j = 1, no_el_in_blk(k) ! loop on elements in elem blk ieb
        jj = jj + 1                         ! sequential elem index
        ne = seq_el_no(jj)                  ! ne = global no of elem
        ! global corner coordinates of elem # ne
        x_g(:) = x(nn(:,ne))
        y_g(:) = y(nn(:,ne))
        ! array of 4 nodal T's for this elem no
        t_gn(:) = t_li(nn(:,ne))
        tbar = pt25*(sum(t_gn(:)))  ! avg of 4 nodal T's for this elem no
        ! get therm prop at avg elem T for elem mat_no_eb
        call get_prop(tbar, mat_no_eb, cond, dk_dt, cp, dcp_dt, int_en) 
        call elem_cond_quad_iso(ne, x_g, y_g, t_gn, cond, dk_dt, el_jac, res)
        call el_matrix_to_global_matrix(ne, nnpe, el_jac, res)

        if(ss == 'yes') cycle       ! this is ss problem
        ! element capacitance
        call elem_cap_quad(ne,rho(mat_no_eb),cp,int_en,el_jac,res)
        call el_matrix_to_global_matrix(ne, nnpe, el_jac, res)

     end do elem_in_blk

  end do elem_blk


end subroutine cond_cap
