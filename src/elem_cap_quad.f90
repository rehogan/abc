
subroutine elem_cap_quad(ne,rho,cp,int_en,el_jac,res)

  ! calculate the element capacitance contribution corresponding to
  ! sub-control volumes for quadralateral element
  ! input:
  ! ne         = element no
  ! xg, yg = global (x,y) coordinates for 4-node quad
  ! rho    = density of elem ne, Kg/m^3
  ! cp     = specific heat of elem ne, J/Kg-K
  ! int_en = internal energy of elem ne, J/Kg

  ! output: 
  ! jac    = Jacobian, dR_i/dT_j
  ! res    = residual, R_i
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
  real(dp), intent(out), dimension(4,4) :: el_jac
  real(dp), intent(out), dimension(4)   :: res
  integer                               :: i, j, nscv
  real(dp)                              :: rho_int_en, coef_1, coef_2
  real(dp),            dimension(4)     :: aa, vol

  ! energy content en_cnt(i,j,k), i=el no, j=time index, k=local node no
  rho_int_en = rho*int_en
! loop on # of sub-control volumes
  vol(:) = v_scv(:,ne)       
  en_cnt(ne,np1,:) = rho_int_en*vol(:) 
!  do nscv = 1,4 ! loop on # of sub-control volumes
!     vol(nscv) = v_scv(nscv,ne)
!     en_cnt(ne,np1,nscv) = rho_int_en*vol(nscv)
!  end do

  ! residual calculation
  coef_1 = one/(theta*dtime(np1))
  res(:) = coef_1*(en_cnt(ne,np1,:) - en_cnt(ne,n,:)) ! residual

  ! Jacobian calculation, each row has identical entries
  coef_2 = coef_1*rho*cp/four
  aa(:) = coef_2*vol(:)
  do i=1,4
     do j=1,4
        el_jac(i,j) = aa(i)    ! elements of Jacobian
     end do
  end do

!  write(6,"('  element Jacobian matrix, ne = ',i3)") ne
!  do i=1,4
!     write(6,"(4es13.3)") (el_jac(i,j), j=1,4)
!  end do

end subroutine elem_cap_quad
