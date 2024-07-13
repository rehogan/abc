
subroutine get_prop(tbar, mat_no_eb, cond, dk_dt, cp, dcp_dt, int_en)
  use mod_precision
  use mod_constants
  use mod_lin_interp
  use mod_mat_prop
  use mod_time_params
  implicit none
  real(dp), intent(inout)  :: tbar
  integer,  intent(in)     :: mat_no_eb
  real(dp), intent(out)    :: cond, dk_dt, cp, dcp_dt, int_en
  real(dp)                 :: cond_1, cond_2, f, onemf, &
                                   t_1, t_2, k_1, k_2
  integer                       :: jlo, m
  interface
     subroutine find (iv, nmax, iv_val, jlo, f)
       use mod_precision
       use mod_constants
       implicit none
       integer,  intent(inout) ::  jlo
       integer,  intent(in)    ::  nmax
       real(dp), intent(in)    :: iv_val
       real(dp), intent(out)   :: f
       real(dp), dimension(*)  :: iv
     end subroutine find

  end interface
  m = mat_no_eb ! shorten variable name, matl no elem block

  ! thermal conductivity interpolation
  call find(TT_cond(1,m),no_cond_ent(m),tbar,jlo,f)
!  write(6,"(' f = ',f6.3, ' k_1 = ',f6.0, ' k_2 = ',f6.0)") &
!      f,t_cond(jlo,m),t_cond(jlo+1,m)
  cond = lin_interp_fun(f,t_cond(jlo,m),t_cond(jlo+1,m))
  dk_dt = (t_cond(jlo+1,m) - t_cond(jlo,m))/(tt_cond(jlo+1,m) - tt_cond(jlo,m))
!  write(6,"(' cond = ',f7.3, ' dk/dT = ',f7.3)") cond, dk_dt
  if(ss == 'yes') return     ! ss does not require heat capacity
  ! heat capacity interpolation
  call find(TT_cp(1,m),no_cp_ent(m),tbar,jlo,f)
!  write(6,"(' f = ',f6.3, ' cp_1 = ',f6.0, ' cp_2 = ',f6.0)") &
!       f,t_cp(jlo,m),t_cp(jlo+1,m)
  cp = lin_interp_fun(f,t_cp(jlo,m),t_cp(jlo+1,m))
  dcp_dt = (t_cp(jlo+1,m) - t_cp(jlo,m))/(tt_cp(jlo+1,m) - tt_cp(jlo,m))
!  write(6,"(' cp = ',f7.1,' dcp/dT = ',f7.1)") cp, dcp_dt
  ! internal energy related stuff
  ! parabolic interpolation for internal energy
  int_en = t_int_en(jlo,m) + t_cp(jlo,m)*(tbar - tt_cp(jlo,m)) + &
       pt5*dcp_dt*(tbar - tt_cp(jlo,m))**2
!  write(6,"(' int_en = ', es12.3)") int_en
end subroutine get_prop
