
subroutine flux_natl_conv(geom, t_w, t_inf, qdot_pp, &
     h_conv, dh_conv_dt_w, dh_conv_dt_inf)
  ! thermophysical properties of air for 225 < T < 325 K
  ! from Gebhart, Heat Conduction and Mass Diffusion, McGraw-Hill, p600
  ! original source was R. C. Weast, ed(1970),
  ! Handbook of Tables for Applied Engineering Science,
  ! CRC Press, Boca Raton, FL
  ! units are SI, but cp, mu, and cond have been scaled
  ! tabular values are cp*e-3 mu*e6 and cond*e3
  ! working units are J, Kg, m, s
  ! the leading "t" in the variable names is suggestive of
  ! "tabular", i.e. tt = tabular temperature
  ! inside this routine, T must be in K

  use mod_precision
  use mod_constants
  use mod_bc
  use mod_air_properties
  use mod_lin_interp
  implicit none
  save
  character(len=7), intent(in) :: geom
  real(dp),         intent(in) :: t_w, t_inf
  real(dp), intent(out)        :: qdot_pp,  h_conv, & 
                                  dh_conv_dt_w, dh_conv_dt_inf
  integer                 :: jlo
  integer                 :: j
  ! rough shop dimensions: Dz = 36', W= 17', H = 13'
  real(dp), parameter     :: rm_h = 13._dp     ! room height, ft
  real(dp), parameter     :: rm_w = 17._dp     ! room width, ft
  real(dp), parameter     :: rm_d = 36._dp     ! room depth, ft, pipes this dir
  real(dp)                :: f, ra_no        ! interp fac, Ray no
  real(dp)                :: beta,  len_eff
  real(dp)                :: del_t
  real(dp)                :: nu_a, alpha_a, cond_a
  real(dp)                :: t_eff
  complex(dp)             :: c_nu_a, c_alpha_a, c_cond_a, c_h_conv, &
                                  c_qdot_pp, c_ra_no, c_t_w, c_t_inf

  real(dp)                :: dqdot_dt_inf, dqdot_dt_w

  interface

     complex function c_prop(jlo,f,del_t,p,tt)
       use mod_precision
       use mod_constants
       implicit none
       integer,       intent(in)               :: jlo
       real(dp), intent(in)               :: f, del_t
       real(dp), intent(in), dimension(*) :: p, tt
     end function c_prop

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

     subroutine natl_conv_corr(geom, len_eff, c_t_w, c_t_inf, c_nu_a, &
          c_alpha_a, c_cond_a, c_ra_no, c_h_conv)
       use mod_precision
       use mod_constants
       implicit none
       save
       character(len=7), intent(in)  :: geom
       real(dp),    intent(in)  :: len_eff
       complex(dp), intent(in)  :: c_t_w, c_t_inf, c_nu_a, c_alpha_a, &
            c_cond_a, c_ra_no
       complex(dp), intent(out) :: c_h_conv
     end subroutine natl_conv_corr

  end interface

  if(first_time /= 'no') then
     jlo = 1
     del_t = 1.e-16_dp
     first_time = 'no'
  end if

  ! scale factors have been removed and everything has
  ! consistent units of J, m, s, K

  ! interpolate in the above air tables


!  t_eff = pt5*(t_inf + t_w) + conv_c_to_k    ! convert to K
  t_eff = pt5*(t_inf + t_w)                  ! consistent K
  call find(tt(1), 21, t_eff, jlo, f)        ! locate place in tables
  beta = one/t_eff
  len_eff = rm_w/(two*(one + rm_w/rm_d))           ! area/perimeter, ft
  len_eff = len_eff*conv_ft_to_m             ! ft -> m

  ! kinematic viscosity (nu) m^/s
   nu_a = lin_interp_fun(f,tnu_a(jlo),tnu_a(jlo+1))
   c_nu_a = c_prop(jlo,f,del_t,tnu_a(1),tt(1))

  ! thermal diffusivity (alpha), m^2/s
  alpha_a = lin_interp_fun(f,talpha_a(jlo),talpha_a(jlo+1))
  c_alpha_a = c_prop(jlo,f,del_t,talpha_a(1),tt(1))

  ! thermal conductivity (k), W/m-K
  cond_a = lin_interp_fun(f,tcond_a(jlo),tcond_a(jlo+1))
  c_cond_a = c_prop(jlo,f,del_t,tcond_a(1),tt(1))

  ! d(qdot")/dT_w
  c_t_w = cmplx(t_w, del_t, dp)
  c_t_inf = cmplx(t_inf, zero, dp)
  c_ra_no = g_si*beta*(c_t_w - c_t_inf)*len_eff**3/(c_nu_a*c_alpha_a)
  call natl_conv_corr(geom, len_eff, c_t_w, c_t_inf, c_nu_a, &
     c_alpha_a, c_cond_a, c_ra_no, c_h_conv)
  c_qdot_pp = c_h_conv*(c_t_w - c_t_inf)
  h_conv = real(c_h_conv)
  dh_conv_dt_w = aimag(c_h_conv)/del_t
  qdot_pp = real(c_qdot_pp)
  dqdot_dt_w = aimag(c_qdot_pp)/del_t
  write(6,"('   geom         h         dh/dT_w      qdot_pp   d(qdot_pp)/dT_w')")
  write(6,"(2x,a7,4es13.5)") geom, h_conv, dh_conv_dt_w, qdot_pp, dqdot_dt_w

  ! d(qdot")/dT_inf
  c_t_w = cmplx(t_w, zero, dp)
  c_t_inf = cmplx(t_inf, del_t, dp)
  c_ra_no = g_si*beta*(c_t_w - c_t_inf)*len_eff**3/(c_nu_a*c_alpha_a)
  call natl_conv_corr(geom, len_eff, c_t_w, c_t_inf, c_nu_a, &
     c_alpha_a, c_cond_a, c_ra_no, c_h_conv)
  c_qdot_pp = c_h_conv*(c_t_w - c_t_inf)
  h_conv = real(c_h_conv)
  dh_conv_dt_inf = aimag(c_h_conv)/del_t
  qdot_pp = real(c_qdot_pp)
  dqdot_dt_inf = aimag(c_qdot_pp)/del_t
  write(6,"('   geom         h         dh/dT_inf    qdot_pp   d(qdot_pp)/dT_inf')")
  write(6,"(2x,a7,4es13.5)") geom, h_conv, dh_conv_dt_inf, qdot_pp, dqdot_dt_inf

end subroutine flux_natl_conv
