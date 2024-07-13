
subroutine natl_conv_corr(geom, len_eff, c_t_w, c_t_inf, c_nu_a, &
     c_alpha_a, c_cond_a, c_ra_no, c_h_conv)
  ! natural convection correlation
  ! evaluate experimental correlations for Nusselt no as function
  ! of Rayleigh no, as taken from Incropera, DeWitt, Bergman, & Lavine,
  ! 6th ed., 2007
  use mod_precision
  use mod_constants
  implicit none
  save
  character(len=7), intent(in)  :: geom
  real(dp),    intent(in)  :: len_eff
  complex(dp), intent(in)  :: c_t_w, c_t_inf, c_nu_a, c_alpha_a, &
                                   c_cond_a, c_ra_no
  complex(dp), intent(out) :: c_h_conv
  ! local variables
  complex(dp)              :: c_a, c_b, c_c
  real(dp)                 :: ra_no
  ! h in W/m^2-K
  ra_no = real(c_ra_no,dp)     ! real part of complex Ra
  if(geom == "horiz") then
     if(ra_no <= 1.e07_dp) then
        c_h_conv = (c_cond_a/len_eff)*0.54_dp*(c_ra_no)**pt25
     else
        c_h_conv = (c_cond_a/len_eff)*0.15_dp*(c_ra_no)**(one/three)
     end if

  else if(geom == "vert") then
     if(ra_no <= 1.e09_dp) then
        c_a = 0.670_dp*c_ra_no**(one/four)
        c_b = (one + 0.492_dp*c_alpha_a/c_nu_a)**(nine/16._dp)
        c_c = c_b**(four/nine)
        c_h_conv = (c_cond_a/len_eff)*(0.68_dp + c_a/c_c)
     else
        c_a = 0.387_dp*c_ra_no**(one/six)
        c_b = (one + 0.492_dp*c_alpha_a/c_nu_a)**(nine/16._dp)
        c_c = c_b**(eight/27._dp)
        c_h_conv = (c_cond_a/len_eff)*(0.825_dp + c_a/c_c)**2
     end if
  else if(geom == "verif") then
     ! dummy ss verification problem for convection on front face
     ! and specified T of backface. compare to analytical solution
     ! of single nonlinear algebraic equation for surface T
     c_h_conv = twelve*(c_t_inf - c_t_w)**pt25
  end if

end subroutine natl_conv_corr
