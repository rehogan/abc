
subroutine air_props
  ! fill the air thermal property arrays using the data in
  ! module mod_air_table
  use mod_constants
  use mod_air_properties
  use mod_air_table
  implicit none
  integer :: j
  tt(:) = air_table(:)%tt                       ! T, K
  trho_a(:) = air_table(:)%trho_a               ! rho, Kg/m^3
  tcp_a(:) = air_table(:)%tcp_a*1.e3_dp         ! c_p, J/Kg-K
  tmu_a(:) = air_table(:)%tmu_a*1.e-6_dp        ! visc mu, N-s/m^2
  tcond_a(:) = air_table(:)%tcond_a*1.e-3_dp    ! cond k, W/m-K
  tnu_a(:) = tmu_a(:)/trho_a(:)                 ! visc nu, m^2/s
  talpha_a(:) = tcond_a(:)/(trho_a(:)*tcp_a(:)) ! diffusivity, m^2/s

  ! print statements below are left over from debugging

  ! print the needed values
!  print *, ' needed values for h calculation'
!  print *, '     TT       Talpha_a     Tnu_a       Tcond_a'
!  write(6,"(4es12.5)") (tt(j), talpha_a(j), tnu_a(j), tcond_a(j),j=1,21)

! print the original table with scaled values
!  print *, 'original table'
!  print *, ' TT     Trho_a   Tcp_a    T_mu_a  Tcond_a'
!  write(6,"(f6.0, 2x, f6.3, 2x, f7.3, 2x, f7.2, 2x, f7.2)") &
!       (tt(j), trho_a(j), tcp_a(j)*1.e-3_dp, tmu_a(j)*1.e6_dp, &
!       tcond_a(j)*1.e03_dp,j=1,21)

end subroutine air_props
