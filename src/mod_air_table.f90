
module mod_air_table
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
  implicit none

!  integer, parameter :: dp = kind (0.0d0)

  type air_table_t
    real(dp) :: tt
    real(dp) :: trho_a
    real(dp) :: tcp_a
    real(dp) :: tmu_a
    real(dp) :: tcond_a
  end type

  !                 K     Kg/m^3    KJ/Kg-K   N-s/m^2*e6 W/m-K*e3
  type(air_table_t), parameter :: air_table(21) = [  &
     air_table_t (225._dp, 1.572_dp,  1.006_dp, 14.67_dp,  20.20_dp), &
     air_table_t (230._dp, 1.537_dp,  1.006_dp, 14.94_dp,  20.62_dp), &
     air_table_t (235._dp, 1.505_dp,  1.006_dp, 15.20_dp,  21.04_dp), &
     air_table_t (240._dp, 1.473_dp,  1.005_dp, 15.47_dp,  21.45_dp), &
     air_table_t (245._dp, 1.443_dp,  1.005_dp, 15.73_dp,  21.86_dp), &
     air_table_t (250._dp, 1.413_dp,  1.005_dp, 15.99_dp,  22.27_dp), &
     air_table_t (255._dp, 1.386_dp,  1.005_dp, 16.25_dp,  22.68_dp), &
     air_table_t (260._dp, 1.359_dp,  1.005_dp, 16.50_dp,  23.08_dp), &
     air_table_t (265._dp, 1.333_dp,  1.005_dp, 16.75_dp,  23.48_dp), &
     air_table_t (270._dp, 1.308_dp,  1.006_dp, 17.00_dp,  23.88_dp), &
     air_table_t (275._dp, 1.284_dp,  1.006_dp, 17.26_dp,  24.28_dp), &
     air_table_t (280._dp, 1.261_dp,  1.006_dp, 17.50_dp,  24.67_dp), &
     air_table_t (285._dp, 1.240_dp,  1.006_dp, 17.74_dp,  25.06_dp), &
     air_table_t (290._dp, 1.218_dp,  1.006_dp, 17.98_dp,  25.47_dp), &
     air_table_t (295._dp, 1.197_dp,  1.006_dp, 18.22_dp,  25.85_dp), &
     air_table_t (300._dp, 1.177_dp,  1.006_dp, 18.46_dp,  26.24_dp), &
     air_table_t (305._dp, 1.158_dp,  1.006_dp, 18.70_dp,  26.63_dp), &
     air_table_t (310._dp, 1.139_dp,  1.007_dp, 18.93_dp,  27.01_dp), &
     air_table_t (315._dp, 1.121_dp,  1.007_dp, 19.15_dp,  27.40_dp), &
     air_table_t (320._dp, 1.103_dp,  1.007_dp, 19.39_dp,  27.78_dp), &
     air_table_t (325._dp, 1.086_dp,  1.008_dp, 19.63_dp,  28.15_dp)  &
]

end module mod_air_table
