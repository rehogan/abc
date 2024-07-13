module mod_rad_encl_lib

contains

subroutine mat_mul_diag_dxa(n,d_vec,a,prod)
  ! multiply two matricees D and A where matrix D has non-zero diagonal
  ! elements and zero everywhere else
  ! this routine should be much more efficient than using intrinsic 
  ! F90 routines to perform the matrix multiplication DxA since most 
  ! elements of D are zero
  ! input:
  !   N     = number of surfaces
  !   D_vec = N-vector that represent diagonal elements of (NxN) matrix
  !           D(1,1), D(2,2), ..., D(N,N) with zero for off diagonal terms
  !   A     = (NxN) matrix

  ! output:
  !   PROD  = (NxN) matrix resulting from DxA

  use mod_precision
  use mod_constants
  implicit none                                        
  integer      , intent(in)                  :: n      ! (# surfaces)
  real(kind=dp), intent(in),  dimension(n,n) :: a      ! matrix a(n,n)
  real(kind=dp), intent(in),  dimension(n)   :: d_vec  ! vector (n)
  real(kind=dp), intent(out), dimension(n,n) :: prod   ! [d][a]

  ! local variables
  INTEGER                                    :: I, J

  ! compute matrix product [D]*[A] where D is a (nxn) matrix with
  ! zeros everywhere except the diagonal and [A] is a full matrix. 
  ! Only the diagonal elements are stored as
  !    D_vec(1) = D(1,1)
  !    D_vec(2) = D(2,2)
  !    ... 
  !    D_vec(N) = D(N,N)

  do i=1,n                           ! loop over rows
     do j=1,n                        ! loop over columns
        prod(i,j) = d_vec(i)*a(i,j)
     end do
  end do

end subroutine mat_mul_diag_dxa

subroutine encl(np,eps,f,U_mat,s,tw,sigma,whichex)
  !    assembly of 4 sided radiation enclosure equations
  !    radiative flux is linear in black body emissive power
  !    input:
  !      np      = no of surfaces in enclosure
  !      eps     = vector of surface emittance eps(n)
  !      f       = matrix of view factors f(n,n)
  !      s       = vector of side lengths (+ diagonals) s(6)
  !      tw      = vector of wall temperatures tw(4)
  !      sigma   = Stefan-Boltzmann constant
  !      whichex = title of example problem

  use mod_precision
  use mod_constants
  implicit none


  integer,  intent(in)                   :: np       ! no surfaces in encl
  real(dp), intent(in), dimension(np)    :: eps, tw
  real(dp), intent(in), dimension(np,np) :: f
  real(dp), intent(in), dimension(np,np) :: U_mat
  real(dp), intent(in), dimension(np+2)  :: s        ! side lengths + 2 diag
  real(dp), intent(in)                   :: sigma
  character(len=42), intent(in)          :: whichex

  real(dp)                   :: row_ex, dum, qxatot, sing_flag
  real(dp), dimension(np)    :: b, eb, qxa, qr, res, sum, t
  real(dp), dimension(np,np) :: d, c, cs
  INTEGER                    :: i, j, k, l, m, n, nsurf 
  integer,  dimension(np)    :: indx 

  character(len=3) :: debug
  !    SI units
  real*8 :: rank = 0.0d0

  debug = 'yes'
  nsurf = np

  do i=1,nsurf
     t(i) = tw(i) + rank   ! just in case input T was F and you want R
     eb(i) = sigma*t(i)*t(i)*t(i)*t(i) ! black body emissive power
  end do

  ! using known surface temperatures from problem specification, compute
  ! radiative fluxes to machine precision
  ! this gives the known flux vector to drive the iterative solution to 
  ! [Jac]{DT} = - {Res}, {Res} = {q_r} - [U]{e_b}

  call mat_vec_prod(4,4,U_mat(1,1),eb(1),qr) ! {q_r} = [U]{e_b}
  call newton(sigma,T,Tw,U_mat,qr,res)

  write(6,"('   I     RES              T         error')")
  do i=1,nsurf
     write(6,"(x,i3,es13.5, 2x, 2es13.5)") i,res(i),T(i),abs(t(i) - tw(i))
  end do

  !    Radiant Flux Results
  write(6,"('**Radiant Flux Results')")
  do i=1,nsurf
     write(6,"(2x,i3,es13.5)") i,qr(i)
  end do

  !    sum the flux times area in enclosure
  qxatot = zero
  do i=1,4
     qxa(i) = s(i)*qr(i)
     qxatot = qxatot + qxa(i)
  end do
  write(6,"('**(Radiant flux) x (area) + sum (qxA)')")
  write(6,"(2x,es12.5)") (qxa(i), i=1,4), qxatot

end subroutine encl

subroutine get_encl_u_mat(U_mat)
  ! {q_r} = [U]{e_b}, [U] = U_mat
  use mod_precision
  use mod_constants
  use mod_rad_encl
  implicit none
  real(dp),  intent(out),  dimension(ns,ns) :: U_mat  ! {q_r} = [U]{e_b}
  ! local variables
  integer :: i, j
  integer,                 dimension(ns)    :: indx   ! vector 
  real(dp),                dimension(ns)    :: b      ! vector b(ns)
  real(dp),                dimension(ns,ns) :: cinv   ! matrix cs(ns,ns)
  real(dp),                dimension(ns,ns) :: cs     ! matrix cs(ns,ns)
  real(dp),                dimension(ns,ns) :: ident_mat ! identity matrix
  real(dp),                dimension(ns,ns) :: ImF
  real(dp),                dimension(ns,ns) :: mat_1
  real(dp)                                  :: row_ex, sing_flag

  ! generate identity matrix
  i_mat(:,:) = zero              ! fill array with zero's
  do i=1,ns
     i_mat(i,i) = one            ! place 1's along diagonal
  end do

  eps_mat(:,:) = zero            ! initialize array with zero's
  ! fill emittance matrix with diagonal entries
  do i=1,ns
     eps_mat(i,i) = eps(i)       ! enter diagonal entries
  end do

  ! reflectance matrix
  rho_mat(:,:) = zero            ! fill array with zero's
  do i=1,ns
     rho_vec(i) = one - eps(i)
     rho_mat(i,i) = rho_vec(i) ! enter diagonal entries
  end do

  ! compute nxn matrix [I] - [F][rho]
  call mat_mul_Axdiag_D(ns,F,rho_vec,fxrho) ! compute product [F][rho]
  imfxrho = i_mat - fxrho                  ! [I] - [F][rho]

  cs = imfxrho            ! save the original matrix [C] as [CS]
  call ludcmp(imfxrho,ns,ns,indx,row_ex,sing_flag)
  ! calculate inverse of matrix [C]
  do i=1,ns
     b(:) = zero
     b(i) = one
     call lubksb(imfxrho,ns,ns,indx,b)
     cinv(:,i) = b(:)
  end do
  ! check cinv*cs
  ident_mat = matmul(cs,cinv)

  write(6,"(' inverse matrix ([I] - [F][eps])^-1  ')")
  do i=1,4
     write(6,"(4es13.5)") (cinv(i,j),j=1,4)
  end do

  write(6,"(' [I] identity matrix check  ')")
  do i=1,4
     write(6,"(4es13.5)") (ident_mat(i,j),j=1,4)
  end do

  ! ([I] - [F]) matrix
  ImF = I_mat - F
  write(6,"(' ([I] - [F]) matrix  ')")
  do i=1,4
     write(6,"(4es13.5)") (ImF(i,j),j=1,4)
  end do

  mat_1 = matmul(cinv,ImF)
  write(6,"(' mat_1 matrix  ')")
  do i=1,4
     write(6,"(4es13.5)") (mat_1(i,j),j=1,4)
  end do

  call mat_mul_diag_dxa(ns,eps,mat_1,U_mat) ! U_mat returned

  write(6,"(' U_mat  ')")
  do i=1,4
     write(6,"(4es13.5)") (U_mat(i,j),j=1,4)
  end do

end subroutine get_encl_u_mat

subroutine lubksb(a,n,np,indx,b)

  !this subroutine performs the Lower-Upper forward and back substitution
  ! see ``numerical recipes'', pg 36, for more details

  use mod_precision
  use mod_constants
  implicit none
  integer,  intent(in)                      :: n, np
  real(dp), intent(in),    dimension(np,np) :: a
  integer,  intent(in),    dimension(n)     :: indx
  real(dp), intent(inout), dimension(n)     :: b
  real(dp)                                  :: sum
  integer                                   :: i, ii, j, ll

  ii = 0
  do 20 i = 1 , n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if (ii .ne. 0) then
        do 10 j = ii , i - 1
           sum = sum - a(i,j) * b(j)
10         continue
           elseif (sum .ne. 0.) then
           ii = i
        endif
        b(i) = sum
20      continue
        do 40 i = n , 1 , -1
           sum = b(i)
           if (i .lt. n) then
              do 30 j = i + 1 , n
                 sum = sum - a(i,j) * b(j)
30               continue
              endif
              b(i) = sum / a(i,i)
40            continue

end subroutine lubksb

subroutine ludcmp(a,n,np,indx,row_ex,flag)

  ! this subroutine performs the lower-upper decompositon of a 
  ! coefficient matrix
  ! see "numerical recipes", pg 35, for more details
  ! converted to f90 by ben blackwell

  use mod_precision
  use mod_constants
  implicit none
  real(dp), intent(inout), dimension(np,np) :: a
  integer,  intent(in)                      :: n, np
  integer,  intent(inout), dimension(n)     :: indx
  real(dp), intent(out)                     :: row_ex, flag                               
  integer                                   :: i, imax, j, k, max_facet
  real(dp)                                  :: aamax, dum, tiny, sum
  parameter (max_facet = 100, tiny = 1.0d-30)
  integer,              dimension(max_facet) :: vv
  if (n .gt. max_facet) then
     write( 6,* ) 'subroutine ludcmp -- redimension max_facet'
     stop
  end if
  row_ex = one
  flag = one
  do 20 i = 1 , n
     aamax = zero
     do 10 j = 1 , n
        if (dabs(a(i,j)) .gt. aamax) aamax = dabs(a(i,j))
10      continue
        if (aamax .eq. 0.) then
           write (*,*) 'singular matrix.'
           flag = -one
           return
        end if
        vv(i) = one/aamax
20      continue
        do 90 j = 1 , n
           if (j .gt. 1) then
              do 40 i = 1 , j - 1
                 sum = a(i,j)
                 if (i .gt. 1) then
                    do 30 k = 1 , i - 1
                       sum = sum - a(i,k) * a(k,j)
30                     continue
                       a(i,j) = sum
                    endif
40                  continue
                 endif
                 aamax = 0.
                 do 60 i = j , n
                    sum = a(i,j)
                    if (j .gt. 1) then
                       do 50 k = 1 , j - 1
                          sum = sum - a(i,k) * a(k,j)
50                        continue
                          a(i,j) = sum
                       endif
                       dum = vv(i) * dabs(sum)
                       if (dum .ge. aamax) then
                          imax = i
                          aamax = dum
                       endif
60                     continue
                       if (j .ne. imax) then
                          do 70 k = 1 , n
                             dum = a(imax,k)
                             a(imax,k) = a(j,k)
                             a(j,k) = dum
70                           continue
                             row_ex = -row_ex
                             vv(imax) = vv(j)
                          endif
                          indx(j) = imax
                          if (j .ne. n) then
                             if (a(j,j) .eq. 0.) a(j,j) = tiny
                             dum = 1. / a(j,j)
                             do 80 i = j + 1 , n
                                a(i,j) = a(i,j) * dum
80                              continue
                             endif
90                           continue
                             if (a(n,n) .eq. 0.) a(n,n) = tiny
                             return
end subroutine ludcmp

subroutine mat_vec_prod(ni,nj,a,b,c)    

  !    multiplication of matrix times vector 
  !    c(ni)=a(ni x nj)*b(nj)        

  !    input:

  !         ni, nj, a, b
  !         matrix a has ni rows and nj columns
  !         matrix b has nj rows

  !    output:

  !         matrix c has ni rows

  use mod_precision
  use mod_constants
  implicit none
  integer,  intent(in)                    :: ni, nj
  real(dp), intent(in),  dimension(ni,nj) :: a
  real(dp), intent(in),  dimension(nj)    :: b
  real(dp), intent(out), dimension(ni)    :: c
  integer                                 :: i, j

  do i = 1 , ni
     c(i) = zero
     do j = 1 , nj
        c(i) = c(i) + a(i,j)*b(j)
     end do
  end do

end subroutine mat_vec_prod

subroutine mat_mul_axdiag_d(n,a,d_vec,prod)
  ! multiply two matricees A and D where matrix D has non-zero diagonal
  ! elements and zero everywhere else
  ! this routine should be much more efficient than using intrinsic 
  ! F90 routines to perform the matrix multiplication AxD since most 
  ! elements of D are zero
  ! input:
  !   N = number of surfaces
  !   A     = (NxN) matrix
  !   D_vec = N-vector that represent diagonal elements of (NxN) matrix
  !       D(1,1), D(2,2), ..., D(N,N) with zero for off diagonal terms

  ! output:
  !   PROD = (NxN) matrix resulting from AxD

  use mod_precision
  use mod_constants
  implicit none                                        
  integer      , intent(in)                  :: n      ! (# surfaces)
  real(kind=dp), intent(in),  dimension(n,n) :: a      ! matrix a(n,n)
  real(kind=dp), intent(in),  dimension(n)   :: d_vec  ! vector (n)
  real(kind=dp), intent(out), dimension(n,n) :: prod   ! 

  ! local variables
  INTEGER                                    :: I, J

  ! compute matrix product [A]*[D] where A is (nxn) and D is a (nxn) matrix with
  ! zeros everywhere except the diagonal. Only the diagonal elements are stored
  ! as D(1), D(2), ..., D(N)

  do i=1,n                           ! loop over rows
     do j=1,n                        ! loop over columns
        prod(i,j) = a(i,j)*d_vec(j)
     end do
  end do

end subroutine mat_mul_axdiag_d

subroutine newton(sigma,T_nup1,T_w,U,qr,res)

  ! newton's method for computing T-vector from q-vector
  ! calculate residual for the enclosure equations
  ! intrinsic f90 matrix functions used in some instances
  ! {R} = {q_r} - [U]*{e_b} 

  USE MOD_PRECISION
  USE MOD_CONSTANTS
  IMPLICIT NONE                                        
  integer i, iter, j, n
  PARAMETER(N=4)
  REAL(DP), INTENT(IN)                    :: SIGMA
  REAL(DP), INTENT(INOUT), DIMENSION(N)   :: T_nup1 ! vector T_nup1(N)
  REAL(DP), INTENT(IN),    DIMENSION(N)   :: T_w    ! vector T_w(N)
  REAL(DP), INTENT(IN),    DIMENSION(N,N) :: U      ! matrix U(N,N)
  REAL(DP), INTENT(IN),    DIMENSION(N,1) :: qr     ! vector QR(N,1)
  REAL(DP), INTENT(OUT),   DIMENSION(N,1) :: res    ! vector RES(N,1)
  ! local variables
  REAL(DP),                DIMENSION(N,1) :: eb     ! vector EB(N,1)
  REAL(DP),                DIMENSION(N)   :: B      ! vector B(N)
  REAL(DP),                DIMENSION(N)   :: T_nu   ! vector T_nu(N)
  REAL(DP),                DIMENSION(N,10):: Ti     ! matrix RES(N,10)
  REAL(DP),                DIMENSION(N,N) :: jac    ! matrix JAC(N,N)
  REAL(DP),                DIMENSION(N)   :: sigt3  ! vector sigma*T**3
  REAL(DP),                DIMENSION(N)   :: dt     ! vector T correction
  REAL(DP),                DIMENSION(N,1) :: negres ! vector -RES(N,1)
  integer,                 DIMENSION(N)   :: indx   ! vector 
  integer, parameter                      :: jmax = 9
  real(dp),                dimension(n,n) :: ident_mat ! identity matrix
  REAL(DP) :: row_ex, sing_flag, res_norm, sol_norm, rel_norm, L_2_norm, sum


  ! perturb "exact" initial guess to check for convergence
  do i=1,n
     T_nup1(i) = T_nup1(i) + six*five**2*(-1)**i
  end do

  write(6,"(' Iter no    L_2 Norm on T')")

  iter = 0
  do j=1,jmax          ! begin iteration loop
     iter = iter + 1
     T_nu = T_nup1     ! set old solution to previously converged solution
     do i=1,n
        sigt3(i) = sigma*t_nu(i)*t_nu(i)*t_nu(i)
        eb(i,1) = sigt3(i)*t_nu(i)        ! vector {e_b}
        sigt3(i) = four*sigt3(i)          ! factor of 4 from deriv
     end do
  !  res = {q_r} - [U]*{e_b}
     negres = -(qr - matmul(U,eb))         ! negative of res vector
     call mat_mul_axdiag_d(n,U,sigt3,jac)
     jac = -jac    ! compute Jacobian
     call ludcmp(jac,n,n,indx,row_ex,sing_flag)
     call lubksb(jac,n,n,indx,negres)
     sum = zero
     do i=1,n
        dt(i) = negres(i,1)                  ! temperature correction
        t_nup1(i) = t_nu(i) + dt(i)          ! T for iter j+1
        sum = sum + (t_w(i) - t_nup1(i))**2
     end do
     L_2_norm = sqrt(sum/float(n))
     write(6,"(i5,6x,es12.5)") j,L_2_norm

     ! compute norms of residuals and solution
     res_norm = sqrt(dot_product(dt,dt)/float(n))
     sol_norm = sqrt(dot_product(t_nup1,t_nup1)/float(n))
     rel_norm = res_norm/sol_norm
     if(rel_norm < v_small) exit
  end do
  write(6,"('  n = ',i3,2x,'rel norm = ',es13.5)") n,rel_norm

end subroutine newton

subroutine vfac(w,h,f,s)
  !    input:
  !      w = width of enclosure
  !      h = height of enclosure
  !    output:
  !      f = view factor array
  !      s = side lengths plus diagonals for 4-sided enclosure

  !    calculate view factors for 4 sided enclosure using 
  !    cross string method 
  !    there are two forms of the equations:
  !    one for adjacent (touching) sides and 
  !    one for (non-touching) opposing sides

  USE MOD_PRECISION
  USE MOD_CONSTANTS
  IMPLICIT NONE                                        
  real(dp), intent(in)                  :: w, h
  real(dp), intent(out), dimension(6)   :: s
  real(dp), intent(out), dimension(4,4) :: f

  real(dp), dimension(4)   :: x, y
  real(dp), dimension(4,4) :: a, c
  real(dp)                 :: sum
  integer                  :: i, j
  x(:) = zero
  y(:) = zero
  x(2) = w
  x(3) = w
  y(3) = h
  y(4) = h
  f(:,:) = zero       ! fill VF array with zeros
  s(1) = len(1,2,x,y) ! side lengths plus diagonals of quad
  s(2) = len(2,3,x,y)
  s(3) = len(3,4,x,y)
  s(4) = len(4,1,x,y)
  s(5) = len(1,3,x,y) ! diagonals
  s(6) = len(2,4,x,y)
  write(6,"('**4 side lengths, two diagonals')")
  write(6,"(6es12.5)") s(1),s(2),s(3),s(4),s(5),s(6)
  write(6,"('**View Factor Matrix, row sum check, 1 - sum')")
  f(1,2) = vfun3(1,2,3)
  f(1,3) = vfun4(1,2,3,4)
  f(1,4) = vfun3(2,1,4)
  !    row 1
  sum = zero
  do i=1,4
     sum = sum + f(1,i)
  end do
  write(6,"(4f12.5,es13.5)") f(1,1), f(1,2), f(1,3), f(1,4), one-sum
  !    row 2
  f(2,3) = vfun3(2,3,4)
  f(2,4) = vfun4(2,3,4,1)
  f(2,1) = vfun3(3,2,1)
  sum = zero
  do i=1,4
     sum = sum + f(2,i)
  end do
  write(6,"(4f12.5,es13.5)") f(2,1), f(2,2), f(2,3), f(2,4), one-sum
  !    row 3
  f(3,4) = vfun3(3,4,1)
  f(3,1) = vfun4(3,4,1,2)
  f(3,2) = vfun3(4,3,2)
  sum = zero
  do i=1,4
     sum = sum + f(3,i)
  end do
  write(6,"(4f12.5,es13.5)") f(3,1), f(3,2), f(3,3), f(3,4), one-sum
  !    row 4
  f(4,1) = vfun3(4,1,2)
  f(4,2) = vfun4(4,1,2,3)
  f(4,3) = vfun3(1,4,3)
  sum = zero
  do i=1,4
     sum = sum + f(4,i)
  end do
  write(6,"(4f12.5,es13.5)") f(4,1), f(4,2), f(4,3), f(4,4), one-sum

  ! temporary hard wire for 1-d radiation gap
  f(1:4,1:4) = zero
  f(1,3) = one
  f(3,1) = one
  write(6,"(' View Factors for 1-d radiation gap')")
  do i=1,4
     write(6,"(4f5.1)") (f(i,j),j=1,4)
  end do

  contains
  function len(i,j,x,y)
    implicit none
    integer, intent(in)   :: i, j
    real(dp), intent(in), dimension(4) :: x, y
    real(dp)              :: len
    len = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2)
  end function len

  real function vfun3(i,j,k)
    integer, intent(in) :: i, j, k
    vfun3 = (len(i,j,x,y) + len(j,k,x,y) - len(k,i,x,y))/(2.*len(i,j,x,y))
  end function vfun3

  real function vfun4(i,j,k,m)
    integer, intent(in) :: i, j, k, m
    vfun4 = (len(i,k,x,y) + len(j,m,x,y) -(len(j,k,x,y) + len(i,m,x,y)))&
         /(2.*len(i,j,x,y))
  end function vfun4

end subroutine vfac

subroutine bc_del_encl_rad(U)
  use mod_constants
  use mod_files
  use mod_global_parameters
  use mod_nodes_elements
  use mod_rad_encl
  implicit none
  real(dp), intent(in), dimension(4,4)   :: U

  ! local variables
  integer                  :: i, j, m, mr
  integer,  dimension(2)   :: nn_b
  real(dp)                 :: t_i, t_j, coef
  real(dp), dimension(4)   :: q_r, tbar, e_b, sij
  real(dp), dimension(4,4) :: dqdotdT
  real(dp), dimension(4,4) :: dTbar4dT
  integer,  dimension(2)   :: surf_ij_m
  real(dp), dimension(2)   :: res
  real(dp), dimension(2,4) :: jac
  interface
     ! required since this function exists outside mod_rad_encl_lib.f90
     function mat_vec_mult(A, v) result (w) ! required since this exists outside mod_rad_encl_lib.f90
       use mod_precision
       implicit none
       real(kind=dp), dimension(:,:), intent(in) :: A
       real(kind=dp), dimension(:),   intent(in) :: v
       real(kind=dp), dimension( SIZE(A,1) )     :: w
     end function mat_vec_mult

  end interface

  mr(i,j,nhbp1) = j + nhbp1 - i   ! local to global no transformation

  ! calculate side temperatures as average of two nodal 
  ! temperatures on each enclosure side
  do m=1,4
!     nn_er(1:4) = nn_encl_rad(m,1:4) ! create 1-d nodal array
     i = surf_ij(m,1)
     j = surf_ij(m,2)
     t_i = t_li(i)
     t_j = t_li(j)
     Tbar(m) = pt5*(t_i + t_j) ! avg T for edge i-j
     e_b(m) = sigma_si*tbar(m)*tbar(m)*tbar(m)*tbar(m) ! emissive power
     sij(m) = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2) ! side length (area)
  end do
  ! calculate heat flux for all enclosure sides
  ! [U] is the matrix defined by {q_r} = [U]{e_b} and was calculated earlier
  q_r = mat_vec_mult(U,e_b)  ! in encl eqn, outflow flux is +
  ! debug information
  write(6,"(' Tbar = ',4es13.5)") (tbar(m),m=1,4)
  write(6,"(' q_r = ',4es13.5)") (q_r(m),m=1,4)

  ! calculate matrix dTbar^4/dT, Tbar is array of enclosure side temps
  !                              T    is array of enclosure node temps
  dTbar4dT = zero    ! initialize array to zeros
  ! fill in non-zero terms, refer to notes for matrix structure
  ! diagonal terms, note sigma included for scaling
  do m=1,4
     dTbar4dT(m,m) = sigma_si*Tbar(m)**3
  end do
  ! super diagonal terms
  do m=1,3
     dTbar4dT(m,m+1) = dTbar4dT(m,m)
  end do
  ! lower left corner term
  dTbar4dT(4,1) = dTbar4dT(4,4)

  ! dqdot"/dT sensitivity of face flux to nodal temperatures
  ! note that sigma is built into dTbar^4/dT
  dqdotdT = two*matmul(U,dTbar4dT)

  write(6,"(4i4)") (nn_encl_rad(j),j=1,4)
  ! loop on no of sides in enclosure
  do m=1,4
!!$     q_r(1) = 14674.24042_dp
!!$     q_r(2) = zero
!!$     q_r(3) = -q_r(1)
!!$     q_r(4) = zero
     surf_ij_m(:) = surf_ij(m,:)
     res = zero
     res(1) = q_r(m)*sij(m)*pt5
     res(2) = res(1)
     jac(1,1:4) = pt5*sij(m)*dqdotdt(m,1:4)
     jac(2,1:4) = jac(1,1:4)
!     jac = zero
     
     ! merge elem bc contribution into global matrix
     write(6,"(3i4)") m, (surf_ij_m(j), j=1,2)
     call el_encl_rad_bc_matrix_to_global_matrix &
          (2, 4, surf_ij_m, nn_encl_rad, jac, res)
  end do

end subroutine bc_del_encl_rad

subroutine el_encl_rad_bc_matrix_to_global_matrix &
     (nnps, nn_s, surf_ij_m, nn_er, jac, res)

  ! nnps      = no nodes per side (2 for now)
  ! nn_s      = no sides per enclosure (4 for now)
  ! surf_ij_m = node pair (i,j) for side m
  ! nn_er     = array of nodes defining enclosure
  ! jac       = Jacobian, (2 x 4 for encl rad bc)
  ! res       = Residual, (2 for encl rad bc)

  use mod_precision
  use mod_constants
  use mod_nodes_elements
  implicit none
  integer,        intent(in)                       :: nnps
  integer,        intent(in)                       :: nn_s
  integer,        intent(in), dimension(nnps)      :: surf_ij_m
  integer,        intent(in), dimension(nn_s)      :: nn_er
  real(kind=dp),  intent(in), dimension(nnps,nn_s) :: jac
  real(kind=dp),  intent(in), dimension(nnps)      :: res

  ! local variables
  integer                                          :: i, ii, j, jj, kk, mr

  mr(i,j,nhbp1) = j + nhbp1 - i   ! local to global no transformation

  ! merge encl rad bc contribution into global matrix
  ! the row index stays the same but the column index is shifted
  ! because of the banded storage structure of the [C] matrix

  do ii = 1,nnps           ! ii = local node no
     i = surf_ij_m(ii)     ! node no of bc node
     b(i) = b(i) - res(ii) ! residual contribution
     do jj = 1,nn_s
        j = nn_er(jj)
        kk = mr(i,j,nhbp1)
        c(i,kk) = c(i,kk) + jac(ii,jj)  ! global matrix contribution
!        write(6,"(' i, j, i, kk, ii, jj',6i4)") i,j,i,kk,ii,jj
     end do
  end do

end subroutine el_encl_rad_bc_matrix_to_global_matrix
  
end module mod_rad_encl_lib

