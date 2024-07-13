

subroutine band(n,nfb,c,b,x)

  !     Ben Blackwell and Roy Hogan adapted this unsymmetric band solver
  !     "band algorithm for unsymmetric matrices in finite
  !     element and semi-discrete methods," t. c. gopalakrishnan and
  !     a. b. palaniappan, international journal for numerical methods
  !     in engineering, vol. 18, pp1197-1211 (1982).
  !     n     = number of unknowns
  !     nhb   = number of bands above (below) diagonal but not including
  !     the diagonal (read as no of half band)
  !     nhbp1 = nhb + 1, no of bands above (below) diagonal plus diagonal
  !     nhbm1 = nhb - 1 (read as nhb  minus 1)
  !     nfb   = 2*nhb + 1, full band width (read as no of full band)
  !     nk .ge. n + nhb and must match dimension in calling routine
  !     nl .ge. nfb and must match dimension in calling routine
  !     nk .ge. n + nhb and must match dimension in calling routine
  !     nl .ge. nhb+1 and must match dimension in calling routine
  !     c     = coefficient matrix with dimensions matching calling
  !             routine
  !     b     = rhs vector with dimensions matching calling routine
  !     x     = solution vector with dimensions matching calling routine
  !     if a(i,j) is a square coefficient matrix that multiplies the
  !     unknown temperature vector, it is stored in the c(i,k) matrix
  !     which is a rectangular array (n by nfb). they are related by
  !     c(i,j-i+nhb+1) = a(i,j); note that the row index is the same
  !     in both the rectangular and band storage schemes. Only the
  !     column number changes. the diagonal element is c(i,nhb+1) = a(i,i).
  !
  ! example of band storage scheme for nhb = 2, nfb = 2*nhb +  1 = 5
  ! row
  !  1     3   4   5   6   7   8   9  10  11  12 
  !  2     2   3   4   5   6   7   8   9  10  11
  !  3     1   2   3   4   5   6   7   8   9  10
  !  4     0   1   2   3   4   5   6   7   8   9
  !  5    -1   0   1   2   3   4   5   6   7   8
  !  6    -2  -1   0   1   2   3   4   5   6   7
  !  7    -3  -2  -1   0   1   2   3   4   5   6
  !  8    -4  -3  -2  -1   0   1   2   3   4   5
  !  9    -5  -4  -3  -2  -1   0   1   2   3   4
  ! 10    -6  -5  -4  -3  -2  -1   0   1   2   3

  use mod_precision
  implicit none
  integer,                         intent(in)    :: n, nfb    !! nk, nl
  real(kind=dp),  dimension(:,:),  intent(inout) :: c
  real(kind=dp),  dimension(:),    intent(out)   :: b, x
  integer                :: nhb, nhbm1, nhbp1
  integer                :: i, im1, ke, l, k, jbm1, jb, ik1, kk, j, m
  real(kind=dp)          :: f, s
  !!
!  write(6,*) size(c,1), size(c,2), size(b), size(x)    !! check size of c, b, and x
  !!
  nhb   = (nfb-1)/2
  nhbm1 = nhb-1
  nhbp1 = nhb+1
  !!
  do i = 2 , n
     im1 = i - 1
     ke = nhb
     !
     !     elimination        ++
     !
     do l = 1 , nhbm1
        if (i .eq. (n - l + 1)) ke = l
     end do
     do k = 1 , ke
        jbm1 = nhbp1 - k
        jb = jbm1 + 1
        ik1 = i + k - 1
        if (c(ik1,jbm1) .ne. 0.0) then
           kk = nfb - k
           f = c(ik1,jbm1) / c(im1,nhbp1)
           b(ik1) = b(ik1) - b(im1) * f
           do j = jb , kk
              c(ik1,j) = c(ik1,j) - c(im1,j + k) * f
           end do
        endif
     end do
  end do
  !
  !     back substitution
  !
  do j = 1 , n
     i = n - j + 1
     s = 0.
     do m = 1 , nhb
        s = s + x(i + m) * c(i,nhb + m + 1)
     end do
     x(i) = (b(i) - s) / c(i,nhbp1)
  end do
  !
  !     all done
  !
end subroutine band
