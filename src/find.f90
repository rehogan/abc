
SUBROUTINE FIND (IV, NMAX, IV_VAL, JLO, F)

!***********************************************************************

!     PURPOSE:

!       This routine will determine the index J such that for the array
!       IV(.), IV(J) < IV_VAL < IV(J+1) with JLO as the initial guess for J.
!       This is a preliminary calculation for doing linear
!       interpolation.

!       Example: The temperature dependence of thermal conductivity is
!       input in a tabular form with INCREASING temperature. IV(.)
!       represents the array of temperatures and IV_VAL represents the
!       value of temperature for which we desire the thermal
!       conductivity.

!     TEMP(1)     COND(1)
!     TEMP(2)     COND(2)
!        .           .
!        .           .
!        .           .
!     TEMP(NMAX)  COND(NMAX)
!


!     INPUTS:

!       IV( )  = tabular array of independent variable
!       NMAX   = no of elements in IV( ) array
!       IV_VAL = value of independent variable for interpolation
!       JLO    = starting guess for index J; JLO = 1 or JLO = N are
!                acceptable if you do not have any idea where to start

!     OUTPUTS:
!
!       JLO   = final value of index J
!       F     = fractional distance of IV_VAL in the interval
!               IV(J) < IV_VAL < IV(J+1); used in subsequent interpolations
!               of the form DV_VAL = DV(J)*(1.0 - F) + DV(J+1)*F
!               where DV() is the dependent variable array and DV_VAL
!               is the dependent variable value corresponding to IV_VAL
!

!
!     DESCRIPTION:

!       This routine was modeled after a more general one presented in
!       Numerical Recipes, Cambridge University Press. The one
!       presented here should run faster than the one presented in
!       Numerical Recipes since it does not have to also worry about
!       a descending ordered array in addition to the ascending array.
!       It is assumed that the table is arranged in ascending
!       order; if this is not the case, then this routine will fail.
!       Precautions are taken for the case when IV_VAL <= IV(1) or
!       IV_VAL >= IV(N). If IV_VAL <= IV(1), then JLO = 0; if IV_VAL >= IV(N),
!       then JLO = N - 1.

!       The first step is to insure that you are inside table limits.
!       The table is divided into an upper and a lower portion
!       with IV(JLO) being the dividing point.

!         Lower portion:  IV(1)   to IV(JLO)
!         Upper portion:  IV(JLO) to IV(N)

!       The next step is to determine if you are in the upper or
!       lower portion of the table. Once this decision has been made,
!       the next step is to determine the region IV(JLO) < IV_VAL < IV(JHI)
!       where JHI is not necessarily JLO + 1. The first try is always
!       IV(JLO) < IV_VAL < IV(JLO+1); if this fails, the increment between
!       JLO and JHI is doubled. The doubling is repeated until
!       IV(JLO) < IV_VAL < IV(JHI). Once the bracketing has been successful,
!       the interval is successively halved until  IV(JLO) < IV_VAL < IV(JLO+1)
!       has been satisfied.


!     AUTHOR:

!       Ben Blackwell
!       Blackwell Consulting
!       Corrales, NM 87048
!       Phone: (505) 897-5090   e-mail: bblackwell13@comcast.net

!     REFERENCES:

!       NONE

!     LOCAL VARIABLES:


!     ROUTINES CALLED:

!       none

!     REVISION HISTORY (MM/DD/YY):
!
!       01/14/97 - BFB - Subroutine initially written
!       11/27/02 - BFB - Converted from f77 to f90
!       08/11/05 - BFB - Modified to perform linear extrapolation
!                        out bottom and top of table
!       08/25/08 - BFB - added independent and dependent variable nomenclature
!
!***********************************************************************
!
  use mod_precision
  USE MOD_CONSTANTS
  IMPLICIT NONE
  INTEGER,  INTENT(INOUT) ::  JLO
  INTEGER,  INTENT(IN)    ::  NMAX
  REAL(DP), INTENT(IN)    :: IV_VAL
  REAL(DP), INTENT(OUT)   :: F
  REAL(DP), DIMENSION(*)  :: IV
! local
  INTEGER :: INC, JHI, JM

!**************** FIRST EXECUTABLE STATEMENT OF FIND ******************

! is first guess below lower table limit or above upper table limit?

  IF(JLO .LE. 0 .OR. JLO .GT. NMAX) THEN ! You are outside of table limits.
     JLO = 1
  ENDIF

  ! Check to see if IV_VAL .LT. IV(1) or IV_VAL .GT. IV(NMAX); if so, extrapolate

  IF (IV_VAL <= IV(1)) THEN ! you are below the bottom table limit
     JLO = 1
     F = (IV_VAL - IV(JLO))/(IV(JLO+1) - IV(JLO)) ! F < 0, extrapolate below
!     NUM_EXTRAP = NUM_EXTRAP + 1                  ! increment extrap counter
     F = ZERO                                     ! extrapolate with zero slope
     RETURN
  ELSE IF (IV_VAL >= IV(NMAX)) THEN ! you are above the upper table limit
     JLO = NMAX - 1
     F = (IV_VAL - IV(JLO))/(IV(JLO+1) - IV(JLO)) ! F > 1, extrapolate above
!     NUM_EXTRAP = NUM_EXTRAP + 1                  ! increment extrap counter
     F = ONE                                      ! extrapolate with zero slope
     RETURN
  ENDIF

  ! You are inside the table limits so locate where we are in the
  ! table.

  INC=1
  IF(IV_VAL .GE. IV(JLO)) THEN ! You are in upper portion of table.
10   CONTINUE

     JHI = MIN(JLO + INC,NMAX)
     IF (IV_VAL >= IV(JHI)) THEN
        JLO = JHI
        INC = INC + INC !     Double the increment
        GO TO 10
     ENDIF
     GO TO 30
  ELSE         ! You are in lower portion of table.

20   CONTINUE
     IF (IV_VAL < IV(JLO)) THEN
        JHI = JLO
        JLO = MAX(JHI - INC,1)
        INC = INC + INC ! Double the increment
        GO TO 20
     ENDIF
  ENDIF

  !     We now have the location bracketed, IV(JLO) <= IV_VAL <= IV(JHI)
  !     but JHI may not be JLO + 1; if necessary, divide
  !     interval into half to home in on precise location.

30 CONTINUE
  IF((JHI - JLO) == 1) THEN ! We are all finished, compute F and return
     F = (IV_VAL - IV(JLO))/(IV(JHI) - IV(JLO))
     RETURN
  ELSE                        ! Cut the interval in half and try again.

     JM = (JHI + JLO)*pt5
     IF(IV_VAL > IV(JM)) THEN
        JLO = JM
     ELSE
        JHI = JM
     ENDIF
     GO TO 30
  ENDIF

END SUBROUTINE FIND
