      SUBROUTINE ESTIM(Q, E, N, L, THETA, TOL1, TOL2, IWARN)
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER N, L, IWARN
      DOUBLE PRECISION THETA, TOL1, TOL2
C     .. Array Arguments ..
      DOUBLE PRECISION Q(N), E(N)
C     .. External Subroutines/Functions ..
      INTEGER NSINGV
      DOUBLE PRECISION DAMIN
      EXTERNAL DAMIN, NSINGV
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX
C     .. Local Scalars ..
      INTEGER NUM, I, NUMZ
      DOUBLE PRECISION TH, H1, H2, Y, Z, SUMZ
C     .. Executable Statements ..
C
      IWARN = 0
      IF ((L .LT. 0) .OR. (L .GT. N)) RETURN
C
C     Step 1: Initialization of THETA.
C             -----------------------
      IF (L .EQ. 0) THETA = 0.0D0
      IF (THETA .LT. 0.0D0) THEN
         IF (L .EQ. 1) THEN
C
C           An upper bound which is close if S(N-1) >> S(N):
C
            THETA = DAMIN(Q, N, 1)
         ELSE
C
C           An experimentally established estimate which is good if
C           S(N-L) >> S(N-L+1):
C
            THETA = ABS(Q(N-L+1))
         END IF
      END IF
C
C     Step 2: Check quality initial estimate THETA.
C             ------------------------------------
      NUM = NSINGV(Q, E, N, THETA, TOL1, TOL2)
      IF (NUM .EQ. L) RETURN
C
C     Step 3: Initialization starting values for bisection method.
C             ---------------------------------------------------
C     Let S(i), i=1,...,N, be the singular values of J in decreasing
C     order. Then, the computed Y and Z will be such that
C        (number of S(i) <= Y) < L < (number of S(i) <= Z).
C
      IF (NUM .LT. L) THEN
         Y = THETA
         TH = ABS(Q(1))
         Z = TH
         NUMZ = N
         DO 10 I = 2, N
            H1 = ABS(E(I))
            H2 = ABS(Q(I))
            SUMZ = MAX(TH + H1, H2 + H1)
            IF (SUMZ .GT. Z) Z = SUMZ
            TH = H2
   10    CONTINUE
      ELSE
         Z = THETA
         Y = 0.0D0
         NUMZ = NUM
      END IF
C
C     Step 4: Bisection method for finding the upper bound on the L
C             smallest singular values of the bidiagonal.
C             ------------------------------------------
C     A sequence of subintervals [Y,Z] is produced such that
C            (number of S(i) <= Y) < L < (number of S(i) <= Z).
C     NUM : number of S(i) <= TH,
C     NUMZ: number of S(i) <= Z.
C
C     WHILE ((NUM .NE. L) .AND. (Z-Y) .GT. TOL1) DO
   20 IF ((NUM .NE. L) .AND. (Z-Y) .GT. TOL1) THEN
         TH = (Y + Z)/2.0D0
         NUM = NSINGV(Q, E, N, TH, TOL1, TOL2)
         IF (NUM .LT. L) THEN
            Y = TH
         ELSE
            Z = TH
            NUMZ = NUM
         END IF
         GOTO 20
      END IF
C     END WHILE 2
C
C     If (Z - Y) <= TOL1, then at least two singular values of J lie in
C     the interval [Y,Z] within a distance < TOL1 from each other.
C     S(N-L) ans S(N-L+1) are then assumed to coincide.
C     L is increased, and a warning is given.
C
      IF (NUM .NE. L) THEN
         THETA = Z
         L = NUMZ
         IWARN = 1
      ELSE
         THETA = TH
      END IF
      RETURN
C *** Last line of ESTIM **********************************************
      END
