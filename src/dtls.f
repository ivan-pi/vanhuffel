      SUBROUTINE DTLS(C, LDC, M, N, L, S, X, LDX, WRK, RANK, TOL1, TOL2,
     *                COMPRT, IERR, IWARN)
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDC, M, N, L, RANK, LDX, IERR, IWARN
      DOUBLE PRECISION TOL1, TOL2
      CHARACTER COMPRT*1
C     .. Array Arguments ..
      DOUBLE PRECISION C(LDC,*), S(N+L), WRK(N+L+M), X(LDX,*)
C     .. External Subroutines/Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DCOPY, DDOT, DAXPY, DQRDC, DSVDC, HOUSH, TR2
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX, MIN, SQRT
C     .. Local Scalars ..
      INTEGER R1, NL, P, N1, NJ, J1, K, I, J, MC, NL1
      DOUBLE PRECISION SMAX, TEMP, SMAX2
      LOGICAL CRANK, CTOL, ZERO
C     .. Local Arrays ..
      INTEGER DUMAR1(1)
      DOUBLE PRECISION DUMAR2(1,1)
C     .. Parameters ..
      DOUBLE PRECISION NUMEPS
      PARAMETER (NUMEPS = 0.1D-15)
C     .. Executable Statements ..
C
C     Determine whether RANK and/or TOL1 is to be computed.
C
      CRANK = .TRUE.
      CTOL = .TRUE.
      IF ((COMPRT .EQ. 'N') .OR. (COMPRT .EQ. 'n')) THEN
         CRANK = .FALSE.
         CTOL = .FALSE.
      ELSE
         IF ((COMPRT .EQ. 'R') .OR. (COMPRT .EQ. 'r')) CTOL = .FALSE.
         IF ((COMPRT .EQ. 'T') .OR. (COMPRT .EQ. 't')) CRANK = .FALSE.
      END IF
C
      IERR = 0
      IWARN = 0
      NL = N + L
      K = MAX(M, NL)
      P = MIN(M, N)
      IF (M .LT. 1) IERR = 1
      IF (N .LT. 1) IERR = 2
      IF (L .LT. 1) IERR = 3
      IF (LDC .LT. K) IERR = 4
      IF (LDX .LT. N) IERR = 5
      IF (.NOT.CRANK) THEN
         IF (RANK .GT. P) IERR = 6
      END IF
      TOL1 = MAX(TOL1, NUMEPS)
      TOL2 = MAX(TOL2, NUMEPS)
      IF (IERR .NE. 0) RETURN
C
C     Initialize the solution matrix X.
C
      DO J = 1, L
         DO I = 1, N
            X(I,J) = 0.0D0
         END DO
      END DO
C
C     Subroutine DTLS solves a set of linear equations by a Total Least
C     Squares Approximation.
C
C     Step 1:
C     1.a. : if M .GE. 5*NL/3, then transform [A   ;B   ] into upper
C                                               M,N  M,L
C            triangular form R by Householder transformations.
C
      MC = M
      IF (3*M .GE. 5*NL) THEN
         CALL DQRDC(C, LDC, M, NL, S, DUMAR1, WRK, 0)
         DO J = 1, NL
            J1 = J + 1
            DO I = J1, NL
               C(I,J) = 0.0D0
            END DO
         END DO
         MC = NL
      END IF
C                                     T
C     1.b. : compute the SVD of  U S V  of [A;B] (or R).
C
      CALL DSVDC(C, LDC, MC, NL, S, WRK, DUMAR2, 1, C, LDC, WRK(NL+1),
     *           1, IERR)
      IF (IERR .NE. 0) THEN
         IERR = IERR + 1000
         RETURN
      END IF
C
C     Step 2: Compute the rank of the approximation [A+DA;B+DB].
C
      SMAX = TOL1
      IF (CTOL) SMAX = SQRT(2.0D0*K) * SMAX
      SMAX2 = SMAX**2
      IF (CRANK) THEN
         RANK = P
C        WHILE (RANK .GT. 0) .AND. (S(RANK) .LE. SMAX) DO
    5    IF (RANK .GT. 0) THEN
            IF (S(RANK) .LE. SMAX) THEN
               RANK = RANK - 1
               GO TO 5
            END IF
         END IF
C        END WHILE
      END IF
C
C     Step 3: Compute the Householder matrix Q and matrices F and Y
C     such that F is nonsingular.
C
C     REPEAT
C
C        Adjust the rank if S(RANK) has multiplicity > 1.
C
   6     R1 = RANK + 1
C        WHILE (RANK .GT. 0) .AND. (S(RANK)**2 - S(R1)**2 .LE. SMAX2) DO
   7     IF (RANK .GT. 0) THEN
            IF (S(RANK)**2 - S(R1)**2 .LE. SMAX2) THEN
               RANK = RANK - 1
               IWARN = 1
               GO TO 7
            END IF
         END IF
C        END WHILE
         IF (RANK .EQ. 0) RETURN
         R1 = RANK + 1
C
C        Compute the Householder matrix Q and matrices F and Y.
C
         NL1 = MAX(N, R1) + 1
         ZERO = .FALSE.
         I = NL
C        WHILE ((.NOT.ZERO) .AND. (I .GE. NL1)) DO
   8     IF ((.NOT.ZERO) .AND. (I .GE. NL1)) THEN
            K = I - RANK
            CALL DCOPY(K, C(I,R1), LDC, WRK, 1)
            CALL HOUSH(WRK, K, K, TOL2, ZERO, TEMP)
            IF (.NOT.ZERO) THEN
               CALL TR2(C, LDC, WRK, TEMP, 1, I, RANK, K)
            END IF
            I = I - 1
            GOTO 8
         END IF
C        END WHILE
         N1 = N + 1
         IF (ZERO .OR. ABS(C(N1,N1)) .LE. TOL2) THEN
            RANK = RANK - 1
            IWARN = 2
            GO TO 6
         END IF
C     UNTIL ((.NOT.ZERO) .AND. (ABS(C(N1,N1) .GT. TOL2))
C
C     Step 4: Solve X F = -Y by forward elimination,
C             (F is upper triangular).
C
      CALL DAXPY(N, -1.0D0/C(N1,N1), C(1,N1), 1, X, 1)
      DO J = 2, L
         NJ = N + J
         TEMP = C(NJ,NJ)
         J1 = J - 1
         DO I = 1, N
            X(I,J) = -(C(I,NJ) + DDOT(J1, C(N1,NJ), 1, X(I,1), LDX))
     *              /TEMP
         END DO
      END DO
      RETURN
C *** Last line of DTLS ***********************************************
      END
