      SUBROUTINE PTLS(C, LDC, M, N, L, RANK, THETA, X, LDX, Q, INUL,
     *                WRK, IWRK, LWRK, TOL1, TOL2, IERR, IWARN)
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDC, M, N, L, RANK, LDX, IERR, IWARN
      DOUBLE PRECISION THETA, TOL1, TOL2
C     .. Array Arguments ..
      INTEGER IWRK(N+L)
      DOUBLE PRECISION C(LDC,*), X(LDX,*), Q(*), WRK(*)
      LOGICAL INUL(N+L), LWRK(N+L)
C     .. External Subroutines/Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DAXPY,DDOT,DCOPY,DQRDC, BIDIAG, CANCEL, HOUSH, INIT, QRQL
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX, MIN
C     .. Local Scalars ..
      INTEGER I, J, P, NL, PP1, J1, MC, NJ, K, II, I1, N1, MNL1
      DOUBLE PRECISION TEMP, HH, INPROD
      LOGICAL ZERO
C     .. Local Arrays ..
      INTEGER DUMAR1(1)
      DOUBLE PRECISION DUMAR2(1,1)
C     .. Executable Statements ..
C
      IERR = 0
      IWARN = 0
      IF (M .LT. 1) IERR = 1
      IF (N .LT. 1) IERR = 2
      IF (L .LT. 1) IERR = 3
      IF (LDC .LT. MAX(M,N+L)) IERR = 4
      IF (LDX .LT. N) IERR = 5
      IF (RANK .GT. MIN(M,N)) IERR = 6
      IF ((RANK .LT. 0) .AND. (THETA .LT. 0.0D0)) IERR = 7
      IF (TOL1 .LT. 0.0D0) IERR = 8
      IF (TOL2 .LT. 0.0D0) IERR = 9
      IF (IERR .NE. 0) RETURN
C
C     Initializations.
C
      NL = N + L
      MNL1 = M + NL + 1
      P = MIN(M,NL)
      DO I = 1, P
         INUL(I) = .FALSE.
         LWRK(I) = .FALSE.
      END DO

      PP1 = P + 1
      DO I = PP1, NL
         INUL(I) = .TRUE.
         LWRK(I) = .FALSE.
      END DO
C
      DO J = 1, L
         DO I = 1, N
            X(I,J) = 0.0D0
         END DO
      END DO
C
C     Subroutine PTLS solves a set of linear equations by a Total Least
C     Squares Approximation, based on the Partial SVD.
C
C     Step 1: Bidiagonalization phase
C             -----------------------
C     1.a): If M >= 5*(N+L)/3, transform C into upper triangular form R
C           by Householder transformations.
C
      MC = M
      IF (3*M .GE. 5*NL) THEN
         CALL DQRDC(C, LDC, M, NL, Q, DUMAR1, WRK, 0)
         DO J = 1, NL
            J1 = J + 1
            DO I = J1, NL
               C(I,J) = 0.0D0
            END DO
         END DO
         MC = NL
      END IF
C
C     1.b): Transform C (or R) into bidiagonal form Q using Householder
C           transformations.
C
      CALL BIDIAG(C, LDC, MC, NL, Q, Q(P+1), WRK)
C
C     Store the Householder transformations performed onto the rows of C
C     in the last storage locations of the work array WRK.
C
      MC = MIN(NL-2,M)
      IF (MC .GT. 0) THEN
         K = MNL1
         DO II = 1, MC
            J = MC - II + 1
            NJ = NL - J
            CALL DCOPY(NJ, C(J,J+1), LDC, WRK(K), 1)
            K = K + NJ
         END DO
      END IF
C
C     1.c): Initialize the right singular base matrix V with the identi-
C           ty matrix (V overwrites C).
C
      CALL INIT(C, LDC, NL, NL)
C
C     1.d): If M < N+L, bring the bidiagonal Q to M by M by cancelling
C           its last superdiagonal element using Givens rotations.
C
      PP1 = P
      IF (M .LT. NL) THEN
         PP1 = P + 1
         CALL CANCEL(DUMAR2, 1, C, LDC, Q, Q(PP1), P, PP1, PP1, PP1,
     *               TOL2, .FALSE., .TRUE.)
      END IF
C
C     REPEAT
C
C        Compute the Householder matrix Q and matrices F and Y such that
C        F is nonsingular.

C        Step 2: Partial diagonalization phase.
C                -----------------------------
C        Diagonalize the bidiagonal Q partially until convergence to
C        the desired right singular subspace.
C
   6     CALL QRQL(DUMAR2, 1, C, LDC, P, PP1, RANK, THETA, Q, Q(P+1),
     *             INUL, TOL1, TOL2, .FALSE., .TRUE., IERR, IWARN)
C
         IF (IERR.NE.0) RETURN
C
C        Step 3: Back transformation phase.
C                -------------------------
C        Apply the Householder transformations (stored in WRK) perfor-
C        med onto the rows of C during the bidiagonalization phase, to
C        the selected base vectors in the right singular base matrix
C        of C.
C
         DO I = 1, NL
            IF (INUL(I) .AND. (.NOT. LWRK(I))) THEN
               K = MNL1
               DO II = 1, MC
                  J = MC - II + 1
                  NJ = NL - J
                  J1 = J + 1
                  IF (ABS(WRK(K)) .GT. TOL2) THEN
                     TEMP = -DDOT(NJ, WRK(K), 1, C(J1,I), 1)/WRK(K)
                     CALL DAXPY(NJ, TEMP, WRK(K), 1, C(J1,I), 1)
                     K = K + NJ
                  END IF
               END DO
               LWRK(I) = .TRUE.
            END IF
         END DO
         IF (RANK.LE.0) RETURN
C
C        Step 4: Compute matrices F and Y using Householder transf. Q.
C                ------------------------
         K = 0
         DO I = 1, NL
            IF (INUL(I)) THEN
               K = K + 1
               IWRK(K) = I
            END IF
         END DO
C
         IF (K .LT. L) THEN
C
C           Rank TLS approximation is larger than min(M,N).
C
            IERR = 3
            RETURN
         END IF

         N1 = N + 1
         ZERO = .FALSE.
         I = NL
C
C        WHILE ((K > 1) .AND. (I > N) .AND. (.NOT.ZERO)) DO
   10    IF ((K.GT.1) .AND. (I.GT.N) .AND. (.NOT.ZERO)) THEN
            DO J = 1, K
               WRK(J) = C(I,IWRK(J))
            END DO
C
C           Compute Householder transformation.
C
            CALL HOUSH(WRK, K, K, TOL2, ZERO, TEMP)
            IF (.NOT.ZERO) THEN
C
C              Apply Householder transformation onto the selected base
C              vectors.
C
               DO I1 = 1, I
                  INPROD = 0.0D0
                  DO J = 1, K
                     INPROD = INPROD + WRK(J) * C(I1,IWRK(J))
                  END DO
                  HH = INPROD * TEMP
                  DO J = 1, K
                     J1 = IWRK(J)
                     C(I1,J1) = C(I1,J1) - WRK(J) * HH
                  END DO
               END DO
C
               K = K - 1
            END IF
            I = I - 1
            GOTO 10
         END IF
C        END WHILE 10
C
         IF (.NOT.ZERO) K = N1 - RANK
C
C        If F singular, lower the rank of the TLS approximation .
C
         IF (ABS(C(N1,IWRK(K))) .LE. TOL2) THEN
            RANK = RANK - 1
            IWARN = 2
            THETA = -1.0D0
            GO TO 6
         END IF
C     UNTIL ((.NOT.ZERO) .AND. (F nonsingular))
C
C     Step 5: Compute TLS solution
C             --------------------
C     Solve X F = -Y  by forward elimination  (F is upper triangular).

      NJ = IWRK(K)
      CALL DAXPY(N, -1.0D0/C(N1,NJ), C(1,NJ), 1, X, 1)

      DO J = 2, L
         J1 = J - 1
         NJ = IWRK(K+J1)
         TEMP = C(N+J,NJ)
         DO I = 1, N
            X(I,J) = -(C(I,NJ) + DDOT(J1,C(N1,NJ),1,X(I,1),LDX))/TEMP
         END DO
      END DO
      RETURN
C *** Last line of PTLS ***********************************************
      END
