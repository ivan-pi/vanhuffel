      SUBROUTINE QRQL(U, LDU, V, LDV, M, N, RANK, THETA, Q, E, INUL,
     *                TOL1, TOL2, WANTU, WANTV, IERR, IWARN)
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDU, LDV, M, N, RANK, IERR, IWARN
      DOUBLE PRECISION THETA, TOL1, TOL2
      LOGICAL WANTU, WANTV
C     .. Array Arguments ..
      DOUBLE PRECISION U(LDU,*), V(LDV,*), Q(*), E(*)
      LOGICAL INUL(*)
C     .. External Subroutines/Functions ..
      INTEGER NSINGV
      EXTERNAL CANCEL, ESTIM, QRSTEP, QLSTEP, NSINGV
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MIN
C     .. Local Scalars ..
      INTEGER I1, K, J, I, NUMEIG, ITER, MAXIT, P, R
      DOUBLE PRECISION X, SHIFT
      LOGICAL NOC12
C     .. Executable Statements ..
C
      IERR = 0
      P = MIN(M,N)
C
C     Estimate THETA (if not fixed by the user), and set R.
C
      IF (RANK .GE. 0) THEN
         J = P - RANK
         CALL ESTIM(Q, E, P, J, THETA, TOL1, TOL2, IWARN)
         IF (J .LE. 0) RETURN
         R = P - J
      ELSE
         R = 0
      END IF
C
      MAXIT = 50
      RANK = P
      DO 1 I = 1, P
         IF (INUL(I)) RANK = RANK - 1
   1  CONTINUE
      E(1) = 0.0D0
C
C     From now K is the smallest known index such that the subbidia-
C     gonals with indices > K belong to C1 or C2.
C     RANK = P - SUM(dimensions of known elements of C2).
C
      K = P
C     WHILE (C3 NOT EMPTY) DO
   2  IF (RANK .GT. R .AND. K .GT. 0) THEN
C        WHILE (K .GT. 0 .AND INUL(K)) DO
C
C        Search for the rightmost index of a subbidiagonal of C1 or C3.
C
   3     IF (K .GT. 0) THEN
            IF (INUL(K)) THEN
               K = K - 1
               GOTO 3
            END IF
         END IF
C        END WHILE 3
C
         IF (K .EQ. 0) RETURN
C
         ITER = 0
         NOC12 = .TRUE.
C        WHILE ((ITER < MAXIT) .AND. (No element of C1 or C2 found)) DO
   4     IF ((ITER .LT. MAXIT) .AND. NOC12) THEN
C
C           Search for negligible Q(I) or E(I).
C
            I = K
            X = ABS(Q(I))
            SHIFT = X
C           WHILE (ABS(Q(I)) > TOL2 .AND. ABS(E(I)) > TOL2) DO
   5        IF ((X .GT. TOL2) .AND. (ABS(E(I)) .GT. TOL2)) THEN
               I = I - 1
               X = ABS(Q(I))
               IF (X .LT. SHIFT) SHIFT = X
               GOTO 5
            END IF
C           END WHILE 5
C
C           Classify the subbidiagonal found.
C
            J = K - I + 1
            IF ((X .LE. TOL2) .OR. (K .EQ. I)) THEN
               NOC12 = .FALSE.
            ELSE
               NUMEIG = NSINGV(Q(I), E(I), J, THETA, TOL1, TOL2)
               IF (NUMEIG .GE. J .OR. NUMEIG .LE. 0) NOC12 = .FALSE.
            END IF
            IF (NOC12) THEN
               IF (SHIFT .GE. THETA) SHIFT = 0.0D0
               IF (ABS(Q(K)) .LE. ABS(Q(I))) THEN
                  CALL QRSTEP(U, LDU, V, LDV, Q, E, M, N, I, K, SHIFT,
     *                        WANTU, WANTV)
               ELSE
                  CALL QLSTEP(U, LDU, V, LDV, Q, E, M, N, I, K, SHIFT,
     *                        WANTU, WANTV)
               END IF
               ITER = ITER + 1
            END IF
            GOTO 4
         END IF
C        END WHILE 4
C
         IF (ITER .EQ. MAXIT) THEN
            IERR = 10
            RETURN
         END IF
C
         IF (X .LE. TOL2) THEN
C
C           Split at negligible diagonal element abs(Q(I)) <= TOL2.
C
            CALL CANCEL(U, LDU, V, LDV, Q, E, M, N, I, K, TOL2,
     *                  WANTU, WANTV)
            INUL(I) = .TRUE.
            RANK = RANK - 1
         ELSE
C
C           A negligible superdiagonal element abs(E(I)) <= TOL2 has
C           been found, the corresponding subbidiagonal belongs to
C           C1 or C2. Treat this subbidiagonal.
C
            IF (J .GE. 2) THEN
               IF (NUMEIG .EQ. J) THEN
                  DO 6 I1 = I, K
                     INUL(I1) = .TRUE.
   6              CONTINUE
                  RANK = RANK - J
                  K = K - J
               ELSE
                  K = I - 1
               END IF
            ELSE
               IF (X .LE. (THETA + TOL1)) THEN
                  INUL(I) = .TRUE.
                  RANK = RANK - 1
               END IF
               K = K - 1
            END IF
         END IF
         GOTO 2
      END IF
C     END WHILE 2
      RETURN
C *** Last line of QRQL ***********************************************
      END
