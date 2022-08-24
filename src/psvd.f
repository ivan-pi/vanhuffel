      SUBROUTINE PSVD(A, LDA, M, N, RANK, THETA, U, LDU, V, LDV, Q,
     *                INUL, WRK, TOL1, TOL2, MODE, IERR, IWARN)
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDA, M, N, RANK, LDU, LDV, MODE, IERR, IWARN
      DOUBLE PRECISION THETA, TOL1, TOL2
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), U(LDU,*), V(LDV,*), Q(*), WRK(*)
      LOGICAL INUL(*)
C     .. External Subroutines/Functions ..
      EXTERNAL DCOPY, DQRDC, BIDIAG, CANCEL, INIT, QRQL, RESTOR
C     .. Intrinsic Functions ..
      INTRINSIC MAX, MIN, MOD
C     .. Local Scalars ..
      INTEGER J, P, K, I, J1, MA, NM1, NJ, PP1
      LOGICAL WANTU, WANTV, QR, ALL
C     .. Local Arrays ..
      INTEGER DUMAR1(1)
C     .. Executable Statements ..
C
      IWARN = 0
      IERR  = 0
      IF (M .LT. 1) IERR = 1
      IF (N .LT. 1) IERR = 2
      IF (LDA .LT. M) IERR = 3
      IF ((MODE .LT. 0) .OR. (MODE .GE. 100)) IERR = 11
      IF (TOL1 .LT. 0.0D0) IERR = 8
      IF (TOL2 .LT. 0.0D0) IERR = 9
C
C     Determine whether U and/or V is to be computed.
C
      ALL   = .FALSE.
      WANTU = .FALSE.
      WANTV = .FALSE.
      J = MOD(MODE,100)/10
      I = MOD(MODE,10)
      IF ((J.EQ.1 .AND. M.GT.N) .OR. (I.EQ.1 .AND. M.LT.N)) ALL = .TRUE.
      IF (J .NE. 0) WANTU = .TRUE.
      IF (I .NE. 0) WANTV = .TRUE.
C
      IF (WANTU .AND. (LDU .LT. M)) IERR = 4
      IF (WANTV .AND. (LDV .LT. N)) IERR = 5
      P = MIN(M,N)
      IF (RANK .GT. P) IERR = 6
      IF ((RANK .LT. 0) .AND. (THETA .LT. 0.0D0)) IERR = 7
      IF (IERR .NE. 0) RETURN
C
C     Initializations.
C
      IF (3*M .GE. 5*N) THEN
         QR = .TRUE.
      ELSE
         QR = .FALSE.
      END IF
C
      K = MAX(M,N)
      DO 1 I = 1, K
         INUL(I) = .FALSE.
   1  CONTINUE
      IF (ALL .AND. (.NOT.QR)) THEN
         PP1 = P + 1
         DO 2 I = PP1, K
            INUL(I) = .TRUE.
   2     CONTINUE
      END IF
C
C     Step 1: Bidiagonalization phase
C             -----------------------
C     1.a.: if M >= 5*N/3, transform A into upper triangular form R by
C           Householder transformations.
C
      MA = M
      IF (QR) THEN
         CALL DQRDC(A, LDA, M, N, Q, DUMAR1, WRK, 0)
C
C        If (WANTU), store the Householder transformations performed on
C        the columns of A in the last N*(N+1)/2 extra storage locations
C        of work array WRK.
C
         K = M + N + 1
         NM1 = N - 1
         DO 4 J = 1, NM1
            J1 = J + 1
            IF (WANTU) THEN
               NJ = N - J
               WRK(K) = Q(J)
               K = K + 1
               CALL DCOPY(NJ, A(J1,J), 1, WRK(K), 1)
               K = K + NJ
            END IF
            DO 3 I = J1, N
               A(I,J) = 0.0D0
   3        CONTINUE
   4     CONTINUE
         IF (WANTU) WRK(K) = Q(N)
         MA = N
      END IF
C
C     1.b.: Transform A (or R) into bidiagonal form Q using Householder
C           transformations.
C
      CALL BIDIAG(A, LDA, MA, N, Q, Q(P+1), WRK)
C
C     1.c.: Initialize U (if WANTU) and V (if WANTV) with the identity
C           matrix.
C
      PP1 = MIN(M+1,N)
      J = P
      IF (WANTU) THEN
         IF (ALL) J = M
         CALL INIT(U, LDU, M, J)
      END IF
      J = PP1
      IF (WANTV) THEN
         IF (ALL) J = N
         CALL INIT(V, LDV, N, J)
      END IF
C
C     1.d.: If M < N, bring the bidiagonal to M by M by cancelling its
C           last superdiagonal element using Givens rotations.
C
      IF (M .LT. N) CALL CANCEL(U, LDU, V, LDV, Q, Q(PP1), P, PP1, PP1,
     *               PP1, TOL2, WANTU, WANTV)
C
C     Step 2: Partial diagonalization phase.
C             -----------------------------
C             Diagonalize the bidiagonal Q partially until convergence
C             to  the desired left and/or right singular subspace.
C
      CALL QRQL(U, LDU, V, LDV, P, PP1, RANK, THETA, Q, Q(P+1), INUL,
     *          TOL1, TOL2, WANTU, WANTV, IERR, IWARN)
      IF (IERR.NE.0) RETURN
C
C     Step 3: Back transformation phase.
C             -------------------------
C     3.a.: Apply the Householder transformations of the bidiagonaliza-
C           tion onto the base vectors associated with the desired sub-
C           bidiagonals.
C
      CALL RESTOR(A, LDA, MA, N, U, LDU, V, LDV, INUL, TOL2, WANTU,
     *            WANTV)
C
C     3.b.: If M >= 5*N/3 and WANTU, apply the Householder trans-
C           formations of the triangularization of A onto the desired
C           base vectors.
C
      IF (QR .AND. WANTU) THEN
         IF (ALL) THEN
            J = P + 1
            DO 5 I = J, M
               INUL(I) = .TRUE.
   5        CONTINUE
         END IF
         K = M + N + 1
         I = N
         DO 6 J = 1, N
            CALL DCOPY(I, WRK(K), 1, A(J,J), 1)
            K = K + I
            I = I - 1
   6     CONTINUE
         CALL RESTOR(A, LDA, M, N, U, LDU, WRK, 1, INUL, TOL2, .TRUE.,
     *               .FALSE.)
      END IF
      RETURN
C *** Last line of PSVD ***********************************************
      END
