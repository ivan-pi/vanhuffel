      PARAMETER (IN = 5, IOUT = 6)
      PARAMETER (LDA = 30, LDU = 30, LDV = 11, LN = 11, LM = 30,
     *           LQ = 22, LW = 107, LINUL = 30)
      INTEGER M, N, RANK, MODE, IERR, IWARN
      DOUBLE PRECISION A(LDA,LN), U(LDU,LM), V(LDV,LN), Q(LQ),
     *       WRK(LW),
     *       TOL1, TOL2, THETA
      INTRINSIC MIN
      LOGICAL INUL(LINUL)
      INTEGER I, J, K
      READ(IN,12) M, N, RANK, THETA
      WRITE(IOUT,13) M, N, RANK, THETA
      WRITE(IOUT,9)
      DO 1 I = 1, M
         READ(IN,5) (A(I,J), J = 1, N)
  1   CONTINUE
      WRITE(IOUT,6)
      DO 2 I = 1, M
         WRITE(IOUT,7) (A(I,J), J = 1, N)
  2   CONTINUE
      MODE = 11
      TOL1 = 1.0D-08
      TOL2 = 1.0D-10
      CALL PSVD(A, LDA, M, N, RANK, THETA, U, LDU, V, LDV, Q, INUL,
     *          WRK, TOL1, TOL2, MODE, IERR, IWARN)
      WRITE(IOUT,8) IERR, IWARN
      WRITE(IOUT,18) RANK, THETA
      K = MIN(M,N)
      WRITE(IOUT,10) (Q(I), I = 1, K)
      WRITE(IOUT,11) (Q(K+I), I = 2, K)
      WRITE(IOUT,14)
      K = 0
      DO 3 I = 1, M
         IF (INUL(I)) THEN
            K = K + 1
            WRITE(IOUT,15) K, I
            WRITE(IOUT,7) (U(J,I), J = 1, M)
         END IF
  3   CONTINUE
      WRITE(IOUT,16)
      K = 0
      DO 4 I = 1, N
         IF (INUL(I)) THEN
            K = K + 1
            WRITE(IOUT,17) K, I
            WRITE(IOUT,7) (V(J,I), J = 1, N)
         END IF
  4   CONTINUE
      STOP
  5   FORMAT(4D15.0)
  6   FORMAT(/' ',5X,'A(*,I)',8X,'A(*,I+1)',7X,'A(*,I+2)',7X,
     *       'A(*,I+3)')
  7   FORMAT(' ',4D15.6)
  8   FORMAT(/' IERR = ',I3,'   IWARN = ',I3)
  9   FORMAT(/' PARTIAL SVD (PSVD) TEST PROGRAM'/' ',32('-')/)
 10   FORMAT(/' DIAGONAL OF THE PARTIALLY DIAGONALIZED BIDIAGONAL='
     *       /(1X,4D15.6))
 11   FORMAT(/' SUPERDIAGONAL OF THE PARTIALLY DIAGONALIZED ',
     *       'BIDIAGONAL=',/(1X,4D15.6))
 12   FORMAT(3I3,D12.5)
 13   FORMAT(/' M =',I3,'  N =',I3,'  RANK =',I3,'  THETA =',D12.5)
 14   FORMAT(/' BASIS OF THE COMPUTED LEFT SINGULAR SUBSPACE :'/1X,
     *       44('-'))
 15   FORMAT(' THE ',I2,'-TH BASE VECTOR U(*,',I2,') =')
 16   FORMAT(/' BASIS OF THE COMPUTED RIGHT SINGULAR SUBSPACE :'/1X,
     *       45('-'))
 17   FORMAT(' THE ',I2,'-TH BASE VECTOR V(*,',I2,') =')
 18   FORMAT(/' COMPUTED RANK =',I3,4X,' COMPUTED BOUND THETA = ',
     *       D12.5)
      END