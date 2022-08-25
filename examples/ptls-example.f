      PARAMETER (IN = 5, IOUT = 6)
      PARAMETER (LDC = 30, LDX = 10)
      INTEGER M, N, L, NL, RANK, IERR, IWARN, IWRK(11)
      LOGICAL INUL(11), LWRK(11)
      DOUBLE PRECISION C(LDC,11), X(LDX,10), Q(22),
     *       WRK(44),
     *       THETA, TOL1, TOL2
      INTEGER I, J, K
      INTRINSIC MIN
      READ(IN,12) M, N, L, RANK, THETA
      WRITE(IOUT,13) M, N, L, RANK, THETA
      NL = N + L
      WRITE(IOUT,9)
      DO 1 I = 1, M
         READ(IN,5) (C(I,J), J=1,NL)
  1   CONTINUE
      WRITE(IOUT,6)
      DO 2 I = 1, M
         WRITE(IOUT,7) (C(I,J), J=1,NL)
  2   CONTINUE
      TOL1 = 1.0D-08
      TOL2 = 1.0D-10
      CALL PTLS(C, LDC, M, N, L, RANK, THETA, X, LDX, Q, INUL, WRK,
     *          IWRK, LWRK, TOL1, TOL2, IERR, IWARN)
      WRITE(IOUT,8) IERR, IWARN
      WRITE(IOUT,18) RANK, THETA
      K = MIN(M,NL)
      WRITE(IOUT,10) (Q(I), I = 1, K)
      WRITE(IOUT,11) (Q(K+I), I = 2, K)
      WRITE(IOUT,14)
      K = 0
      DO 3 I = 1, NL
         IF (INUL(I)) THEN
            K = K + 1
            WRITE(IOUT,15) K, I
            WRITE(IOUT,7) (C(J,I), J = 1, NL)
         END IF
  3   CONTINUE
      WRITE(IOUT,16)
      DO 4 I = 1, L
         WRITE(IOUT,17) I, (X(J,I), J = 1, N)
  4   CONTINUE
      STOP
  5   FORMAT(4D15.0)
  6   FORMAT(/' ',5X,'C(*,I)',8X,'C(*,I+1)',7X,'C(*,I+2)',7X,
     *       'C(*,I+3)')
  7   FORMAT(' ',4D15.6)
  8   FORMAT(/' IERR = ',I3,'   IWARN = ',I3)
  9   FORMAT(/' PARTIAL TOTAL LEAST SQUARES TEST PROGRAM'/1X,40('-')/)
 10   FORMAT(/' DIAGONAL OF THE PARTIALLY DIAGONALIZED BIDIAGONAL='
     *       /(1X,4D15.6))
 11   FORMAT(/' SUPERDIAGONAL OF THE PARTIALLY DIAGONALIZED ',
     *       'BIDIAGONAL='/(1X,4D15.6))
 12   FORMAT(4I3,D12.5)
 13   FORMAT(/' M =',I3,'  N =',I3,'  L =',I3, '  RANK =',I3,
     *        '  THETA =',D12.5)
 14   FORMAT(/' BASIS OF THE COMPUTED RIGHT SINGULAR SUBSPACE :'/1X,
     *       45('-'))
 15   FORMAT(' THE ',I2,'-TH BASE VECTOR V(*,',I2,') =')
 16   FORMAT(/' TLS SOLUTION :'/1X,12('*'))
 17   FORMAT(/' X(*,',I2,') = '/(' ',4D15.6))
 18   FORMAT(/' COMPUTED RANK =',I3,4X,' COMPUTED BOUND THETA = ',
     *       D12.5)
      END
      