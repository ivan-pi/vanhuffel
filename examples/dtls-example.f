
      PARAMETER (IN=5, IOUT=6)
      PARAMETER (LC=14, LM=10, LN=10, LL=4)
      INTEGER LDC, M, N, L, LDX, IERR, NL, RANK, IWARN
      DOUBLE PRECISION C(LC,LN+LL), X(LN,LL), WRK(LM+LN+LL), S(LN+LL),
     *                 SDEV, TOL2
      INTEGER I, J, K
      CHARACTER MODE*1
      LDC = LC
      LDX = LN
      READ(IN, 11) M, N, L
      WRITE(IOUT, 12) M, N, L
      NL = N + L
      WRITE(IOUT, 8)
      DO 1 I = 1, M
        READ(IN, 4) (C(I,J), J=1,NL)
    1 CONTINUE
      WRITE(IOUT, 5)
      DO 2 I = 1, M
         WRITE (IOUT, 6) (C(I,J), J=1,NL)
    2 CONTINUE
      MODE = 'X'
      SDEV = 1.0D-04
      TOL2  = 1.0D-04
      CALL DTLS(C, LDC, M, N, L, S, X, LDX, WRK, RANK, SDEV, TOL2,
     *          MODE, IERR, IWARN)
      WRITE(IOUT, 7) IERR, IWARN
      WRITE(IOUT,14) RANK
      WRITE(IOUT, 9) (S(I), I=1,NL)
      WRITE(IOUT,13)
      DO 3 I = 1, L
         WRITE (IOUT, 10) I, (X(J,I), J=1,N)
    3 CONTINUE
      STOP
    4 FORMAT(4D15.0)
    5 FORMAT(/' ',5X,'C(*,I)',8X,'C(*,I+1)',7X,'C(*,I+2)',7X,
     *        'C(*,I+3)')
    6 FORMAT(' ',4D15.6)
    7 FORMAT(/' IERR = ',I3,'   IWARN = ',I3)
    8 FORMAT(/' TOTAL LEAST SQUARES TEST PROGRAM'/1X,32('-')/)
    9 FORMAT(/' SINGULAR VALUES S(I) OF C = '/(' ',4D15.6))
   10 FORMAT(/' X(*,'I2,') = '/(' ',4D15.6))
   11 FORMAT(3I3)
   12 FORMAT(/' M = ',I3,3X,'N = ',I3,3X,'L = ',I3)
   13 FORMAT(/' TLS SOLUTION :'/1X,12('*'))
   14 FORMAT(/' COMPUTED RANK =',I3)
      END
       