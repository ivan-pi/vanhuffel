      SUBROUTINE BIDIAG(X, LDX, N, P, Q, E, WRK)
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDX, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION X(LDX,*), Q(*), E(*), WRK(N+P)
C     .. External Subroutines/Functions ..
      DOUBLE PRECISION DNRM2, DDOT
      EXTERNAL DSCAL, DNRM2, DDOT, DAXPY
C     .. Intrinsic Functions ..
      INTRINSIC MAX, MIN, SIGN
C     .. Local Scalars ..
      INTEGER NCT, NRT, LU, L, LP1, J, I, NPP, PP1, M
      DOUBLE PRECISION T
C     .. Executable Statements ..
C
C     Reduce X to bidiagonal form, storing the diagonal elements in Q
C     and the superdiagonal elements in E.
C
      NPP = N + P
      PP1 = P + 1
      NCT = MIN(N-1,P)
      NRT = MAX(0,MIN(P-2,N))
      LU = MAX(NCT,NRT)
      DO 7 L = 1, LU
         LP1 = L + 1
         IF (L .LE. NCT) THEN
C
C           Compute the transformation for the L-th column and place the
C           L-th diagonal in Q(L).
C
            Q(L) = DNRM2(N-L+1, X(L,L), 1)
            IF (Q(L) .NE. 0.0D0) THEN
               IF (X(L,L) .NE. 0.0D0) Q(L) = SIGN(Q(L), X(L,L))
               CALL DSCAL(N-L+1, 1.0D0/Q(L), X(L,L), 1)
               X(L,L) = 1.0D0 + X(L,L)
            END IF
            Q(L) = -Q(L)
         END IF
C
         DO 2 J = LP1, P
            IF (L .GT. NCT) GO TO 1
            IF (Q(L) .EQ. 0.0D0) GO TO 1
C
C              Apply the transformation.
C
               T = -DDOT(N-L+1, X(L,L), 1, X(L,J), 1)/X(L,L)
               CALL DAXPY(N-L+1, T, X(L,L), 1, X(L,J), 1)
    1       CONTINUE
C
C           Place the L-th row of X into WRK for the subsequent
C           calculation of the row transformation.
C
            WRK(J) = X(L,J)
    2    CONTINUE
         IF (L .LE. NRT) THEN
C
C           Compute the L-th row transformation and place the L-th
C           superdiagonal in E(L).
C
            WRK(L) = DNRM2(P-L, WRK(LP1),1)
            IF (WRK(L) .NE. 0.0D0) THEN
               IF (WRK(LP1) .NE. 0.0D0) WRK(L) = SIGN(WRK(L), WRK(LP1))
               CALL DSCAL(P-L, 1.0D0/WRK(L), WRK(LP1),1)
               WRK(LP1) = 1.0D0 + WRK(LP1)
            END IF
            WRK(L) = -WRK(L)
            E(LP1) = WRK(L)
            IF (LP1 .LE. N .AND. WRK(L) .NE. 0.0D0) THEN
C
C              Apply the transformation.
C
               DO 3 I = PP1, NPP
                  WRK(I) = 0.0D0
    3          CONTINUE
               DO 4 J = LP1, P
                  CALL DAXPY(N-L, WRK(J), X(LP1,J), 1, WRK(PP1),1)
    4          CONTINUE
               DO 5 J = LP1, P
                  CALL DAXPY(N-L,-WRK(J)/WRK(LP1),WRK(PP1),1,X(LP1,J),1)
    5          CONTINUE
            END IF
            DO 6 I = LP1, P
               X(L,I) = WRK(I)
    6       CONTINUE
         END IF
    7 CONTINUE
C
C     Set up the final bidiagonal matrix elements, if necessary.
C
      E(1) = 0.0D0
      M = MIN(P-1,N)
      IF (NCT .LT. P) Q(N) = X(N,N)
      IF (NRT .LT. M) E(P) = X(P-1,P)
      RETURN
C *** Last line of BIDIAG *********************************************
      END
