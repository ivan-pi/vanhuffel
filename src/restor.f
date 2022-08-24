      SUBROUTINE RESTOR(X, LDX, M, N, U, LDU, V, LDV, INUL, TOL,
     *                  WANTU, WANTV)

C
C     PURPOSE:
C
C     The subroutine RESTOR applies the Householder transformations Pj
C     stored in factored form into the columns of the array X, to the
C     desired columns of the matrix U by premultiplication, and/or
C     the Householder transformations Qj stored in factored form into
C     the rows of the array X, to the desired columns of the matrix V by
C     premultiplication.
C
C     ARGUMENT LIST:
C
C     X - DOUBLE PRECISION array of DIMENSION (LDX,N).
C         The leading M x N part contains in the columns of its lower
C         triangle the Householder transformations Pj, and in the rows
C         of its upper triangle the Householder transformations Qj in
C         factored form
C     LDX - INTEGER.
C         LDX is the leading dimension of the array X (LDX >= M).
C     M - INTEGER.
C         M is the number of rows of the matrix X.
C     N - INTEGER.
C         N is the number of columns of the matrix X.
C     U - DOUBLE PRECISION array of DIMENSION (LDU,M).
C         On entry, U contains the M by M matrix U.
C         On return, the Householder transformations Pj have been
C         applied to each column i of U corresponding to a parameter
C         INUL(i) = .TRUE.
C         NOTE: U is not referenced if WANTU = .FALSE.
C     LDU - INTEGER.
C         LDU is the leading dimension of the array U (LDU >= M).
C     V - DOUBLE PRECISION of DIMENSION (LDV,N).
C         On entry, U contains the N by N matrix V.
C         On return, the Householder transformations Qj have been
C         applied to each column i of V corresponding to a parameter
C         INUL(i) = .TRUE.
C         NOTE: V is not referenced if WANTV = .FALSE.
C     LDV - INTEGER.
C         LDV is the leading dimension of the array V (LDV >= N).
C     INUL - LOGICAL array of DIMENSION (max(M,N)).
C         INUL(i) = .TRUE. if the i-th column of U and/or V is to be
C         transformed. (1 <= i <= max(M,N)).
C     TOL - DOUBLE PRECISION.
C         Specifies that a Householder transformation Pj (resp. Qj) is
C         considered to be an identity transformation when the corres-
C         ponding element X(j,j) for Pj (resp. X(j,j+1) for Qj) is
C         <= TOL in absolute value.
C     WANTU - LOGICAL.
C         WANTU = .TRUE. if columns in U have to be transformed by the
C         routine RESTOR, else .FALSE.
C     WANTV - LOGICAL.
C         WANTV = .TRUE. if columns in V have to be transformed by the
C         routine RESTOR, else .FALSE.
C
C     EXTERNAL SUBROUTINES AND FUNCTIONS:
C
C     DDOT, DAXPY from BLAS.
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDX, M, N, LDU, LDV
      DOUBLE PRECISION TOL
      LOGICAL WANTU, WANTV
C     .. Array Arguments ..
      DOUBLE PRECISION X(LDX,*), U(LDU,*), V(LDV,*)
      LOGICAL INUL(*)
C     .. External Subroutines/Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT, DAXPY
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MIN
C     .. Local Scalars ..
      INTEGER I, LL, L, NL1, ML1, L1, IM
      DOUBLE PRECISION T
C     .. Executable Statements ..
C
C     Apply the Householder transformations Pj onto the desired columns
C     of U.
C
      IM = MIN(M-1,N)
      IF (WANTU .AND. (IM .GT. 0)) THEN
         DO 20 I = 1, M
            IF (INUL(I)) THEN
               DO 10 LL = 1, IM
                  L = IM - LL + 1
                  IF (ABS(X(L,L)) .GT. TOL) THEN
                     ML1 = M - L + 1
                     T = -DDOT(ML1, X(L,L), 1, U(L,I), 1)/X(L,L)
                     CALL DAXPY(ML1, T, X(L,L), 1, U(L,I), 1)
                  END IF
   10          CONTINUE
            END IF
   20    CONTINUE
      END IF
C
C     Apply the Householder transformations Qj onto the desired columns
C     of V.
C
      IM = MIN(N-2,M)
      IF (WANTV .AND. (IM .GT. 0)) THEN
         DO 40 I = 1, N
            IF (INUL(I)) THEN
               DO 30 LL = 1, IM
                  L = IM - LL + 1
                  L1 = L + 1
                  IF (ABS(X(L,L1)) .GT. TOL) THEN
                     NL1 = N - L
                     T = -DDOT(NL1, X(L,L1), LDX, V(L1,I), 1)/X(L,L1)
                     CALL DAXPY(NL1, T, X(L,L1), LDX, V(L1,I), 1)
                  END IF
   30          CONTINUE
            END IF
   40    CONTINUE
      END IF
      RETURN
C *** Last line of RESTOR *********************************************
      END
