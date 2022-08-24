      SUBROUTINE INIT(X, LDX, M, N)
C
C     PURPOSE:
C
C     The subroutine INIT initializes an M by N matrix X with the M by N
C     identity matrix, characterized by unit diagonal entries and zero
C     off-diagonal elements.
C
C     ARGUMENT LIST:
C
C     X - DOUBLE PRECISION array of DIMENSION (LDX,N)
C         On return, X contains the M by N identity matrix.
C     LDX - INTEGER
C         LDX is the leading dimension of the array X (LDX >= M).
C     M - INTEGER
C         M is the number of rows of the matrix X.
C     N - INTEGER
C         N is the number of columns of the matrix X.
C
C     CONTRIBUTOR: S. Van Huffel, (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION X(LDX,*)
C     .. Local Scalars ..
      INTEGER I, J
C     .. Executable Statements ..
C
      DO 20 J = 1, N
         DO 10 I = 1, M
            X(I,J) = 0.0D0
   10    CONTINUE
         X(J,J) = 1.0D0
   20 CONTINUE
      RETURN
C *** Last line of INIT ***********************************************
      END
