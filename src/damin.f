      DOUBLE PRECISION FUNCTION DAMIN(X, NX, INCX)
C
C     PURPOSE:
C
C     The function DAMIN computes the absolute minimal value of NX
C     elements in an array.
C     The function returns the value zero if NX < 1.
C
C     ARGUMENT LIST:
C
C     X - DOUBLE PRECISION array of DIMENSION (NX x INCX).
C         X is the one-dimensional array of which the absolute minimal
C         value of the elements is to be computed.
C     NX - INTEGER.
C         NX is the number of elements in X to be examined.
C     INCX - INTEGER.
C         INCX is the increment to be taken in the array X, defining
C         the distance between two consecutive elements.
C         INCX = 1, if all elements are contiguous in memory.
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER NX, INCX
C     .. Array Arguments ..
      DOUBLE PRECISION X(NX*INCX)
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     .. Local Scalars ..
      INTEGER I, IX, NINCX
      DOUBLE PRECISION DX
C     .. Executable Statements ..
C
      DAMIN = 0.0D0
      IF (NX .LT. 1) RETURN
C
      DAMIN = ABS(X(1))
      IF (NX .EQ. 1) RETURN
C
      IX = 1 + INCX
      NINCX = NX*INCX
      DO I = IX, NINCX, INCX
         DX = ABS(X(I))
         IF (DX .LT. DAMIN) DAMIN = DX
      END DO
      RETURN
C *** Last line of DAMIN **********************************************
      END
