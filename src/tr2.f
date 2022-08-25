      SUBROUTINE TR2(A, LDA, U, S, I1, I2, J1, J2)
C
C     PURPOSE:
C     The subroutine TR2 performs the Householder transformation
C     H = I - S x UU' on the columns J1+1 to J1+J2  of matrix A, this
C     from rows I1 to I2.
C
C     PARAMETERS:
C
C     A - DOUBLE PRECISION array of DIMENSION (LDA,J1+J2).
C         Contains the submatrix of A onto which the Householder
C         transformation H is applied.
C     LDA - INTEGER.
C         The leading dimension of the array A (LDA >= I2).
C     U - DOUBLE PRECISION array of DIMENSION (J2).
C         Contains the transformation vector of the transformation
C         matrix H.
C     S - DOUBLE PRECISION.
C         Contains the scalar S of the transformation matrix H.
C     I1 - INTEGER.
C         Contains the first row index of A (see purpose).
C     I2 - INTEGER.
C         Contains the last row index of A (see purpose, I2 >= I1).
C     J1 - INTEGER.
C         Contains the first column index of A (see purpose).
C     J2 - INTEGER.
C         Contains the last column index of A (see purpose).
C
C     CONTRIBUTOR: P. Van Dooren (Philips Res. Laboratory, Brussels).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER I1, I2, J1, J2, LDA
      DOUBLE PRECISION S
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), U(J2)
C     .. Local Scalars ..
      INTEGER I, J
      DOUBLE PRECISION INPROD, Y
C     .. Executable Statements ..
C
      DO I = I1, I2
         INPROD = 0.0D0
         DO J = 1, J2
            INPROD = INPROD + U(J) * A(I,J1+J)
         END DO
         Y = INPROD * S
         DO J = 1, J2
            A(I,J1+J) = A(I,J1+J) - U(J) * Y
         END DO
      END DO
      RETURN
C *** Last line of TR2 ************************************************
      END
