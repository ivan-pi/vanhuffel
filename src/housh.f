      SUBROUTINE HOUSH(DUMMY, K, J, TOL, ZERO, S)
C
C     PURPOSE:
C
C     The subroutine HOUSH computes a Householder transformation H,
C     H = I - S x UU', that 'mirrors' a vector DUMMY(1,...,K) to the
C     J-th unit vector.
C
C     ARGUMENT LIST :
C
C     DUMMY - DOUBLE PRECISION array of DIMENSION (K).
C         A row or column vector of a matrix that has to be mirrored
C         to the corresponding unit vector EJ = (0,...,1,0,...,0).
C         On return, DUMMY contains the U-vector of the transformation
C         matrix H = I - S x UU'.
C     K - INTEGER.
C         The dimension of DUMMY.
C     J - INTEGER.
C         The transformation preserves the J-th element of DUMMY to
C         become zero. All the other elements are transformed to zero.
C     TOL - DOUBLE PRECISION.
C         If on entry, norm(DUMMY) < TOL, ZERO is put equal to .TRUE.
C     ZERO - LOGICAL.
C         See the description of TOL.
C     S - DOUBLE PRECISION.
C         On return, S contains the scalar S of the transformation
C         matrix H.
C
C     REFERENCES:
C
C     [1] A. Emami-Naeine and P. Van Dooren, Computation of Zeros of
C         Linear Multivariable Systems.
C         Automatica, 18,No.4 (1982), 415-430.
C
C     CONTRIBUTOR: P. Van Dooren (Philips Res. Laboratory, Brussels).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER K, J
      DOUBLE PRECISION TOL, S
      LOGICAL ZERO
C     .. Array Arguments ..
      DOUBLE PRECISION DUMMY(K)
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION ALFA, DUM1
C     .. Executable Statements ..
C
      ZERO = .TRUE.
      S = 0.0D0
      DO 10 I = 1, K
         S = S + DUMMY(I)**2
   10 CONTINUE
      ALFA = SQRT(S)
      IF (ALFA .LE. TOL) RETURN
      ZERO = .FALSE.
      DUM1 = DUMMY(J)
      IF (DUM1 .GT. 0.0D0) ALFA = -ALFA
      DUMMY(J) = DUM1 - ALFA
      S = 1.0D0/(S - ALFA * DUM1)
      RETURN
C *** Last line of HOUSH **********************************************
      END
