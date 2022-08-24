      INTEGER FUNCTION NSINGV(Q, E, K, THETA, TOL1, TOL2)
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER K
      DOUBLE PRECISION THETA, TOL1, TOL2
C     .. Array Arguments ..
      DOUBLE PRECISION Q(K), E(K)
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     .. Local Scalars ..
      INTEGER J, NUMEIG
      DOUBLE PRECISION R, T
C     .. Executable Statements ..
C
      IF (THETA .LT. 0.0D0) THEN
         NSINGV = 0
         RETURN
      END IF
      T = -THETA - TOL1
      NUMEIG = K
      IF (ABS(Q(1)) .LE. TOL2) THEN
         R = T
      ELSE
         R = T - Q(1) * (Q(1)/T)
         IF (R .GT. 0.0D0) NUMEIG = NUMEIG - 1
      END IF
C
      DO 1 J = 2, K
         IF (ABS(E(J)) .LE. TOL2) THEN
            R = T
         ELSE
            R = T - E(J) * (E(J)/R)
            IF (R .GT. 0.0D0) NUMEIG = NUMEIG - 1
         END IF
         IF (ABS(Q(J)) .LE. TOL2) THEN
            R = T
         ELSE
            R = T - Q(J) * (Q(J)/R)
            IF (R .GT. 0.0D0) NUMEIG = NUMEIG - 1
         END IF
    1 CONTINUE
C
      NSINGV = NUMEIG
      RETURN
C *** Last line of NSINGV *********************************************
      END
