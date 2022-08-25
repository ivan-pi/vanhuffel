      SUBROUTINE QLSTEP(U, LDU, V, LDV, Q, E, M, N, I, K, SHIFT,
     *                  WANTU, WANTV)
C
C     PURPOSE:
C
C     The subroutine QLSTEP performs one QL iteration step onto the
C     unreduced subbidiagonal Jk:
C
C              !Q(i) E(i+1)  0  ...    0  !
C              ! 0   Q(i+1) E(i+2)     .  !
C         Jk = ! .                     .  !
C              ! .                        !
C              ! .                    E(k)!
C              ! 0   ...              Q(k)!
C
C     with k <= p and i >= 1, p = min(M,N), of the bidiagonal J:
C
C              !Q(1) E(2)  0    ...   0  !
C              ! 0   Q(2) E(3)        .  !
C          J = ! .                    .  !
C              ! .                   E(p)!
C              ! 0   ...             Q(p)!
C
C     Hereby, Jk is transformed to  S'Jk T with S and T products of
C     Givens rotations. These Givens rotations S (resp.,T) will be post-
C     multiplied into U (resp.,V), if WANTU (resp.,WANTV) = .TRUE.
C
C     ARGUMENT LIST:
C
C     U - DOUBLE PRECISION array of DIMENSION (LDU,min(M,N)).
C         On entry, U may contain the M by p (p=min(M,N)) left transfor-
C         mation matrix.
C         On return, if WANTU = .TRUE., the Givens rotations S on the
C         left have been postmultiplied into U.
C         NOTE: U is not referenced if WANTU = .FALSE.
C     LDU - INTEGER.
C         LDU is the leading dimension of the array U (LDU >= M).
C     V - DOUBLE PRECISION array of DIMENSION (LDV,min(M,N)).
C         On entry, V may contain the N by p (p=min(M,N)) right trans-
C         formation matrix.
C         On return, if WANTV = .TRUE., the Givens rotations T on the
C         right have been postmultiplied into V.
C         NOTE: V is not referenced if WANTV is .false.
C     LDV - INTEGER.
C         LDV is the leading dimension of the array V (LDV >= N).
C     Q - DOUBLE PRECISION array of DIMENSION (min(M,N)).
C         On entry, Q contains the diagonal entries of the bidiagonal J.
C         On return, Q contains the diagonal entries of the transformed
C         matrix S' J T.
C     E - DOUBLE PRECISION array of DIMENSION (min(M,N)).
C         On entry, E contains the superdiagonal entries of J.
C         On return, E contains the superdiagonal entries of the trans-
C         formed matrix S' J T. E(k+1) = 0 if k < min(M,N).
C     M - INTEGER.
C         M is the number of rows of the matrix U.
C     N - INTEGER.
C         N is the number of rows of the matrix V.
C     I - INTEGER.
C         I is the index of the first diagonal entry of the considered
C         unreduced subbidiagonal Jk of J.
C     K - INTEGER.
C         K is the index of the last diagonal entry of the considered
C         unreduced subbidiagonal Jk of J.
C     SHIFT - DOUBLE PRECISION.
C         Value of the shift used in the QL iteration step.
C     WANTU - LOGICAL.
C         WANTU = .TRUE. if the Givens rotations S must be postmulti-
C         plied on the left into U, else .FALSE.
C     WANTV - LOGICAL.
C         WANTV = .TRUE. if the Givens rotations T must be postmulti-
C         plied  on the left into V, else .FALSE.
C
C     EXTERNAL SUBROUTINES AND FUNCTIONS:
C
C     DROT from BLAS.
C
C     METHOD DESCRIPTION:
C
C     QL iterations diagonalize the bidiagonal by zeroing the super-
C     diagonal elements of Jk from top to bottom.
C     The routine QLSTEP overwrites Jk with the bidiagonal matrix
C     S' Jk T where S and T are Givens rotations.
C     T is essentially the orthogonal matrix that would be obtained by
C     applying one implicit symmetric shift QL step onto the matrix
C     Jk'Jk. This step factors the matrix (Jk'Jk - shift*I) into a
C     product of an orthogonal matrix T and a lower triangular matrix.
C     See [1,Sec.8.2-8.3] and [2] for more details.
C
C     REFERENCES:
C     [1] G.H. Golub and C.F. Van Loan, Matrix Computations. The Johns
C         Hopkins University Press, Baltimore,Maryland (1983).
C     [2] H. Bowdler, R.S. Martin and J.H. Wilkinson, The QR and QL
C         algorithms for symmetric matrices. Numer. Math., 11 (1968),
C         293-306.
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDU, LDV, M, N, I, K
      DOUBLE PRECISION SHIFT
      LOGICAL WANTU, WANTV
C     .. Array Arguments ..
      DOUBLE PRECISION U(LDU,*), V(LDV,*), Q(*), E(*)
C     .. External Subroutines/Functions ..
      EXTERNAL DROT
C     .. Intrinsic Functions ..
      INTRINSIC MIN, SQRT
C     .. Local Scalars ..
      INTEGER IK, J, LL, JJ
      DOUBLE PRECISION F, G, H, C, S, X, Y, Z
C     .. Executable Statements ..
C
      X = Q(K)
      G = SHIFT
      F = (X - G) * (X + G)/X
      C = 1.0D0
      S = 1.0D0
      IK = K - I
      DO JJ = 1, IK
         J = K - JJ
         LL= J + 1
         G = E(LL)
         Y = Q(J)
         H = S * G
         G = C * G
         Z = SQRT(F**2 + H**2)
         IF (JJ. GT. 1) E(LL+1) = Z
         C = F/Z
         S = H/Z
         F = X * C + G * S
         G = -X * S + G * C
         H = Y * S
         Y = Y * C
         IF (WANTU) CALL DROT(M, U(1,LL), 1, U(1,J), 1, C, S)
         Z = SQRT(F**2 + H**2)
         Q(LL) = Z
         C = F/Z
         S = H/Z
         F = C * G + S * Y
         X = -S * G + C * Y
         IF (WANTV) CALL DROT(N, V(1,LL), 1, V(1,J), 1, C, S)
      END DO
      E(I+1) = F
      Q(I) = X
      IF (K .LT. MIN(M,N)) E(K+1) = 0.0D0
      RETURN
C Last line of QLSTEP *************************************************
      END
