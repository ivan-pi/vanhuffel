      SUBROUTINE CANCEL(U, LDU, V, LDV, Q, E, M, N, I, K, TOL,
     *                  WANTU, WANTV)
C
C     PURPOSE:
C
C     Either, subroutine CANCEL separates a zero singular value of a
C     subbidiagonal matrix of order k, k <= p, of the bidiagonal
C
C               !Q(1) E(2)  0    ...   0  !
C               ! 0   Q(2) E(3)        .  !
C           J = ! .                    .  !
C               ! .                   E(p)!
C               ! 0   ...             Q(p)!
C
C     with p = min(M,N), by annihilating one or two superdiagonal
C     elements E(i) and/or E(i+1).
C     Or, CANCEL annihilates the element E(M+1) of the bidiagonal matrix
C
C               !Q(1) E(2)  0    ...   0     0   !
C               ! 0   Q(2) E(3)        .     .   !
C           J = ! .                    .     .   !
C               ! .                   E(M)   .   !
C               ! 0   ...             Q(M) E(M+1)!
C
C     ARGUMENT LIST:
C
C     U - DOUBLE PRECISION array of DIMENSION (LDU,p).
C         On entry, U contains the M by p (p = min(M,N)) left trans-
C         formation matrix.
C         On return, the Givens rotations S on the left, annihilating
C         E(i+1), have been postmultiplied into U.
C         NOTE: U is not referenced if WANTU = .FALSE. .
C     LDU - INTEGER.
C         LDU is the leading dimension of the array U (LDU >= M).
C     V - DOUBLE PRECISION array of DIMENSION (LDV,s).
C         On entry, V contains the N by s (s = min(M+1,N)) right trans-
C         formation matrix.
C         On return, the Givens rotations T on the right, annihilating
C         E(i), have been postmultiplied into V.
C         NOTE: V is not referenced if WANTV = .FALSE. .
C     LDV - INTEGER.
C         LDV is the leading dimension of the array V (LDV >= N).
C     Q - DOUBLE PRECISION array of DIMENSION (p).
C         On entry, Q contains the diagonal entries of the bidiagonal J.
C         p = min(M,N).
C         On return, Q contains the diagonal elements of the transformed
C         bidiagonal S' J T.
C     E - DOUBLE PRECISION array of DIMENSION (s).
C         On entry, E(i), i=2,...,s, contain the superdiagonal entries
C         of the bidiagonal J. s = min(M+1,N), E(1) = 0.0D0.
C         On return, E contains the superdiagonal elements of the trans-
C         formed bidiagonal S' J T.
C     M - INTEGER.
C         M is the number of rows of the matrix U.
C     N - INTEGER.
C         N is the number of rows of the matrix V.
C     I - INTEGER.
C         Either, I is the index of the negligible diagonal entry Q(I)
C         of the bidiagonal J, i,e. abs(Q(I)) <= TOL, I <= p.
C         Or, I = M + 1 if E(M+1) is to be annihilated.
C     K - INTEGER.
C         Either, K is the index of the last diagonal entry of the con-
C         sidered subbidiagonal of J, i.e. abs(E(K+1)) <= TOL, K <= p.
C         Or, K = M + 1 if E(M+1) is to be annihilated.
C     TOL - DOUBLE PRECISION.
C         Specifies that matrix elements Q(i), which are <= TOL in
C         absolute value, are considered to be zero.
C     WANTU - LOGICAL.
C         Logical indicating the need for postmultiplying the Givens
C         rotations S on the left into U.
C     WANTV - LOGICAL.
C         Logical indicating the need for postmultiplying the Givens
C         rotations T on the right into V.
C
C     EXTERNAL SUBROUTINES and FUNCTIONS:
C
C     DROT from BLAS.
C
C     METHOD DESCRIPTION:
C
C     Let the considered subbidiagonal be
C
C               !Q(1) E(2)  0                    ...   0  !
C               ! 0   Q(2) E(3)                  ...      !
C               ! .                              ...      !
C               !             Q(i-1) E(i)              .  !
C          Jk = !                    Q(i) E(i+1)       .  !
C               !                         Q(i+1) .        !
C               ! .                              ..       !
C               ! .                                   E(k)!
C               ! 0    ...                       ...  Q(k)!
C
C     A zero singular value of Jk manifests itself by a zero diagonal
C     entry Q(i) or in practice, a negligible value of Q(i).
C     We call Q(i) negligible if abs(Q(i)) <= TOL.
C     When such a negligible diagonal element Q(i) in Jk is present,
C     the subbidiagonal Jk is splitted by the routine CANCEL into 2 or
C     3 unreduced subbidiagonals by annihilating E(i+1) (if i<k) using
C     Givens rotations S on the left and by annihilating E(i) (if i>1)
C     using Givens rotations T on the right until Jk is reduced to the
C     form :
C
C               !Q(1) E(2)  0                ...   0  !
C               ! 0         .                ...      !
C               ! .                          ...      !
C               !         Q(i-1) 0                 .  !
C     S' Jk T = !                0   0             .  !
C               !                   Q(i+1)   .        !
C               ! .                          ..       !
C               ! .                               E(k)!
C               ! 0    ...                   ...  Q(k)!
C
C     For more details, see [1, pp.11.12-11.14].
C     The case of the annihilation of E(M+1) can be treated by the same
C     process. This may be seen by augmenting the matrix J with an extra
C     row of zeros, i.e. by introducing Q(M+1) = 0.
C
C     REFERENCES:
C
C     [1] J.J. Dongarra, J.R. Bunch, C.B. Moler and G.W. Stewart,
C         LINPACK User's Guide. SIAM, Philadelphia (1979).
C
C     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).
C
C     REVISIONS: 1988, February 15.
C
C     .. Scalar Arguments ..
      INTEGER LDU, LDV, M, N, I, K
      DOUBLE PRECISION TOL
      LOGICAL WANTU, WANTV
C     .. Array Arguments ..
      DOUBLE PRECISION U(LDU,*), V(LDV,*), Q(*), E(*)
C     .. External Subroutines/Functions ..
      EXTERNAL DROT
C     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
C     .. Local Scalars ..
      INTEGER I1, L, L1
      DOUBLE PRECISION C, S, F, G, H
C     .. Executable Statements ..
C
      IF (I .LE. M) Q(I) = 0.0D0
C
C     Annihilate E(I+1) (if I < K).
C
      IF (I .LT. K) THEN
         C = 0.0D0
         S = 1.0D0
         I1 = I + 1
         DO L = I1, K
            G = E(L)
            F = S * G
            E(L) = C * G
            IF (ABS(F) .LE. TOL) GO TO 2
            G = Q(L)
            H = SQRT(F**2 + G**2)
            Q(L) = H
            C = G/H
            S = -F/H
            IF (WANTU) CALL DROT(M, U(1,I), 1, U(1,L), 1, C, S)
         END DO
      END IF
C
C     Annihilate E(I) (if I > 1).
C
    2 IF (I .GT. 1) THEN
         I1 = I - 1
         F = E(I)
         E(I) = 0.0D0
         DO L1 = 1, I1
            IF (ABS(F) .LE. TOL) RETURN
            L = I - L1
            G = Q(L)
            IF (ABS(G) .LE. TOL) THEN
               G = 0.0D0
               H = ABS(F)
            ELSE
               H = SQRT(F**2 + G**2)
            END IF
            Q(L) = H
            C = G/H
            S = -F/H
            G = E(L)
            F = S * G
            E(L) = C * G
            IF (WANTV) CALL DROT(N, V(1,I), 1, V(1,L), 1, C, S)
         END DO
         E(1) = 0.0D0
      END IF
      RETURN
C *** Last line of CANCEL *********************************************
      END
