These codes were written in connection with the Ph.D. thesis
Sabine Van Huffel, 1987, Katholieke Universiteit Leuven.  A revised and
expanded version of the thesis is in SIAM Frontiers book series.
This particular version was received on 21 March 1988.

The primary routines to be called by the user are dtls, psvd, and ptls.
A complete list is:
bidiag	auxiliary routines called by dtls, psvd, ptls
dtls	solves Ax=b, where A as well as b may be inaccurate
estim	singular values of a bidiagonal by bisection
nsingv	number of singular values left of a bound, of a bidiagonal matrix
psvd	basis for the left/right space corresponding to small singular values
ptls	like dtls, but using the faster Partial SVD
qrql	diagonalizes a bidiagonal matrix

For each name x in this list, there is a file x.f containing a
Fortran subroutine and a file x-doc containing documentaion.
Some of the documentation files also contain sample driver programs.

By default, when asking for a top-level routine like dtls you will also
get any necessary auxiliary routines. In keeping with netlib policy,
routines from linpack and the blas have been omitted here since they
are available directly.  To include them in your return mail, say
    send dtls from vanhuffel linpack blas.
    send dtls-doc from vanhuffel

The routines required from BLAS include: DSCAL, DNRM2, DDOT, DAXPY, DROT, DCOPY. 

The routines required from LINPACK are DQRDC, DSVDC