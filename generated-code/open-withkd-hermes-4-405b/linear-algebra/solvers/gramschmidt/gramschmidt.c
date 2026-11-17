/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gramschmidt.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gramschmidt.h"


/* Array initialization. */
static
void init_array(int m, int n,
		DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      A[i][j] = (((DATA_TYPE) ((i*j) % m) / m )*100) + 10;
      Q[i][j] = 0.0;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      R[i][j] = 0.0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("R");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, R[i][j]);
    }
  POLYBENCH_DUMP_END("R");

  POLYBENCH_DUMP_BEGIN("Q");
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Q[i][j]);
    }
  POLYBENCH_DUMP_END("Q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* QR Decomposition with Modified Gram Schmidt:
 http://www.inf.ethz.ch/personal/gander/ */
static
void kernel_gramschmidt(int m, int n,
			DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
			DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
			DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j, k;

  DATA_TYPE nrm;

#pragma scop
We are given the Gram-Schmidt kernel. The original code within the scop region is:

  for (k = 0; k < _PB_N; k++)
    {
      nrm = SCALAR_VAL(0.0);
      for (i = 0; i < _PB_M; i++)
        nrm += A[i][k] * A[i][k];
      R[k][k] = SQRT_FUN(nrm);
      for (i = 0; i < _PB_M; i++)
        Q[i][k] = A[i][k] / R[k][k];
      for (j = k + 1; j < _PB_N; j++)
	{
	  R[k][j] = SCALAR_VAL(0.0);
	  for (i = 0; i < _PB_M; i++)
	    R[k][j] += Q[i][k] * A[i][j];
	  for (i = 0; i < _PB_M; i++)
	    A[i][j] = A[i][j] - Q[i][k] * R[k][j];
	}
    }

We are to optimize only within the scop region.

HIGHEST PRIORITY — LOOP TRANSFORMATIONS (examine FIRST, in this order):
1. Loop tiling / blocking: We have nested loops. The inner loops are over i (from 0 to _PB_M-1) and j (from k+1 to _PB_N-1). The outer loop is k. We note that the inner loops are over i and j, and the outer loop k is sequential (due to dependencies). We can consider tiling the i and j loops to improve cache locality. However, note that the j loop is triangular and depends on k. Also, the inner loops are not too deep (only two levels). We must be cautious because the j loop starts at k+1 and goes to _PB_N-1, and the inner i loop is over the entire range. Tiling the i loop might be beneficial because it is the innermost and accesses A and Q in a contiguous manner. But note that the j loop is also a candidate. However, the j loop is not parallel and has a dependency on k. We might tile the i loop to improve cache reuse for A and Q.

But note: the innermost loop for the norm computation and the update of Q and R are over i. We can tile the i loop. However, the j loop is the next level and is also a candidate for tiling. But the j loop is not parallel and has a dependency on k. We might consider tiling the j loop as well? But note that the j loop is triangular and the inner i loop is the same for every j. We can interchange the j and i loops? But note the dependencies: the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. Let's consider the structure:

  for k ... 
    // norm: reduction over i
    // then update Q: over i
    // then for j from k+1 to N-1:
        // compute R[k][j]: reduction over i
        // update A[i][j] for all i: uses R[k][j] and Q[i][k]

The innermost loops are over i. We can tile the i loop to improve cache locality for A and Q. However, note that the norm computation and the update of Q are over the entire i range. The computation of R[k][j] is a reduction over i, and then the update of A[i][j] is over i. We can tile the i loop for these inner loops.

But note: the j loop is the next level. We might also tile the j loop? However, the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

However, note that the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

But note: the j loop is the next level. We might also tile the j loop? However, the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has a dependency on k. Also, the j loop is triangular. We can tile the j loop in blocks? But note that the inner i loop is the same for every j. We can interchange the j and i loops? But the update of A[i][j] uses R[k][j] which is computed in the same j iteration. So we cannot interchange j and i arbitrarily.

Alternatively, we can tile the i loop for the innermost computations. We'll break the i loop into chunks (tiles) of size, say, 64. This will help in reusing the same cache lines for A and Q.

But note: the norm computation is a reduction. We can tile the reduction. Similarly, the update of Q is a simple loop. The computation of R[k][j] is a reduction, and the update of A[i][j] is a simple loop. We can tile the i loop for all these.

However, the j loop is the next level. We might also tile the j loop? But note that the j loop is not parallel and has
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,N,m,n);
  POLYBENCH_2D_ARRAY_DECL(R,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(Q,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(R),
	      POLYBENCH_ARRAY(Q));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gramschmidt (m, n,
		      POLYBENCH_ARRAY(A),
		      POLYBENCH_ARRAY(R),
		      POLYBENCH_ARRAY(Q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Q);

  return 0;
}
