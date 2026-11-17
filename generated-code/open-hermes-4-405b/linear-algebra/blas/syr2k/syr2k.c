We are given a code for the syr2k kernel. The scop region is:

  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j <= i; j++)
      C[i][j] *= beta;
    for (k = 0; k < _PB_M; k++)
      for (j = 0; j <= i; j++)
        {
          C[i][j] += A[j][k]*alpha*B[i][k] + B[j][k]*alpha*A[i][k];
        }
  }

We are to optimize only within the scop region.

First, we look at the loop structure. There are two inner loops: one that multiplies the lower triangle of C by beta, and then a nested loop (over k and j) that updates the same lower triangle of C.

We note that the two inner loops (the first j-loop and the k-j loop) are adjacent and have the same iteration space for j (0 to i). However, the first loop is a single statement and the second is a double loop. We can consider loop fusion? But note: the first loop is a single j-loop and the second is a k-loop containing a j-loop. They are not at the same level. So fusion in the typical sense (merging two adjacent loops at the same level) is not directly applicable.

Alternatively, we can restructure the code to have a single j-loop that does both operations? But note the k-loop is nested inside the i-loop and then the j-loop. The current structure is:

  for i
    for j (0 to i): C[i][j] *= beta
    for k
      for j (0 to i): update C[i][j] with A and B

We can change the order? But note that the k-loop is independent of the first j-loop. However, the first j-loop must complete before the k-loop because the k-loop uses the updated C[i][j] (after multiplying by beta). So we cannot interchange the k-loop and the first j-loop.

But we can fuse the two j-loops? Actually, we can restructure to have:

  for i
    for j (0 to i):
        C[i][j] *= beta;
        for k: ... update C[i][j] ...

But wait, the k-loop is currently outside the j-loop. We can move the k-loop inside the j-loop? That would change the order of operations. Currently, for a fixed i, we do:
  Step 1: for each j in [0, i]: multiply by beta.
  Step 2: for each k, and for each j in [0, i]: update.

If we move the k-loop inside the j-loop, then for a fixed i and j, we do:
  Step 1: multiply by beta.
  Step 2: for each k: update.

This is equivalent because the update for a fixed (i,j) is independent of other j's? Actually, the update for (i,j) uses A[j][k] and B[i][k] and A[i][k] and B[j][k]. There is no dependence between different j's in the inner k-loop. So we can interchange the k and j loops? But note: the current order is k then j. We want to change to j then k? That would be:

  for i
    for j (0 to i):
        C[i][j] *= beta;
        for k: ... update ...

But note: the update for (i,j) uses A[j][k] and B[i][k] and A[i][k] and B[j][k]. The indices j and i are fixed in the inner k-loop. So it is safe.

However, the original code does the entire j-loop for the beta multiplication first, then the entire k-loop (which then does a j-loop). The restructured version does for each j: beta multiplication and then the k-loop for that j. This is a loop interchange (interchanging the k and j loops) and then fusing the two j-loops? Actually, we are fusing the two operations (beta multiplication and the update) for the same (i,j) and then changing the order of the k and j loops.

But note: the original k-loop is over k from 0 to _PB_M-1, and the j-loop is over j from 0 to i. The restructured version would have:

  for i
    for j from 0 to i:
        C[i][j] *= beta;
        for k from 0 to _PB_M-1:
            C[i][j] += ... 

This is valid because the update for (i,j) is independent of other j's and the k-loop does not depend on j (except through the fixed j in the inner loop). Also, the beta multiplication is done once per (i,j) and then the k-loop updates the same (i,j) multiple times. This is the same as the original.

But note: the original code does the beta multiplication for all j in [0,i] first, then the k-loop which updates all j in [0,i] for each k. The restructured code does for each j: beta multiplication and then the k-loop for that j. The order of updates for different j's is changed? Actually, the original updates for a fixed k and varying j are done in increasing j order. The restructured version does for a fixed j, all k. Then moves to the next j. This is a loop interchange (interchanging the k and j loops) and then fusing the two j-loops (the one for beta and the one for the update) into one.

This transformation is beneficial because it improves locality: for a fixed (i,j), we load C[i][j] once, multiply by beta, and then update it for all k. This keeps C[i][j] in a register. Also, we can then consider tiling the k-loop? But note the k-loop is now the innermost.

But wait: the original code has two nested loops: the outer i-loop, then inside the i-loop we have two parts: a j-loop and a k-loop that contains a j-loop. The restructured code has:

  for i
    for j (0 to i)
        ... 
        for k

This is a nested loop of depth 3: i, j, k.

Now, we can consider tiling? The innermost loop is k, which is over _PB_M. We can tile the k-loop? But note that the j-loop is triangular (j from 0 to i). Tiling a triangular loop is tricky. Alternatively, we can tile the i and j loops? But the j-loop is not rectangular.

Alternatively, we can tile the k-loop? The k-loop is a simple loop from 0 to _PB_M-1. We can break it into tiles of size TB. Then we would have:

  for i
    for j (0 to i)
        C[i][j] *= beta;
        for (kk = 0; kk < _PB_M; kk += TB)
            for (k = kk; k < kk+TB; k++)
                C[i][j] += ... 

This might help if _PB_M is large, by reducing the number of times we have to load the same C[i][j] and also by reusing A and B in cache? But note that A and B are accessed as A[j][k] and B[i][k] and A[i][k] and B[j][k]. For fixed i and j, as k varies, we are accessing the k-th column of A and B. So tiling k would help in reusing the same cache line for consecutive k? But the inner loop is over k, so without tiling we are already accessing consecutive k. Tiling might not help much for A and B because they are accessed by row? Actually, in C, arrays are row-major. So A[j][k] is the j-th row, k-th column. As k increases, we are moving along the row. So consecutive k are consecutive in memory. Therefore, the inner loop over k is already good for spatial locality.

But note: we are accessing two different rows: A[j][k] and A[i][k] (and similarly for B). So for a fixed j and i, we are accessing two rows of A and two rows of B. The cache lines for these rows might be reused as k increases? But without tiling, we are streaming through the entire row. Tiling the k-loop would break the row into chunks. This might cause more cache misses because we break the sequential access? So tiling the k-loop might not be beneficial.

Alternatively, we can consider tiling the i and j loops? But the j-loop is triangular. We can tile the i-loop and then for each tile in i, we do a rectangular block of j? But note the j-loop goes from 0 to i, so it's a triangle. We can use a technique for tiling triangular loops? For example, we can break the i-loop into blocks of size TI, and then for each block, we do j from 0 to i (which is the current block's upper bound). But then we have to handle the diagonal. This is complex and might not be necessary.

Given the complexity and the fact that the inner k-loop is already memory friendly, we might skip tiling.

Now, let's consider loop fusion: we have already fused the two operations (beta multiplication and the update) for the same (i,j) by restructuring. So we have achieved fusion.

Now, what about loop skewing? There doesn't seem to be a carried dependence that would benefit from skewing.

So we move to OpenMP.

We have a nested loop: i, then j (0 to i), then k. We want to parallelize. The outermost loop is i. We can parallelize the i-loop? But note: the j-loop goes from 0 to i, so the amount of work per i is increasing. This might lead to load imbalance. Also, the inner k-loop is the same for every i. So the total work for iteration i is (i+1) * _PB_M. This is increasing with i. So static scheduling might be bad.

We can use dynamic scheduling? But the overhead might be high. Alternatively, we can collapse the i and j loops? But note: the j-loop is triangular. Collapsing a triangular loop is not directly supported by OpenMP. We can transform the triangular loop into a rectangular one? We can change the loop to:

  for i
    for j=0 to i

into:

  for jj=0 to _PB_N-1
    i = ...; j = ...;

But note: the original j goes from 0 to i, and i from 0 to _PB_N-1. We can set:

  for (i=0; i<_PB_N; i++)
    for (j=0; j<=i; j++)

We can reindex: let jj = j, and then i must be at least j. We can do:

  for (jj=0; jj<_PB_N; jj++)
    for (i=jj; i<_PB_N; i++)

But then the inner loop is over i, and we have to do the k-loop for each (i, jj). This changes the order: originally, for a fixed i, we do j from 0 to i. Now, for a fixed jj, we do i from jj to _PB_N-1. This is a different order. But note: the update for (i,jj) is independent of other (i,jj) pairs? Actually, the update for (i,jj) uses A[jj][k] and B[i][k] and A[i][k] and B[jj][k]. There is no dependence between different (i,jj) pairs? So we can change the loop order.

But note: the original code does the beta multiplication and then the update for (i,j) in the same j-loop. In the restructured version we are considering, we would have:

  for (jj=0; jj<_PB_N; jj++)
    for (i=jj; i<_PB_N; i++) {
        // Here, j = jj
        C[i][jj] *= beta;
        for (k=0; k<_PB_M; k++)
            C[i][jj] += A[jj][k]*alpha*B[i][k] + B[jj][k]*alpha*A[i][k];
    }

This is valid. Now, the outer loop is jj, and the inner loop is i. The work per jj is (_PB_N - jj) * _PB_M. This is decreasing with jj. So the first iterations (small jj) have more work. This is still imbalanced.

We can collapse the two loops? The total number of iterations is the number of pairs (i,jj) with jj<=i<_PB_N. This is _PB_N*(_PB_N+1)/2. We can collapse the two loops and then use a static schedule? But note: the work per iteration is _PB_M (the k-loop). So if we collapse, we have a single loop of _PB_N*(_PB_N+1)/2 iterations, each doing _PB_M work. This is balanced? Yes, because each iteration does the same amount of work. So we can do:

  #pragma omp parallel for collapse(2) schedule(static)
  for (jj=0; jj<_PB_N; jj++)
    for (i=jj; i<_PB_N; i++) {
        ... 
    }

But note: the k-loop is inside. We are not parallelizing the k-loop. The k-loop is the innermost and we are leaving it sequential. This is acceptable because the k-loop is over _PB_M, which might be large, but we are parallelizing the outer two loops.

Alternatively, we can collapse the i and j loops without reindexing? The original loops are:

  for (i=0; i<_PB_N; i++)
    for (j=0; j<=i; j++)

We can collapse these two? But OpenMP collapse requires rectangular loops. The j-loop is triangular. So we cannot collapse without reindexing.

So we choose the reindexed version.

But note: the original code uses _PB_N and _PB_M. We must use the same.

Now, let's write the transformed code:

  #pragma omp parallel for collapse(2) schedule(static) default(none) \
    private(i, j, k, jj) shared(C, A, B, alpha, beta, _PB_N, _PB_M)
  for (jj=0; jj<_PB_N; jj++) {
    for (i=jj; i<_PB_N; i++) {
        C[i][jj] *= beta;
        for (k=0; k<_PB_M; k++) {
            C[i][jj] += A[jj][k] * alpha * B[i][k] + B[jj][k] * alpha * A[i][k];
        }
    }
  }

But note: we have to be cautious about variable sharing. We list the shared variables: C, A, B, alpha, beta, _PB_N, _PB_M. The private variables: i, j, k, jj. But note: j is not used in the inner loop? We removed j and replaced by jj. So we don't need j. Also, k is private by default? But we list it to be safe.

But wait: the original code had a variable j that was used in the inner loop. We replaced it by jj. So we don't have j anymore.

Alternatively, we can avoid reindexing and use the original i and j, but then we cannot collapse. We can parallelize the i-loop and then use dynamic scheduling to account for the increasing work? 

  #pragma omp parallel for schedule(dynamic) default(none) \
    private(i, j, k) shared(C, A, B, alpha, beta, _PB_N, _PB_M)
  for (i=0; i<_PB_N; i++) {
    for (j=0; j<=i; j++) {
        C[i][j] *= beta;
        for (k=0; k<_PB_M; k++) {
            C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
        }
    }
  }

This avoids reindexing. The dynamic schedule will assign chunks of i to threads. Since the work per i is (i+1)*_PB_M, which increases with i, dynamic scheduling can balance the load. The overhead of dynamic scheduling might be acceptable if _PB_N is large.

But note: the inner k-loop is heavy (_PB_M iterations). So the granularity per i is (i+1)*_PB_M, which is large for big i. So dynamic scheduling might be efficient.

We have two options. We choose the dynamic scheduling without reindexing because it is simpler and avoids changing the loop structure.

But note: the original code had the k-loop outside the j-loop? We restructured to have the k-loop inside the j-loop. So we are already changing the structure. We did that for locality. So we are going to use:

  for (i=0; i<_PB_N; i++) {
    for (j=0; j<=i; j++) {
        C[i][j] *= beta;
        for (k=0; k<_PB_M; k++) {
            C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
        }
    }
  }

and then parallelize the i-loop with dynamic scheduling.

But note: the inner k-loop is the same for every j in the same i. We can interchange the j and k loops? That would make the k-loop outer and j-loop inner. Then we can collapse the i and k loops? But the k-loop is over _PB_M, which might be large. And the j-loop is triangular. So we would have:

  for (i=0; i<_PB_N; i++) {
    for (k=0; k<_PB_M; k++) {
        for (j=0; j<=i; j++) {
            // But note: the beta multiplication must be done only once per (i,j). 
            // We cannot move the beta multiplication inside the k-loop because it would be done _PB_M times.
        }
    }
  }

So we cannot do that.

Therefore, we stick with the restructured version (with k inside j) and parallelize the i-loop with dynamic scheduling.

Now, about data sharing: 
  - C, A, B, alpha, beta, _PB_N, _PB_M are shared.
  - i, j, k are private.

We can use 'default(none)' to force explicit specification.

Also, note that the inner k-loop is heavy, so we might consider SIMD? We can add 'simd' to the inner k-loop? But note: the inner k-loop is the one we are parallelizing? Actually, we are parallelizing the i-loop. The inner k-loop is sequential. We can add:

  #pragma omp simd
  for (k=0; k<_PB_M; k++)

But we must check for dependencies: the k-loop is a reduction on C[i][j]. So it is safe.

We can do:

  for (i=0; i<_PB_N; i++) {
    for (j=0; j<=i; j++) {
        C[i][j] *= beta;
        #pragma omp simd
        for (k=0; k<_PB_M; k++) {
            C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
        }
    }
  }

But note: the entire i-loop is parallelized. The inner j-loop is sequential. The k-loop is vectorized.

Alternatively, we can use 'parallel for simd' on the k-loop? But that would be nested parallelism, which is not allowed. We are already parallelizing the i-loop.

So we do:

  #pragma omp parallel for schedule(dynamic) default(none) \
    private(i, j, k) shared(C, A, B, alpha, beta, _PB_N, _PB_M)
  for (i=0; i<_PB_N; i++) {
    for (j=0; j<=i; j++) {
        C[i][j] *= beta;
        #pragma omp simd
        for (k=0; k<_PB_M; k++) {
            C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
        }
    }
  }

But note: the original code used _PB_N and _PB_M. We must use the same.

Also, we must preserve the scop markers.

Let's write the entire scop region:

#pragma scop
  #pragma omp parallel for schedule(dynamic) default(none) \
    private(i, j, k) shared(C, A, B, alpha, beta, _PB_N, _PB_M)
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j <= i; j++) {
        C[i][j] *= beta;
        #pragma omp simd
        for (k = 0; k < _PB_M; k++) {
            C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
        }
    }
  }
#pragma endscop

But note: the original code had the k-loop outside the j-loop. We changed the order. We believe it is safe.

We must check: the original code did the beta multiplication for all j in [0,i] first, then the k-loop which updated all j in [0,i] for each k. The new code does for each j: beta multiplication and then the k-loop for that j. This is equivalent because the update for (i,j) is independent of other j's and the k-loop for a fixed (i,j) is independent of the order of j's.

Therefore, we output the code with the scop region replaced by the above.

But note: the original code used variables _PB_N and _PB_M. We use the same.

Also, we must not change anything outside the scop.

We write the entire kernel_syr2k function with the scop region replaced.

Note: the original code had:

  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j <= i; j++)
      C[i][j] *= beta;
    for (k = 0; k < _PB_M; k++)
      for (j = 0; j <= i; j++)
        {
          C[i][j] += A[j][k]*alpha*B[i][k] + B[j][k]*alpha*A[i][k];
        }
  }

We replace it with the restructured and parallelized version.

Final code for the scop region:

#pragma scop
  #pragma omp parallel for schedule(dynamic) default(none) \
    private(i, j, k) shared(C, A, B, alpha, beta, _PB_N, _PB_M)
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j <= i; j++) {
        C[i][j] *= beta;
        #pragma omp simd
        for (k = 0; k < _PB_M; k++) {
            C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
        }
    }
  }
#pragma endscop

We output the entire code with this change.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "syr2k.h"


/* Array initialization. */
static
void init_array(int n, int m,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		DATA_TYPE POLYBENCH_2D(B,N,M,n,m))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      A[i][j] = (DATA_TYPE) ((i*j+1)%n) / n;
      B[i][j] = (DATA_TYPE) ((i*j+2)%m) / m;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      C[i][j] = (DATA_TYPE) ((i*j+3)%n) / m;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_syr2k(int n, int m,
		  DATA_TYPE alpha,
		  DATA_TYPE beta,
		  DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		  DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		  DATA_TYPE POLYBENCH_2D(B,N,M,n,m))
{
  int i, j, k;

//BLAS PARAMS
//UPLO  = 'L'
//TRANS = 'N'
//A is NxM
//B is NxM
//C is NxN
#pragma scop
  #pragma omp parallel for schedule(dynamic) default(none) \
    private(i, j, k) shared(C, A, B, alpha, beta, _PB_N, _PB_M)
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j <= i; j++) {
        C[i][j] *= beta;
        #pragma omp simd
        for (k = 0; k < _PB_M; k++) {
            C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
        }
    }
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,N,M,n,m);

  /* Initialize array(s). */
  init_array (n, m, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_syr2k (n, m,
		alpha, beta,
		POLYBENCH_ARRAY(C),
		POLYBENCH_ARRAY(A),
		POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}