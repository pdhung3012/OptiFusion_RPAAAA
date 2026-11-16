We are given a code for LU decomposition (ludcmp) and are to optimize the region between #pragma scop and #pragma endscop.
 The code has three main parts:
  1. LU decomposition of matrix A (the first big loop over i, j, k)
  2. Forward substitution to solve L y = b (the second loop over i, j)
  3. Backward substitution to solve U x = y (the third loop over i, j)

 We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

 Step 1: Loop Tiling
   We look for nested loops that can be tiled to improve cache locality. The LU decomposition part has three nested loops (i, j, k). 
   However, note that the LU decomposition is done in two parts: the lower triangle (j < i) and the upper triangle (j >= i). 
   The inner loop (k) is over a triangular region. Tiling such loops is complex and may require careful handling of boundaries.

   But note: the problem size is _PB_N, which is typically large. We can try to tile the i and j loops. However, the dependencies in LU decomposition are such that:
     - For the lower triangle: A[i][j] depends on A[i][k] and A[k][j] for k<j, and then we write A[i][j] = ... / A[j][j].
     - For the upper triangle: A[i][j] depends on A[i][k] and A[k][j] for k<i, and then we write A[i][j] = ... .

   The standard LU decomposition has a wavefront pattern. Tiling is possible but non-trivial. We might consider tiling the i-loop and then the j-loop, but note that the inner loop (k) is triangular and the bounds depend on j (for the lower triangle) and i (for the upper triangle).

   Given the complexity and the fact that the inner loop is already over a small set (k from 0 to j-1 or i-1), we might not get much benefit from tiling the inner loop. Instead, we can try to tile the outer loops (i and j) to improve locality in the access to A.

   However, note that the two inner loops (for the lower and upper triangle) are adjacent and have the same outer loop (i). We might consider fusing them? But they are already in the same i-loop. Actually, the j-loop is split into two: one for j from 0 to i-1 and one for j from i to _PB_N-1.

   Alternatively, we can try to tile the i-loop and then the j-loop. But note that the j-loop in the lower triangle goes from 0 to i-1, which is triangular. We can use a rectangular tile and then adjust the bounds.

   However, after careful thought, the LU decomposition is often optimized by blocking (tiling) the entire algorithm. But the provided code is the naive triple loop. We are allowed to restructure the loops. We can try to tile the i-loop and then the j-loop, but note that the inner k-loop is over a small set and might not benefit from tiling.

   Considering the priority, we will attempt to tile the outer two loops (i and j) for the LU decomposition part. We choose a tile size, say 64, for the i-loop and then for the j-loop. But note that the j-loop in the lower triangle is triangular. We can do:

     for (ii = 0; ii < _PB_N; ii += TILE_I)
        for (i = ii; i < min(ii+TILE_I, _PB_N); i++)
          for (j = 0; j < i; j++)   // lower triangle
             ... 

     and similarly for the upper triangle? But note the upper triangle j-loop goes from i to _PB_N-1.

   However, the two j-loops are adjacent and we can try to fuse them? Actually, they are already in the same i-loop. We can try to tile the i-loop and then the j-loop for both parts together? But the j-loop for the lower triangle is from 0 to i-1 and for the upper triangle from i to _PB_N-1. We can combine them into one j-loop from 0 to _PB_N-1 and then conditionally do the lower or upper part? That would be:

        for (i = 0; i < _PB_N; i++)
          for (j = 0; j < _PB_N; j++)
             if (j < i) { ... } else { ... }

   But note: the original code does the lower triangle first and then the upper triangle. The condition j < i is the lower triangle and j>=i is the upper. However, the inner k-loop for the lower triangle goes to j-1 and for the upper triangle to i-1. So we can write:

        for (i = 0; i < _PB_N; i++) {
          for (j = 0; j < _PB_N; j++) {
             w = A[i][j];
             if (j < i) {
                 for (k = 0; k < j; k++)
                    w -= A[i][k] * A[k][j];
                 A[i][j] = w / A[j][j];
             } else {
                 for (k = 0; k < i; k++)
                    w -= A[i][k] * A[k][j];
                 A[i][j] = w;
             }
          }
        }

   This is a fusion of the two j-loops. But note: the original code does the lower triangle for j from 0 to i-1 and then the upper for j from i to _PB_N-1. The fused version does the same but in one j-loop. This is legal because the two j-loops are independent (they write to different parts of A[i][j]). However, the inner k-loop for the lower triangle uses j as the upper bound and for the upper triangle uses i. This is acceptable.

   Now, we can tile the i and j loops. We choose a tile size, say 64. We do:

        #define TILE_I 64
        #define TILE_J 64

        for (ii = 0; ii < _PB_N; ii += TILE_I)
          for (i = ii; i < min(ii+TILE_I, _PB_N); i++)
            for (jj = 0; jj < _PB_N; jj += TILE_J)
              for (j = jj; j < min(jj+TILE_J, _PB_N); j++) {
                 ... // the same body as above

   But note: the inner k-loop for the lower triangle (j < i) goes to j-1, which is within the tile? Not necessarily, because k goes from 0 to j-1, which might be in a previous tile. Similarly, for the upper triangle, k goes to i-1, which might be in a previous tile. So we are not improving the locality of the inner loop. The inner loop is over k, which is not tiled. Therefore, tiling the i and j loops might not help the inner k-loop because the k-loop accesses A[i][k] and A[k][j] for k from 0 to j-1 (or i-1). The access to A[i][k] is row-wise and contiguous, but the access to A[k][j] is column-wise and non-contiguous. Tiling the i and j loops might help in keeping a tile of A in cache, but the inner k-loop might still cause cache misses for A[k][j] because j is fixed and k varies.

   Alternatively, we can try to interchange the loops? But the dependencies are such that we cannot interchange arbitrarily.

   Given the complexity and the fact that the inner k-loop is the most time-consuming, we might consider tiling the k-loop? But the k-loop is triangular and the bounds depend on j or i. This is difficult.

   After careful thought, we decide that tiling the entire LU decomposition is beyond the scope of a simple transformation. We will leave the LU decomposition as is and focus on the forward and backward substitutions, which are simpler.

   However, note that the forward and backward substitutions are also nested loops. We can try to tile them? But they are not heavy enough to justify tiling? The inner loop is over j, which is triangular and the bounds are small (j from 0 to i-1 for forward, and j from i+1 to _PB_N-1 for backward). So we skip tiling for the entire scop.

 Step 2: Loop Fusion
   We have three separate loops: LU decomposition, forward substitution, and backward substitution. They are independent and cannot be fused because they are in sequence and the second and third depend on the first. So no fusion.

 Step 3: Loop Skewing
   We look for carried dependencies that can be skewed to expose parallelism. The LU decomposition has a carried dependence: the computation of A[i][j] depends on A[k][j] for k<j (in the lower triangle) and on A[i][k] and A[k][j] for k<i (in the upper triangle). The forward substitution: y[i] depends on y[j] for j<i. The backward substitution: x[i] depends on x[j] for j>i. These are all carried dependencies that prevent parallelization of the outer loop.

   However, note that the LU decomposition has a wavefront pattern. We can skew the loops to expose parallelism? For example, we can change the order of computation to a wavefront. But the standard LU decomposition is often parallelized by wavefront or by using a block cyclic approach. However, the problem does not require changing the algorithm to a block version. We are to use loop transformations.

   We can try to skew the i and j loops? For the lower triangle, we have:

        for (i=0; i<N; i++)
          for (j=0; j<i; j++)
             A[i][j] = ... 

   The dependence: A[i][j] depends on A[i][k] for k<j (which is in the same row, so no dependence on previous i) and on A[k][j] for k<j (which is in the same column, so we have a dependence on previous rows). Specifically, the value A[k][j] is computed in a previous iteration (k<j, and k<i because j<i). So the dependence is along the k dimension. We can try to skew the i and j loops? 

   We can introduce a new index: s = i + j. Then we can iterate over s and then over j? But note that the inner k-loop is over k from 0 to j-1, which is independent of i. This might not help.

   Alternatively, we can try to parallelize the inner loops? The inner k-loop is a reduction? Actually, it is a linear combination. We can try to parallelize the k-loop? But the k-loop is not parallel because each iteration updates w. We can use a reduction? But note that the k-loop is a simple linear update and the order does not matter. We can change the k-loop to a parallel reduction? But the k-loop is small (j-1 or i-1) and the overhead of parallelization might be too high.

   Given the complexity, we decide not to apply skewing.

 Step 4: OpenMP Directives

   We now apply OpenMP directives to the three main loops.

   First loop: LU decomposition (the big i-loop). We note that the outer i-loop cannot be parallelized because of the carried dependencies (each row depends on previous rows). However, we can parallelize the inner j-loop? But note that the j-loop for the lower triangle and the j-loop for the upper triangle are independent for fixed i? Actually, the lower triangle j-loop (j from 0 to i-1) and the upper triangle j-loop (j from i to _PB_N-1) are independent. But the lower triangle j-loop has a dependence on previous rows (via A[k][j] for k<j) and the upper triangle j-loop also has a dependence on previous rows (via A[k][j] for k<i). However, for a fixed i, the j-loop iterations are independent? 

   Let's check the lower triangle: for fixed i, the j-loop iterations (j from 0 to i-1) are independent? The inner k-loop for a given j uses A[i][k] and A[k][j]. The A[i][k] is from the same row and is read-only for the inner k-loop? Actually, we are writing A[i][j] at the end. But note: the inner k-loop for j1 and j2 (with j1 != j2) do not interfere because they write to different A[i][j]. However, they both read A[i][k] for k in [0, j-1] (for j1) and [0, j2-1] (for j2). Since j1 and j2 are different, the sets of k might overlap? But the reads are from the same row and are read-only. So the j-loop for the lower triangle is parallel for fixed i? 

   Similarly, the upper triangle j-loop: for fixed i, the j-loop iterations (j from i to _PB_N-1) are independent? The inner k-loop uses A[i][k] and A[k][j]. Again, the writes are to different A[i][j] and the reads from A[i][k] are from the same row and read-only. The reads from A[k][j] are from different columns and are read-only. So the j-loop for the upper triangle is also parallel for fixed i.

   Therefore, we can parallelize the j-loop for each i. But note: the j-loop is split into two. We can combine them and then parallelize the entire j-loop for fixed i? We did the fusion above. So we can do:

        for (i = 0; i < _PB_N; i++) {
          #pragma omp parallel for private(w, k) // we can use 'parallel for' for the j-loop
          for (j = 0; j < _PB_N; j++) {
             ... // the fused body
          }
        }

   However, note that the inner k-loop is over k, and we must declare k as private. Also, w is private to each j iteration.

   But wait: the original code does the lower triangle first and then the upper triangle. The fused version does the same but in one j-loop. The order of j is from 0 to _PB_N-1. This is acceptable because the lower triangle j (j<i) and the upper triangle j (j>=i) are independent.

   However, note that the inner k-loop for the lower triangle uses j as the upper bound and for the upper triangle uses i. This is still correct.

   But note: the original code for the lower triangle does j from 0 to i-1 and then the upper triangle does j from i to _PB_N-1. The fused version does j from 0 to _PB_N-1. The condition j < i is the lower triangle and j>=i is the upper. This is the same.

   Therefore, we can parallelize the j-loop for each i.

   Second loop: forward substitution (the i-loop for y). This loop has a carried dependence: y[i] depends on y[j] for j<i. So we cannot parallelize the i-loop. But we can parallelize the inner j-loop? The inner j-loop is a reduction on w. We can use a reduction? But note that the inner j-loop is:

        w = b[i];
        for (j=0; j<i; j++)
           w -= A[i][j] * y[j];

   This is a reduction on w. We can parallelize the j-loop with a reduction on w? But note that the j-loop is triangular and the number of iterations is i, which is small for small i and large for large i. The overhead of parallelization might be too high for small i. We can use a schedule that uses static scheduling and hope that for large i the parallelization is beneficial? But note that the j-loop is not heavy enough to justify parallelization? The inner loop body is a single statement. We might leave it sequential.

   Third loop: backward substitution (the i-loop for x). This loop runs from _PB_N-1 down to 0. It has a carried dependence: x[i] depends on x[j] for j>i. So we cannot parallelize the i-loop. The inner j-loop is:

        w = y[i];
        for (j=i+1; j<_PB_N; j++)
           w -= A[i][j] * x[j];

   This is a reduction on w. Similarly, we can parallelize the j-loop with a reduction? But again, the j-loop is triangular and the number of iterations is _PB_N - i - 1, which is small for large i and large for small i. We might leave it sequential.

   However, note that the forward and backward substitutions are not the most time-consuming part. The LU decomposition is the heavy part.

   Therefore, we will only parallelize the j-loop in the LU decomposition.

   We also note that the LU decomposition j-loop is over _PB_N iterations, which is large. We can use a static schedule.

   We must be cautious: the inner k-loop is over k, and we must declare k as private. Also, w is private to each j iteration.

   We write:

        for (i = 0; i < _PB_N; i++) {
          #pragma omp parallel for private(j, w, k) // j is the loop variable, so we don't need to declare it private? Actually, the loop variable is private by default. But we have k and w.
          for (j = 0; j < _PB_N; j++) {
             w = A[i][j];
             if (j < i) {
                 for (k = 0; k < j; k++)
                    w -= A[i][k] * A[k][j];
                 A[i][j] = w / A[j][j];
             } else {
                 for (k = 0; k < i; k++)
                    w -= A[i][k] * A[k][j];
                 A[i][j] = w;
             }
          }
        }

   But note: the original code had two separate j-loops. We fused them. This is a transformation we did for the purpose of parallelization. We must ensure that the fused loop is correct.

   However, we must check the data sharing: 
     - A is shared (read-write). We are writing to A[i][j] in different j iterations. Since each j iteration writes to a different element, it is safe.
     - b, x, y are not used in this loop.

   We also note that the inner k-loop for the lower triangle uses A[k][j] for k<j. Since j is the current iteration and k<j, and we are processing j in parallel, we must ensure that the value of A[k][j] is not being written by another thread? But note: for a fixed i, we are only writing to A[i][j] (the current row). The A[k][j] for k<j are from previous rows and are not being written in this i iteration. So it is safe.

   Similarly, for the upper triangle, we use A[k][j] for k<i. These are from previous rows and are not being written in this i iteration. So it is safe.

   Therefore, we can parallelize the j-loop for each i.

   We also note that the entire LU decomposition is now:

        for (i=0; i<_PB_N; i++) {
          #pragma omp parallel for private(w, k)
          for (j=0; j<_PB_N; j++) {
             ... 
          }
        }

   We can also consider collapsing the i and j loops? But note that the i-loop is sequential (because of the carried dependence on i) and the j-loop is parallel. We cannot collapse because the i-loop is not parallel.

   We can also consider using SIMD for the inner k-loop? The inner k-loop is a reduction. We can use #pragma omp simd reduction(+:w) for the k-loop? But note that the k-loop is small and the compiler might vectorize it automatically. We can try:

        for (k = 0; k < j; k++)
           w -= A[i][k] * A[k][j];

        becomes:

        #pragma omp simd reduction(+:w)
        for (k = 0; k < j; k++)
           w -= A[i][k] * A[k][j];

   Similarly for the other k-loop.

   We can do this for both the lower and upper triangle.

   For the forward substitution, we can also try to vectorize the inner j-loop? It is a reduction. We can do:

        w = b[i];
        #pragma omp simd reduction(+:w)
        for (j=0; j<i; j++)
           w -= A[i][j] * y[j];

   Similarly for the backward substitution:

        w = y[i];
        #pragma omp simd reduction(+:w)
        for (j=i+1; j<_PB_N; j++)
           w -= A[i][j] * x[j];

   We must be cautious: the j-loop in the forward substitution goes from 0 to i-1, and in the backward from i+1 to _PB_N-1. The bounds are triangular. The compiler might handle it.

   We also note that the forward and backward substitutions are not parallelized at the outer level, but we can try to use SIMD for the inner loop.

   We will add SIMD directives to the inner k-loops in the LU decomposition and to the inner j-loops in the forward and backward substitutions.

   We also note that the data might be aligned? We don't know. We can add an aligned clause if we know the alignment, but we don't. We skip.

   We also note that the reduction in the inner loops is a subtraction. We use reduction(+:w) because the operation is subtraction? Actually, the reduction is for the entire expression: w = w - (A[i][k] * A[k][j]). This is equivalent to w += (- (A[i][k] * A[k][j])). So we can use reduction(+:w). The compiler will handle it.

   Therefore, we will add:

        In the LU decomposition, for the k-loop in the lower triangle:
          #pragma omp simd reduction(+:w)
          for (k = 0; k < j; k++)

        and similarly for the upper triangle.

        In the forward substitution:
          #pragma omp simd reduction(+:w)
          for (j=0; j<i; j++)

        In the backward substitution:
          #pragma omp simd reduction(+:w)
          for (j=i+1; j<_PB_N; j++)

   We also note that the outer i-loop in the LU decomposition is sequential. We can try to use multiple threads for the j-loop, which we are already doing.

   We also note that the forward and backward substitutions are sequential in the outer loop, but we are vectorizing the inner loop.

   We must be cautious: the inner loop in the forward substitution for small i might not be vectorized efficiently. But the compiler will handle it.

   We also note that the variable w in the inner loops must be private to the SIMD lane. The reduction clause takes care of that.

   We also note that the data types are DATA_TYPE, which is typically double. We assume the compiler can vectorize double precision.

   We also note that the original code uses _PB_N, which is a macro for the problem size.

   We will now write the optimized code for the scop region.

   Summary of changes:

     - Fused the two j-loops in the LU decomposition into one j-loop from 0 to _PB_N-1, with a condition on j.
     - Added OpenMP parallel for to the j-loop in the LU decomposition (for each i).
     - Added OpenMP simd with reduction(+:w) to the inner k-loops in the LU decomposition and to the inner j-loops in the forward and backward substitutions.

   We do not change the forward and backward substitutions outer loops.

   We also note that the original code uses variables i, j, k, w. We must declare them appropriately. The original code declares them at the beginning of the function. We are not changing the function signature. We are only changing the scop region.

   We must be cautious: the variable w is declared outside the loops. We are using it inside the parallel region. We must make it private in the parallel for. We are already doing that.

   We also note that the inner k-loop in the LU decomposition uses k, which we declare as private in the parallel for.

   We also note that the forward and backward substitutions use w, which is declared at the beginning. We are using it in the inner loop and then assigning to y[i] or x[i]. We are not changing the scope of w.

   We write the code accordingly.

   Note: We must preserve the original variable declarations. We are not changing them.

   We also note that the original code uses _PB_N, which is a macro. We leave it as is.

   We also note that the scop region includes the entire kernel_ludcmp function body? Actually, the scop region is marked and includes the three loops.

   We output the entire code with the scop region replaced.

   However, note that the original code has:

        for (i = 0; i < _PB_N; i++) {
          for (j = 0; j <i; j++) { ... }
          for (j = i; j < _PB_N; j++) { ... }
        }

   We are replacing it with:

        for (i = 0; i < _PB_N; i++) {
          #pragma omp parallel for private(w, k)
          for (j = 0; j < _PB_N; j++) {
             w = A[i][j];
             if (j < i) {
                 #pragma omp simd reduction(+:w)
                 for (k = 0; k < j; k++) {
                    w -= A[i][k] * A[k][j];
                 }
                 A[i][j] = w / A[j][j];
             } else {
                 #pragma omp simd reduction(+:w)
                 for (k = 0; k < i; k++) {
                    w -= A[i][k] * A[k][j];
                 }
                 A[i][j] = w;
             }
          }
        }

   And then the forward substitution:

        for (i = 0; i < _PB_N; i++) {
           w = b[i];
           #pragma omp simd reduction(+:w)
           for (j = 0; j < i; j++)
              w -= A[i][j] * y[j];
           y[i] = w;
        }

   And the backward substitution:

        for (i = _PB_N-1; i >=0; i--) {
           w = y[i];
           #pragma omp simd reduction(+:w)
           for (j = i+1; j < _PB_N; j++)
              w -= A[i][j] * x[j];
           x[i] = w / A[i][i];
        }

   We also note that the original code uses the same variable w for all. We are reusing it. This is safe because the scop region is the entire function body? Actually, the scop region is the three loops. We are not changing the variable w outside the loops.

   We also note that the original code has the scop markers. We must preserve them.

   We output the entire code with the scop region replaced by the above.

   However, note that the original code uses _PB_N, which is a macro. We leave it as is.

   We also note that the function parameters are the same.

   We also note that the code uses POLYBENCH_2D and POLYBENCH_1D macros. We leave them as is.

   We also note that the code uses DATA_TYPE. We leave it as is.

   We also note that the code uses the variable names as in the original.

   We also note that the original code has the scop markers exactly as:

        #pragma scop
        ... code ...
        #pragma endscop

   We must preserve these lines.

   Therefore, we write the optimized code accordingly.

   Let's write the entire kernel_ludcmp function with the scop region replaced.

   Note: We are not changing the function signature.

   We also note that the original code has the scop region including the three loops. We are replacing the entire scop region.

   We write:
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "ludcmp.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_1D(b,N,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n))
{
  int i, j;
  DATA_TYPE fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
    {
      x[i] = 0;
      y[i] = 0;
      b[i] = (i+1)/fn/2.0 + 4;
    }

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
	A[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      for (j = i+1; j < n; j++) {
	A[i][j] = 0;
      }
      A[i][i] = 1;
    }

  /* Make the matrix positive semi-definite. */
  /* not necessary for LU, but using same code as cholesky */
  int r,s,t;
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      (POLYBENCH_ARRAY(B))[r][s] = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	(POLYBENCH_ARRAY(B))[r][s] += A[r][t] * A[s][t];
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	A[r][s] = (POLYBENCH_ARRAY(B))[r][s];
  POLYBENCH_FREE_ARRAY(B);

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(x,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("x");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, x[i]);
  }
  POLYBENCH_DUMP_END("x");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_ludcmp(int n,
		   DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		   DATA_TYPE POLYBENCH_1D(b,N,n),
		   DATA_TYPE POLYBENCH_1D(x,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n))
{
  int i, j, k;

  DATA_TYPE w;

#pragma scop
  for (i = 0; i < _PB_N; i++) {
    #pragma omp parallel for private(j, w, k) default(none) shared(A, i, _PB_N)
    for (j = 0; j < _PB_N; j++) {
      w = A[i][j];
      if (j < i) {
        #pragma omp simd reduction(+:w) safelen(64)
        for (k = 0; k < j; k++) {
          w -= A[i][k] * A[k][j];
        }
        A[i][j] = w / A[j][j];
      } else {
        #pragma omp simd reduction(+:w) safelen(64)
        for (k = 0; k < i; k++) {
          w -= A[i][k] * A[k][j];
        }
        A[i][j] = w;
      }
    }
  }

  for (i = 0; i < _PB_N; i++) {
    w = b[i];
    #pragma omp simd reduction(+:w) safelen(64)
    for (j = 0; j < i; j++) {
      w -= A[i][j] * y[j];
    }
    y[i] = w;
  }

  for (i = _PB_N-1; i >= 0; i--) {
    w = y[i];
    #pragma omp simd reduction(+:w) safelen(64)
    for (j = i+1; j < _PB_N; j++) {
      w -= A[i][j] * x[j];
    }
    x[i] = w / A[i][i];
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(b, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(b),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(y));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_ludcmp (n,
		 POLYBENCH_ARRAY(A),
		 POLYBENCH_ARRAY(b),
		 POLYBENCH_ARRAY(x),
		 POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(x)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(b);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}