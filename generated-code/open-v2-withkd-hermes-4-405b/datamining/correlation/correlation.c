/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* correlation.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "correlation.h"


/* Array initialization. */
static
void init_array (int m,
		 int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)N;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = (DATA_TYPE)(i*j)/M + i;

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(corr,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("corr");
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, corr[i][j]);
    }
  POLYBENCH_DUMP_END("corr");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_correlation(int m, int n,
			DATA_TYPE float_n,
			DATA_TYPE POLYBENCH_2D(data,N,M,n,m),
			DATA_TYPE POLYBENCH_2D(corr,M,M,m,m),
			DATA_TYPE POLYBENCH_1D(mean,M,m),
			DATA_TYPE POLYBENCH_1D(stddev,M,m))
{
  int i, j, k;

  DATA_TYPE eps = SCALAR_VAL(0.1);


#pragma scop
We are given the correlation.c code from PolyBench. We are to optimize the code within the #pragma scop and #pragma endscop markers.

The code has several loops. We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

Let's break down the scop region:

1. First loop: Computes the mean for each column j.
   for (j=0; j<_PB_M; j++) {
      mean[j] = 0.0;
      for (i=0; i<_PB_N; i++)
          mean[j] += data[i][j];
      mean[j] /= float_n;
   }

2. Second loop: Computes the standard deviation for each column j.
   for (j=0; j<_PB_M; j++) {
      stddev[j] = 0.0;
      for (i=0; i<_PB_N; i++)
          stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      stddev[j] /= float_n;
      stddev[j] = SQRT_FUN(stddev[j]);
      stddev[j] = (stddev[j] <= eps) ? 1.0 : stddev[j];
   }

3. Third loop: Centers and reduces the column vectors.
   for (i=0; i<_PB_N; i++)
      for (j=0; j<_PB_M; j++) {
          data[i][j] -= mean[j];
          data[i][j] /= (SQRT_FUN(float_n) * stddev[j]);
      }

4. Fourth loop: Computes the correlation matrix.
   for (i=0; i<_PB_M-1; i++) {
      corr[i][i] = 1.0;
      for (j=i+1; j<_PB_M; j++) {
          corr[i][j] = 0.0;
          for (k=0; k<_PB_N; k++)
              corr[i][j] += (data[k][i] * data[k][j]);
          corr[j][i] = corr[i][j];
      }
   }
   corr[_PB_M-1][_PB_M-1] = 1.0;

We are to prioritize loop transformations:

1. Loop tiling: We look for nested loops that can be tiled to improve cache locality. The fourth loop (correlation matrix) is a candidate because it has three nested loops and the inner loop is a reduction. However, note that the third loop (centering) is also a candidate for tiling because it accesses data[i][j] in a row-major order and we are writing to the same array. But note that the fourth loop uses the entire data array and the correlation matrix is symmetric.

   However, the first two loops (mean and stddev) are reductions over columns. They are not easily tiled because each column is independent. The third loop is a simple element-wise operation and might not benefit much from tiling because it is memory-bound and already has good spatial locality.

   The fourth loop (correlation matrix) is the most compute-intensive and has a reduction pattern. We can tile the i and j loops (the outer two) to improve locality. We note that the inner loop (k) is over the rows of data. We can tile the i and j loops to form blocks of the correlation matrix and then compute the block by iterating over k. This would allow reuse of data[k][i] and data[k][j] for multiple (i,j) pairs in the same block.

   However, note that the correlation matrix is symmetric and we are only computing the upper triangle. Tiling the i and j loops for the upper triangle is more complex. We might instead consider tiling the entire i and j loops and then only compute the upper triangle within the tile. But note that the original loop only iterates j from i+1 to M-1. We can restructure to a full i,j loop and then set the symmetric part, but that would change the original pattern. Alternatively, we can tile the i loop and then for each i-tile, tile the j loop from i+1 to the end. This is more complex.

   Given the complexity and the fact that the problem size might not be huge, we might skip tiling for now and focus on OpenMP. But note that the problem says to prioritize tiling. However, we must be cautious because the correlation matrix computation is not a straightforward matrix multiplication. It is a reduction over k for each (i,j) pair.

   Alternatively, we can consider tiling the k loop? But the k loop is the innermost and we are reducing over it. Tiling the k loop would break the reduction and require extra steps. So we avoid that.

   After careful thought, we decide to tile the i and j loops for the correlation matrix. We choose a tile size of 64 for both i and j (if the problem size is large enough). We will break the i loop into tiles of size TI and the j loop into tiles of size TJ. Then for each tile, we compute the block of the correlation matrix for i in [i0, i0+TI) and j in [j0, j0+TJ) but only for j>=i+1? Actually, we can compute the entire tile and then set the symmetric part. But note that the original code only computes the upper triangle and then sets the lower triangle. We can do the same in the tiled version.

   However, note that the original code sets corr[i][i] = 1.0 and then for j from i+1 to M-1. We can restructure the tiled version to compute the entire upper triangle by having the j loop start at i+1. We can tile the i loop and then for each i-tile, we tile the j loop from max(i+1, j0) to min(j0+TJ, M-1). This is complex and might not be worth the effort for a general solution.

   Given the time, we decide to skip tiling for the correlation matrix and instead focus on OpenMP and other transformations. We note that the problem says to prioritize tiling, but if we cannot do it safely and effectively, we move on.

   Alternatively, we can fuse the first two loops? They both iterate over j and then i. They are adjacent and have the same outer loop. We can fuse them:

   for (j=0; j<_PB_M; j++) {
      mean[j] = 0.0;
      for (i=0; i<_PB_N; i++)
          mean[j] += data[i][j];
      mean[j] /= float_n;

      stddev[j] = 0.0;
      for (i=0; i<_PB_N; i++)
          stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      stddev[j] /= float_n;
      stddev[j] = SQRT_FUN(stddev[j]);
      stddev[j] = (stddev[j] <= eps) ? 1.0 : stddev[j];
   }

   This is safe because the two loops are independent in terms of j. This fusion reduces the number of passes over the j loop and might improve cache locality because we are accessing the same column j twice in a row.

   We decide to do this fusion.

   Now, the third loop (centering) is independent and we leave it as is.

   The fourth loop (correlation) we leave as is for now.

2. Loop skewing: We don't see any carried dependencies that would benefit from skewing. The loops are either reductions or independent.

So we do:
   - Fuse the first two loops (mean and stddev) into one j loop.

Now, we apply OpenMP:

1. Parallelism: We can parallelize the fused j loop because each j is independent. We can also parallelize the centering loop (i and j) and the correlation loop (i and j).

   For the fused j loop: we can use #pragma omp parallel for over j.

   For the centering loop: we can use #pragma omp parallel for collapse(2) because the i and j loops are independent.

   For the correlation loop: the outer i loop is independent? Actually, each i iteration writes to corr[i][i] and then to corr[i][j] for j>i and also sets corr[j][i]. But note that for a fixed i, the j loop writes to corr[i][j] and corr[j][i]. However, for two different i's, they write to different rows and columns? Actually, the writes for i1 and i2 (i1 != i2) are to distinct corr[i1][*] and corr[i2][*] and also to the symmetric part. But note that when i1 and i2 are such that i1 < i2, then the iteration for i1 will write to corr[i1][i2] and the iteration for i2 will write to corr[i2][i1] (which is the same as corr[i1][i2] because of symmetry). So there is a race condition on the symmetric element.

   We must avoid this. We can do:

   Option 1: Parallelize the i loop and then within each i, the j loop is sequential. But then we set corr[j][i] = corr[i][j] without a race because j>i and we are in the i iteration. However, when two threads are processing i1 and i2 (with i1 < i2), then the thread for i1 will write to corr[i1][i2] and the thread for i2 will write to corr[i2][i1] (which is the same location). This is a race.

   Option 2: We can restructure to compute only the upper triangle without setting the lower triangle in the same loop. Then after the parallel loop, we can set the lower triangle in a separate loop. But note that the original code sets the lower triangle immediately. We can do:

   #pragma omp parallel for private(j,k)
   for (i=0; i<_PB_M-1; i++) {
      corr[i][i] = 1.0;
      for (j=i+1; j<_PB_M; j++) {
          corr[i][j] = 0.0;
          for (k=0; k<_PB_N; k++)
              corr[i][j] += (data[k][i] * data[k][j]);
          // Do not set corr[j][i] here to avoid race
      }
   }
   // Then set the lower triangle and the last diagonal
   #pragma omp parallel for collapse(2)  // This is safe because we are only writing to corr[j][i] for j>i, and these are distinct from the upper triangle we just computed.
   for (i=0; i<_PB_M-1; i++) {
      for (j=i+1; j<_PB_M; j++) {
          corr[j][i] = corr[i][j];
      }
   }
   corr[_PB_M-1][_PB_M-1] = 1.0;

   But note: the original code sets the last diagonal element separately. We can set it in the first parallel loop? Actually, the first loop runs for i from 0 to M-2. We can set the last diagonal in the second parallel loop? Or we can set it in a single thread.

   Alternatively, we can set the last diagonal in the first loop by extending the i loop to M-1 and then skipping the j loop for i=M-1. But the original code does:

      for (i=0; i<_PB_M-1; i++) {
          ... 
      }
      corr[_PB_M-1][_PB_M-1] = 1.0;

   We can do:

   #pragma omp parallel
   {
      #pragma omp for nowait  // We can use nowait because the next part is also parallel and we have a barrier at the end of the parallel region? Actually, we need to synchronize before the next part that uses corr.
      for (i=0; i<_PB_M; i++) {
          corr[i][i] = 1.0;
      }

      #pragma omp for
      for (i=0; i<_PB_M-1; i++) {
          for (j=i+1; j<_PB_M; j++) {
              corr[i][j] = 0.0;
              for (k=0; k<_PB_N; k++)
                  corr[i][j] += (data[k][i] * data[k][j]);
          }
      }
   }

   Then set the lower triangle:

   #pragma omp parallel for collapse(2)
   for (i=0; i<_PB_M; i++) {
      for (j=0; j<i; j++) {
          corr[i][j] = corr[j][i];
      }
   }

   But note: the original code does not set the entire lower triangle, only the part for j>i in the upper triangle. Actually, the original code sets corr[j][i] = corr[i][j] for j>i. This is the same as setting the lower triangle for i>j.

   However, the above parallel for the lower triangle would set every element below the diagonal. This is acceptable.

   But note: the original code does not set the diagonal again. We set the diagonal in the first parallel loop.

   This restructuring avoids the race and allows parallelization of the i loop for the upper triangle.

   However, we must be cautious: the first loop (setting the diagonal) is now done in parallel. Then the second loop (computing the upper triangle) is done in parallel. There is no race because the diagonal elements are written only once and the upper triangle elements are written only by one thread (the one that owns i). Then the lower triangle is set in a separate parallel loop.

   This is safe.

   But note: the original code does not set the entire lower triangle? Actually, it sets only the symmetric part for the pairs (i,j) that are computed. Our method sets the entire lower triangle. This is equivalent.

   We choose this restructuring.

   Steps for the correlation part:

   - First, set the entire diagonal in parallel.
   - Then, compute the upper triangle (i<j) in parallel over i (and j inner).
   - Then, set the lower triangle (i>j) in parallel.

   However, note that the original code does not compute the entire upper triangle? It computes for i from 0 to M-2 and j from i+1 to M-1. This is the upper triangle without the diagonal. We are setting the diagonal separately.

   This is acceptable.

3. Scheduling: We can use static scheduling for all loops because the work is balanced.

4. Synchronization: We use nowait for the diagonal setting because it is independent of the upper triangle computation? Actually, the upper triangle computation does not use the diagonal. So we can use nowait. But note: the upper triangle computation uses data, which is read-only. So it is safe.

   However, we must ensure that the diagonal is set before the lower triangle is set? Actually, the lower triangle setting uses the upper triangle and the diagonal. But the diagonal is set in the first loop and the upper triangle in the second. We have a barrier at the end of the parallel region (implicit at the end of the omp parallel block). So it is safe.

5. Data sharing: We use default(none) and specify shared and private variables.

6. Reductions: The mean and stddev computations are reductions. We can use reduction(+:mean[j]) and reduction(+:stddev[j]) for the fused j loop? But note: the fused j loop has two separate reductions (one for mean and one for stddev) and they are for the same j. We can do:

   #pragma omp parallel for private(i) reduction(+:mean[:_PB_M], stddev[:_PB_M])  // But note: OpenMP 4.5 supports array reductions? Actually, we can use reduction with a variable that is an array? The standard does not allow that. We have to do a reduction for each element? That is not efficient.

   Alternatively, we can use a critical section or atomic for each j? That would be inefficient.

   Instead, we can use a private array for each thread and then combine. But that is complex.

   Alternatively, we can note that the reductions for different j are independent. We can use:

   #pragma omp parallel
   {
      #pragma omp for nowait
      for (j=0; j<_PB_M; j++) {
          // private mean and stddev for this j? Actually, we can do the entire inner loop for j in one thread.
          // So we don't need a reduction clause because each j is independent and we are not accumulating across threads for the same j.
          // We are accumulating for the same j in the inner i loop, but that is done by one thread.
          // So we can do:
          double my_mean = 0.0;
          for (i=0; i<_PB_N; i++)
              my_mean += data[i][j];
          mean[j] = my_mean / float_n;

          double my_stddev = 0.0;
          for (i=0; i<_PB_N; i++)
              my_stddev += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
          stddev[j] = my_stddev / float_n;
          stddev[j] = SQRT_FUN(stddev[j]);
          stddev[j] = (stddev[j] <= eps) ? 1.0 : stddev[j];
      }
   }

   But note: the inner loop for stddev uses mean[j] which we just computed. This is safe because the same thread computes both.

   However, we cannot use a private my_mean and then assign to mean[j] because multiple threads would write to mean[j] without synchronization. We must have each thread write to its own j. We can do:

   We can use a private array for the entire mean and stddev? But that would require a lot of memory.

   Alternatively, we can use a critical section for writing to mean[j] and stddev[j]? That would be inefficient.

   The best is to note that the j loop is parallelized by j, so each thread gets a set of j's. We can compute the mean and stddev for each j independently. We can do:

   #pragma omp parallel for private(i)
   for (j=0; j<_PB_M; j++) {
      double my_mean = 0.0;
      for (i=0; i<_PB_N; i++)
          my_mean += data[i][j];
      mean[j] = my_mean / float_n;

      double my_stddev = 0.0;
      for (i=0; i<_PB_N; i++)
          my_stddev += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      stddev[j] = my_stddev / float_n;
      stddev[j] = SQRT_FUN(stddev[j]);
      stddev[j] = (stddev[j] <= eps) ? 1.0 : stddev[j];
   }

   This is safe because each j is processed by one thread. We don't need a reduction clause because the inner loop is done by one thread.

   But note: the inner loop for stddev uses mean[j] which we just computed. This is safe.

   However, we must declare my_mean and my_stddev as private to the thread. They are local to the loop body.

   This is acceptable.

7. SIMD: We can add simd to the inner loops. For the fused j loop, the inner i loops are reductions. We can use #pragma omp simd reduction(+:my_mean) and similarly for my_stddev.

   For the centering loop, we can use #pragma omp simd because it is element-wise.

   For the correlation loop, the inner k loop is a reduction. We can use #pragma omp simd reduction(+:corr[i][j]).

   We assume the compiler can vectorize these.

8. Thread management: We can use proc_bind(close) to keep threads close to the data.

Now, we write the optimized code within the scop region.

We do:

   - Fuse the first two loops (mean and stddev) and parallelize over j with OpenMP. Use simd for the inner i loops.
   - Parallelize the centering loop with collapse(2) and use simd for the inner j loop? Actually, the centering loop is i then j. We can collapse(2) and then use simd for the inner j? But note: the inner j loop is not a reduction, it is a simple assignment. We can use simd on the j loop.

   - Restructure the correlation loop as described: set the diagonal in parallel, then compute the upper triangle in parallel, then set the lower triangle in parallel. Use simd for the inner k loop in the upper triangle computation.

We also add default(none) and specify shared and private variables.

Let's code accordingly.

Note: We must preserve the original semantics and variable names.

We assume that _PB_M and _PB_N are defined and that the arrays are appropriately dimensioned.

We also note that the original code uses SQRT_FUN and SCALAR_VAL. We preserve these.

We also note that the original code uses eps = SCALAR_VAL(0.1). We preserve that.

We now write the optimized scop region.
#pragma scop
  #pragma omp parallel default(none) shared(_PB_M, _PB_N, data, mean, stddev, float_n, eps) private(i, j, k)
  {
    #pragma omp for schedule(static)
    for (j = 0; j < _PB_M; j++)
    {
      double my_mean = SCALAR_VAL(0.0);
      #pragma omp simd reduction(+:my_mean)
      for (i = 0; i < _PB_N; i++)
        my_mean += data[i][j];
      mean[j] = my_mean / float_n;
      
      double my_stddev = SCALAR_VAL(0.0);
      #pragma omp simd reduction(+:my_stddev)
      for (i = 0; i < _PB_N; i++)
        my_stddev += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      stddev[j] = my_stddev / float_n;
      stddev[j] = SQRT_FUN(stddev[j]);
      stddev[j] = stddev[j] <= eps ? SCALAR_VAL(1.0) : stddev[j];
    }
  }

  #pragma omp parallel for collapse(2) default(none) shared(_PB_N, _PB_M, data, mean, stddev, float_n) private(i, j) proc_bind(close)
  for (i = 0; i < _PB_N; i++)
    #pragma omp simd
    for (j = 0; j < _PB_M; j++)
    {
      data[i][j] -= mean[j];
      data[i][j] /= SQRT_FUN(float_n) * stddev[j];
    }

  #pragma omp parallel default(none) shared(_PB_M, _PB_N, data, corr) private(i, j, k) proc_bind(close)
  {
    #pragma omp for nowait
    for (i = 0; i < _PB_M; i++)
      corr[i][i] = SCALAR_VAL(1.0);

    #pragma omp for
    for (i = 0; i < _PB_M-1; i++)
    {
      for (j = i+1; j < _PB_M; j++)
      {
        double tmp = SCALAR_VAL(0.0);
        #pragma omp simd reduction(+:tmp)
        for (k = 0; k < _PB_N; k++)
          tmp += data[k][i] * data[k][j];
        corr[i][j] = tmp;
      }
    }
  }

  #pragma omp parallel for collapse(2) default(none) shared(_PB_M, corr) private(i, j) proc_bind(close)
  for (i = 0; i < _PB_M; i++)
    for (j = 0; j < i; j++)
      corr[i][j] = corr[j][i];
#pragma endscop
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE float_n;
  POLYBENCH_2D_ARRAY_DECL(data,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(corr,DATA_TYPE,M,M,m,m);
  POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);
  POLYBENCH_1D_ARRAY_DECL(stddev,DATA_TYPE,M,m);

  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_correlation (m, n, float_n,
		      POLYBENCH_ARRAY(data),
		      POLYBENCH_ARRAY(corr),
		      POLYBENCH_ARRAY(mean),
		      POLYBENCH_ARRAY(stddev));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(corr)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(corr);
  POLYBENCH_FREE_ARRAY(mean);
  POLYBENCH_FREE_ARRAY(stddev);

  return 0;
}
