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
We are given a code for the PolyBench correlation kernel. We are to optimize the code within the '#pragma scop' and '#pragma endscop' markers.

The original code has several loops:

1. First loop: Computes the mean for each column j.
   for (j = 0; j < _PB_M; j++)
   {
      mean[j] = 0.0;
      for (i = 0; i < _PB_N; i++)
          mean[j] += data[i][j];
      mean[j] /= float_n;
   }

2. Second loop: Computes the standard deviation for each column j.
   for (j = 0; j < _PB_M; j++)
   {
      stddev[j] = 0.0;
      for (i = 0; i < _PB_N; i++)
          stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      stddev[j] /= float_n;
      stddev[j] = SQRT_FUN(stddev[j]);
      stddev[j] = stddev[j] <= eps ? 1.0 : stddev[j];
   }

3. Third loop: Centers and reduces the column vectors.
   for (i = 0; i < _PB_N; i++)
      for (j = 0; j < _PB_M; j++)
      {
          data[i][j] -= mean[j];
          data[i][j] /= SQRT_FUN(float_n) * stddev[j];
      }

4. Fourth loop: Computes the correlation matrix.
   for (i = 0; i < _PB_M-1; i++)
   {
      corr[i][i] = 1.0;
      for (j = i+1; j < _PB_M; j++)
      {
          corr[i][j] = 0.0;
          for (k = 0; k < _PB_N; k++)
              corr[i][j] += (data[k][i] * data[k][j]);
          corr[j][i] = corr[i][j];
      }
   }
   corr[_PB_M-1][_PB_M-1] = 1.0;

We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

Analysis:

1. Loop Tiling:
   - The first two loops (mean and stddev) are over columns and then rows. They are independent per column. We can consider tiling the inner loop (over i) to improve cache locality? But note that the inner loop is a reduction and the data is accessed by column. Since the data is stored in row-major order (as per the declaration: POLYBENCH_2D(data,N,M,n,m) -> data[i][j] with i being the row index and j the column index), accessing by column (i.e., data[i][j] for fixed j and varying i) is strided. Tiling the inner loop might not help because the stride is large. Instead, we might consider fusing the first two loops? But note that the second loop depends on the mean computed in the first loop. So they cannot be fused arbitrarily.

   - The third loop (centering) is over rows and then columns. This loop accesses data by row (contiguous) and then by column. We can tile the inner loop (over j) to improve locality? But note that the inner loop is over columns and the access is contiguous per row. Tiling the outer loop (over i) might be more beneficial? However, the outer loop is over rows and each row is contiguous. We can also consider tiling both loops to form blocks.

   - The fourth loop (correlation) is a nested loop with a reduction. The inner loop (over k) accesses two columns of the data matrix. We can tile the k loop to improve cache locality? Also, the outer loops (i and j) are over the correlation matrix. We can tile the i and j loops? But note that the correlation matrix is symmetric and we are only computing the upper triangle. Tiling might be complex.

   Given the complexity and the fact that the first two loops are column-wise and the third is row-wise, we might not get a big win from tiling without a more sophisticated approach. However, we note that the third loop (centering) is a good candidate for tiling because it is a row-wise traversal and we can form blocks of rows and columns.

   But note: the problem sizes (N and M) are not known. We must choose a tile size that is legal and beneficial. We'll choose a tile size of 64 for the centering loop? But we must be cautious because the inner loop is over columns and the tile size should be chosen to fit in cache.

   However, after careful thought, we decide to focus on the centering loop for tiling. We can tile the i (row) loop and the j (column) loop. We'll use a tile size of 64 for both dimensions? But note that the inner loop is over j and we are accessing data[i][j] which is contiguous in j. So tiling the j loop might not be necessary. Instead, we can tile the i loop to form blocks of rows. This would allow us to reuse the mean and stddev arrays for a block of rows? Actually, the mean and stddev arrays are indexed by j (column) and are the same for all rows. So we don't get reuse of mean and stddev by tiling the i loop. However, we can tile the j loop to get better cache locality for the data array? But note that the data array is being written and we are accessing it row by row. The access pattern is already contiguous per row. So tiling might not be critical.

   Alternatively, we can consider fusing the first two loops? But they are separated by a division and a square root. And the second loop depends on the entire first loop. So fusion is not possible.

   Therefore, we decide to skip tiling for now and focus on OpenMP.

2. Loop Fusion: 
   - The first two loops are over the same j loop. However, the second loop depends on the entire first loop (because it uses mean[j] for every j). So we cannot fuse them arbitrarily. We can only fuse the inner loops? But the inner loops are over i and are independent per j. However, the first loop must complete for all j before the second loop starts. So fusion is not possible.

   - The third loop is independent of the first two? Actually, it uses mean and stddev which are computed in the first two loops. So it must come after. And the fourth loop uses the centered data, so it must come after the third.

   Therefore, we cannot fuse any of these loops.

3. Loop Skewing: 
   - There are no carried dependencies that would benefit from skewing. The loops are either reductions (first two) or independent (third and fourth). So we skip.

Now, we move to OpenMP:

We will parallelize the outer loops where possible.

1. First loop (mean): 
   The outer loop is over j (columns). Each j is independent. We can parallelize the j loop.

   We note that the inner loop is a reduction. We can use a reduction clause? But note: the reduction is over the inner loop and the outer loop is parallel. We can do:

   #pragma omp parallel for private(i) reduction(+:mean) ... 

   However, note that mean is an array and we are updating one element per j. Actually, each j iteration updates only mean[j]. So we don't need a reduction over the entire array. We can simply make mean shared and have each thread update its own j. But note: the inner loop is a reduction for a single element. We can do:

   #pragma omp parallel for private(i)
   for (j=0; j<_PB_M; j++) {
      mean[j] = 0.0;
      for (i=0; i<_PB_N; i++) {
          mean[j] += data[i][j];
      }
      mean[j] /= float_n;
   }

   This is safe because each thread writes to a distinct element of mean.

2. Second loop (stddev): 
   Similarly, we can parallelize the j loop. Each j is independent.

   #pragma omp parallel for private(i)
   for (j=0; j<_PB_M; j++) {
      stddev[j] = 0.0;
      for (i=0; i<_PB_N; i++) {
          stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      }
      stddev[j] /= float_n;
      stddev[j] = SQRT_FUN(stddev[j]);
      stddev[j] = stddev[j] <= eps ? 1.0 : stddev[j];
   }

   Note: mean[j] is read and is shared. But it is only read and not written, so it is safe.

3. Third loop (centering):
   This loop is over i (rows) and then j (columns). We can parallelize the outer i loop? But note: the inner loop is over j and we are writing to data[i][j]. Since each i iteration writes to a distinct row, we can parallelize the i loop.

   #pragma omp parallel for private(j)
   for (i=0; i<_PB_N; i++) {
      for (j=0; j<_PB_M; j++) {
          data[i][j] -= mean[j];
          data[i][j] /= SQRT_FUN(float_n) * stddev[j];
      }
   }

   However, note that the inner loop is over j and the access to data is contiguous. We can also consider collapsing the two loops? But the total number of iterations is N*M, which might be too large for the OpenMP scheduler. We can use collapse(2) to expose more parallelism? But note: the loop is already over rows and then columns. Collapsing might help if the number of rows is small. We'll use collapse(2) to let the runtime decide.

   Alternatively, we can use:
   #pragma omp parallel for collapse(2) private(i,j) 
   But note: the original loop order is i then j. We can collapse the two loops? Yes, because the iterations are independent.

   However, we must be cautious: the inner loop is over j and the outer over i. Collapsing will create a single loop of size N*M. This might be acceptable.

   But note: the data array is stored in row-major order. The original loop order (i then j) is good for spatial locality. Collapsing the loops might break the contiguous access pattern? Actually, the OpenMP collapse does not change the order of iterations, it just flattens the two loops into one. The access pattern remains the same. So it is safe.

   We'll use collapse(2) to expose more parallelism.

4. Fourth loop (correlation):
   This loop is more complex. The outer loop is over i (from 0 to M-2) and then j (from i+1 to M-1). The inner loop is over k (from 0 to N-1). We are computing the upper triangle of the correlation matrix.

   We can parallelize the outer i loop? But note: the j loop depends on i. We can parallelize the i loop because each i iteration is independent? Actually, the i iteration writes to corr[i][j] and corr[j][i] for j>i. And different i iterations write to different rows and columns? Not exactly: when i=0, we write to corr[0][1..M-1] and corr[1..M-1][0]. When i=1, we write to corr[1][2..M-1] and corr[2..M-1][1]. So the writes are to distinct elements? Yes, because the correlation matrix is symmetric and we are only writing the upper triangle and then setting the lower triangle. So each (i,j) pair is written exactly once (for j>i) and then the symmetric element is set. Therefore, we can parallelize the i loop.

   However, note that the inner k loop is a reduction. We can do:

   #pragma omp parallel for private(j,k)
   for (i=0; i<_PB_M-1; i++) {
      corr[i][i] = 1.0;
      for (j=i+1; j<_PB_M; j++) {
          corr[i][j] = 0.0;
          for (k=0; k<_PB_N; k++) {
              corr[i][j] += (data[k][i] * data[k][j]);
          }
          corr[j][i] = corr[i][j];
      }
   }

   But note: the inner k loop is a reduction. We can use a reduction clause? However, the reduction is for a single element (corr[i][j]) and we are writing to a distinct element per (i,j). So we don't need a reduction clause. We can simply make corr shared and have each thread update its own (i,j) element.

   However, we must be cautious: the inner loop is over k and the data is accessed by k. We are accessing two columns of data: data[k][i] and data[k][j]. The access pattern is by column, which is strided. We might consider tiling the k loop? But we are not doing that now.

   Also, note that the number of iterations for the i loop is M-1, which might be small. We can collapse the i and j loops? But the j loop depends on i. We cannot collapse because the inner loop bounds depend on the outer loop variable.

   Alternatively, we can restructure the loop to have a single loop over the upper triangle? But that would require a more complex indexing.

   We decide to parallelize the i loop and leave the j and k loops as they are.

   Also, note the last line: corr[_PB_M-1][_PB_M-1] = 1.0; This is a single assignment and can be done outside the parallel region.

Now, we also consider scheduling and other clauses.

- We'll use static scheduling for all loops because the work is uniform.

- We'll use 'default(none)' and specify the shared and private variables.

- We note that the first two loops (mean and stddev) are over j and the inner loop is over i. We can use 'schedule(static)' and let the chunk size be determined by the runtime.

- For the centering loop, we are using collapse(2). We can use 'schedule(static)'.

- For the correlation loop, we are parallelizing the i loop. We can use 'schedule(static)'.

We also note that the first two loops are independent and can be run in parallel. But they are separated by a barrier? Actually, the second loop depends on the first. So we must run them sequentially.

We'll add OpenMP directives accordingly.

Let's write the optimized code within the scop region.

Note: We must preserve the original semantics and variable scopes.

We'll do:

#pragma scop

// First loop: mean
#pragma omp parallel for default(none) shared(_PB_M, _PB_N, data, mean, float_n) private(i,j)
for (j = 0; j < _PB_M; j++)
{
   mean[j] = SCALAR_VAL(0.0);
   for (i = 0; i < _PB_N; i++)
      mean[j] += data[i][j];
   mean[j] /= float_n;
}

// Second loop: stddev
#pragma omp parallel for default(none) shared(_PB_M, _PB_N, data, mean, stddev, float_n, eps) private(i,j)
for (j = 0; j < _PB_M; j++)
{
   stddev[j] = SCALAR_VAL(0.0);
   for (i = 0; i < _PB_N; i++)
      stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
   stddev[j] /= float_n;
   stddev[j] = SQRT_FUN(stddev[j]);
   stddev[j] = stddev[j] <= eps ? SCALAR_VAL(1.0) : stddev[j];
}

// Third loop: centering
#pragma omp parallel for collapse(2) default(none) shared(_PB_N, _PB_M, data, mean, stddev, float_n) private(i,j)
for (i = 0; i < _PB_N; i++)
   for (j = 0; j < _PB_M; j++)
   {
      data[i][j] -= mean[j];
      data[i][j] /= SQRT_FUN(float_n) * stddev[j];
   }

// Fourth loop: correlation
#pragma omp parallel for default(none) shared(_PB_M, _PB_N, data, corr) private(i,j,k)
for (i = 0; i < _PB_M-1; i++)
{
   corr[i][i] = SCALAR_VAL(1.0);
   for (j = i+1; j < _PB_M; j++)
   {
      corr[i][j] = SCALAR_VAL(0.0);
      for (k = 0; k < _PB_N; k++)
          corr[i][j] += (data[k][i] * data[k][j]);
      corr[j][i] = corr[i][j];
   }
}
corr[_PB_M-1][_PB_M-1] = SCALAR_VAL(1.0);

#pragma endscop

But note: in the first two loops, we have a private variable i. We must declare it as private. Also, in the centering loop, we have i and j as private. In the correlation loop, we have i, j, k as private.

We also note that the first two loops are over j and the inner loop is over i. We have declared j as the loop variable and i as private. This is correct.

However, in the first two loops, we have:

   for (j=0; j<_PB_M; j++) {
      ... 
      for (i=0; i<_PB_N; i++) ...

We have declared j as the loop variable and i as private. But note: the OpenMP standard requires that the loop variable of a for loop in a parallel for must be private. So j is implicitly private. We don't need to list it in the private clause? Actually, we have listed j in the private clause? We have:

   private(i,j)

But j is the loop variable and is already private. We can remove j from the private list? Actually, it is harmless to list it. But to be precise, we can remove j because it is the loop variable and is private by default.

Similarly, in the centering loop, we have:

   for (i=0; ...) 
      for (j=0; ...)

We have collapse(2) and the loop variables i and j are private by default. So we don't need to list them? But we have listed them. It is acceptable.

We'll leave as is.

Also, note that the first two loops have a reduction in the inner loop? But we are not using a reduction clause because the reduction is for a single element per j. And we are updating distinct elements. So it is safe.

We also note that the centering loop uses collapse(2). This is valid because the two loops are perfectly nested and the number of iterations is known.

Finally, we note that the correlation loop has a triangular structure. We are parallelizing the outer i loop. The inner j loop runs from i+1 to M-1. This is acceptable.

We'll output the code accordingly.

But note: the original code uses _PB_M and _PB_N. We must use the same.

Let's write the code.
#pragma scop
#pragma omp parallel for default(none) shared(_PB_M, _PB_N, data, mean, float_n) private(i)
for (j = 0; j < _PB_M; j++)
{
  mean[j] = SCALAR_VAL(0.0);
  for (i = 0; i < _PB_N; i++)
    mean[j] += data[i][j];
  mean[j] /= float_n;
}

#pragma omp parallel for default(none) shared(_PB_M, _PB_N, data, mean, stddev, float_n, eps) private(i)
for (j = 0; j < _PB_M; j++)
{
  stddev[j] = SCALAR_VAL(0.0);
  for (i = 0; i < _PB_N; i++)
    stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
  stddev[j] /= float_n;
  stddev[j] = SQRT_FUN(stddev[j]);
  stddev[j] = stddev[j] <= eps ? SCALAR_VAL(1.0) : stddev[j];
}

#pragma omp parallel for collapse(2) default(none) shared(_PB_N, _PB_M, data, mean, stddev, float_n)
for (i = 0; i < _PB_N; i++)
  for (j = 0; j < _PB_M; j++)
  {
    data[i][j] -= mean[j];
    data[i][j] /= SQRT_FUN(float_n) * stddev[j];
  }

#pragma omp parallel for default(none) shared(_PB_M, _PB_N, data, corr) private(j,k)
for (i = 0; i < _PB_M-1; i++)
{
  corr[i][i] = SCALAR_VAL(1.0);
  for (j = i+1; j < _PB_M; j++)
  {
    corr[i][j] = SCALAR_VAL(0.0);
    for (k = 0; k < _PB_N; k++)
      corr[i][j] += (data[k][i] * data[k][j]);
    corr[j][i] = corr[i][j];
  }
}
corr[_PB_M-1][_PB_M-1] = SCALAR_VAL(1.0);
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
