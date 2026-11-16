We are given the bicg.c code and are to optimize the region between #pragma scop and #pragma endscop.
 The original code in the scop region:
   for (i = 0; i < _PB_M; i++)
     s[i] = 0;
   for (i = 0; i < _PB_N; i++)
     {
       q[i] = SCALAR_VAL(0.0);
       for (j = 0; j < _PB_M; j++)
         {
           s[j] = s[j] + r[i] * A[i][j];
           q[i] = q[i] + A[i][j] * p[j];
         }
     }

 We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

 Step 1: Loop Transformations
   - Tiling: We have two nested loops. The inner loop is over j, and the outer over i. The inner loop has two statements that access A[i][j], and also s[j] and p[j]. 
     However, note that the outer loop initializes q[i] and then the inner loop updates s[j] (which is an array of size M) and q[i] (which is an array of size N). 
     The inner loop has a reduction on s[j] (across i) and a reduction on q[i] (across j). 

     We cannot tile the inner loop arbitrarily because of the reduction on s[j] (which is accumulated over i). However, we can tile the outer loop (i) and the inner loop (j) to improve locality for A[i][j] and p[j]. 
     But note: the inner loop writes to s[j] for every i, so we cannot change the order of i without breaking the reduction. 

     Alternatively, we can try to tile the j loop? But the j loop is the inner loop and we are accumulating s[j] over i. We cannot tile the j loop without breaking the reduction.

     However, observe that the inner loop has two independent reductions: one for s[j] (which is a reduction over i) and one for q[i] (which is a reduction over j). 
     We can try to interchange the loops? But note: the outer loop is over i and the inner over j. The access to A[i][j] is row-major. 

     Let's consider the two statements separately:

        s[j] = s[j] + r[i] * A[i][j];   // This is: for each j, we are accumulating over i.

        q[i] = q[i] + A[i][j] * p[j];   // This is: for each i, we are accumulating over j.

     We can split the loop into two separate loops? But note that the original code does both in the same inner loop. 

     However, we can try loop fusion? Actually, the two statements are already fused. 

     Alternatively, we can try to tile the i loop? But the i loop is the outer loop and we are accumulating s[j] over i. We cannot tile the i loop without breaking the reduction for s[j] because the reduction is over the entire i range.

     But note: the reduction for s[j] is over i, and the reduction for q[i] is over j. We can try to tile the j loop? The j loop is the inner loop and we are accumulating q[i] over j. We can tile the j loop without breaking the reduction for q[i] because the reduction is over j. However, the s[j] update is also in the same inner loop and we are accumulating over i. 

     We can try to tile the j loop and then the i loop? But the i loop is the outer loop. 

     Alternatively, we can try to interchange the loops? 

        Original: i (outer) -> j (inner)

        After interchange: j (outer) -> i (inner)

        Then the inner loop would be over i. Then we can tile the j loop? 

        But note: the reduction for s[j] becomes: for a fixed j, we are accumulating over i -> that is the same as before. 
        The reduction for q[i] becomes: for a fixed i, we are accumulating over j -> but now the inner loop is over i, so we cannot accumulate over j in the inner loop? 

        Actually, if we interchange, we get:

          for (j=0; j<_PB_M; j++) {
            for (i=0; i<_PB_N; i++) {
              s[j] = s[j] + r[i] * A[i][j];
              q[i] = q[i] + A[i][j] * p[j];
            }
          }

        Now, the inner loop is over i. The update for s[j] is still accumulating over i (which is now the inner loop) -> that is fine. 
        But the update for q[i] is accumulating over j? But note: the inner loop is over i, and j is fixed in the outer loop. So we are accumulating q[i] for a fixed j? That is not the same. 

        Actually, the original reduction for q[i] is over j. In the interchanged version, we are accumulating q[i] for each j in the outer loop? That would be:

          q[i] = ... for j=0, then j=1, ... so we are accumulating over j? But note: the inner loop is over i, so for a fixed j, we are updating q[i] for every i. Then we move to the next j and update q[i] again? 

        This is the same as the original: the original also updated q[i] for every j. So the interchanged version is semantically equivalent? 

        However, the access to A[i][j] in the interchanged version is column-major, which is bad for locality. 

     Therefore, we avoid interchange.

     We can try to tile the inner j loop? We can break the j loop into tiles of size TB. Then we would have:

        for (i=0; i<_PB_N; i++) {
          q[i] = 0;
          for (jj=0; jj<_PB_M; jj+=TB) {
            for (j=jj; j<min(jj+TB, _PB_M); j++) {
              s[j] = s[j] + r[i] * A[i][j];
              q[i] = q[i] + A[i][j] * p[j];
            }
          }
        }

     This tiling of j does not break the reduction for q[i] because we are still accumulating over j. And for s[j], we are still accumulating over i. 

     But note: the inner loop is now over j in tiles. This might improve locality for p[j] and s[j] because we are accessing contiguous j in a tile. Also, A[i][j] is accessed row-wise, so contiguous in j. 

     However, the reduction for s[j] is over i, and we are not changing the order of i. So it is safe.

     But note: the initialization of s[i] is done in a separate loop. We are not changing that.

     So we can do:

        // Initialize s
        for (i = 0; i < _PB_M; i++)
          s[i] = 0;

        // Then the main loop
        for (i = 0; i < _PB_N; i++) {
          q[i] = SCALAR_VAL(0.0);
          for (jj = 0; jj < _PB_M; jj += TB) {
            int j_end = min(jj+TB, _PB_M);
            for (j = jj; j < j_end; j++) {
              s[j] = s[j] + r[i] * A[i][j];
              q[i] = q[i] + A[i][j] * p[j];
            }
          }
        }

     We choose a tile size TB. Since we are targeting cache, we can choose TB to be a multiple of the cache line size and also to fit the working set. 
     The working set in the inner tile: 
        A[i][jj:j_end] -> TB elements
        p[jj:j_end] -> TB elements
        s[jj:j_end] -> TB elements (but note s is both read and written, so we need to consider it as part of the working set)
        r[i] is a scalar for the outer i, and q[i] is a scalar for the outer i.

     We can choose TB = 256? But note that the inner loop is over j, and we are accumulating q[i] which is a scalar. The main issue is that we are accessing A[i][j] and p[j] and s[j] in the inner loop. 

     However, note that the inner loop is over j, and we are accessing contiguous j. So without tiling, we already have good spatial locality. Tiling might not help much? 

     But observe: without tiling, the inner loop is over j from 0 to M-1. The entire j loop is done for one i. Then we move to the next i. The next i will access the same j indices? But note that A is accessed by row, so the next i will access the next row. The p and s arrays are accessed by j, so they are reused for every i. 

     Therefore, without tiling, the entire p and s arrays are accessed for every i. If M is large, p and s might not fit in cache. Tiling the j loop would break the j loop into chunks, so that we only need a tile of p and s in cache at a time. 

     So tiling the j loop is beneficial for reducing the working set of p and s.

     We choose TB = 256 as a starting point? But note that the problem size might be small. We can use a smaller tile, say 64 or 128. However, the code uses _PB_M and _PB_N which are defined by the benchmark. We can use a fixed tile size of 256, but we must handle the boundary.

     Alternatively, we can use a tile size that is a multiple of the cache line (e.g., 64 bytes -> 16 doubles). But we are using DATA_TYPE which might be float or double. We'll assume double for safety. Then 16 doubles = 128 bytes. We can choose TB=256 (which is 16 cache lines) or 128 (8 cache lines). 

     However, the benchmark might run with small M. We'll use TB=256 and use min(jj+TB, _PB_M) for the upper bound.

 Step 2: OpenMP Directives
   - We can parallelize the outer i loop? The outer i loop runs from 0 to _PB_N-1. Each iteration i writes to q[i] and reads r[i] and the entire row A[i][0..M-1]. Also, it updates s[0..M-1] (which is shared). 

     The update to s[j] is a reduction: s[j] += ... and it is accumulated over i. So we have a reduction on s. 

     We can use:

        #pragma omp parallel for reduction(+:s[:_PB_M]) ... 

     But note: the reduction on an array of unknown size at compile time? We can use:

        #pragma omp parallel for reduction(+:s[0:_PB_M])

     However, the OpenMP standard allows array reduction only in OpenMP 5.0 and later. Alternatively, we can use a critical section or atomic, but that would be inefficient.

     Alternatively, we can restructure: since the reduction for s is over i, and we are parallelizing i, we can use a private array for s and then do a reduction at the end. But note that s is also read in the update? Actually, we are updating s[j] for each i. 

     We can do:

        Create a private array for s for each thread, then after the parallel loop, combine them.

     But note: the original code initializes s to zero. Then we are accumulating. We can do:

        double *s_private = (double*) calloc(_PB_M, sizeof(double));
        #pragma omp parallel
        {
          double *s_local = (double*) calloc(_PB_M, sizeof(double)); // or use stack if small
          #pragma omp for nowait
          for (i=0; i<_PB_N; i++) {
            // ... use s_local instead of s
          }
          #pragma omp critical
          {
            for (j=0; j<_PB_M; j++)
              s[j] += s_local[j];
          }
          free(s_local);
        }
        free(s_private);

     But this is expensive because of the allocation and the critical section.

     Alternatively, we can use an array reduction if the compiler supports it. Since the problem says we can use reduction, and we are targeting a modern compiler, we can use:

        #pragma omp parallel for reduction(+:s[:_PB_M])

     But note: the syntax for array reduction in OpenMP is:

        reduction(+:s[0:_PB_M])

     However, the benchmark uses _PB_M which is a macro. We can use:

        #pragma omp parallel for reduction(+:s[0:_PB_M])

     But note: the original code initializes s to zero. The reduction clause will initialize the private s to zero and then combine by summing. So it is safe.

     Also, we have the q array: each thread writes to its own q[i] (since i is the loop index). So q is private to the thread? Actually, no: the loop is over i, and each thread gets a set of i's. So q[i] is written only by the thread that handles that i. So we don't need a reduction for q. But we must declare q as shared? Actually, we are writing to q[i] and each i is unique. So we can use:

        #pragma omp parallel for reduction(+:s[0:_PB_M]) private(j,jj,j_end) 

     But note: the inner loops use j, jj, j_end. We must declare them private.

     Alternatively, we can use:

        #pragma omp parallel for reduction(+:s[0:_PB_M]) firstprivate(_PB_M) private(j,jj,j_end)

     But _PB_M is a macro, so it is a constant. We don't need to pass it.

     However, the tile size TB is a constant we choose. We can define it as a constant.

     Also, note that the initialization of q[i] is inside the loop. We can do:

        q[i] = 0.0;   // This is done by each thread for its own i.

     So we don't need to do anything special for q.

     But note: the original code initializes s to zero. We are using a reduction on s, so the reduction will take care of initializing the private copy to zero and then accumulating. However, the original s is already initialized to zero. The reduction will add the private copies to the original s (which is zero) so it is correct.

     However, the reduction(+:s[0:_PB_M]) will create a private copy of the entire s array for each thread, which is of size _PB_M. This might be expensive in memory if we have many threads and _PB_M is large.

     Alternatively, we can avoid the array reduction by restructuring the loop. We can split the loop into two: one for the s update and one for the q update? But note they are fused. 

     We can do:

        // Initialize s to zero (as in the original)
        // Then parallelize the i loop without reduction, but then we have a race on s.

     So we stick with the array reduction.

     We also note that the inner j loop is tiled. We can add a collapse? The outer loop is i, and we are tiling j. We cannot collapse because the j loop is inside and we have a tile loop. 

     We can try to use schedule(static) to distribute the i iterations.

     Also, we can consider SIMD: the inner j loop is contiguous and we are doing:

        s[j] = s[j] + r[i] * A[i][j];
        q[i] = q[i] + A[i][j] * p[j];

     We can use #pragma omp simd on the inner j loop? But note: the inner j loop is the one we are tiling. We can put the simd on the innermost loop.

     We can do:

        for (j = jj; j < j_end; j++) {
          s[j] = s[j] + r[i] * A[i][j];
          q[i] = q[i] + A[i][j] * p[j];
        }

     We can add:

        #pragma omp simd

     But note: the inner loop is over j, and we are writing to s[j] and reading s[j]. This is a reduction. We can use:

        #pragma omp simd reduction(+:s[j])   // But this is not standard.

     Alternatively, we can rely on the compiler to vectorize the inner loop. We can use:

        #pragma omp simd

     without reduction because the reduction is over i, and the inner loop is over j. The inner loop does not have a reduction over j? Actually, the inner loop does not have a reduction over j for s[j]? It is updating s[j] for a fixed j? No, j is the loop index. So for a fixed j, we are updating s[j] only once per i? Actually, no: the inner loop is over j, so for a fixed i, we update s[j] for every j. And we are accumulating over i. So the inner loop does not have a reduction over j. 

     Therefore, the inner loop is a simple loop without reduction. We can vectorize it.

     But note: the update to q[i] is a reduction over j. The inner loop is accumulating q[i] for a fixed i. So we have a reduction on q[i] over j. We can use:

        #pragma omp simd reduction(+:q[i])

     However, the OpenMP simd reduction clause for a scalar is supported. But note: q[i] is a scalar for the inner loop? Yes, because i is fixed in the outer loop. 

     So we can do:

        #pragma omp simd reduction(+:q_i) 
        for (j = jj; j < j_end; j++) {
          s[j] = s[j] + r_i * A[i][j];
          q_i = q_i + A[i][j] * p[j];
        }
        q[i] = q_i;   // But note: we are accumulating q[i] over the entire j loop, so we cannot break it into tiles without accumulating.

     Actually, we are already accumulating q[i] over the entire j loop. We can do:

        double q_i = 0.0;
        for (jj=0; jj<_PB_M; jj+=TB) {
          int j_end = min(jj+TB, _PB_M);
          #pragma omp simd
          for (j=jj; j<j_end; j++) {
            s[j] = s[j] + r[i] * A[i][j];
            q_i = q_i + A[i][j] * p[j];
          }
        }
        q[i] = q_i;

     This avoids the reduction on q[i] in the inner loop and instead uses a scalar that is accumulated over the entire j loop. Then we assign to q[i] at the end. This is safe and allows vectorization without reduction.

     But note: the original code did the accumulation in the inner loop. This is equivalent.

     So we restructure the inner part:

        double q_i = 0.0;
        for (jj=0; jj<_PB_M; jj+=TB) {
          int j_end = min(jj+TB, _PB_M);
          #pragma omp simd
          for (j=jj; j<j_end; j++) {
            s[j] = s[j] + r[i] * A[i][j];
            q_i = q_i + A[i][j] * p[j];
          }
        }
        q[i] = q_i;

     This also avoids the need for a reduction on q[i] in the inner loop.

     Now, the update to s[j] is still a reduction over i, but we are handling that with the outer OpenMP reduction.

     We can also use vectorization for the s update? The s update is:

        s[j] = s[j] + r[i] * A[i][j];

     This is a vector load of s[j], a vector load of A[i][j], a vector load of r[i] (broadcasted), then a multiply and add, then store s[j]. 

     The compiler should be able to vectorize this.

     We'll add #pragma omp simd to the inner j loop.

 Step 3: Code Structure

   We define a tile size, say TB = 256.

   We restructure the scop region as:

      // Initialize s
      for (i = 0; i < _PB_M; i++)
        s[i] = 0;

      // Main loop
      #pragma omp parallel for reduction(+:s[0:_PB_M]) private(j,jj,j_end,q_i) schedule(static)
      for (i = 0; i < _PB_N; i++) {
        double q_i = 0.0;
        for (jj = 0; jj < _PB_M; jj += TB) {
          j_end = (jj + TB < _PB_M) ? (jj + TB) : _PB_M;
          #pragma omp simd
          for (j = jj; j < j_end; j++) {
            s[j] = s[j] + r[i] * A[i][j];
            q_i = q_i + A[i][j] * p[j];
          }
        }
        q[i] = q_i;
      }

   We must be cautious: the array reduction on s[0:_PB_M] is supported in OpenMP 5.0 and later. If the compiler does not support it, we might need an alternative. But the problem says we can use reduction.

   Also, note that the original code uses SCALAR_VAL(0.0) for q[i]. We use 0.0.

   We assume that the compiler supports array sections in reduction.

   We also note that the tile size TB is a constant. We can define it as:

        #define TB 256

   But we must be cautious: if _PB_M is less than TB, then the inner loop runs from 0 to _PB_M.

   We use:

        j_end = (jj + TB < _PB_M) ? (jj + TB) : _PB_M;

   Alternatively, we can use:

        j_end = jj + TB;
        if (j_end > _PB_M) j_end = _PB_M;

   But the ternary operator is fine.

   We also note that the inner loop is vectorized. We assume the compiler can vectorize the inner loop.

   We also add schedule(static) to the parallel for.

   We also note that the variables j, jj, j_end, and q_i are private to each thread.

   We do not change the initialization of s.

   We do not change the function signature or anything outside the scop.

   We preserve the #pragma scop and #pragma endscop.

   We output the entire code with the scop region replaced.

   However, note: the original code uses _PB_M and _PB_N. We use the same.

   We also note that the original code uses DATA_TYPE. We use double for the scalar q_i? But the arrays are DATA_TYPE. We should use DATA_TYPE for q_i.

   We change:

        DATA_TYPE q_i = SCALAR_VAL(0.0);

   But note: SCALAR_VAL(0.0) is used in the original. We can use:

        DATA_TYPE q_i = SCALAR_VAL(0.0);

   However, the original uses SCALAR_VAL(0.0) for q[i]. We do the same.

   Also, note that the original code uses:

        q[i] = SCALAR_VAL(0.0);

   We are replacing that with:

        DATA_TYPE q_i = SCALAR_VAL(0.0);
        ... then assign q[i] = q_i;

   This is equivalent.

   We also note that the inner loop uses r[i] and p[j]. We are not changing that.

   We also note that the array A is accessed as A[i][j]. We are not changing that.

   We also note that the tile size TB is fixed. We choose 256. We can make it a define at the top of the scop region? But the problem says we can adjust associated C statements. We can define it inside the scop region.

   We do:

        #define TB 256

   But note: we must not leave the define visible outside the scop. We can undef it at the end? Or we can use a constant. We'll use a constant and hope the compiler optimizes.

   Alternatively, we can use:

        const int TB = 256;

   inside the scop region, before the parallel loop.

   We do:

        const int TB = 256;

   Then use it.

   We also note that the min function: we can use a ternary operator.

   We write:

        int j_end = (jj + TB < _PB_M) ? (jj + TB) : _PB_M;

   But note: _PB_M is a macro, so it is a constant.

   We are ready to write the code.

   However, note: the original code uses:

        for (i = 0; i < _PB_M; i++)
          s[i] = 0;

   We leave that as is.

   Then the parallel loop.

   We also note: the original code uses:

        for (i = 0; i < _PB_N; i++)

   We use the same.

   We also note: the inner loop variable j is an integer. We declare it inside the scop region? Actually, the original code declares i and j at the beginning of the function. We are allowed to adjust associated C statements. We can redeclare j inside the scop region? But the original j is declared at the top of the function. We can use a new variable? We are using j in the inner loop. We can use the same j? But note: the outer loop uses i, and we are using j in the inner loop. We are also using jj and j_end. We declare them inside the scop region? We can declare them at the beginning of the scop region? Or inside the parallel region? We declare them as private in the OpenMP directive, so we can declare them inside the parallel loop.

   We do:

        #pragma omp parallel for ... private(j,jj,j_end,q_i)

   Then inside the loop, we declare:

        DATA_TYPE q_i = SCALAR_VAL(0.0);
        for (jj=0; jj<_PB_M; jj+=TB) {
          int j_end = (jj+TB < _PB_M) ? (jj+TB) : _PB_M;
          #pragma omp simd
          for (j=jj; j<j_end; j++) {
            ...
          }
        }

   We also note: the variable i is the loop variable of the parallel for, so it is private.

   We are done.

   We output the entire code with the scop region replaced.

   Note: We must preserve the #pragma scop and #pragma endscop.

   We write the code accordingly.

   However, note: the original code uses _PB_M and _PB_N. We use the same.

   We also note: the original code uses SCALAR_VAL(0.0). We use the same.

   We also note: the original code uses DATA_TYPE. We use the same.

   We also note: the original code uses POLYBENCH_2D and POLYBENCH_1D. We do not change that.

   We only change the code between #pragma scop and #pragma endscop.

   We write the optimized scop region.

   Let's code it.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "bicg.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		 DATA_TYPE POLYBENCH_1D(r,N,n),
		 DATA_TYPE POLYBENCH_1D(p,M,m))
{
  int i, j;

  for (i = 0; i < m; i++)
    p[i] = (DATA_TYPE)(i % m) / m;
  for (i = 0; i < n; i++) {
    r[i] = (DATA_TYPE)(i % n) / n;
    for (j = 0; j < m; j++)
      A[i][j] = (DATA_TYPE) (i*(j+1) % n)/n;
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_1D(s,M,m),
		 DATA_TYPE POLYBENCH_1D(q,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("s");
  for (i = 0; i < m; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, s[i]);
  }
  POLYBENCH_DUMP_END("s");
  POLYBENCH_DUMP_BEGIN("q");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, q[i]);
  }
  POLYBENCH_DUMP_END("q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_bicg(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		 DATA_TYPE POLYBENCH_1D(s,M,m),
		 DATA_TYPE POLYBENCH_1D(q,N,n),
		 DATA_TYPE POLYBENCH_1D(p,M,m),
		 DATA_TYPE POLYBENCH_1D(r,N,n))
{
  int i, j;

#pragma scop
  for (i = 0; i < _PB_M; i++)
    s[i] = 0;
  
  const int TB = 256;
  #pragma omp parallel for reduction(+:s[0:_PB_M]) private(j,jj,j_end) schedule(static)
  for (i = 0; i < _PB_N; i++) {
    DATA_TYPE q_i = SCALAR_VAL(0.0);
    for (jj = 0; jj < _PB_M; jj += TB) {
      int j_end = (jj + TB < _PB_M) ? (jj + TB) : _PB_M;
      #pragma omp simd
      for (j = jj; j < j_end; j++) {
        s[j] = s[j] + r[i] * A[i][j];
        q_i = q_i + A[i][j] * p[j];
      }
    }
    q[i] = q_i;
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, M, n, m);
  POLYBENCH_1D_ARRAY_DECL(s, DATA_TYPE, M, m);
  POLYBENCH_1D_ARRAY_DECL(q, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(p, DATA_TYPE, M, m);
  POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(r),
	      POLYBENCH_ARRAY(p));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_bicg (m, n,
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(s),
	       POLYBENCH_ARRAY(q),
	       POLYBENCH_ARRAY(p),
	       POLYBENCH_ARRAY(r));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(s), POLYBENCH_ARRAY(q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(s);
  POLYBENCH_FREE_ARRAY(q);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(r);

  return 0;
}