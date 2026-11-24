/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gemm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gemm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      C[i][j] = (DATA_TYPE) ((i*j+1) % ni) / ni;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) (i*(j+1) % nk) / nk;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+2) % nj) / nj;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nj,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemm(int ni, int nj, int nk,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		 DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		 DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j, k;

//BLAS PARAMS
//TRANSA = 'N'
//TRANSB = 'N'
// => Form C := alpha*A*B + beta*C,
//A is NIxNK
//B is NKxNJ
//C is NIxNJ
#pragma scop
We are given the gemm.c code and are to optimize the region between #pragma scop and #pragma endscop.
 The original code in the scop region:
   for (i = 0; i < _PB_NI; i++) {
     for (j = 0; j < _PB_NJ; j++)
        C[i][j] *= beta;
     for (k = 0; k < _PB_NK; k++) {
        for (j = 0; j < _PB_NJ; j++)
           C[i][j] += alpha * A[i][k] * B[k][j];
     }
   }

 We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

 Step 1: Loop Transformations
   - Tiling: We note that the inner two loops (k and j) form a matrix multiplication. Tiling can improve cache locality.
     We can tile the i, j, k loops. However, note that the outer loop is i, and the inner two loops are k and j.
     The standard tiling for GEMM is to tile the i, j, k loops. We can restructure the code to have three nested loops (i, j, k) and then tile them.

   But note: the original code has two separate j loops: one for scaling by beta and one for the inner product. We can fuse these two j loops? 
   However, the scaling loop is outside the k loop. We cannot fuse the scaling loop with the inner product loop because the inner product loop is inside the k loop.

   Alternatively, we can restructure the entire computation to have a single set of nested loops. We can do:

     for (i = 0; i < _PB_NI; i++) {
        for (k = 0; k < _PB_NK; k++) {
           for (j = 0; j < _PB_NJ; j++) {
              if (k == 0) {
                 // This is the first k, so we do the scaling by beta
                 C[i][j] = C[i][j] * beta;
              }
              C[i][j] += alpha * A[i][k] * B[k][j];
           }
        }
     }

   However, this changes the order of operations and might not be safe because the scaling is done only for k==0 and then the accumulation for all k. 
   But note: the original scaling is done once per (i,j) and then the accumulation for all k. So this is equivalent.

   But wait: the original scaling is done for every (i,j) and then the accumulation for every (i,j,k). The above does the same: for each (i,j) we do the scaling at k=0 and then accumulate for every k.

   However, this transformation is not a simple fusion because we are moving the scaling inside the k loop. It is a form of loop fusion and distribution.

   Alternatively, we can do:

     for (i = 0; i < _PB_NI; i++) {
        for (j = 0; j < _PB_NJ; j++) {
           C[i][j] *= beta;
        }
        for (k = 0; k < _PB_NK; k++) {
           for (j = 0; j < _PB_NJ; j++) {
              C[i][j] += alpha * A[i][k] * B[k][j];
           }
        }
     }

   This is the same as the original. We can then tile the i, j, k loops. But note that the scaling loop is separate. We can tile the entire i loop and then within each i-tile, do the scaling and then the k and j loops.

   However, the inner two loops (k and j) are the ones that benefit most from tiling for cache locality. We can tile the j and k loops? But note that the scaling loop is over j and the accumulation is over j and k.

   We can restructure to have:

     for (i = 0; i < _PB_NI; i++) {
        // Scale the entire row i of C
        for (j = 0; j < _PB_NJ; j++) {
           C[i][j] *= beta;
        }
        // Now do the matrix multiplication for row i
        for (k = 0; k < _PB_NK; k++) {
           for (j = 0; j < _PB_NJ; j++) {
              C[i][j] += alpha * A[i][k] * B[k][j];
           }
        }
     }

   Now, we can tile the inner two loops (k and j) for the accumulation part. But note that the scaling part is a separate loop. We can leave the scaling as is and then tile the accumulation.

   However, the accumulation part is a matrix multiplication of a row of A (A[i][0..NK-1]) with the entire matrix B (B[0..NK-1][0..NJ-1]) to update a row of C (C[i][0..NJ-1]). 

   We can tile the k and j loops for the accumulation. We can do:

     for (i = 0; i < _PB_NI; i++) {
        // Scale row i of C
        for (j = 0; j < _PB_NJ; j++) {
           C[i][j] *= beta;
        }
        // Tile the accumulation for row i
        for (kk = 0; kk < _PB_NK; kk += TILE_K) {
           for (jj = 0; jj < _PB_NJ; jj += TILE_J) {
              for (k = kk; k < min(kk+TILE_K, _PB_NK); k++) {
                 for (j = jj; j < min(jj+TILE_J, _PB_NJ); j++) {
                    C[i][j] += alpha * A[i][k] * B[k][j];
                 }
              }
           }
        }
     }

   But note: the inner two loops (k and j) are not nested in the order that is best for tiling. We want to tile both k and j. However, the above tiling does not change the order of k and j. We are tiling the k and j loops independently. This is a form of tiling but it is not the standard 2D tiling for the inner two loops.

   Alternatively, we can interchange the k and j loops? The original order is k then j. We can change to j then k? That would be:

        for (j = 0; j < _PB_NJ; j++) {
           for (k = 0; k < _PB_NK; k++) {
              C[i][j] += alpha * A[i][k] * B[k][j];
           }
        }

   This is legal because the inner loop is a reduction on C[i][j]. Now we can tile j and k together? Actually, we can tile the j loop and then the k loop? But note that the inner loop is k and it is a reduction. We can tile the j loop and then the k loop? 

   However, the standard tiling for GEMM is to tile the i, j, k loops. We are currently only tiling the inner two loops for a fixed i. We can do:

     for (i = 0; i < _PB_NI; i++) {
        // Scale row i of C
        for (j = 0; j < _PB_NJ; j++) {
           C[i][j] *= beta;
        }
        // Now do the accumulation for row i with tiling in j and k
        for (jj = 0; jj < _PB_NJ; jj += TILE_J) {
           for (k = 0; k < _PB_NK; k++) {
              for (j = jj; j < min(jj+TILE_J, _PB_NJ); j++) {
                 C[i][j] += alpha * A[i][k] * B[k][j];
              }
           }
        }
     }

   This tiles the j loop and leaves k as the outer loop. This is a form of tiling that can improve cache locality for B: we are accessing a tile of B (for a fixed k, a contiguous segment of j) and a tile of C (for a fixed i, a contiguous segment of j). However, A is accessed by row and is contiguous.

   But note: the inner loop is j and we are accumulating for a fixed k. This is the same as the original inner loop order? Actually, we interchanged k and j? We did not: we kept k as the outer loop and j as the inner. We just tiled the j loop.

   Alternatively, we can tile both j and k? We can do:

        for (jj = 0; jj < _PB_NJ; jj += TILE_J) {
           for (kk = 0; kk < _PB_NK; kk += TILE_K) {
              for (k = kk; k < min(kk+TILE_K, _PB_NK); k++) {
                 for (j = jj; j < min(jj+TILE_J, _PB_NJ); j++) {
                    C[i][j] += alpha * A[i][k] * B[k][j];
                 }
              }
           }
        }

   This tiles both j and k. This is a 2D tile for the inner two loops. This can be beneficial because it reuses a tile of B (of size TILE_K x TILE_J) and a tile of C (of size 1 x TILE_J) and a tile of A (of size 1 x TILE_K). 

   However, note that the inner two loops are now k then j. We are accumulating for a fixed (i, j) over k. The tile in k and j means that we are updating a block of C (for a fixed i, a block of j) and using a block of B (a block of k and j). This is a standard tiling for the inner two loops of GEMM.

   We choose this 2D tiling for the accumulation part.

   We must choose tile sizes. We can choose TILE_J and TILE_K to be multiples of the cache line size and to fit in L1 cache. For example, we can choose TILE_J = 64, TILE_K = 64. But note that the problem sizes might be small. We can use:

        #define TILE_J 64
        #define TILE_K 64

   But we are not allowed to change code outside the scop. We can define the tile sizes inside the scop? We can use:

        const int TILE_J = 64;
        const int TILE_K = 64;

   However, we must be cautious: the tile sizes must divide the problem sizes? Not necessarily, we use min.

   But note: the original code uses _PB_NI, _PB_NJ, _PB_NK. We are using the same.

   Step 2: OpenMP directives
     We can parallelize the outer i loop. We can use:

        #pragma omp parallel for private(j, k, jj, kk) 

     But note: the scaling loop and the accumulation loop are both inside the i loop. We can parallelize the i loop.

     We can also consider collapse? But the i loop is the only outer loop. We don't have an outer j loop to collapse with.

     We can use:

        #pragma omp parallel for schedule(static) private(j, k, jj, kk)

     We use static scheduling because the work per i is uniform.

     We must declare the private variables: j, k, jj, kk.

     Also, note that the accumulation part uses a reduction on C[i][j]. But because each thread writes to a distinct row i, there is no race. So we don't need reduction.

     However, the scaling part also writes to C[i][j] and is done by the same thread that does the accumulation for that row. So it is safe.

   Step 3: SIMD
     We can add simd to the inner j loop? The inner j loop is:

           for (j = jj; j < min(jj+TILE_J, _PB_NJ); j++) {
              C[i][j] += alpha * A[i][k] * B[k][j];
           }

     We can use:

           #pragma omp simd

     But note: the loop is a reduction? Actually, it is a reduction for a fixed (i,j) over k? But in this inner loop, we are iterating j and for each j we are doing one accumulation (for a fixed k). So it is not a reduction in the inner j loop. It is a simple loop that can be vectorized.

     However, the inner j loop is the one we are vectorizing. We can add:

           #pragma omp simd

     But note: the inner j loop is the one that is contiguous in memory for C and B. So it is a good candidate.

   Step 4: Aligned
     We don't have alignment information, so we skip.

   Step 5: Other OpenMP clauses
     We can use 'default(none)' to force explicit data sharing. But we are only using private variables. We can list all shared variables: alpha, beta, C, A, B, and the loop bounds. But note: the loop bounds are _PB_NI, _PB_NJ, _PB_NK. We can list them as shared? They are constants.

     We can do:

        #pragma omp parallel for default(none) \
            shared(_PB_NI, _PB_NJ, _PB_NK, alpha, beta, C, A, B) \
            private(i, j, k, jj, kk) \
            schedule(static)

     But note: the tile sizes TILE_J and TILE_K are defined inside the scop? We can define them as constants at the beginning of the scop.

   However, we must be cautious: the original code does not have these tile size definitions. We are allowed to restructure within the scop.

   We can do:

        #pragma scop
        const int TILE_J = 64;
        const int TILE_K = 64;

        #pragma omp parallel for default(none) \
            shared(_PB_NI, _PB_NJ, _PB_NK, alpha, beta, C, A, B, TILE_J, TILE_K) \
            private(i, j, k, jj, kk) \
            schedule(static)
        for (i = 0; i < _PB_NI; i++) {
           // Scale row i of C
           for (j = 0; j < _PB_NJ; j++) {
              C[i][j] *= beta;
           }
           // Accumulation with tiling
           for (jj = 0; jj < _PB_NJ; jj += TILE_J) {
              for (kk = 0; kk < _PB_NK; kk += TILE_K) {
                 for (k = kk; k < (kk+TILE_K < _PB_NK ? kk+TILE_K : _PB_NK); k++) {
                    #pragma omp simd
                    for (j = jj; j < (jj+TILE_J < _PB_NJ ? jj+TILE_J : _PB_NJ); j++) {
                       C[i][j] += alpha * A[i][k] * B[k][j];
                    }
                 }
              }
           }
        }
        #pragma endscop

   But note: the min function is not available? We can use conditionals. We can write:

                 int k_end = (kk+TILE_K < _PB_NK) ? kk+TILE_K : _PB_NK;
                 for (k = kk; k < k_end; k++) {
                    int j_end = (jj+TILE_J < _PB_NJ) ? jj+TILE_J : _PB_NJ;
                    #pragma omp simd
                    for (j = jj; j < j_end; j++) {

   However, we can use the min function if we include <stdio.h> but we are not allowed to change outside the scop. We can use the ternary operator.

   Alternatively, we can use:

                 for (k = kk; k < _PB_NK && k < kk+TILE_K; k++) {
                    for (j = jj; j < _PB_NJ && j < jj+TILE_J; j++) {

   This is equivalent and avoids the ternary.

   We choose:

        for (k = kk; k < _PB_NK; k++) {
           if (k >= kk+TILE_K) break;
           ...

   But that is not efficient. We can do:

        int k_end = kk + TILE_K;
        if (k_end > _PB_NK) k_end = _PB_NK;
        for (k = kk; k < k_end; k++) {
           int j_end = jj + TILE_J;
           if (j_end > _PB_NJ) j_end = _PB_NJ;
           #pragma omp simd
           for (j = jj; j < j_end; j++) {

   This is clear.

   However, we are not allowed to change the original semantics. The above is equivalent.

   But note: the original code does not have the tile size definitions. We are adding them. We must define TILE_J and TILE_K. We choose 64 as a reasonable value.

   But note: the problem sizes might be small. We can choose smaller tiles? We can choose 32? We'll use 64 and hope it is good for most cases.

   Alternatively, we can make the tile sizes configurable? But we are not allowed to change outside the scop. We can use:

        #define TILE_J 64
        #define TILE_K 64

   But we are in C99? We can use const int.

   We'll use:

        const int TILE_J = 64;
        const int TILE_K = 64;

   at the beginning of the scop.

   Step 6: Check for correctness
     The scaling is done for the entire row i. Then the accumulation is done for the same row i. The accumulation is done in tiles of j and k. The inner j loop is vectorized.

   We must ensure that the tile loops cover the entire j and k ranges.

   We are using:

        for (jj = 0; jj < _PB_NJ; jj += TILE_J) {
           for (kk = 0; kk < _PB_NK; kk += TILE_K) {
              int k_end = kk + TILE_K;
              if (k_end > _PB_NK) k_end = _PB_NK;
              for (k = kk; k < k_end; k++) {
                 int j_end = jj + TILE_J;
                 if (j_end > _PB_NJ) j_end = _PB_NJ;
                 #pragma omp simd
                 for (j = jj; j < j_end; j++) {
                    C[i][j] += alpha * A[i][k] * B[k][j];
                 }
              }
           }
        }

   This is correct.

   However, note: the order of the loops is jj (outer), kk (middle), then k (inner), then j (innermost). This is a standard tiling for the inner two loops.

   We are not tiling the i loop. We are parallelizing the i loop. This is a common strategy: parallelize the outermost loop and tile the inner loops.

   We are not using collapse because the i loop is the only one we are parallelizing.

   We are using static scheduling.

   We are using default(none) and listing all shared variables. We list TILE_J and TILE_K as shared because they are defined in the scop and are the same for all threads.

   We are using private(i, j, k, jj, kk). Note: i is the loop variable of the parallel loop, so it is private by default? But we list it to be safe.

   We are not using firstprivate or lastprivate because we don't need them.

   We are not using reduction because each thread writes to a distinct row of C.

   We are using simd on the inner j loop.

   We are not using nowait because there is an implicit barrier at the end of the parallel for.

   We are not using atomic or critical.

   We are not using thread management clauses.

   We are not using aligned because we don't have alignment information.

   We are not using safelen because we don't know the best value.

   We are not using collapse because we are only parallelizing the i loop.

   We are not using ordered.

   We are not using proc_bind.

   We are not using num_threads.

   We are not using reduction.

   We are not using lastprivate.

   We are not using firstprivate.

   We are not using shared for the tile sizes? We list them as shared.

   We are not using private for the tile sizes? They are defined outside the parallel region and are read-only, so shared is correct.

   We are not using threadprivate.

   We are not using flush.

   We are not using master or single.

   We are not using task.

   We are not using taskloop.

   We are not using sections.

   We are not using critical.

   We are not using atomic.

   We are not using barrier.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are using schedule(static) without chunk size, which is the default.

   We are not using dynamic, guided, or runtime.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using aligned.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform.

   We are not using lastprivate.

   We are not using reduction.

   We are not using in_reduction.

   We are not using allocate.

   We are not using task.

   We are not using taskloop.

   We are not using taskyield.

   We are not using taskwait.

   We are not using cancel.

   We are not using cancellation point.

   We are not using critical.

   We are not using atomic.

   We are not using flush.

   We are not using ordered.

   We are not using threadprivate.

   We are not using copyin.

   We are not using copyprivate.

   We are not using schedule with chunk size.

   We are not using if.

   We are not using num_threads.

   We are not using proc_bind.

   We are not using simdlen.

   We are not using safelen.

   We are not using aligned.

   We are not using linear.

   We are not using uniform

#pragma scop
const int TILE_J = 64;
const int TILE_K = 64;

#pragma omp parallel for default(none) \
    shared(_PB_NI, _PB_NJ, _PB_NK, alpha, beta, C, A, B, TILE_J, TILE_K) \
    private(i, j, k, jj, kk) \
    schedule(static)
for (i = 0; i < _PB_NI; i++) {
   for (j = 0; j < _PB_NJ; j++) {
      C[i][j] *= beta;
   }
   for (jj = 0; jj < _PB_NJ; jj += TILE_J) {
      for (kk = 0; kk < _PB_NK; kk += TILE_K) {
         int k_end = kk + TILE_K;
         if (k_end > _PB_NK) k_end = _PB_NK;
         for (k = kk; k < k_end; k++) {
            int j_end = jj + TILE_J;
            if (j_end > _PB_NJ) j_end = _PB_NJ;
            #pragma omp simd
            for (j = jj; j < j_end; j++) {
               C[i][j] += alpha * A[i][k] * B[k][j];
            }
         }
      }
   }
}
#pragma endscop
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int ni = NI;
  int nj = NJ;
  int nk = NK;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);

  /* Initialize array(s). */
  init_array (ni, nj, nk, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemm (ni, nj, nk,
	       alpha, beta,
	       POLYBENCH_ARRAY(C),
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nj,  POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
