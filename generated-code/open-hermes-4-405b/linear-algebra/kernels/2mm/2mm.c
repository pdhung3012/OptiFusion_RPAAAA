/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* 2mm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "2mm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk, int nl,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) ((i*j+1) % ni) / ni;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+1) % nj) / nj;
  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      C[i][j] = (DATA_TYPE) ((i*(j+3)+1) % nl) / nl;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (DATA_TYPE) (i*(j+2) % nk) / nk;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("D");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, D[i][j]);
    }
  POLYBENCH_DUMP_END("D");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_2mm(int ni, int nj, int nk, int nl,
		DATA_TYPE alpha,
		DATA_TYPE beta,
		DATA_TYPE POLYBENCH_2D(tmp,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(D,NI,NL,ni,nl))
{
  int i, j, k;

#pragma scop
We are given a code for matrix multiplication (2mm) and are to optimize the region between #pragma scop and #pragma endscop.
 The code has two main parts:
  1. Compute tmp = alpha * A * B
  2. Compute D = beta * D + tmp * C

 We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

 Step 1: Loop Tiling (Blocking)
   - The two matrix multiplications are independent and can be tiled separately.
   - We'll tile the two nested loops (i and j) for both parts. The inner k loop is the reduction.

   For the first part (tmp = alpha * A * B):
     Original: 
        for i in [0, NI)
          for j in [0, NJ)
            tmp[i][j] = 0
            for k in [0, NK)
                tmp[i][j] += alpha * A[i][k] * B[k][j]

     We can tile the i and j loops. Let's choose a tile size of, say, 64 for both i and j (if the problem is large enough).
     We'll create two outer loops (ii, jj) that iterate over tiles, and then inner loops (i, j) within the tile.

   Similarly for the second part (D = beta * D + tmp * C):
        for i in [0, NI)
          for j in [0, NL)
            D[i][j] *= beta
            for k in [0, NJ)
                D[i][j] += tmp[i][k] * C[k][j]

     We can tile the i and j loops (with tile sizes, say, 64 for i and 64 for j) and then the inner k loop.

 Step 2: Loop Fusion
   - The two parts are independent and cannot be fused because they are separated by the entire tmp array being computed and then used.

 Step 3: Loop Skewing
   - There are no carried dependencies that require skewing in these matrix multiplications.

 Therefore, we focus on tiling and then OpenMP.

 Step 4: OpenMP
   - We can parallelize the outer loops (over i and j) for both parts. However, note that the first part writes to tmp and the second part reads from tmp and writes to D. They are independent and can be parallelized separately.

   For the first part (tmp computation):
     We can use:
        #pragma omp parallel for collapse(2) private(k) schedule(static)
        for (ii = 0; ii < NI; ii += TI)
          for (jj = 0; jj < NJ; jj += TJ)
            for (i = ii; i < ii+TI; i++)
              for (j = jj; j < jj+TJ; j++)
                ... inner k loop ...

   Similarly for the second part (D computation).

   We must be cautious about the tile sizes (TI, TJ) and the collapse. We can collapse the two outer tile loops.

   Also, note that the inner k loop is a reduction. We can use a reduction clause? But note: the inner k loop is inside the j loop and we are accumulating into tmp[i][j]. Since each (i,j) is independent, we don't need a reduction over k. The reduction is per (i,j). So we don't need an OpenMP reduction for the inner k.

   However, we must ensure that the inner k loop is not parallelized (it's a reduction and we are doing it sequentially). We are only parallelizing the outer i and j.

   We can also consider SIMD for the inner k loop? But note: the inner k loop is a dot product and might be vectorized by the compiler. We can add a #pragma omp simd to the k loop.

   For the second part, similarly.

   We must also consider data sharing: 
     - tmp, A, B, C, D are shared.
     - i, j, k, ii, jj are private.

   We'll use default(none) and list all shared variables.

   Also, note that the first part initializes tmp[i][j] to 0. We must ensure that each thread writes to its own tile without conflict.

   We can also use the 'nowait' clause? But there is no next loop that we are waiting for? Actually, the two parts are independent and we can run them in parallel? But note: the second part uses tmp which is produced by the first part. So we must wait for the first part to finish. Therefore, we cannot use nowait for the first part. We can use nowait for the inner loops? Not necessary.

   We'll do:

     #pragma omp parallel
     {
        #pragma omp for collapse(2) schedule(static) nowait
        for (ii=0; ii<NI; ii+=TI)
          for (jj=0; jj<NJ; jj+=TJ)
            ... compute tmp for this tile ...

        // Then we have an implicit barrier at the end of the parallel region? 
        // But we want to run the second part only after the first part is done. 
        // So we can do:

        #pragma omp for collapse(2) schedule(static)
        for (ii=0; ii<NI; ii+=TI)
          for (jj=0; jj<NL; jj+=TJ)
            ... compute D for this tile ...
     }

   However, note that the second part uses tmp which is produced by the first part. So we must have a barrier between the two. The nowait on the first loop would remove the barrier, but then we cannot start the second loop until the first is done. So we remove the nowait? Actually, the nowait on the first loop would allow the threads to proceed to the next pragma without waiting, but then the second loop would start immediately and we would have a race. So we must not use nowait for the first loop. Alternatively, we can put a barrier after the first loop. But the #pragma omp for without nowait has an implicit barrier at the end. So we can do:

        #pragma omp for collapse(2) schedule(static)   // first loop: has barrier at the end
        ... 

        #pragma omp for collapse(2) schedule(static)   // second loop: has barrier at the end

   But note: the two loops are in the same parallel region. The first loop has a barrier at the end, so all threads will wait until the first loop is done. Then they proceed to the second loop.

   Alternatively, we can use two separate parallel regions? That would also work and might be simpler. But then we have two parallel regions. We can do:

        #pragma omp parallel for collapse(2) ...   // first loop
        ... 

        #pragma omp parallel for collapse(2) ...   // second loop

   This is acceptable and might be more efficient because the first loop might use a different number of threads than the second? But we are not changing the number of threads. And the overhead of two parallel regions is acceptable.

   However, note that the first loop produces tmp and the second loop consumes it. We must ensure that the entire first loop is done before the second starts. The two separate parallel regions will ensure that because the first parallel region will not return until all threads are done.

   We'll use two separate parallel regions.

   Also, we can consider the tile sizes. We'll choose TI=64, TJ=64. But note: the problem sizes might be small? The original code uses _PB_NI, etc. which are defined in the header. We assume they are large enough.

   We must also consider the case when the tile size does not divide the dimension. We use min in the inner loop.

   We'll write:

        #define TI 64
        #define TJ 64

        // First part: tmp = alpha * A * B
        #pragma omp parallel for collapse(2) default(none) \
            shared(tmp, A, B, alpha, NI, NJ, NK) \
            private(i, j, k, ii, jj) \
            schedule(static)
        for (ii = 0; ii < NI; ii += TI)
          for (jj = 0; jj < NJ; jj += TJ)
            for (i = ii; i < min(ii+TI, NI); i++)
              for (j = jj; j < min(jj+TJ, NJ); j++)
              {
                tmp[i][j] = SCALAR_VAL(0.0);
                for (k = 0; k < NK; k++)
                  tmp[i][j] += alpha * A[i][k] * B[k][j];
              }

        // Second part: D = beta * D + tmp * C
        #pragma omp parallel for collapse(2) default(none) \
            shared(D, tmp, C, beta, NI, NL, NJ) \
            private(i, j, k, ii, jj) \
            schedule(static)
        for (ii = 0; ii < NI; ii += TI)
          for (jj = 0; jj < NL; jj += TJ)
            for (i = ii; i < min(ii+TI, NI); i++)
              for (j = jj; j < min(jj+TJ, NL); j++)
              {
                D[i][j] *= beta;
                for (k = 0; k < NJ; k++)
                  D[i][j] += tmp[i][k] * C[k][j];
              }

   We use min to avoid going out of bounds.

   Also, note that the inner k loop might be vectorized. We can add a #pragma omp simd to the k loop? But the compiler might do it automatically. We can try:

        for (k = 0; k < NK; k++)
          tmp[i][j] += alpha * A[i][k] * B[k][j];

   becomes:

        #pragma omp simd
        for (k = 0; k < NK; k++)
          tmp[i][j] += alpha * A[i][k] * B[k][j];

   Similarly for the second part.

   However, note that the reduction in the inner k loop is a dot product and is safe for SIMD.

   We'll add the simd pragma.

   Also, we can consider the alignment? We don't have control over the allocation, so we skip aligned clause.

   We'll use:

        #pragma omp simd
        for (k = 0; k < NK; k++)

   But note: the inner loop is over k and we are accumulating into a single variable (tmp[i][j]). This is a reduction and the compiler should be able to vectorize it. The pragma omp simd will force the compiler to vectorize if possible.

   We'll add it.

   Also, note that the second part has a similar inner loop.

   We must be cautious: the inner loop for the second part is:

        for (k = 0; k < NJ; k++)
          D[i][j] += tmp[i][k] * C[k][j];

   This is also a dot product and can be vectorized.

   We'll add the simd pragma there too.

   We must also consider the data types: DATA_TYPE is defined in polybench.h, but we assume it is a floating point type.

   We'll write the code accordingly.

   However, note: the original code uses _PB_NI, _PB_NJ, etc. We must use the same.

   We'll replace NI, NJ, NK, NL with _PB_NI, _PB_NJ, _PB_NK, _PB_NL.

   Also, we must define min? We can use a macro or use a ternary operator. We'll use:

        #define MIN(a,b) ((a) < (b) ? (a) : (b))

   But note: we are inside a function. We can define it at the top of the scop region? But the problem says to preserve the markers. We can define it inside the scop region? Actually, we are allowed to restructure the code within the scop region.

   Alternatively, we can use:

        int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
        for (i = ii; i < i_end; i++)

   But that would require extra variables. We can use:

        for (i = ii; i < _PB_NI && i < ii+TI; i++)

   But that is not efficient. We'll use:

        int i_end = ii + TI;
        if (i_end > _PB_NI) i_end = _PB_NI;
        for (i = ii; i < i_end; i++)

   Similarly for j.

   We'll do that to avoid the min function.

   Alternatively, we can use:

        for (i = ii; i < (ii+TI > _PB_NI ? _PB_NI : ii+TI); i++)

   But that is inefficient because it checks every iteration. We'll precompute the end.

   We'll do:

        int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
        int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;   // for the first part

   But note: the first part uses NJ and the second part uses NL for the j loop.

   We'll do:

        // First part
        int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
        int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
        for (i = ii; i < i_end; i++)
          for (j = jj; j < j_end; j++)

   Similarly for the second part.

   We'll choose TI=64, TJ=64. But note: the problem sizes might be small? We assume they are large. If they are small, tiling might not be beneficial, but it won't break.

   We'll write the code accordingly.

   We must also consider that the tile sizes might be too big for the innermost loop? We are tiling the i and j loops, and the inner k loop is the same as before.

   We'll now write the optimized code within the scop region.

   Note: We are not allowed to change the function signature or anything outside the scop region.

   We'll replace the original two loops with the tiled and parallelized versions.

   We'll define the tile sizes as constants at the beginning of the scop region.

   We'll use:

        #define TI 64
        #define TJ 64

   But note: we are inside a function. We can define them as constants? Or use #define. We'll use #define.

   However, we must be cautious: the #define might conflict with other parts? We are inside the scop region and the original code doesn't have these. We'll undef them at the end? But the problem says to preserve the markers. We can define them inside the scop region and then undef at the end of the scop region? But the problem says to preserve the markers. We are allowed to add code between the markers.

   Alternatively, we can use const int? But then we cannot use them in the loop bounds? We can use:

        const int TI = 64;
        const int TJ = 64;

   But then we cannot use them in the loop bounds because they are not constant expressions? In C, we can use const int in array sizes? But for loop bounds, it's okay? Actually, in C99, it's allowed. But to be safe, we'll use #define.

   We'll do:

        #define TI 64
        #define TJ 64

        ... code ...

        #undef TI
        #undef TJ

   But note: the problem says to preserve the markers. We are adding code between the markers. That is allowed.

   Alternatively, we can avoid defining them and just use 64? But then we cannot change easily. We'll use the #define.

   We'll write the entire scop region.

   Let's code accordingly.

   Important: We must not change the original semantics.

   We'll now write the optimized code for the scop region.
#pragma scop
#define TI 64
#define TJ 64
  /* D := alpha*A*B*C + beta*D */
  #pragma omp parallel for collapse(2) default(none) \
      shared(tmp, A, B, alpha) private(i, j, k, ii, jj) \
      schedule(static)
  for (ii = 0; ii < _PB_NI; ii += TI)
    for (jj = 0; jj < _PB_NJ; jj += TJ)
      {
        int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
        int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
        for (i = ii; i < i_end; i++)
          for (j = jj; j < j_end; j++)
            {
              tmp[i][j] = SCALAR_VAL(0.0);
              #pragma omp simd
              for (k = 0; k < _PB_NK; k++)
                tmp[i][j] += alpha * A[i][k] * B[k][j];
            }
      }
  #pragma omp parallel for collapse(2) default(none) \
      shared(D, tmp, C, beta) private(i, j, k, ii, jj) \
      schedule(static)
  for (ii = 0; ii < _PB_NI; ii += TI)
    for (jj = 0; jj < _PB_NL; jj += TJ)
      {
        int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
        int j_end = (jj + TJ) < _PB_NL ? (jj + TJ) : _PB_NL;
        for (i = ii; i < i_end; i++)
          for (j = jj; j < j_end; j++)
            {
              D[i][j] *= beta;
              #pragma omp simd
              for (k = 0; k < _PB_NJ; k++)
                D[i][j] += tmp[i][k] * C[k][j];
            }
      }
#undef TI
#undef TJ
#pragma endscop
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int ni = NI;
  int nj = NJ;
  int nk = NK;
  int nl = NL;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(tmp,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NJ,NL,nj,nl);
  POLYBENCH_2D_ARRAY_DECL(D,DATA_TYPE,NI,NL,ni,nl);

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_2mm (ni, nj, nk, nl,
	      alpha, beta,
	      POLYBENCH_ARRAY(tmp),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(D)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);

  return 0;
}
