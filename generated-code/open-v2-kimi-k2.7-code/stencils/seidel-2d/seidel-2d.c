/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* seidel-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "seidel-2d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_seidel_2d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int t, i, j;

#pragma scop
{
  /* Block size for spatial tiling; chosen to keep a block plus halo
     in cache while exposing enough independent blocks per wavefront. */
  const int B = 64;

  const int M = n - 2;
  const int nb_i = (M + B - 1) / B;
  const int nb_j = (M + B - 1) / B;
  const int qmax = 2 * (nb_i - 1) + (nb_j - 1);

  #pragma omp parallel default(none) \
                       shared(A, n, tsteps, B, nb_i, nb_j, M, qmax) \
                       private(t)
  {
    for (t = 0; t <= tsteps - 1; t++) {
      /* Iterate block wavefronts sequentially.  Within a wavefront q,
         all blocks are independent and can be updated in parallel. */
      for (int q = 0; q <= qmax; q++) {
        int a = q - (nb_j - 1);
        int bi_min = (a > 0) ? (a + 1) / 2 : a / 2;
        if (bi_min < 0)
          bi_min = 0;
        int bi_max = (nb_i - 1 < q / 2) ? nb_i - 1 : q / 2;

        #pragma omp for schedule(dynamic, 1)
        for (int bi_idx = bi_min; bi_idx <= bi_max; bi_idx++) {
          int bj_idx = q - 2 * bi_idx;

          int bi = 1 + bi_idx * B;
          int bj = 1 + bj_idx * B;

          int ii_max = bi + B - 1;
          if (ii_max > n - 2)
            ii_max = n - 2;

          int jj_max = bj + B - 1;
          if (jj_max > n - 2)
            jj_max = n - 2;

          for (int ii = bi; ii <= ii_max; ii++)
            for (int jj = bj; jj <= jj_max; jj++)
              A[ii][jj] = (A[ii-1][jj-1] + A[ii-1][jj] + A[ii-1][jj+1]
                           + A[ii][jj-1] + A[ii][jj] + A[ii][jj+1]
                           + A[ii+1][jj-1] + A[ii+1][jj] + A[ii+1][jj+1])
                          / SCALAR_VAL(9.0);
        }
      }
    }
  }
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_seidel_2d (tsteps, n, POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}