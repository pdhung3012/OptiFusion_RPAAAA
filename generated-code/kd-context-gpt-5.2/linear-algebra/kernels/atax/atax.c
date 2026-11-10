/*******************************************************
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* atax.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "atax.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n))
{
  int i, j;
  DATA_TYPE fn;
  fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
      x[i] = 1 + (i / fn);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A[i][j] = (DATA_TYPE) ((i+j) % n) / (5*m);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("y");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, y[i]);
  }
  POLYBENCH_DUMP_END("y");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_atax(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n),
		 DATA_TYPE POLYBENCH_1D(tmp,M,m))
{
  int i, j;

#pragma scop
  /* Block over j for the A^T * tmp accumulation to improve locality.
     Parallelize over j-tiles; each thread produces a private y-tile to
     avoid atomics on y. Keep i sequential inside each tile for correctness. */
  {
    const int N = _PB_N;
    const int M = _PB_M;
    const int TJ = 256; /* cache-friendly starting point */

#pragma omp parallel default(none) shared(A, x, y, tmp) firstprivate(N, M) private(i, j)
    {
      /* Initialize y in parallel. */
#pragma omp for schedule(static)
      for (j = 0; j < N; j++)
        y[j] = 0;

      /* Compute tmp in parallel over i (independent rows). */
#pragma omp for schedule(static)
      for (i = 0; i < M; i++)
      {
        DATA_TYPE sum = SCALAR_VAL(0.0);
#pragma omp simd reduction(+:sum)
        for (j = 0; j < N; j++)
          sum += A[i][j] * x[j];
        tmp[i] = sum;
      }

      /* Accumulate y = A^T * tmp using j-tiling and per-thread private tiles. */
#pragma omp for schedule(static)
      for (int jj = 0; jj < N; jj += TJ)
      {
        const int jend = (jj + TJ < N) ? (jj + TJ) : N;

        /* Private accumulator for this tile. */
        DATA_TYPE yloc[TJ];
        const int tlen = jend - jj;

        for (int t = 0; t < tlen; t++)
          yloc[t] = SCALAR_VAL(0.0);

        for (i = 0; i < M; i++)
        {
          const DATA_TYPE ti = tmp[i];
#pragma omp simd
          for (j = jj; j < jend; j++)
            yloc[j - jj] += A[i][j] * ti;
        }

        for (j = jj; j < jend; j++)
          y[j] += yloc[j - jj];
      }
    }
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, M, N, m, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, M, m);

  /* Initialize array(s). */
  init_array (m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(x));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_atax (m, n,
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(x),
	       POLYBENCH_ARRAY(y),
	       POLYBENCH_ARRAY(tmp));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(tmp);

  return 0;
}