/**
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
  {
    int ii, jj;
    const int BM = 32;
    const int BN = 256;

    #pragma omp parallel default(none) shared(A,x,y,tmp) private(i,j,ii,jj)
    {
      /* y[:] = 0 */
      #pragma omp for schedule(static)
      for (j = 0; j < _PB_N; j++)
        y[j] = 0;

      /* tmp[i] = sum_j A[i][j]*x[j] */
      #pragma omp for schedule(static)
      for (ii = 0; ii < _PB_M; ii += BM)
        for (i = ii; i < _PB_M && i < ii + BM; i++) {
          DATA_TYPE sum = SCALAR_VAL(0.0);
          #pragma omp simd reduction(+:sum)
          for (j = 0; j < _PB_N; j++)
            sum += A[i][j] * x[j];
          tmp[i] = sum;
        }

      /* Ensure all tmp is computed before using it */
      #pragma omp barrier

      /* y[j] += sum_i A[i][j]*tmp[i] with j-tiling for locality; reduction via thread-local buffer */
      {
        DATA_TYPE* restrict ypriv = (DATA_TYPE*) polybench_alloc_data(_PB_N, sizeof(DATA_TYPE));
        #pragma omp for schedule(static)
        for (j = 0; j < _PB_N; j++)
          ypriv[j] = SCALAR_VAL(0.0);

        #pragma omp for schedule(static)
        for (jj = 0; jj < _PB_N; jj += BN) {
          int jend = jj + BN;
          if (jend > _PB_N) jend = _PB_N;
          for (i = 0; i < _PB_M; i++) {
            const DATA_TYPE t = tmp[i];
            #pragma omp simd
            for (j = jj; j < jend; j++)
              ypriv[j] += A[i][j] * t;
          }
        }

        #pragma omp for schedule(static)
        for (j = 0; j < _PB_N; j++)
          y[j] += ypriv[j];

        free(ypriv);
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