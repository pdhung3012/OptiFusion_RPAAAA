/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* ludcmp.c: this file is part of PolyBench/C */

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
{
  int TS = 32;
  int NT = (_PB_N + TS - 1) / TS;
  int t, ii, jj, ii_min, ii_max, i0, i1, j0, j1;
  DATA_TYPE sum;

  #pragma omp parallel private(t, ii, jj, i, j, k, w, ii_min, ii_max, i0, i1, j0, j1) \
                       shared(A, _PB_N, TS, NT)
  {
    for (t = 0; t <= 2*(NT-1); t++) {
      ii_min = (t > NT-1) ? t - (NT-1) : 0;
      ii_max = (t < NT-1) ? t : NT-1;

      #pragma omp for schedule(dynamic,1)
      for (ii = ii_min; ii <= ii_max; ii++) {
        jj = t - ii;
        i0 = ii * TS;
        i1 = i0 + TS - 1;
        if (i1 >= _PB_N) i1 = _PB_N - 1;
        j0 = jj * TS;
        j1 = j0 + TS - 1;
        if (j1 >= _PB_N) j1 = _PB_N - 1;

        for (i = i0; i <= i1; i++) {
          for (j = j0; j <= j1; j++) {
            if (j < i) {
              w = A[i][j];
              for (k = 0; k < j; k++) {
                w -= A[i][k] * A[k][j];
              }
              A[i][j] = w / A[j][j];
            } else {
              w = A[i][j];
              for (k = 0; k < i; k++) {
                w -= A[i][k] * A[k][j];
              }
              A[i][j] = w;
            }
          }
        }
      }
    }
  }

  for (i = 0; i < _PB_N; i++) {
    sum = 0.0;
    #pragma omp simd reduction(+:sum)
    for (j = 0; j < i; j++) {
      sum += A[i][j] * y[j];
    }
    y[i] = b[i] - sum;
  }

  for (i = _PB_N-1; i >= 0; i--) {
    sum = 0.0;
    #pragma omp simd reduction(+:sum)
    for (j = i+1; j < _PB_N; j++) {
      sum += A[i][j] * x[j];
    }
    x[i] = (y[i] - sum) / A[i][i];
  }
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