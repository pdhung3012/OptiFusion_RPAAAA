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
  int d, ii, jj, kk, ti, i0, i1, tj, jend;
  const int BS = 32;
  int imin, imax, len;

  /* LU factorization: wavefront skewing of the i/j loops plus
     fusion of the lower/upper updates into one anti-diagonal sweep. */
  #pragma omp parallel default(none) shared(A) \
                       private(d, imin, imax, len, ti, i0, i1, ii, jj, kk, w) \
                       proc_bind(close)
  {
    for (d = 0; d < 2 * _PB_N - 1; d++) {
        imin = (d < _PB_N) ? 0 : d - _PB_N + 1;
        imax = (d < _PB_N) ? d : _PB_N - 1;
        len = imax - imin + 1;

        #pragma omp for schedule(dynamic, 1)
        for (ti = 0; ti < (len + BS - 1) / BS; ti++) {
            i0 = imin + ti * BS;
            i1 = i0 + BS;
            if (i1 > imax + 1) i1 = imax + 1;

            for (ii = i0; ii < i1; ii++) {
                jj = d - ii;
                if (jj < ii) {
                    w = A[ii][jj];
                    for (kk = 0; kk < jj; kk++)
                        w -= A[ii][kk] * A[kk][jj];
                    A[ii][jj] = w / A[jj][jj];
                } else {
                    w = A[ii][jj];
                    for (kk = 0; kk < ii; kk++)
                        w -= A[ii][kk] * A[kk][jj];
                    A[ii][jj] = w;
                }
            }
        }
    }
  }

  /* Forward substitution with tiled j-reduction and SIMD. */
  for (ii = 0; ii < _PB_N; ii++) {
      w = b[ii];
      for (tj = 0; tj < ii; tj += BS) {
          jend = tj + BS;
          if (jend > ii) jend = ii;
          #pragma omp simd reduction(-:w)
          for (jj = tj; jj < jend; jj++)
              w -= A[ii][jj] * y[jj];
      }
      y[ii] = w;
  }

  /* Back substitution with tiled j-reduction and SIMD. */
  for (ii = _PB_N - 1; ii >= 0; ii--) {
      w = y[ii];
      for (tj = ii + 1; tj < _PB_N; tj += BS) {
          jend = tj + BS;
          if (jend > _PB_N) jend = _PB_N;
          #pragma omp simd reduction(-:w)
          for (jj = tj; jj < jend; jj++)
              w -= A[ii][jj] * x[jj];
      }
      x[ii] = w / A[ii][ii];
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