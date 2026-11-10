/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* lu.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "lu.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

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
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_lu(int n,
	       DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j, k;

#pragma scop
  const int B = 32;
  #pragma omp parallel private(i, j, k)
  {
    for (int b = 0; b < _PB_N; b += B) {
      int b_end = (b + B < _PB_N) ? (b + B) : _PB_N;

      /* 1. Diagonal block factorization (and panel below it) */
      #pragma omp single
      {
        for (k = b; k < b_end; k++) {
          DATA_TYPE inv_diag = (DATA_TYPE)1.0 / A[k][k];
          for (i = k + 1; i < b_end; i++) {
            A[i][k] *= inv_diag;
            DATA_TYPE aik = A[i][k];
            #pragma omp simd
            for (j = k + 1; j < b_end; j++) {
              A[i][j] -= aik * A[k][j];
            }
          }
        }
      }

      /* 2. Update upper block-row (U12) and lower block-column (L21) concurrently */
      #pragma omp for schedule(static) nowait
      for (int jj = b_end; jj < _PB_N; jj += B) {
        int jj_end = (jj + B < _PB_N) ? (jj + B) : _PB_N;
        for (i = b; i < b_end; i++) {
          for (k = b; k < i; k++) {
            DATA_TYPE aik = A[i][k];
            #pragma omp simd
            for (j = jj; j < jj_end; j++) {
              A[i][j] -= aik * A[k][j];
            }
          }
        }
      }

      #pragma omp for schedule(static)
      for (i = b_end; i < _PB_N; i++) {
        for (j = b; j < b_end; j++) {
          DATA_TYPE sum = A[i][j];
          for (k = b; k < j; k++) {
            sum -= A[i][k] * A[k][j];
          }
          A[i][j] = sum / A[j][j];
        }
      }

      /* 3. Trailing submatrix update (A22 -= L21 * U12) */
      #pragma omp for schedule(static)
      for (int ii = b_end; ii < _PB_N; ii += B) {
        int ii_end = (ii + B < _PB_N) ? (ii + B) : _PB_N;
        for (int jj = b_end; jj < _PB_N; jj += B) {
          int jj_end = (jj + B < _PB_N) ? (jj + B) : _PB_N;
          for (i = ii; i < ii_end; i++) {
            for (k = b; k < b_end; k++) {
              DATA_TYPE aik = A[i][k];
              #pragma omp simd
              for (j = jj; j < jj_end; j++) {
                A[i][j] -= aik * A[k][j];
              }
            }
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

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_lu (n, POLYBENCH_ARRAY(A));

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