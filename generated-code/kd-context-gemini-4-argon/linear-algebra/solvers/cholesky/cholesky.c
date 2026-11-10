/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* cholesky.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "cholesky.h"


/* Array initialization. */
static
void init_array(int n,
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
    for (j = 0; j <= i; j++) {
    if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
  }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_cholesky(int n,
		     DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j, k;


#pragma scop
  const int B = 32;
  #pragma omp parallel private(i, j, k)
  {
    for (int jb = 0; jb < _PB_N; jb += B) {
      int j_end = (jb + B < _PB_N) ? (jb + B) : _PB_N;

      /* 1. Diagonal block factorization (left-looking with column-wise trailing rows parallelization) */
      for (j = jb; j < j_end; j++) {
        #pragma omp single
        {
          DATA_TYPE diag = A[j][j];
          for (k = 0; k < j; k++) {
            diag -= A[j][k] * A[j][k];
          }
          A[j][j] = SQRT_FUN(diag);
        }

        DATA_TYPE inv_diag = (DATA_TYPE)1.0 / A[j][j];
        #pragma omp for schedule(static)
        for (i = j + 1; i < j_end; i++) {
          DATA_TYPE sum = A[i][j];
          for (k = 0; k < j; k++) {
            sum -= A[i][k] * A[j][k];
          }
          A[i][j] = sum * inv_diag;
        }
      }

      /* 2. Panel update for rows i >= j_end across columns jb <= j < j_end */
      #pragma omp for schedule(static)
      for (int ib = j_end; ib < _PB_N; ib += B) {
        int i_end = (ib + B < _PB_N) ? (ib + B) : _PB_N;

        /* Tile over k < jb to maximize cache reuse of panels */
        for (int kb = 0; kb < jb; kb += B) {
          int k_end = (kb + B < jb) ? (kb + B) : jb;
          for (i = ib; i < i_end; i++) {
            for (j = jb; j < j_end; j++) {
              DATA_TYPE sum = A[i][j];
              for (k = kb; k < k_end; k++) {
                sum -= A[i][k] * A[j][k];
              }
              A[i][j] = sum;
            }
          }
        }

        /* Finish remaining k in [jb, j) and divide by diagonal */
        for (i = ib; i < i_end; i++) {
          for (j = jb; j < j_end; j++) {
            DATA_TYPE sum = A[i][j];
            for (k = jb; k < j; k++) {
              sum -= A[i][k] * A[j][k];
            }
            A[i][j] = sum / A[j][j];
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
  kernel_cholesky (n, POLYBENCH_ARRAY(A));

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