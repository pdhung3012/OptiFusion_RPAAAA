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
  {
    const int bs = 64;
#pragma omp parallel private(i, j, k)
    {
      int kb;
      for (kb = 0; kb < _PB_N; kb += bs) {
        int ke = (kb + bs < _PB_N) ? kb + bs : _PB_N;

        /* Factor the diagonal block (sequential, small). */
#pragma omp single
        {
          for (k = kb; k < ke; k++) {
            for (i = k + 1; i < ke; i++) {
              A[i][k] /= A[k][k];
              DATA_TYPE aik = A[i][k];
              for (j = k + 1; j < ke; j++)
                A[i][j] -= aik * A[k][j];
            }
          }
        }

        /* Row panel (U part): rows in block, columns >= ke. */
#pragma omp for schedule(static) nowait
        for (int jb = ke; jb < _PB_N; jb += bs) {
          int je = (jb + bs < _PB_N) ? jb + bs : _PB_N;
          for (int ii = kb + 1; ii < ke; ii++) {
            for (int kk = kb; kk < ii; kk++) {
              DATA_TYPE aik = A[ii][kk];
              for (int jj = jb; jj < je; jj++)
                A[ii][jj] -= aik * A[kk][jj];
            }
          }
        }

        /* Column panel (L part): rows >= ke, columns in block. */
#pragma omp for schedule(static)
        for (i = ke; i < _PB_N; i++) {
          for (int kk = kb; kk < ke; kk++) {
            A[i][kk] /= A[kk][kk];
            DATA_TYPE aik = A[i][kk];
            for (int jj = kk + 1; jj < ke; jj++)
              A[i][jj] -= aik * A[kk][jj];
          }
        }

        /* Trailing submatrix update, tiled. */
#pragma omp for collapse(2) schedule(static)
        for (int ib = ke; ib < _PB_N; ib += bs) {
          for (int jb = ke; jb < _PB_N; jb += bs) {
            int ie = (ib + bs < _PB_N) ? ib + bs : _PB_N;
            int je = (jb + bs < _PB_N) ? jb + bs : _PB_N;
            for (int ii = ib; ii < ie; ii++) {
              for (int kk = kb; kk < ke; kk++) {
                DATA_TYPE aik = A[ii][kk];
                for (int jj = jb; jj < je; jj++)
                  A[ii][jj] -= aik * A[kk][jj];
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