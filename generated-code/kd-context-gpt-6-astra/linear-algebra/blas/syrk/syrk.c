/*************************************************************
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* syrk.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "syrk.h"


/* Array initialization. */
static
void init_array(int n, int m,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(A,N,M,n,m))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      A[i][j] = (DATA_TYPE) ((i*j+1)%n) / n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      C[i][j] = (DATA_TYPE) ((i*j+2)%m) / m;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_syrk(int n, int m,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m))
{
  int i, j, k;

//BLAS PARAMS
//TRANS = 'N'
//UPLO  = 'L'
// =>  Form  C := alpha*A*A**T + beta*C.
//A is NxM
//C is NxN
#pragma scop
  {
    /* Triangular-domain tiling for locality and better parallel granularity. */
    const int TI = 32;
    const int TJ = 32;
    const int TK = 32;

    int ii, jj, kk;

    #pragma omp parallel for schedule(static) private(ii,jj,kk,i,j,k) shared(C,A,alpha,beta)
    for (ii = 0; ii < _PB_N; ii += TI) {
      int i_end = ii + TI;
      if (i_end > _PB_N) i_end = _PB_N;

      for (jj = 0; jj <= ii; jj += TJ) {
        int j_end = jj + TJ;
        if (j_end > _PB_N) j_end = _PB_N;

        /* Scale C by beta for this (i,j) tile (lower triangle only). */
        for (i = ii; i < i_end; ++i) {
          int jlim = i + 1;
          int jmax = j_end < jlim ? j_end : jlim;
          for (j = jj; j < jmax; ++j) {
            C[i][j] *= beta;
          }
        }

        /* Block reduction over k to improve reuse of A blocks. */
        for (kk = 0; kk < _PB_M; kk += TK) {
          int k_end = kk + TK;
          if (k_end > _PB_M) k_end = _PB_M;

          for (i = ii; i < i_end; ++i) {
            int jlim = i + 1;
            int jmax = j_end < jlim ? j_end : jlim;
            for (k = kk; k < k_end; ++k) {
              DATA_TYPE aik = A[i][k];
              #pragma omp simd
              for (j = jj; j < jmax; ++j) {
                C[i][j] += alpha * aik * A[j][k];
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
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,N,M,n,m);

  /* Initialize array(s). */
  init_array (n, m, &alpha, &beta, POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_syrk (n, m, alpha, beta, POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}