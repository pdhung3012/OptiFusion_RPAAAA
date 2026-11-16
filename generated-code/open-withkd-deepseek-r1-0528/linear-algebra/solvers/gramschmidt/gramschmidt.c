/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gramschmidt.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gramschmidt.h"


/* Array initialization. */
static
void init_array(int m, int n,
		DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      A[i][j] = (((DATA_TYPE) ((i*j) % m) / m )*100) + 10;
      Q[i][j] = 0.0;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      R[i][j] = 0.0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("R");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, R[i][j]);
    }
  POLYBENCH_DUMP_END("R");

  POLYBENCH_DUMP_BEGIN("Q");
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Q[i][j]);
    }
  POLYBENCH_DUMP_END("Q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* QR Decomposition with Modified Gram Schmidt:
 http://www.inf.ethz.ch/personal/gander/ */
static
void kernel_gramschmidt(int m, int n,
			DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
			DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
			DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j, k;

  DATA_TYPE nrm;

#pragma scop
  const int TI = 64;
  const int TJ = 64;
  int ii, jj, i, j_inner, ii2;
  for (k = 0; k < _PB_N; k++) {
    nrm = SCALAR_VAL(0.0);
    
    #pragma omp parallel for reduction(+:nrm) schedule(static) private(ii, i)
    for (ii = 0; ii < _PB_M; ii += TI) {
      for (i = ii; i < (ii + TI < _PB_M ? ii + TI : _PB_M); i++) {
        nrm += A[i][k] * A[i][k];
      }
    }
    
    R[k][k] = SQRT_FUN(nrm);
    
    #pragma omp parallel for schedule(static) private(ii, i)
    for (ii = 0; ii < _PB_M; ii += TI) {
      for (i = ii; i < (ii + TI < _PB_M ? ii + TI : _PB_M); i++) {
        Q[i][k] = A[i][k] / R[k][k];
      }
    }
    
    #pragma omp parallel for schedule(static) private(j_inner, ii2, i)
    for (jj = k + 1; jj < _PB_N; jj += TJ) {
      for (j_inner = jj; j_inner < (jj + TJ < _PB_N ? jj + TJ : _PB_N); j_inner++) {
        R[k][j_inner] = SCALAR_VAL(0.0);
      }
      for (ii2 = 0; ii2 < _PB_M; ii2 += TI) {
        for (i = ii2; i < (ii2 + TI < _PB_M ? ii2 + TI : _PB_M); i++) {
          #pragma omp simd
          for (j_inner = jj; j_inner < (jj + TJ < _PB_N ? jj + TJ : _PB_N); j_inner++) {
            R[k][j_inner] += Q[i][k] * A[i][j_inner];
          }
        }
      }
      for (ii2 = 0; ii2 < _PB_M; ii2 += TI) {
        for (i = ii2; i < (ii2 + TI < _PB_M ? ii2 + TI : _PB_M); i++) {
          #pragma omp simd
          for (j_inner = jj; j_inner < (jj + TJ < _PB_N ? jj + TJ : _PB_N); j_inner++) {
            A[i][j_inner] = A[i][j_inner] - Q[i][k] * R[k][j_inner];
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
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,N,m,n);
  POLYBENCH_2D_ARRAY_DECL(R,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(Q,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(R),
	      POLYBENCH_ARRAY(Q));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gramschmidt (m, n,
		      POLYBENCH_ARRAY(A),
		      POLYBENCH_ARRAY(R),
		      POLYBENCH_ARRAY(Q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Q);

  return 0;
}
