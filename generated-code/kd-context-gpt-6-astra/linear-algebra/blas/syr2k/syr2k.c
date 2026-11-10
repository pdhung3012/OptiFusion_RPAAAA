/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* syr2k.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "syr2k.h"


/* Array initialization. */
static
void init_array(int n, int m,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		DATA_TYPE POLYBENCH_2D(B,N,M,n,m))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      A[i][j] = (DATA_TYPE) ((i*j+1)%n) / n;
      B[i][j] = (DATA_TYPE) ((i*j+2)%m) / m;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      C[i][j] = (DATA_TYPE) ((i*j+3)%n) / m;
    }
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
void kernel_syr2k(int n, int m,
		  DATA_TYPE alpha,
		  DATA_TYPE beta,
		  DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		  DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		  DATA_TYPE POLYBENCH_2D(B,N,M,n,m))
{
  int i, j, k;

//BLAS PARAMS
//UPLO  = 'L'
//TRANS = 'N'
//A is NxM
//B is NxM
//C is NxN
#pragma scop
  {
    const int Nn = _PB_N;
    const int Mm = _PB_M;

    /* Tile sizes chosen to improve cache locality while keeping overhead low. */
    const int TI = 32;
    const int TJ = 32;
    const int TK = 32;

    int ii, jj, kk;

    #pragma omp parallel default(none) shared(C,A,B,alpha,beta,Nn,Mm) private(ii,jj,kk,i,j,k)
    {
      /* Phase 1: C = beta * C over lower triangle */
      #pragma omp for schedule(static)
      for (ii = 0; ii < Nn; ii += TI) {
        int i_end = ii + TI; if (i_end > Nn) i_end = Nn;
        for (i = ii; i < i_end; ++i) {
          int j_end = i + 1;
          for (jj = 0; jj < j_end; jj += TJ) {
            int jb_end = jj + TJ; if (jb_end > j_end) jb_end = j_end;
            #pragma omp simd
            for (j = jj; j < jb_end; ++j)
              C[i][j] *= beta;
          }
        }
      }

      /* Phase 2: rank-2k update, tiled in k and j for locality */
      #pragma omp for schedule(static)
      for (ii = 0; ii < Nn; ii += TI) {
        int i_end = ii + TI; if (i_end > Nn) i_end = Nn;
        for (kk = 0; kk < Mm; kk += TK) {
          int k_end = kk + TK; if (k_end > Mm) k_end = Mm;
          for (i = ii; i < i_end; ++i) {
            DATA_TYPE a_ik, b_ik;
            for (jj = 0; jj <= i; jj += TJ) {
              int jb_end = jj + TJ; if (jb_end > i + 1) jb_end = i + 1;
              for (k = kk; k < k_end; ++k) {
                b_ik = B[i][k] * alpha;
                a_ik = A[i][k] * alpha;
                #pragma omp simd
                for (j = jj; j < jb_end; ++j) {
                  C[i][j] += A[j][k] * b_ik + B[j][k] * a_ik;
                }
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
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,N,M,n,m);

  /* Initialize array(s). */
  init_array (n, m, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_syr2k (n, m,
		alpha, beta,
		POLYBENCH_ARRAY(C),
		POLYBENCH_ARRAY(A),
		POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}