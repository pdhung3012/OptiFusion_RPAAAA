/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* heat-3d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "heat-3d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		 DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        A[i][j][k] = B[i][j][k] = (DATA_TYPE) (i + j + (n-k))* 10 / (n);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n))

{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {
         if ((i * n * n + j * n + k) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
         fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
      }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_heat_3d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		      DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int t, i, j, k;

#pragma scop
    const int TILE_SIZE = 32;
    for (t = 1; t <= TSTEPS; t++) {
        #pragma omp parallel
        {
            #pragma omp for collapse(2) schedule(static)
            for (int i = 1; i < _PB_N-1; i += TILE_SIZE) {
                for (int j = 1; j < _PB_N-1; j += TILE_SIZE) {
                    for (int k = 1; k < _PB_N-1; k += TILE_SIZE) {
                        int i_end = (i + TILE_SIZE) < (_PB_N-1) ? (i + TILE_SIZE) : (_PB_N-1);
                        int j_end = (j + TILE_SIZE) < (_PB_N-1) ? (j + TILE_SIZE) : (_PB_N-1);
                        int k_end = (k + TILE_SIZE) < (_PB_N-1) ? (k + TILE_SIZE) : (_PB_N-1);
                        for (int ii = i; ii < i_end; ii++) {
                            for (int jj = j; jj < j_end; jj++) {
                                #pragma omp simd
                                for (int kk = k; kk < k_end; kk++) {
                                    B[ii][jj][kk] =   SCALAR_VAL(0.125) * (A[ii+1][jj][kk] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii-1][jj][kk])
                                                     + SCALAR_VAL(0.125) * (A[ii][jj+1][kk] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii][jj-1][kk])
                                                     + SCALAR_VAL(0.125) * (A[ii][jj][kk+1] - SCALAR_VAL(2.0) * A[ii][jj][kk] + A[ii][jj][kk-1])
                                                     + A[ii][jj][kk];
                                }
                            }
                        }
                    }
                }
            }
            
            #pragma omp for collapse(2) schedule(static)
            for (int i = 1; i < _PB_N-1; i += TILE_SIZE) {
                for (int j = 1; j < _PB_N-1; j += TILE_SIZE) {
                    for (int k = 1; k < _PB_N-1; k += TILE_SIZE) {
                        int i_end = (i + TILE_SIZE) < (_PB_N-1) ? (i + TILE_SIZE) : (_PB_N-1);
                        int j_end = (j + TILE_SIZE) < (_PB_N-1) ? (j + TILE_SIZE) : (_PB_N-1);
                        int k_end = (k + TILE_SIZE) < (_PB_N-1) ? (k + TILE_SIZE) : (_PB_N-1);
                        for (int ii = i; ii < i_end; ii++) {
                            for (int jj = j; jj < j_end; jj++) {
                                #pragma omp simd
                                for (int kk = k; kk < k_end; kk++) {
                                    A[ii][jj][kk] =   SCALAR_VAL(0.125) * (B[ii+1][jj][kk] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii-1][jj][kk])
                                                     + SCALAR_VAL(0.125) * (B[ii][jj+1][kk] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii][jj-1][kk])
                                                     + SCALAR_VAL(0.125) * (B[ii][jj][kk+1] - SCALAR_VAL(2.0) * B[ii][jj][kk] + B[ii][jj][kk-1])
                                                     + B[ii][jj][kk];
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
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_3D_ARRAY_DECL(A, DATA_TYPE, N, N, N, n, n, n);
  POLYBENCH_3D_ARRAY_DECL(B, DATA_TYPE, N, N, N, n, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_heat_3d (tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

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
