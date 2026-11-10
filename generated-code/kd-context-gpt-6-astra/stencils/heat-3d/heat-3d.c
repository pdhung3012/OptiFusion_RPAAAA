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
    {
      const int n0 = _PB_N;
      const int i_lb = 1, i_ub = n0 - 2;
      const int j_lb = 1, j_ub = n0 - 2;
      const int k_lb = 1, k_ub = n0 - 2;

      /* Spatial tiling to improve cache locality; keep time steps sequential. */
      const int TI = 16;
      const int TJ = 16;
      const int TK = 64;

      for (t = 1; t <= TSTEPS; t++) {

        /* Use one persistent team per time step; barrier between phases is required. */
#pragma omp parallel default(none) shared(A, B) firstprivate(i_lb, i_ub, j_lb, j_ub, k_lb, k_ub)
        {
          int ii, jj, kk, i2, j2, k2;

#pragma omp for collapse(2) schedule(static)
          for (ii = i_lb; ii <= i_ub; ii += TI) {
            for (jj = j_lb; jj <= j_ub; jj += TJ) {
              for (kk = k_lb; kk <= k_ub; kk += TK) {
                const int iend = (ii + TI - 1 <= i_ub) ? (ii + TI - 1) : i_ub;
                const int jend = (jj + TJ - 1 <= j_ub) ? (jj + TJ - 1) : j_ub;
                const int kend = (kk + TK - 1 <= k_ub) ? (kk + TK - 1) : k_ub;

                for (i2 = ii; i2 <= iend; i2++) {
                  for (j2 = jj; j2 <= jend; j2++) {
#pragma omp simd
                    for (k2 = kk; k2 <= kend; k2++) {
                      DATA_TYPE aijk = A[i2][j2][k2];
                      B[i2][j2][k2] =
                          SCALAR_VAL(0.125) * (A[i2+1][j2][k2] - SCALAR_VAL(2.0) * aijk + A[i2-1][j2][k2])
                        + SCALAR_VAL(0.125) * (A[i2][j2+1][k2] - SCALAR_VAL(2.0) * aijk + A[i2][j2-1][k2])
                        + SCALAR_VAL(0.125) * (A[i2][j2][k2+1] - SCALAR_VAL(2.0) * aijk + A[i2][j2][k2-1])
                        + aijk;
                    }
                  }
                }
              }
            }
          }

          /* Implicit barrier here is required before reading B for the A update. */

#pragma omp for collapse(2) schedule(static)
          for (ii = i_lb; ii <= i_ub; ii += TI) {
            for (jj = j_lb; jj <= j_ub; jj += TJ) {
              for (kk = k_lb; kk <= k_ub; kk += TK) {
                const int iend = (ii + TI - 1 <= i_ub) ? (ii + TI - 1) : i_ub;
                const int jend = (jj + TJ - 1 <= j_ub) ? (jj + TJ - 1) : j_ub;
                const int kend = (kk + TK - 1 <= k_ub) ? (kk + TK - 1) : k_ub;

                for (i2 = ii; i2 <= iend; i2++) {
                  for (j2 = jj; j2 <= jend; j2++) {
#pragma omp simd
                    for (k2 = kk; k2 <= kend; k2++) {
                      DATA_TYPE bijk = B[i2][j2][k2];
                      A[i2][j2][k2] =
                          SCALAR_VAL(0.125) * (B[i2+1][j2][k2] - SCALAR_VAL(2.0) * bijk + B[i2-1][j2][k2])
                        + SCALAR_VAL(0.125) * (B[i2][j2+1][k2] - SCALAR_VAL(2.0) * bijk + B[i2][j2-1][k2])
                        + SCALAR_VAL(0.125) * (B[i2][j2][k2+1] - SCALAR_VAL(2.0) * bijk + B[i2][j2][k2-1])
                        + bijk;
                    }
                  }
                }
              }
            }
          }
        } /* omp parallel */
      } /* t */
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