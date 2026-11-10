/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* seidel-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "seidel-2d.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
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
      if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
    }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_seidel_2d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int t, i, j;

#pragma scop
  if (_PB_TSTEPS >= 1 && _PB_N >= 3)
  {
    /* Skewed space: T = t, I = t + i, J = 2t + i + j.
       All dependence distances become non-negative -> rectangular tiling
       legal; tiles on the same wavefront (tt+ii+jj = w) are independent. */
    const int BT = 16, BI = 32, BJ = 32;
    const int TS = _PB_TSTEPS;
    const int NN = _PB_N;
    const int ntt = (TS - 1) / BT;
    const int nii = (TS + NN - 3) / BI;
    const int njj = (2 * TS + 2 * NN - 6) / BJ;
    const int nw = ntt + nii + njj;

#pragma omp parallel private(t, i, j)
    {
      int w;
      for (w = 0; w <= nw; w++) {
        int ttlo = w - nii - njj; if (ttlo < 0) ttlo = 0;
        int tthi = w < ntt ? w : ntt;
        int tt;
#pragma omp for schedule(dynamic, 1)
        for (tt = ttlo; tt <= tthi; tt++) {
          int iilo = w - tt - njj; if (iilo < 0) iilo = 0;
          int iihi = (w - tt) < nii ? (w - tt) : nii;
          int ii;
          for (ii = iilo; ii <= iihi; ii++) {
            int jj = w - tt - ii;
            int tlo = tt * BT;
            int thi = tlo + BT - 1; if (thi > TS - 1) thi = TS - 1;
            for (t = tlo; t <= thi; t++) {
              int Ilo = ii * BI; if (Ilo < t + 1) Ilo = t + 1;
              int Ihi = ii * BI + BI - 1; if (Ihi > t + NN - 2) Ihi = t + NN - 2;
              int I;
              for (I = Ilo; I <= Ihi; I++) {
                i = I - t;
                int jlo = jj * BJ - 2 * t - i; if (jlo < 1) jlo = 1;
                int jhi = jj * BJ + BJ - 1 - 2 * t - i; if (jhi > NN - 2) jhi = NN - 2;
                for (j = jlo; j <= jhi; j++)
                  A[i][j] = (A[i-1][j-1] + A[i-1][j] + A[i-1][j+1]
                             + A[i][j-1] + A[i][j] + A[i][j+1]
                             + A[i+1][j-1] + A[i+1][j] + A[i+1][j+1])/SCALAR_VAL(9.0);
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
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_seidel_2d (tsteps, n, POLYBENCH_ARRAY(A));

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