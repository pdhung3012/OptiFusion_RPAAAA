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
  /* Wavefront (skewed) parallelization.
     Original nest: t, then i (carried dep A[i-1][*]), then j (carried dep A[i][j-1]).
     Skew the time loop against i via the affine map s = i + 2*t:
     every dependence of (t,i,j) maps to a strictly smaller s, so all
     iterations sharing the same s (i.e., pairs (t, i=s-2t)) are
     independent and can run in parallel; j stays sequential per row.
     Reads of row i+1 (wavefront s+1, not yet executed) and rows of
     earlier wavefronts match exactly the in-place values seen by the
     original sequential code. */
  {
    const int _T = _PB_TSTEPS;
    const int _N = _PB_N;
    int _s, _t, _i, _j, _lo, _hi;

#pragma omp parallel private(_s,_t,_i,_j,_lo,_hi)
    for (_s = 1; _s <= (_N - 2) + 2 * (_T - 1); _s++) {
      /* valid t: 1 <= s - 2t <= _N-2  and  0 <= t <= _T-1 */
      _lo = (_s - (_N - 2) + 1) / 2;   /* ceil((s - (N-2)) / 2) */
      if (_lo < 0)
	_lo = 0;
      _hi = (_s - 1) / 2;
      if (_hi > _T - 1)
	_hi = _T - 1;

#pragma omp for schedule(static)
      for (_t = _lo; _t <= _hi; _t++) {
	_i = _s - 2 * _t;
	for (_j = 1; _j <= _N - 2; _j++)
	  A[_i][_j] = (A[_i-1][_j-1] + A[_i-1][_j] + A[_i-1][_j+1]
		       + A[_i][_j-1] + A[_i][_j] + A[_i][_j+1]
		       + A[_i+1][_j-1] + A[_i+1][_j] + A[_i+1][_j+1])
	    / SCALAR_VAL(9.0);
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