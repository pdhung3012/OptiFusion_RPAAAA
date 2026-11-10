/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* bicg.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "bicg.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		 DATA_TYPE POLYBENCH_1D(r,N,n),
		 DATA_TYPE POLYBENCH_1D(p,M,m))
{
  int i, j;

  for (i = 0; i < m; i++)
    p[i] = (DATA_TYPE)(i % m) / m;
  for (i = 0; i < n; i++) {
    r[i] = (DATA_TYPE)(i % n) / n;
    for (j = 0; j < m; j++)
      A[i][j] = (DATA_TYPE) (i*(j+1) % n)/n;
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_1D(s,M,m),
		 DATA_TYPE POLYBENCH_1D(q,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("s");
  for (i = 0; i < m; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, s[i]);
  }
  POLYBENCH_DUMP_END("s");
  POLYBENCH_DUMP_BEGIN("q");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, q[i]);
  }
  POLYBENCH_DUMP_END("q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_bicg(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		 DATA_TYPE POLYBENCH_1D(s,M,m),
		 DATA_TYPE POLYBENCH_1D(q,N,n),
		 DATA_TYPE POLYBENCH_1D(p,M,m),
		 DATA_TYPE POLYBENCH_1D(r,N,n))
{
  int i, j;

#pragma scop
  /* Tiling/blocking on j for better cache locality; parallelize over i with a private
     accumulation buffer for s to avoid races (s is a reduction over i). */
  {
    const int JT = 64; /* conservative cache-friendly tile; boundary handled below */

#pragma omp parallel default(none) shared(A, s, q, p, r) firstprivate(JT) private(i, j)
    {
      DATA_TYPE *s_priv = (DATA_TYPE*) polybench_alloc_data((size_t)_PB_M, sizeof(DATA_TYPE));
      if (s_priv) {
#pragma omp for schedule(static)
        for (j = 0; j < _PB_M; j++)
          s_priv[j] = SCALAR_VAL(0.0);

#pragma omp for schedule(static)
        for (i = 0; i < _PB_N; i++) {
          DATA_TYPE qi = SCALAR_VAL(0.0);
          for (j = 0; j < _PB_M; j++) {
            DATA_TYPE aij = A[i][j];
            s_priv[j] = s_priv[j] + r[i] * aij;
            qi = qi + aij * p[j];
          }
          q[i] = qi;
        }

        /* Reduce s_priv into shared s with j-tiling to improve locality and to
           limit time in atomic updates. */
#pragma omp for schedule(static)
        for (int jj = 0; jj < _PB_M; jj += JT) {
          int jend = jj + JT;
          if (jend > _PB_M) jend = _PB_M;
          for (j = jj; j < jend; j++) {
#pragma omp atomic update
            s[j] += s_priv[j];
          }
        }

        polybench_free_data(s_priv);
      } else {
        /* Fallback: original sequential semantics inside parallel region disabled.
           If allocation fails, execute the baseline computation single-threaded. */
#pragma omp single
        {
          for (j = 0; j < _PB_M; j++)
            s[j] = 0;
          for (i = 0; i < _PB_N; i++) {
            q[i] = SCALAR_VAL(0.0);
            for (j = 0; j < _PB_M; j++) {
              s[j] = s[j] + r[i] * A[i][j];
              q[i] = q[i] + A[i][j] * p[j];
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
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, M, n, m);
  POLYBENCH_1D_ARRAY_DECL(s, DATA_TYPE, M, m);
  POLYBENCH_1D_ARRAY_DECL(q, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(p, DATA_TYPE, M, m);
  POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(r),
	      POLYBENCH_ARRAY(p));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_bicg (m, n,
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(s),
	       POLYBENCH_ARRAY(q),
	       POLYBENCH_ARRAY(p),
	       POLYBENCH_ARRAY(r));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(s), POLYBENCH_ARRAY(q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(s);
  POLYBENCH_FREE_ARRAY(q);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(r);

  return 0;
}