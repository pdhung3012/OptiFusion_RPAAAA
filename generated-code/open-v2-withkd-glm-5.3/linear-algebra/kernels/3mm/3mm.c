/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* 3mm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "3mm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk, int nl, int nm,
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(C,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(D,NM,NL,nm,nl))
{
  int i, j;

  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) ((i*j+1) % ni) / (5*ni);
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) ((i*(j+1)+2) % nj) / (5*nj);
  for (i = 0; i < nj; i++)
    for (j = 0; j < nm; j++)
      C[i][j] = (DATA_TYPE) (i*(j+3) % nl) / (5*nl);
  for (i = 0; i < nm; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = (DATA_TYPE) ((i*(j+2)+2) % nk) / (5*nk);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nl,
		 DATA_TYPE POLYBENCH_2D(G,NI,NL,ni,nl))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("G");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, G[i][j]);
    }
  POLYBENCH_DUMP_END("G");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_3mm(int ni, int nj, int nk, int nl, int nm,
		DATA_TYPE POLYBENCH_2D(E,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
		DATA_TYPE POLYBENCH_2D(F,NJ,NL,nj,nl),
		DATA_TYPE POLYBENCH_2D(C,NJ,NM,nj,nm),
		DATA_TYPE POLYBENCH_2D(D,NM,NL,nm,nl),
		DATA_TYPE POLYBENCH_2D(G,NI,NL,ni,nl))
{
  int i, j, k;
  int ii, jj, kk;
  int it, jt, kt;

#pragma scop
#define TILE 32
  /* E := A*B, F := C*D, G := E*F : cache-blocked, single parallel team */
#pragma omp parallel private(i, j, k, ii, jj, kk, it, jt, kt)
  {
    /* E := A*B */
#pragma omp for collapse(2) schedule(static)
    for (it = 0; it < _PB_NI; it += TILE)
      for (jt = 0; jt < _PB_NJ; jt += TILE) {
	int i_max = (it + TILE < _PB_NI) ? it + TILE : _PB_NI;
	int j_max = (jt + TILE < _PB_NJ) ? jt + TILE : _PB_NJ;
	for (i = it; i < i_max; i++)
	  for (j = jt; j < j_max; j++) {
	    DATA_TYPE sum = E[i][j];
	    sum = SCALAR_VAL(0.0);
	    for (k = 0; k < _PB_NK; ++k)
	      sum += A[i][k] * B[k][j];
	    E[i][j] = sum;
	  }
      }
    /* F := C*D (independent of E, but keep barrier before G) */
#pragma omp for collapse(2) schedule(static) nowait
    for (it = 0; it < _PB_NJ; it += TILE)
      for (jt = 0; jt < _PB_NL; jt += TILE) {
	int i_max = (it + TILE < _PB_NJ) ? it + TILE : _PB_NJ;
	int j_max = (jt + TILE < _PB_NL) ? jt + TILE : _PB_NL;
	for (i = it; i < i_max; i++)
	  for (j = jt; j < j_max; j++) {
	    DATA_TYPE sum = F[i][j];
	    sum = SCALAR_VAL(0.0);
	    for (k = 0; k < _PB_NM; ++k)
	      sum += C[i][k] * D[k][j];
	    F[i][j] = sum;
	  }
      }
    /* Barrier needed: G consumes E and F */
#pragma omp barrier
    /* G := E*F */
#pragma omp for collapse(2) schedule(static)
    for (it = 0; it < _PB_NI; it += TILE)
      for (jt = 0; jt < _PB_NL; jt += TILE) {
	int i_max = (it + TILE < _PB_NI) ? it + TILE : _PB_NI;
	int j_max = (jt + TILE < _PB_NL) ? jt + TILE : _PB_NL;
	for (i = it; i < i_max; i++)
	  for (j = jt; j < j_max; j++) {
	    DATA_TYPE sum = G[i][j];
	    sum = SCALAR_VAL(0.0);
	    for (k = 0; k < _PB_NJ; ++k)
	      sum += E[i][k] * F[k][j];
	    G[i][j] = sum;
	  }
      }
  }
#undef TILE
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int ni = NI;
  int nj = NJ;
  int nk = NK;
  int nl = NL;
  int nm = NM;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(E, DATA_TYPE, NI, NJ, ni, nj);
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, NI, NK, ni, nk);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, NK, NJ, nk, nj);
  POLYBENCH_2D_ARRAY_DECL(F, DATA_TYPE, NJ, NL, nj, nl);
  POLYBENCH_2D_ARRAY_DECL(C, DATA_TYPE, NJ, NM, nj, nm);
  POLYBENCH_2D_ARRAY_DECL(D, DATA_TYPE, NM, NL, nm, nl);
  POLYBENCH_2D_ARRAY_DECL(G, DATA_TYPE, NI, NL, ni, nl);

  /* Initialize array(s). */
  init_array (ni, nj, nk, nl, nm,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_3mm (ni, nj, nk, nl, nm,
	      POLYBENCH_ARRAY(E),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(F),
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(D),
	      POLYBENCH_ARRAY(G));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nl,  POLYBENCH_ARRAY(G)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(E);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(F);
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(D);
  POLYBENCH_FREE_ARRAY(G);

  return 0;
}