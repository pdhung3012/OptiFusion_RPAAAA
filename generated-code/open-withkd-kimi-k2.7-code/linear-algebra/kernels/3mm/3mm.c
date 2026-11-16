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

#pragma scop
  {
#define TS 32
#define MIN2(a,b) ((a)<(b)?(a):(b))
    int ii, jj, kk;

#pragma omp parallel private(i,j,k,ii,jj,kk)
    {
      /* E := A*B */
#pragma omp for collapse(2) schedule(static) nowait
      for (i = 0; i < _PB_NI; i += TS)
	for (j = 0; j < _PB_NJ; j += TS)
	{
	  int i_end = MIN2(i + TS, _PB_NI);
	  int j_end = MIN2(j + TS, _PB_NJ);
	  for (ii = i; ii < i_end; ++ii)
	    for (jj = j; jj < j_end; ++jj)
	      E[ii][jj] = SCALAR_VAL(0.0);
	  for (k = 0; k < _PB_NK; k += TS)
	  {
	    int k_end = MIN2(k + TS, _PB_NK);
	    for (ii = i; ii < i_end; ++ii)
	      for (kk = k; kk < k_end; ++kk)
#pragma omp simd
		for (jj = j; jj < j_end; ++jj)
		  E[ii][jj] += A[ii][kk] * B[kk][jj];
	  }
	}

      /* F := C*D */
#pragma omp for collapse(2) schedule(static)
      for (i = 0; i < _PB_NJ; i += TS)
	for (j = 0; j < _PB_NL; j += TS)
	{
	  int i_end = MIN2(i + TS, _PB_NJ);
	  int j_end = MIN2(j + TS, _PB_NL);
	  for (ii = i; ii < i_end; ++ii)
	    for (jj = j; jj < j_end; ++jj)
	      F[ii][jj] = SCALAR_VAL(0.0);
	  for (k = 0; k < _PB_NM; k += TS)
	  {
	    int k_end = MIN2(k + TS, _PB_NM);
	    for (ii = i; ii < i_end; ++ii)
	      for (kk = k; kk < k_end; ++kk)
#pragma omp simd
		for (jj = j; jj < j_end; ++jj)
		  F[ii][jj] += C[ii][kk] * D[kk][jj];
	  }
	}

      /* G := E*F */
#pragma omp for collapse(2) schedule(static)
      for (i = 0; i < _PB_NI; i += TS)
	for (j = 0; j < _PB_NL; j += TS)
	{
	  int i_end = MIN2(i + TS, _PB_NI);
	  int j_end = MIN2(j + TS, _PB_NL);
	  for (ii = i; ii < i_end; ++ii)
	    for (jj = j; jj < j_end; ++jj)
	      G[ii][jj] = SCALAR_VAL(0.0);
	  for (k = 0; k < _PB_NJ; k += TS)
	  {
	    int k_end = MIN2(k + TS, _PB_NJ);
	    for (ii = i; ii < i_end; ++ii)
	      for (kk = k; kk < k_end; ++kk)
#pragma omp simd
		for (jj = j; jj < j_end; ++jj)
		  G[ii][jj] += E[ii][kk] * F[kk][jj];
	  }
	}
    }
#undef TS
#undef MIN2
  }
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