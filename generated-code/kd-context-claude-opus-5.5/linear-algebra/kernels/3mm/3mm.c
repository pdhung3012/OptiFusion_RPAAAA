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
    const int TI = 32;
    const int TJ = 256;
    const int TK = 32;
    (void)i; (void)j; (void)k;

#pragma omp parallel
    {
      /* E := A*B  (independent of F, so no barrier afterwards) */
#pragma omp for collapse(2) schedule(static) nowait
      for (int ii = 0; ii < _PB_NI; ii += TI)
	for (int jj = 0; jj < _PB_NJ; jj += TJ)
	  {
	    int ie = (ii + TI < _PB_NI) ? ii + TI : _PB_NI;
	    int je = (jj + TJ < _PB_NJ) ? jj + TJ : _PB_NJ;
	    for (int i2 = ii; i2 < ie; i2++)
#pragma omp simd
	      for (int j2 = jj; j2 < je; j2++)
		E[i2][j2] = SCALAR_VAL(0.0);
	    for (int kk = 0; kk < _PB_NK; kk += TK)
	      {
		int ke = (kk + TK < _PB_NK) ? kk + TK : _PB_NK;
		for (int i2 = ii; i2 < ie; i2++)
		  for (int k2 = kk; k2 < ke; k2++)
		    {
		      DATA_TYPE a = A[i2][k2];
#pragma omp simd
		      for (int j2 = jj; j2 < je; j2++)
			E[i2][j2] += a * B[k2][j2];
		    }
	      }
	  }

      /* F := C*D  (barrier at end: G needs both E and F) */
#pragma omp for collapse(2) schedule(static)
      for (int ii = 0; ii < _PB_NJ; ii += TI)
	for (int jj = 0; jj < _PB_NL; jj += TJ)
	  {
	    int ie = (ii + TI < _PB_NJ) ? ii + TI : _PB_NJ;
	    int je = (jj + TJ < _PB_NL) ? jj + TJ : _PB_NL;
	    for (int i2 = ii; i2 < ie; i2++)
#pragma omp simd
	      for (int j2 = jj; j2 < je; j2++)
		F[i2][j2] = SCALAR_VAL(0.0);
	    for (int kk = 0; kk < _PB_NM; kk += TK)
	      {
		int ke = (kk + TK < _PB_NM) ? kk + TK : _PB_NM;
		for (int i2 = ii; i2 < ie; i2++)
		  for (int k2 = kk; k2 < ke; k2++)
		    {
		      DATA_TYPE c = C[i2][k2];
#pragma omp simd
		      for (int j2 = jj; j2 < je; j2++)
			F[i2][j2] += c * D[k2][j2];
		    }
	      }
	  }

      /* G := E*F */
#pragma omp for collapse(2) schedule(static)
      for (int ii = 0; ii < _PB_NI; ii += TI)
	for (int jj = 0; jj < _PB_NL; jj += TJ)
	  {
	    int ie = (ii + TI < _PB_NI) ? ii + TI : _PB_NI;
	    int je = (jj + TJ < _PB_NL) ? jj + TJ : _PB_NL;
	    for (int i2 = ii; i2 < ie; i2++)
#pragma omp simd
	      for (int j2 = jj; j2 < je; j2++)
		G[i2][j2] = SCALAR_VAL(0.0);
	    for (int kk = 0; kk < _PB_NJ; kk += TK)
	      {
		int ke = (kk + TK < _PB_NJ) ? kk + TK : _PB_NJ;
		for (int i2 = ii; i2 < ie; i2++)
		  for (int k2 = kk; k2 < ke; k2++)
		    {
		      DATA_TYPE e = E[i2][k2];
#pragma omp simd
		      for (int j2 = jj; j2 < je; j2++)
			G[i2][j2] += e * F[k2][j2];
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