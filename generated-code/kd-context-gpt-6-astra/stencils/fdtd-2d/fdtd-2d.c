/*******************************************************
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* fdtd-2d.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "fdtd-2d.h"


/* Array initialization. */
static
void init_array (int tmax,
		 int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int i, j;

  for (i = 0; i < tmax; i++)
    _fict_[i] = (DATA_TYPE) i;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      {
	ex[i][j] = ((DATA_TYPE) i*(j+1)) / nx;
	ey[i][j] = ((DATA_TYPE) i*(j+2)) / ny;
	hz[i][j] = ((DATA_TYPE) i*(j+3)) / nx;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nx,
		 int ny,
		 DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("ex");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ex[i][j]);
    }
  POLYBENCH_DUMP_END("ex");
  POLYBENCH_DUMP_FINISH;

  POLYBENCH_DUMP_BEGIN("ey");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, ey[i][j]);
    }
  POLYBENCH_DUMP_END("ey");

  POLYBENCH_DUMP_BEGIN("hz");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++) {
      if ((i * nx + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, hz[i][j]);
    }
  POLYBENCH_DUMP_END("hz");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_fdtd_2d(int tmax,
		    int nx,
		    int ny,
		    DATA_TYPE POLYBENCH_2D(ex,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(ey,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_2D(hz,NX,NY,nx,ny),
		    DATA_TYPE POLYBENCH_1D(_fict_,TMAX,tmax))
{
  int t, i, j;

#pragma scop

  /* Spatial tiling for locality; single persistent parallel region per time step. */
  for (t = 0; t < _PB_TMAX; t++)
    {
      const int nx_ = _PB_NX;
      const int ny_ = _PB_NY;
      const int nxm1 = nx_ - 1;
      const int nym1 = ny_ - 1;

      /* Tile sizes chosen to improve cache locality while keeping overhead low. */
      const int TI = 32;
      const int TJ = 128;

      int ii, jj;

#pragma omp parallel default(none) shared(ex,ey,hz,_fict_) firstprivate(t,nx_,ny_,nxm1,nym1,TI,TJ) private(i,j,ii,jj)
      {
        /* Set ey[0][j] = _fict_[t] */
#pragma omp for schedule(static)
        for (j = 0; j < ny_; j++)
          ey[0][j] = _fict_[t];

        /* Update ey[i][j], i=1..nx-1, j=0..ny-1 */
#pragma omp for schedule(static) collapse(2)
        for (ii = 1; ii < nx_; ii += TI)
          for (jj = 0; jj < ny_; jj += TJ)
            {
              int iend = ii + TI;
              int jend = jj + TJ;
              if (iend > nx_) iend = nx_;
              if (jend > ny_) jend = ny_;
              for (i = ii; i < iend; i++)
                {
#pragma omp simd
                  for (j = jj; j < jend; j++)
                    ey[i][j] = ey[i][j] - SCALAR_VAL(0.5) * (hz[i][j] - hz[i-1][j]);
                }
            }

        /* Update ex[i][j], i=0..nx-1, j=1..ny-1 */
#pragma omp for schedule(static) collapse(2)
        for (ii = 0; ii < nx_; ii += TI)
          for (jj = 1; jj < ny_; jj += TJ)
            {
              int iend = ii + TI;
              int jend = jj + TJ;
              if (iend > nx_) iend = nx_;
              if (jend > ny_) jend = ny_;
              for (i = ii; i < iend; i++)
                {
#pragma omp simd
                  for (j = jj; j < jend; j++)
                    ex[i][j] = ex[i][j] - SCALAR_VAL(0.5) * (hz[i][j] - hz[i][j-1]);
                }
            }

        /* Update hz[i][j], i=0..nx-2, j=0..ny-2 (depends on updated ex/ey) */
#pragma omp for schedule(static) collapse(2)
        for (ii = 0; ii < nxm1; ii += TI)
          for (jj = 0; jj < nym1; jj += TJ)
            {
              int iend = ii + TI;
              int jend = jj + TJ;
              if (iend > nxm1) iend = nxm1;
              if (jend > nym1) jend = nym1;
              for (i = ii; i < iend; i++)
                {
#pragma omp simd
                  for (j = jj; j < jend; j++)
                    hz[i][j] = hz[i][j] - SCALAR_VAL(0.7) * (ex[i][j+1] - ex[i][j] +
                                                            ey[i+1][j] - ey[i][j]);
                }
            }
      } /* end parallel */
    }

#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int tmax = TMAX;
  int nx = NX;
  int ny = NY;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(ex,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(ey,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_2D_ARRAY_DECL(hz,DATA_TYPE,NX,NY,nx,ny);
  POLYBENCH_1D_ARRAY_DECL(_fict_,DATA_TYPE,TMAX,tmax);

  /* Initialize array(s). */
  init_array (tmax, nx, ny,
	      POLYBENCH_ARRAY(ex),
	      POLYBENCH_ARRAY(ey),
	      POLYBENCH_ARRAY(hz),
	      POLYBENCH_ARRAY(_fict_));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_fdtd_2d (tmax, nx, ny,
		  POLYBENCH_ARRAY(ex),
		  POLYBENCH_ARRAY(ey),
		  POLYBENCH_ARRAY(hz),
		  POLYBENCH_ARRAY(_fict_));


  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nx, ny, POLYBENCH_ARRAY(ex),
				    POLYBENCH_ARRAY(ey),
				    POLYBENCH_ARRAY(hz)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(ex);
  POLYBENCH_FREE_ARRAY(ey);
  POLYBENCH_FREE_ARRAY(hz);
  POLYBENCH_FREE_ARRAY(_fict_);

  return 0;
}