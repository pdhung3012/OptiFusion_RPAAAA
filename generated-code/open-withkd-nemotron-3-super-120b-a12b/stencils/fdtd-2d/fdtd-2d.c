/**
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
{
    int t, ii, jj, i, j, i_end, j_end;
    for(t = 0; t < _PB_TMAX; t++)
    {
        #pragma omp parallel private(ii, jj, i, j, i_end, j_end)
        {
            // First kernel: set ey[0][j]
            #pragma omp for schedule(static)
            for (j = 0; j < _PB_NY; j)
                ey[0][j] = _fict_[t];

            // Second kernel: update ey for i>=1
            #pragma omp for collapse(2) schedule(static)
            for (ii = 1; ii < _PB_NX; ii += 32)
            {
                for (jj = 0; jj < _PB_NY; jj += 32)
                {
                    i_end = ii + 32;
                    if (i_end > _PB_NX) i_end = _PB_NX;
                    j_end = jj + 32;
                    if (j_end > _PB_NY) j_end = _PB_NY;
                    for (i = ii; i < i_end; i++)
                    {
                        for (j = jj; j < j_end; j++)
                        {
                            ey[i][j] = ey[i][j] - SCALAR_VAL(0.5) * (hz[i][j] - hz[i-1][j]);
                        }
                    }
                }
            }

            // Third kernel: update ex for j>=1
            #pragma omp for collapse(2) schedule(static)
            for (ii = 0; ii < _PB_NX; ii += 32)
            {
                for (jj = 1; jj < _PB_NY; jj += 32)
                {
                    i_end = ii + 32;
                    if (i_end > _PB_NX) i_end = _PB_NX;
                    j_end = jj + 32;
                    if (j_end > _PB_NY) j_end = _PB_NY;
                    for (i = ii; i < i_end; i++)
                    {
                        for (j = jj; j < j_end; j++)
                        {
                            ex[i][j] = ex[i][j] - SCALAR_VAL(0.5) * (hz[i][j] - hz[i][j-1]);
                        }
                    }
                }
            }

            // Fourth kernel: update hz for i<_PB_NX-1 and j<_PB_NY-1
            #pragma omp for collapse(2) schedule(static)
            for (ii = 0; ii < _PB_NX-1; ii += 32)
            {
                for (jj = 0; jj < _PB_NY-1; jj += 32)
                {
                    i_end = ii + 32;
                    if (i_end > _PB_NX-1) i_end = _PB_NX-1;
                    j_end = jj + 32;
                    if (j_end > _PB_NY-1) j_end = _PB_NY-1;
                    for (i = ii; i < i_end; i++)
                    {
                        for (j = jj; j < j_end; j++)
                        {
                            hz[i][j] = hz[i][j] - SCALAR_VAL(0.7) * (ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
                        }
                    }
                }
            }
        } // end parallel
    } // end t loop
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
