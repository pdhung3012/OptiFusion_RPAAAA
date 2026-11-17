/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* deriche.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "deriche.h"


/* Array initialization. */
static
void init_array (int w, int h, DATA_TYPE* alpha,
		 DATA_TYPE POLYBENCH_2D(imgIn,W,H,w,h),
		 DATA_TYPE POLYBENCH_2D(imgOut,W,H,w,h))
{
  int i, j;

  *alpha=0.25; //parameter of the filter

  //input should be between 0 and 1 (grayscale image pixel)
  for (i = 0; i < w; i++)
     for (j = 0; j < h; j++)
	imgIn[i][j] = (DATA_TYPE) ((313*i+991*j)%65536) / 65535.0f;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int w, int h,
		 DATA_TYPE POLYBENCH_2D(imgOut,W,H,w,h))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("imgOut");
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++) {
      if ((i * h + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, imgOut[i][j]);
    }
  POLYBENCH_DUMP_END("imgOut");
  POLYBENCH_DUMP_FINISH;
}



/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Original code provided by Gael Deest */
static
void kernel_deriche(int w, int h, DATA_TYPE alpha,
       DATA_TYPE POLYBENCH_2D(imgIn, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(y1, W, H, w, h),
       DATA_TYPE POLYBENCH_2D(y2, W, H, w, h)) {
    int i,j;
    DATA_TYPE xm1, tm1, ym1, ym2;
    DATA_TYPE xp1, xp2;
    DATA_TYPE tp1, tp2;
    DATA_TYPE yp1, yp2;

    DATA_TYPE k;
    DATA_TYPE a1, a2, a3, a4, a5, a6, a7, a8;
    DATA_TYPE b1, b2, c1, c2;

#pragma scop
   k = (SCALAR_VAL(1.0)-EXP_FUN(-alpha))*(SCALAR_VAL(1.0)-EXP_FUN(-alpha))/(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*alpha*EXP_FUN(-alpha)-EXP_FUN(SCALAR_VAL(2.0)*alpha));
   a1 = a5 = k;
   a2 = a6 = k*EXP_FUN(-alpha)*(alpha-SCALAR_VAL(1.0));
   a3 = a7 = k*EXP_FUN(-alpha)*(alpha+SCALAR_VAL(1.0));
   a4 = a8 = -k*EXP_FUN(SCALAR_VAL(-2.0)*alpha);
   b1 =  POW_FUN(SCALAR_VAL(2.0),-alpha);
   b2 = -EXP_FUN(SCALAR_VAL(-2.0)*alpha);
   c1 = c2 = 1;

   /* Fused horizontal passes: each row is independent, so the three
      row-wise loops (forward y1, backward y2, and combine into imgOut)
      are fused and parallelized over rows. */
#pragma omp parallel for private(i,j) schedule(static)
   for (i=0; i<_PB_W; i++) {
        DATA_TYPE lxm1, lym1, lym2;
        DATA_TYPE lxp1, lxp2;
        DATA_TYPE lyp1, lyp2;

        lym1 = SCALAR_VAL(0.0);
        lym2 = SCALAR_VAL(0.0);
        lxm1 = SCALAR_VAL(0.0);
        for (j=0; j<_PB_H; j++) {
            y1[i][j] = a1*imgIn[i][j] + a2*lxm1 + b1*lym1 + b2*lym2;
            lxm1 = imgIn[i][j];
            lym2 = lym1;
            lym1 = y1[i][j];
        }

        lyp1 = SCALAR_VAL(0.0);
        lyp2 = SCALAR_VAL(0.0);
        lxp1 = SCALAR_VAL(0.0);
        lxp2 = SCALAR_VAL(0.0);
        for (j=_PB_H-1; j>=0; j--) {
            y2[i][j] = a3*lxp1 + a4*lxp2 + b1*lyp1 + b2*lyp2;
            lxp2 = lxp1;
            lxp1 = imgIn[i][j];
            lyp2 = lyp1;
            lyp1 = y2[i][j];
        }

        for (j=0; j<_PB_H; j++)
            imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);
    }

   /* Fused vertical passes: each column is independent, so the three
      column-wise loops (forward y1, backward y2, and combine into
      imgOut) are fused and parallelized over columns. */
#pragma omp parallel for private(i,j) schedule(static)
   for (j=0; j<_PB_H; j++) {
        DATA_TYPE ltm1, lym1, lym2;
        DATA_TYPE ltp1, ltp2;
        DATA_TYPE lyp1, lyp2;

        ltm1 = SCALAR_VAL(0.0);
        lym1 = SCALAR_VAL(0.0);
        lym2 = SCALAR_VAL(0.0);
        for (i=0; i<_PB_W; i++) {
            y1[i][j] = a5*imgOut[i][j] + a6*ltm1 + b1*lym1 + b2*lym2;
            ltm1 = imgOut[i][j];
            lym2 = lym1;
            lym1 = y1[i][j];
        }

        ltp1 = SCALAR_VAL(0.0);
        ltp2 = SCALAR_VAL(0.0);
        lyp1 = SCALAR_VAL(0.0);
        lyp2 = SCALAR_VAL(0.0);
        for (i=_PB_W-1; i>=0; i--) {
            y2[i][j] = a7*ltp1 + a8*ltp2 + b1*lyp1 + b2*lyp2;
            ltp2 = ltp1;
            ltp1 = imgOut[i][j];
            lyp2 = lyp1;
            lyp1 = y2[i][j];
        }

        for (i=0; i<_PB_W; i++)
            imgOut[i][j] = c2*(y1[i][j] + y2[i][j]);
    }

#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int w = W;
  int h = H;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  POLYBENCH_2D_ARRAY_DECL(imgIn, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(imgOut, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(y1, DATA_TYPE, W, H, w, h);
  POLYBENCH_2D_ARRAY_DECL(y2, DATA_TYPE, W, H, w, h);


  /* Initialize array(s). */
  init_array (w, h, &alpha, POLYBENCH_ARRAY(imgIn), POLYBENCH_ARRAY(imgOut));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_deriche (w, h, alpha, POLYBENCH_ARRAY(imgIn), POLYBENCH_ARRAY(imgOut), POLYBENCH_ARRAY(y1), POLYBENCH_ARRAY(y2));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(w, h, POLYBENCH_ARRAY(imgOut)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(imgIn);
  POLYBENCH_FREE_ARRAY(imgOut);
  POLYBENCH_FREE_ARRAY(y1);
  POLYBENCH_FREE_ARRAY(y2);

  return 0;
}