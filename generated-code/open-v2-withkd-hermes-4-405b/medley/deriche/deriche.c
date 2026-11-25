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
We are given the deriche.c code and are to optimize the region between #pragma scop and #pragma endscop.
 The code has several nested loops. We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

 Steps:
 1. Loop tiling: Look for nested loops that can be tiled to improve cache locality. The deriche kernel has multiple loops over the same arrays (imgIn, imgOut, y1, y2). 
    The loops are:
      - Two forward passes (over j) for y1 and y2 (one in increasing j, one in decreasing j) for each i.
      - Then a fusion step for imgOut from y1 and y2.
      - Then two more passes (over i) for y1 and y2 (one in increasing i, one in decreasing i) for each j.
      - Then a final fusion step.

    However, note that the loops are not easily tileable in the same way as a matrix multiplication because of the dependencies (each step uses the previous element in the sequence).

    But observe: the first two loops (for y1 and y2) are independent per i. Similarly, the next two loops (for y1 and y2) are independent per j. 
    The fusion steps are also independent per (i,j).

    We can consider tiling the i and j loops. However, the inner loops are sequential (due to the recurrence). Therefore, we can tile the outer loops (i and j) to improve locality for the inner sequential loops.

    Specifically, we can tile the i-loop for the first two passes and the j-loop for the next two passes. The fusion steps can also be tiled.

    However, note that the first two passes are over j for fixed i, and the next two passes are over i for fixed j. This suggests that we might want to tile in blocks of (i,j) and then process the entire block.

    But the dependencies are along the inner dimension. So we can do:

        for i0 in 0..w by TI:
          for j0 in 0..h by TJ:
            // Process the block [i0:i0+TI-1, j0:j0+TJ-1] for the first two passes and the first fusion?
            // But note: the first two passes are row-wise (fixed i, varying j) and the next two are column-wise (fixed j, varying i).

    This is complex because the two sets of passes are orthogonal. We might instead tile the row-wise passes by i and the column-wise passes by j.

    Alternatively, we can tile the row-wise passes (the first two) by i and j? Actually, the inner loop is j and is sequential. So we can only tile the outer i loop for the row-wise passes? But then we are processing a set of consecutive rows. This might help because the rows are stored consecutively? Actually, the arrays are 2D and stored in row-major order. So consecutive rows are consecutive in memory. Therefore, tiling the i-loop (by rows) for the row-wise passes would help because we are processing a block of consecutive rows.

    Similarly, for the column-wise passes, we can tile the j-loop (by columns) because consecutive columns are stored with a stride of w. But note: the inner loop is i (which is the row index) and consecutive i's are consecutive in memory. So if we tile the j-loop, we are processing a set of consecutive columns. However, the inner loop over i is sequential and we are accessing the same column j for consecutive i's. This is contiguous in memory? Actually, for a fixed j, the elements for consecutive i are at addresses: base + i * (sizeof(DATA_TYPE)) for a column? But wait, the array is row-major: 
        imgOut[i][j] is at base + i * (h * sizeof(DATA_TYPE)) + j * sizeof(DATA_TYPE)
    So consecutive i's for fixed j are separated by h elements. This is not contiguous. Therefore, tiling the j-loop for the column-wise passes might not be as beneficial.

    Given the complexity and the fact that the inner loops are sequential, we might focus on tiling the outer loops that are parallel and then use OpenMP to parallelize the outer tiles.

    However, the problem says: "investigate tiling of suitable nested loops". The row-wise passes have an outer i loop and inner j loop. We can tile the i loop? But the inner j loop is sequential and we cannot change the order. So we can do:

        for (i0 = 0; i0 < _PB_W; i0 += TI) {
            for (i = i0; i < min(i0+TI, _PB_W); i++) {
                ... original inner j loop ...
            }
        }

    This would process TI consecutive rows at a time. This might improve locality because we are reusing the same rows of imgIn and y1 for the entire block? Actually, the inner j loop for a fixed i uses the entire row i of imgIn and produces the entire row i of y1. Then the next row i+1 uses the next row of imgIn. So if we process TI consecutive rows, we are accessing TI consecutive rows of imgIn and y1. This is good because they are consecutive in memory.

    Similarly, for the column-wise passes, we can tile the j loop? But note the inner loop is over i (which is the row index) and we are accessing the same column j for consecutive i's. This is not contiguous. However, if we tile the j loop, we are processing TJ consecutive columns. For each column j, we access the entire column (all rows) of imgOut and y1. The entire column is stored with a stride of w (the row length). So if we process TJ consecutive columns, we are accessing TJ columns that are stored with a stride of w. This might not be as cache friendly.

    Therefore, we decide to tile the row-wise passes by the outer i loop and the column-wise passes by the outer j loop? But note the column-wise passes have an outer j loop and inner i loop. We can tile the j loop? 

    However, the problem also says: "choose legal tile sizes". We can choose TI and TJ to be 64 or 32? But we must be cautious because the inner loops are sequential and we cannot break the dependence.

    Alternatively, we can use a different approach: tile the entire (i,j) space for the fusion steps? But the fusion steps are independent per (i,j). So we can do:

        for (i0=0; i0<_PB_W; i0+=TI)
          for (j0=0; j0<_PB_H; j0+=TJ)
            for (i=i0; i<min(i0+TI, _PB_W); i++)
              for (j=j0; j<min(j0+TJ, _PB_H); j++)
                imgOut[i][j] = ...;

    This would improve locality for the fusion steps. But note that the fusion steps are not the most expensive part. The row-wise and column-wise passes are the expensive ones.

    Given the time, we decide to tile the row-wise passes (the first two) by the outer i loop and the column-wise passes (the next two) by the outer j loop. We leave the fusion steps as they are (or we can tile them too, but they are simple).

    However, note that the first two passes are:

        for (i=0; i<_PB_W; i++) {
            ... // initialize
            for (j=0; j<_PB_H; j++) {
                ... // recurrence in j
            }
        }

        for (i=0; i<_PB_W; i++) {
            ... // initialize
            for (j=_PB_H-1; j>=0; j--) {
                ... // recurrence in j (backward)
            }
        }

    We can tile the outer i loop for both. Similarly, the next two passes:

        for (j=0; j<_PB_H; j++) {
            ... // initialize
            for (i=0; i<_PB_W; i++) {
                ... // recurrence in i
            }
        }

        for (j=0; j<_PB_H; j++) {
            ... // initialize
            for (i=_PB_W-1; i>=0; i--) {
                ... // recurrence in i (backward)
            }
        }

    We can tile the outer j loop for both.

    We choose tile sizes: let's use 64 for both TI and TJ? But we must make sure that the tile size is at least 1 and at most the dimension. We'll use:

        #define TI 64
        #define TJ 64

    But note: the problem says "choose legal tile sizes", so we can choose any that are safe. We'll use 64 as a reasonable starting point.

 2. Loop fusion: We have two fusion steps. They are already separate. We cannot fuse them with the passes because the passes produce y1 and y2 that are then used in the fusion. Also, the two fusion steps are separated by the two column-wise passes. So fusion is not straightforward.

 3. Loop skewing: The inner loops have carried dependencies (each iteration depends on the immediate previous one). Skewing might break the dependency? We are not going to skew because the dependencies are along the inner loop and we cannot parallelize the inner loop.

 Therefore, we focus on tiling the outer loops of the four main passes.

 Implementation:

   We will tile the first two passes (row-wise) by the outer i loop:

        #define TI 64
        for (i0 = 0; i0 < _PB_W; i0 += TI) {
            int i_end = min(i0+TI, _PB_W);
            for (i = i0; i < i_end; i++) {
                ... original inner j loop for the first pass ...
            }
        }

        Similarly for the second pass (backward j) and then for the two column-wise passes, we tile the outer j loop:

        #define TJ 64
        for (j0 = 0; j0 < _PB_H; j0 += TJ) {
            int j_end = min(j0+TJ, _PB_H);
            for (j = j0; j < j_end; j++) {
                ... original inner i loop for the third pass ...
            }
        }

        and similarly for the fourth pass.

   We leave the fusion steps as they are? Or we can tile them too? The fusion steps are:

        for (i=0; i<_PB_W; i++)
            for (j=0; j<_PB_H; j++)
                imgOut[i][j] = ...;

        and

        for (i=0; i<_PB_W; i++)
            for (j=0; j<_PB_H; j++)
                imgOut[i][j] = ...;

   We can tile these as well to improve locality. We'll use the same TI and TJ.

   So we do:

        for (i0=0; i0<_PB_W; i0+=TI)
          for (j0=0; j0<_PB_H; j0+=TJ)
            for (i=i0; i<min(i0+TI, _PB_W); i++)
              for (j=j0; j<min(j0+TJ, _PB_H); j++)
                imgOut[i][j] = ...;

   But note: the first fusion step uses y1 and y2 that were computed in the row-wise passes. The row-wise passes were tiled by rows. So the rows in a tile are consecutive. The fusion step for a tile of rows and columns would access a block of y1 and y2. This is good.

   Similarly, the second fusion step uses y1 and y2 from the column-wise passes. The column-wise passes were tiled by columns. So the columns in a tile are consecutive. The fusion step for a tile of rows and columns would access a block of y1 and y2. This is also good.

   Therefore, we tile all the loops.

 4. OpenMP: We can parallelize the outer tile loops. For the row-wise passes, we parallelize the i0 loop. For the column-wise passes, we parallelize the j0 loop. For the fusion steps, we parallelize the i0 loop (or we can collapse the i0 and j0 loops?).

    We'll use:

        #pragma omp parallel for
        for (i0=0; i0<_PB_W; i0+=TI) ...

    and similarly for the others.

    We must be cautious: the original code has dependencies between the passes. The first two passes (row-wise) are independent per i, so we can parallelize the i0 loop. Similarly, the next two passes (column-wise) are independent per j, so we can parallelize the j0 loop. The fusion steps are independent per (i,j) so we can parallelize the outer tile loops.

    We'll use static scheduling.

    Also, we can use collapse for the fusion steps? The fusion steps have two outer tile loops (i0 and j0). We can do:

        #pragma omp parallel for collapse(2)
        for (i0=0; i0<_PB_W; i0+=TI)
          for (j0=0; j0<_PB_H; j0+=TJ)
            ...

    This might give better load balancing.

    But note: the tile sizes TI and TJ are fixed. The number of tiles in i0 is ceil(_PB_W/TI) and in j0 is ceil(_PB_H/TJ). The total number of tiles is the product. This might be large enough to keep all threads busy.

    We'll use collapse(2) for the fusion steps.

    For the row-wise passes, we have only one outer loop (i0) so we don't collapse.

    Similarly, for the column-wise passes, we have only one outer loop (j0) so we don't collapse.

 5. We must define the tile sizes. We'll use:

        #define TI 64
        #define TJ 64

    But note: the problem says "choose legal tile sizes". We assume 64 is safe.

 6. We must be cautious about the backward loops. The backward loops run from high to low. We are tiling the outer loop (which is forward) and then the inner loop is backward. This is acceptable.

 7. We must also consider the variables that are declared inside the original loops. We move the declarations of the temporary variables (like ym1, ym2, xm1, etc.) inside the inner loop? Actually, they are declared inside the i loop. We are now having an outer i0 loop and then an inner i loop. We can leave them inside the inner i loop.

 Implementation:

   We restructure the code as follows:

   Step 1: Tile the first row-wise pass (forward j) by i0.

   Step 2: Tile the second row-wise pass (backward j) by i0.

   Step 3: Tile the first fusion step by i0 and j0 (with collapse(2)).

   Step 4: Tile the first column-wise pass (forward i) by j0.

   Step 5: Tile the second column-wise pass (backward i) by j0.

   Step 6: Tile the second fusion step by i0 and j0 (with collapse(2)).

 8. We must be cautious about the min function. We can use:

        int i_end = (i0 + TI) < _PB_W ? (i0 + TI) : _PB_W;

    But we can write: min(i0+TI, _PB_W) -> we can use a macro or a function? We can use:

        #define MIN(a,b) ((a)<(b)?(a):(b))

    But the original code uses _PB_W and _PB_H which are constants? Actually, they are defined in the header. We can use:

        int i_end = (i0 + TI) < _PB_W ? (i0 + TI) : _PB_W;

    However, we are in C99 and we can use the ternary operator.

 9. We must preserve the original semantics. The tile loops must cover the entire range.

 10. We must not change the code outside the scop region.

 11. We must not add data races. The parallel regions are over independent tiles.

 12. We must use OpenMP directives appropriately.

 Let's write the code.

 Note: We are only allowed to change the code between #pragma scop and #pragma endscop.

 We'll define the tile sizes at the beginning of the scop region? Or we can define them as constants. We'll define:

        #define TI 64
        #define TJ 64

    But note: the problem says "choose legal tile sizes", so we can choose any. We choose 64.

 However, the problem says: "Do not add dependencies, data races, invalid OpenMP directives, or unsupported transformations."

 We are confident that the tiling and OpenMP parallelization are safe.

 Let's code accordingly.
#pragma scop  
#define TI 64  
#define TJ 64  
#define MIN(a,b) ((a)<(b)?(a):(b))  

   k = (SCALAR_VAL(1.0)-EXP_FUN(-alpha))*(SCALAR_VAL(1.0)-EXP_FUN(-alpha))/(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*alpha*EXP_FUN(-alpha)-EXP_FUN(SCALAR_VAL(2.0)*alpha));  
   a1 = a5 = k;  
   a2 = a6 = k*EXP_FUN(-alpha)*(alpha-SCALAR_VAL(1.0));  
   a3 = a7 = k*EXP_FUN(-alpha)*(alpha+SCALAR_VAL(1.0));  
   a4 = a8 = -k*EXP_FUN(SCALAR_VAL(-2.0)*alpha);  
   b1 =  POW_FUN(SCALAR_VAL(2.0),-alpha);  
   b2 = -EXP_FUN(SCALAR_VAL(-2.0)*alpha);  
   c1 = c2 = 1;  

   #pragma omp parallel for  
   for (int i0 = 0; i0 < _PB_W; i0 += TI) {  
        int i_end = MIN(i0 + TI, _PB_W);  
        for (int i = i0; i < i_end; i++) {  
            DATA_TYPE ym1 = SCALAR_VAL(0.0);  
            DATA_TYPE ym2 = SCALAR_VAL(0.0);  
            DATA_TYPE xm1 = SCALAR_VAL(0.0);  
            for (int j = 0; j < _PB_H; j++) {  
                y1[i][j] = a1 * imgIn[i][j] + a2 * xm1 + b1 * ym1 + b2 * ym2;  
                xm1 = imgIn[i][j];  
                ym2 = ym1;  
                ym1 = y1[i][j];  
            }  
        }  
    }  

    #pragma omp parallel for  
    for (int i0 = 0; i0 < _PB_W; i0 += TI) {  
        int i_end = MIN(i0 + TI, _PB_W);  
        for (int i = i0; i < i_end; i++) {  
            DATA_TYPE yp1 = SCALAR_VAL(0.0);  
            DATA_TYPE yp2 = SCALAR_VAL(0.0);  
            DATA_TYPE xp1 = SCALAR_VAL(0.0);  
            DATA_TYPE xp2 = SCALAR_VAL(0.0);  
            for (int j = _PB_H - 1; j >= 0; j--) {  
                y2[i][j] = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;  
                xp2 = xp1;  
                xp1 = imgIn[i][j];  
                yp2 = yp1;  
                yp1 = y2[i][j];  
            }  
        }  
    }  

    #pragma omp parallel for collapse(2)  
    for (int i0 = 0; i0 < _PB_W; i0 += TI) {  
        for (int j0 = 0; j0 < _PB_H; j0 += TJ) {  
            int i_end = MIN(i0 + TI, _PB_W);  
            int j_end = MIN(j0 + TJ, _PB_H);  
            for (int i = i0; i < i_end; i++) {  
                for (int j = j0; j < j_end; j++) {  
                    imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);  
                }  
            }  
        }  
    }  

    #pragma omp parallel for  
    for (int j0 = 0; j0 < _PB_H; j0 += TJ) {  
        int j_end = MIN(j0 + TJ, _PB_H);  
        for (int j = j0; j < j_end; j++) {  
            DATA_TYPE tm1 = SCALAR_VAL(0.0);  
            DATA_TYPE ym1 = SCALAR_VAL(0.0);  
            DATA_TYPE ym2 = SCALAR_VAL(0.0);  
            for (int i = 0; i < _PB_W; i++) {  
                y1[i][j] = a5 * imgOut[i][j] + a6 * tm1 + b1 * ym1 + b2 * ym2;  
                tm1 = imgOut[i][j];  
                ym2 = ym1;  
                ym1 = y1[i][j];  
            }  
        }  
    }  

    #pragma omp parallel for  
    for (int j0 = 0; j0 < _PB_H; j0 += TJ) {  
        int j_end = MIN(j0 + TJ, _PB_H);  
        for (int j = j0; j < j_end; j++) {  
            DATA_TYPE tp1 = SCALAR_VAL(0.0);  
            DATA_TYPE tp2 = SCALAR_VAL(0.0);  
            DATA_TYPE yp1 = SCALAR_VAL(0.0);  
            DATA_TYPE yp2 = SCALAR_VAL(0.0);  
            for (int i = _PB_W - 1; i >= 0; i--) {  
                y2[i][j] = a7 * tp1 + a8 * tp2 + b1 * yp1 + b2 * yp2;  
                tp2 = tp1;  
                tp1 = imgOut[i][j];  
                yp2 = yp1;  
                yp1 = y2[i][j];  
            }  
        }  
    }  

    #pragma omp parallel for collapse(2)  
    for (int i0 = 0; i0 < _PB_W; i0 += TI) {  
        for (int j0 = 0; j0 < _PB_H; j0 += TJ) {  
            int i_end = MIN(i0 + TI, _PB_W);  
            int j_end = MIN(j0 + TJ, _PB_H);  
            for (int i = i0; i < i_end; i++) {  
                for (int j = j0; j < j_end; j++) {  
                    imgOut[i][j] = c2 * (y1[i][j] + y2[i][j]);  
                }  
            }  
        }  
    }  

#pragma endscop
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
