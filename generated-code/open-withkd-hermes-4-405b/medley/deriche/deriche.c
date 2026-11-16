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
We are given the Deriche benchmark code. The task is to optimize the code within the '#pragma scop' and '#pragma endscop' markers, focusing on loop transformations (tiling, fusion, skewing) and then OpenMP directives.

The original code has several nested loops. Let's break down the scop region:

1. Initialization of constants (k, a1..a8, b1, b2, c1, c2) - these are scalar assignments and not loop-based, so we leave them as is.

2. First loop nest (i, j) for computing y1:
   for (i=0; i<_PB_W; i++) {
        ym1 = 0.0; ym2 = 0.0; xm1 = 0.0;
        for (j=0; j<_PB_H; j++) {
            y1[i][j] = a1*imgIn[i][j] + a2*xm1 + b1*ym1 + b2*ym2;
            xm1 = imgIn[i][j];
            ym2 = ym1;
            ym1 = y1[i][j];
        }
   }

   This is a 2D loop. The inner loop (j) has a recurrence: ym1 and ym2 are carried along j. So the inner loop is sequential. The outer loop (i) is independent.

3. Second loop nest (i, j) for computing y2 (backward in j):
   for (i=0; i<_PB_W; i++) {
        yp1 = 0.0; yp2 = 0.0; xp1 = 0.0; xp2 = 0.0;
        for (j=_PB_H-1; j>=0; j--) {
            y2[i][j] = a3*xp1 + a4*xp2 + b1*yp1 + b2*yp2;
            xp2 = xp1;
            xp1 = imgIn[i][j];
            yp2 = yp1;
            yp1 = y2[i][j];
        }
   }

   Similarly, the inner loop (j) is sequential (backward) and the outer loop (i) is independent.

4. Third loop nest (i, j) for computing imgOut from y1 and y2:
   for (i=0; i<_PB_W; i++)
        for (j=0; j<_PB_H; j++) {
            imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);
        }

   This is a simple element-wise operation. Both loops are independent.

5. Fourth loop nest (j, i) for computing y1 again (using imgOut as input, and updating y1):
   for (j=0; j<_PB_H; j++) {
        tm1 = 0.0; ym1 = 0.0; ym2 = 0.0;
        for (i=0; i<_PB_W; i++) {
            y1[i][j] = a5*imgOut[i][j] + a6*tm1 + b1*ym1 + b2*ym2;
            tm1 = imgOut[i][j];
            ym2 = ym1;
            ym1 = y1 [i][j];
        }
   }

   Now the outer loop is j and inner loop is i. The inner loop (i) has a recurrence (carried along i). The outer loop (j) is independent.

6. Fifth loop nest (j, i) for computing y2 again (backward in i):
   for (j=0; j<_PB_H; j++) {
        tp1 = 0.0; tp2 = 0.0; yp1 = 0.0; yp2 = 0.0;
        for (i=_PB_W-1; i>=0; i--) {
            y2[i][j] = a7*tp1 + a8*tp2 + b1*yp1 + b2*yp2;
            tp2 = tp1;
            tp1 = imgOut[i][j];
            yp2 = yp1;
            yp1 = y2[i][j];
        }
   }

   The inner loop (i) is sequential (backward) and the outer loop (j) is independent.

7. Sixth loop nest (i, j) for the final imgOut:
   for (i=0; i<_PB_W; i++)
        for (j=0; j<_PB_H; j++)
            imgOut[i][j] = c2*(y1[i][j] + y2[i][j]);

   Again, independent.

Loop Transformation Analysis:

1. Tiling: 
   - The first two loop nests (for y1 and y2) are over (i, j) and have inner sequential loops. Tiling the outer loop (i) might be beneficial for cache locality? But note that the inner loop is sequential and the recurrence is along j. Tiling the j dimension would break the recurrence. So we cannot tile the inner loop. However, we can tile the outer loop (i) to improve locality for the imgIn array? But note that the inner loop accesses imgIn[i][j] in order, so it is already contiguous. Tiling the outer loop might not help much because each i is independent and the inner loop is sequential. Similarly for the other loops.

   - The third and sixth loops (element-wise) are perfect for tiling. We can tile both i and j to improve cache locality.

   - The fourth and fifth loops: the outer loop is j and inner loop is i. The inner loop has a recurrence along i. We cannot tile the inner loop. We can tile the outer loop (j) to improve locality? But note that the inner loop accesses imgOut[i][j] in order (contiguous in i). Tiling j might help if we are accessing multiple j's and the same i's? But the recurrence is along i, so we cannot change the order of i. Tiling j might be beneficial for the imgOut array? However, the inner loop is sequential and we are writing to y1[i][j] and y2[i][j] which are 2D arrays. Tiling j might help by reusing the same j tile for multiple i's? But note that the inner loop is over i and we are writing to a fixed j. So tiling j would break the inner loop into tiles of j, but the inner loop is over i. This might not be straightforward.

   Given the complexity and the fact that the inner loops are sequential, we might avoid tiling for the first, second, fourth, and fifth loops. We will focus on tiling the third and sixth loops.

2. Loop Fusion:
   - The third loop (imgOut = c1*(y1+y2)) and the sixth loop (imgOut = c2*(y1+y2)) are similar. However, they are separated by the fourth and fifth loops that update y1 and y2. So we cannot fuse the third and sixth.

   - The first and second loops: they both use imgIn and produce y1 and y2. They are independent and can be fused? But note they are over the same i and j, but the second loop runs j backwards. We cannot fuse them because the j order is opposite.

   - The fourth and fifth loops: they both use imgOut and update y1 and y2. They are independent and can be fused? But again, the fifth loop runs i backwards. So we cannot fuse.

   Therefore, fusion opportunities are limited.

3. Loop Skewing:
   - The first loop: the inner loop has a recurrence: y1[i][j] depends on ym1 and ym2 which are from the previous j. This is a simple recurrence and skewing might not help because the inner loop is already sequential and we cannot parallelize it.

   - Similarly for the other loops with recurrences.

   So skewing is not beneficial.

OpenMP Directives:

We will focus on parallelizing the outer loops that are independent.

- For the first loop: outer loop i is independent -> parallel for i.
- For the second loop: outer loop i is independent -> parallel for i.
- For the third loop: both i and j are independent -> we can parallelize i and j. We can use collapse(2) to parallelize both loops together.
- For the fourth loop: outer loop j is independent -> parallel for j.
- For the fifth loop: outer loop j is independent -> parallel for j.
- For the sixth loop: both i and j are independent -> collapse(2).

Additionally, we can consider:
- Using a single parallel region for multiple loops to avoid fork/join overhead? But note that the loops are separated by dependencies (the third loop uses y1 and y2 from the first two, then the fourth and fifth use imgOut from the third, and the sixth uses y1 and y2 from the fourth and fifth). So we must have barriers between these phases. We can use one parallel region and then use multiple omp for with implicit barriers (or nowait if safe, but here we have dependencies so we need barriers).

We'll use:
#pragma omp parallel
{
   #pragma omp for
   for (i=0; i<_PB_W; i++) { ... }   // first loop

   #pragma omp for
   for (i=0; i<_PB_W; i++) { ... }   // second loop

   #pragma omp for collapse(2)
   for (i=0; i<_PB_W; i++)
      for (j=0; j<_PB_H; j++) { ... }   // third loop

   #pragma omp for
   for (j=0; j<_PB_H; j++) { ... }   // fourth loop

   #pragma omp for
   for (j=0; j<_PB_H; j++) { ... }   // fifth loop

   #pragma omp for collapse(2)
   for (i=0; i<_PB_W; i++)
      for (j=0; j<_PB_H; j++) { ... }   // sixth loop
}

But note: the first two loops are over i, and the next two over j. We can also consider using collapse(2) for the first two? But the inner loop is sequential so we cannot collapse. Similarly for the fourth and fifth.

We'll also consider scheduling: static is usually good for regular work.

For the element-wise loops (third and sixth) we can use SIMD. We can add simd to the inner loop? But note that the inner loop is j and we are accessing contiguous j. We can use #pragma omp simd on the inner loop? Or we can use #pragma omp for simd on the outer loop? Actually, we are using collapse(2) so the entire 2D loop is parallelized. We can add simd to the inner loop? But the compiler might vectorize automatically. We can also use #pragma omp simd on the inner loop of the third and sixth loops? However, the third and sixth loops are already being parallelized by the outer loop. We can try to use SIMD for the inner loop of the third and sixth loops? But note that the inner loop is j and we are accessing contiguous j, so it is vector-friendly.

Alternatively, we can use:
   #pragma omp for collapse(2) simd
But note: the collapse(2) already breaks the loop into chunks. We can try to use SIMD on the inner loop without collapse? But we are using collapse(2) to parallelize both loops. We can also use:
   #pragma omp parallel for collapse(2) simd
But that might be too aggressive? We'll try to use SIMD for the inner loop of the third and sixth loops by adding a simd directive to the inner loop? However, the inner loop is the j loop and we are already parallelizing the i loop. We can do:

   #pragma omp for collapse(2)
   for (i=0; i<_PB_W; i++)
      #pragma omp simd
      for (j=0; j<_PB_H; j++) { ... }

But note: the collapse(2) already distributes the (i,j) pairs. The inner loop is then a chunk of j's. We can let the compiler vectorize the inner loop without an explicit simd? Or we can add the simd to the inner loop to ensure vectorization.

We'll add simd to the inner loops of the third and sixth loops.

Also, we can consider tiling the third and sixth loops? We decided to tile them. We can do:

   for (i0=0; i0<_PB_W; i0+=TI)
      for (j0=0; j0<_PB_H; j0+=TJ)
         for (i=i0; i<min(i0+TI, _PB_W); i++)
            for (j=j0; j<min(j0+TJ, _PB_H); j++)
               imgOut[i][j] = ...;

But note: we are going to parallelize the outer loops. We can tile and then parallelize the outer tile loops? We can do:

   #pragma omp for collapse(2)
   for (i0=0; i0<_PB_W; i0+=TI)
      for (j0=0; j0<_PB_H; j0+=TJ)
         for (i=i0; i<min(i0+TI, _PB_W); i++)
            for (j=j0; j<min(j0+TJ, _PB_H); j++)
               imgOut[i][j] = ...;

But note: the inner two loops are over a tile. We can also use SIMD on the inner j loop? We can add #pragma omp simd to the inner j loop.

We'll choose a tile size. Since we don't know the exact cache, we can choose 64x64? But note that the data type is DATA_TYPE (float or double). We'll use 64 for both dimensions.

But note: the third and sixth loops are independent and we are tiling them. We'll do the same for both.

However, the problem says: "Do not change observable behavior". Tiling with min/max is safe.

But note: the original code does not have tiling. We are adding it. We must be cautious.

Alternatively, we can avoid tiling and rely on the collapse(2) and SIMD? But the problem says to prioritize tiling. So we will tile the third and sixth loops.

But note: the third loop is:
   for (i=0; i<_PB_W; i++)
      for (j=0; j<_PB_H; j++)
         imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);

We can tile it as:

   #define TI 64
   #define TJ 64

   for (i0=0; i0<_PB_W; i0+=TI)
      for (j0=0; j0<_PB_H; j0+=TJ)
         for (i=i0; i<min(i0+TI, _PB_W); i++)
            #pragma omp simd
            for (j=j0; j<min(j0+TJ, _PB_H); j++)
               imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);

But note: we are going to parallelize the outer two loops (i0, j0) with collapse(2). So:

   #pragma omp for collapse(2)
   for (i0=0; i0<_PB_W; i0+=TI)
      for (j0=0; j0<_PB_H; j0+=TJ)
         for (i=i0; i<min(i0+TI, _PB_W); i++)
            #pragma omp simd
            for (j=j0; j<min(j0+TJ, _PB_H); j++)
               imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);

Similarly for the sixth loop.

But note: the min function might be expensive? We can precompute the bounds? But the compiler might optimize it.

We'll do the same for the sixth loop.

Now, let's write the optimized code within the scop region.

We'll define the tile sizes. We choose 64 for both dimensions.

Note: We must not change the original semantics. The tiling with min is safe.

We'll also add the OpenMP directives as described.

We'll use a single parallel region for the entire scop to avoid multiple fork/join.

We'll use default(none) and specify the shared variables? But note: the variables are shared by default. We can use default(none) and then list shared and private variables.

But note: the variables in the scop are:
   k, a1..a8, b1, b2, c1, c2: these are scalars and shared.
   The arrays: imgIn, imgOut, y1, y2: shared.
   The loop variables: i, j, and the temporaries (ym1, etc.) are private.

We'll use:
   #pragma omp parallel default(none) \
        shared(imgIn, imgOut, y1, y2, _PB_W, _PB_H, a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, c1, c2) \
        private(i, j, ym1, ym2, xm1, yp1, yp2, xp1, xp2, tm1, tp1, tp2)

But note: the temporaries in the first loop: ym1, ym2, xm1 are declared inside the i loop? Actually, they are declared at the beginning of the i loop. We can make them private per thread? Yes, because each thread will have its own i.

Similarly for the other temporaries.

We'll do:

   #pragma omp parallel default(none) \
        shared(imgIn, imgOut, y1, y2, _PB_W, _PB_H, a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, c1, c2)
   {
      // First loop
      #pragma omp for schedule(static)
      for (i=0; i<_PB_W; i++) {
          // ... same as original, but note: the temporaries are private to the thread because they are declared inside the loop? Actually, we moved the declaration of temporaries to the scop? No, they are declared at the beginning of the scop? Actually, they are declared at the beginning of the function. We must make them private.

We have a problem: the temporaries (ym1, etc.) are declared at the beginning of the scop. We cannot make them private per thread because they are shared by default? We must declare them as private in the OpenMP directive? But note: they are used in different loops. We can declare them as private in the parallel directive? But then they are private for the entire parallel region. We can do that.

Alternatively, we can declare them inside the parallel region? But that would change the scope. We can move the declarations of the temporaries to inside the parallel region? But note: the original code declares them at the beginning of the scop. We are allowed to restructure within the scop? The problem says: "You MAY restructure loops and adjust associated C statements within the marked regions when needed for the priority loop transformations below."

We can move the declarations of the temporaries to inside the parallel region? But note: the temporaries are used in multiple loops. We can declare them at the beginning of the parallel block? Then they are private to each thread.

We'll do:

   #pragma omp parallel default(none) \
        shared(imgIn, imgOut, y1, y2, _PB_W, _PB_H, a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, c1, c2) \
        private(i, j)
   {
      // Declare temporaries here so they are private to each thread
      DATA_TYPE ym1, ym2, xm1, tm1, ym1_2, ym2_2, xm1_2; // for the first loop? Actually, we have multiple sets of temporaries for different loops.

But note: the temporaries are used in different loops and have different names? Actually, the same names are reused? In the original code, the same variable names (ym1, etc.) are reused in different loops. We can declare them at the beginning of the parallel block and then use them in the loops? But note: the loops are executed in sequence and the temporaries are reset at the beginning of each loop. So it is safe.

We'll declare all the temporaries that are used in the loops at the beginning of the parallel block. They are: 
   ym1, ym2, xm1, tm1, ym1_2, ym2_2, xm1_2? Actually, the original code uses:
   In the first loop: ym1, ym2, xm1.
   In the second loop: yp1, yp2, xp1, xp2.
   In the fourth loop: tm1, ym1, ym2.
   In the fifth loop: tp1, tp2, yp1, yp2.

We can declare:
   DATA_TYPE ym1, ym2, xm1, yp1, yp2, xp1, xp2, tm1, tp1, tp2;

But note: the same name 'ym1' is used in the first and fourth loops? Yes, but they are in different scopes? Actually, in the original code, they are in the same function scope. We are moving them to the parallel block. We can use the same names? Yes, because the loops are executed in sequence and the variables are reset at the beginning of each loop.

So we'll declare all these temporaries at the beginning of the parallel block.

Now, let's write the code.

We'll define the tile sizes as macros? But the problem says not to change code outside the scop. We can define them inside the scop? We'll define:

   #define TI 64
   #define TJ 64

But note: we are only allowed to change within the scop. We can define them at the beginning of the scop? But the scop starts at '#pragma scop'. We can put the #define after the #pragma scop? Yes.

Alternatively, we can use constants? We'll use:

   const int TI = 64;
   const int TJ = 64;

But note: we are in C, and we are using min. We can use:

   int i0, j0;
   for (i0=0; i0<_PB_W; i0+=TI)
      for (j0=0; j0<_PB_H; j0+=TJ)
         for (i=i0; i< (i0+TI < _PB_W ? i0+TI : _PB_W); i++)
            for (j=j0; j< (j0+TJ < _PB_H ? j0+TJ : _PB_H); j++)

But we can use a macro for min? Or we can use:

   int i_end = (i0+TI < _PB_W) ? i0+TI : _PB_W;
   for (i=i0; i<i_end; i++)

But that would require extra variables. We'll use the ternary operator in the loop condition? It might be inefficient? We can rely on the compiler to optimize.

Alternatively, we can use:

   for (i=i0; i<min(i0+TI, _PB_W); i++)

But min is not standard in C for integers? We can use a macro? Or we can use:

   #include <stdio.h> // already included, but we don't have stdlib for min? Actually, we can use a macro.

We'll define:

   #define MIN(a,b) ((a)<(b)?(a):(b))

But note: we are only allowed to change within the scop. We can define MIN at the beginning of the scop? Yes.

So:

   #pragma scop
   #define MIN(a,b) ((a)<(b)?(a):(b))
   ... rest of the code ...

But note: the original code does not have this macro. We are adding it. It is safe.

Alternatively, we can avoid the macro and use the ternary operator in the loop condition? We'll do that to avoid a macro.

We'll write:

   for (i = i0; i < (i0+TI < _PB_W ? i0+TI : _PB_W); i++)

But that is a bit messy. We'll use:

   int i_end = i0 + TI;
   if (i_end > _PB_W) i_end = _PB_W;
   for (i = i0; i < i_end; i++)

But that adds extra variables. We'll do it for clarity and hope the compiler optimizes.

We'll do the same for j.

Now, let's write the entire scop region.

Note: We must preserve the original variable names and semantics.

We'll do:

   #pragma scop
   // Define tile sizes
   const int TI = 64;
   const int TJ = 64;

   // The scalar assignments remain the same
   k = ...;
   a1 = ...;
   ... up to c2 = 1;

   #pragma omp parallel default(none) \
        shared(imgIn, imgOut, y1, y2, _PB_W, _PB_H, a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, c1, c2)
   {
      // Declare private temporaries for the entire parallel region
      DATA_TYPE ym1, ym2, xm1, yp1, yp2, xp1, xp2, tm1, tp1, tp2;

      // First loop: compute y1
      #pragma omp for schedule(static)
      for (i=0; i<_PB_W; i++) {
          ym1 = SCALAR_VAL(0.0);
          ym2 = SCALAR_VAL(0.0);
          xm1 = SCALAR_VAL(0.0);
          for (j=0; j<_PB_H; j++) {
              y1[i][j] = a1*imgIn[i][j] + a2*xm1 + b1*ym1 + b2*ym2;
              xm1 = imgIn[i][j];
              ym2 = ym1;
              ym1 = y1[i][j];
          }
      }

      // Second loop: compute y2
      #pragma omp for schedule(static)
      for (i=0; i<_PB_W; i++) {
          yp1 = SCALAR_VAL(0.0);
          yp2 = SCALAR_VAL(0.0);
          xp1 = SCALAR_VAL(0.0);
          xp2 = SCALAR_VAL(0.0);
          for (j=_PB_H-1; j>=0; j--) {
              y2[i][j] = a3*xp1 + a4*xp2 + b1*yp1 + b2*yp2;
              xp2 = xp1;
              xp1 = imgIn[i][j];
              yp2 = yp1;
              yp1 = y2[i][j];
          }
      }

      // Third loop: compute imgOut from y1 and y2 (with tiling and SIMD)
      #pragma omp for collapse(2) schedule(static)
      for (int i0=0; i0<_PB_W; i0+=TI) {
          int i_end = i0 + TI;
          if (i_end > _PB_W) i_end = _PB_W;
          for (int j0=0; j0<_PB_H; j0+=TJ) {
              int j_end = j0 + TJ;
              if (j_end > _PB_H) j_end = _PB_H;
              for (i=i0; i<i_end; i++) {
                  #pragma omp simd
                  for (j=j0; j<j_end; j++) {
                      imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);
                  }
              }
          }
      }

      // Fourth loop: compute y1 again (using imgOut)
      #pragma omp for schedule(static)
      for (j=0; j<_PB_H; j++) {
          tm1 = SCALAR_VAL(0.0);
          ym1 = SCALAR_VAL(0.0);
          ym2 = SCALAR_VAL(0.0);
          for (i=0; i<_PB_W; i++) {
              y1[i][j] = a5*imgOut[i][j] + a6*tm1 + b1*ym1 + b2*ym2;
              tm1 = imgOut[i][j];
              ym2 = ym1;
              ym1 = y1[i][j];
          }
      }

      // Fifth loop: compute y2 again (using imgOut)
      #pragma omp for schedule(static)
      for (j=0; j<_PB_H; j++) {
          tp1 = SCALAR_VAL(0.0);
          tp2 = SCALAR_VAL(0.0);
          yp1 = SCALAR_VAL(0.0);
          yp2 = SCALAR_VAL(0.0);
          for (i=_PB_W-1; i>=0; i--) {
              y2[i][j] = a7*tp1 + a8*tp2 + b1*yp1 + b2*yp2;
              tp2 = tp1;
              tp1 = imgOut[i][j];
              yp2 = yp1;
              yp1 = y2[i][j];
          }
      }

      // Sixth loop: final imgOut (with tiling and SIMD)
      #pragma omp for collapse(2) schedule(static)
      for (int i0=0; i0<_PB_W; i0+=TI) {
          int i_end = i0 + TI;
          if (i_end > _PB_W) i_end = _PB_W;
          for (int j0=0; j0<_PB_H; j0+=TJ) {
              int j_end = j0 + TJ;
              if (j_end > _PB_H) j_end = _PB_H;
              for (i=i0; i<i_end; i++) {
                  #pragma omp simd
                  for (j=j0; j<j_end; j++) {
                      imgOut[i][j] = c2*(y1[i][j] + y2[i][j]);
                  }
              }
          }
      }
   } // end parallel

   #pragma endscop

But note: the original code uses _PB_W and _PB_H. We are using them.

We must be cautious: the variables i and j are declared in the scop? Actually, they are declared at the beginning of the function. We are using them in the parallel region. We declared them as private in the parallel directive? We did: private(i, j). But note: we are also using i0, j0, i_end, j_end. We must declare them as private? We can declare them inside the parallel block? We did: we declared them in the for loops? Actually, we declared i0, j0, i_end, j_end inside the parallel block? We can declare them at the beginning of the parallel block? Or we can let them be private by default? Since we used default(none), we must declare them as private? We can add them to the private clause? Or we can declare them inside the parallel block? We'll declare them inside the parallel block? Actually, we are declaring them in the for loops? That is inside the parallel block. So they are private to the thread.

But note: the variables i and j are used in multiple loops. We declared them as private in the parallel directive? Yes.

We'll write the code accordingly.

One more thing: the original code uses SCALAR_VAL(0.0) and EXP_FUN, etc. We are preserving that.

Let's write the final code.
#pragma scop
  k = (SCALAR_VAL(1.0)-EXP_FUN(-alpha))*(SCALAR_VAL(1.0)-EXP_FUN(-alpha))/(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*alpha*EXP_FUN(-alpha)-EXP_FUN(SCALAR_VAL(2.0)*alpha));
  a1 = a5 = k;
  a2 = a6 = k*EXP_FUN(-alpha)*(alpha-SCALAR_VAL(1.0));
  a3 = a7 = k*EXP_FUN(-alpha)*(alpha+SCALAR_VAL(1.0));
  a4 = a8 = -k*EXP_FUN(SCALAR_VAL(-2.0)*alpha);
  b1 =  POW_FUN(SCALAR_VAL(2.0),-alpha);
  b2 = -EXP_FUN(SCALAR_VAL(-2.0)*alpha);
  c1 = c2 = 1;

  const int TI = 64;
  const int TJ = 64;
  
  #pragma omp parallel default(none) \
    shared(imgIn, imgOut, y1, y2, _PB_W, _PB_H, a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, c1, c2) \
    private(i, j, ym1, ym2, xm1, yp1, yp2, xp1, xp2, tm1, tp1, tp2)
  {
    DATA_TYPE ym1, ym2, xm1, yp1, yp2, xp1, xp2, tm1, tp1, tp2;
    
    #pragma omp for schedule(static)
    for (i=0; i<_PB_W; i++) {
      ym1 = SCALAR_VAL(0.0);
      ym2 = SCALAR_VAL(0.0);
      xm1 = SCALAR_VAL(0.0);
      for (j=0; j<_PB_H; j++) {
        y1[i][j] = a1*imgIn[i][j] + a2*xm1 + b1*ym1 + b2*ym2;
        xm1 = imgIn[i][j];
        ym2 = ym1;
        ym1 = y1[i][j];
      }
    }
    
    #pragma omp for schedule(static)
    for (i=0; i<_PB_W; i++) {
      yp1 = SCALAR_VAL(0.0);
      yp2 = SCALAR_VAL(0.0);
      xp1 = SCALAR_VAL(0.0);
      xp2 = SCALAR_VAL(0.0);
      for (j=_PB_H-1; j>=0; j--) {
        y2[i][j] = a3*xp1 + a4*xp2 + b1*yp1 + b2*yp2;
        xp2 = xp1;
        xp1 = imgIn[i][j];
        yp2 = yp1;
        yp1 = y2[i][j];
      }
    }
    
    #pragma omp for collapse(2) schedule(static)
    for (int i0=0; i0<_PB_W; i0+=TI) {
      int i_end = (i0 + TI < _PB_W) ? i0 + TI : _PB_W;
      for (int j0=0; j0<_PB_H; j0+=TJ) {
        int j_end = (j0 + TJ < _PB_H) ? j0 + TJ : _PB_H;
        for (i=i0; i<i_end; i++) {
          #pragma omp simd
          for (j=j0; j<j_end; j++) {
            imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);
          }
        }
      }
    }
    
    #pragma omp for schedule(static)
    for (j=0; j<_PB_H; j++) {
      tm1 = SCALAR_VAL(0.0);
      ym1 = SCALAR_VAL(0.0);
      ym2 = SCALAR_VAL(0.0);
      for (i=0; i<_PB_W; i++) {
        y1[i][j] = a5*imgOut[i][j] + a6*tm1 + b1*ym1 + b2*ym2;
        tm1 = imgOut[i][j];
        ym2 = ym1;
        ym1 = y1[i][j];
      }
    }
    
    #pragma omp for schedule(static)
    for (j=0; j<_PB_H; j++) {
      tp1 = SCALAR_VAL(0.0);
      tp2 = SCALAR_VAL(0.0);
      yp1 = SCALAR_VAL(0.0);
      yp2 = SCALAR_VAL(0.0);
      for (i=_PB_W-1; i>=0; i--) {
        y2[i][j] = a7*tp1 + a8*tp2 + b1*yp1 + b2*yp2;
        tp2 = tp1;
        tp1 = imgOut[i][j];
        yp2 = yp1;
        yp1 = y2[i][j];
      }
    }
    
    #pragma omp for collapse(2) schedule(static)
    for (int i0=0; i0<_PB_W; i0+=TI) {
      int i_end = (i0 + TI < _PB_W) ? i0 + TI : _PB_W;
      for (int j0=0; j0<_PB_H; j0+=TJ) {
        int j_end = (j0 + TJ < _PB_H) ? j0 + TJ : _PB_H;
        for (i=i0; i<i_end; i++) {
          #pragma omp simd
          for (j=j0; j<j_end; j++) {
            imgOut[i][j] = c2*(y1[i][j] + y2[i][j]);
          }
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
