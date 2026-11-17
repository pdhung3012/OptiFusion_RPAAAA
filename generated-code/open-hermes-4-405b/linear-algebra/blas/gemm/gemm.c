We are given a matrix multiplication kernel (gemm) that computes C = alpha * A * B + beta * C.
 The original code has two nested loops: 
   First, it scales C by beta: for each i, for each j, C[i][j] *= beta.
   Then, for each i, for each k, for each j, it does: C[i][j] += alpha * A[i][k] * B[k][j].

 We are to optimize the code within the '#pragma scop' and '#pragma endscop' markers.

 Step 1: Loop Transformations (highest priority)
   a) Loop Tiling: We can tile the i, j, k loops to improve cache locality. Since the innermost loop is j, and we are accessing A[i][k] and B[k][j], we can tile the i and j loops to keep a block of C in cache. Also, we can tile k to keep a block of A and B in cache.

   b) Loop Fusion: The original code has two separate loops: one for scaling by beta and one for the multiplication. We can fuse these two loops to avoid writing and then reading C twice. However, note that the scaling loop is over i and j, and the multiplication loop is over i, k, j. We can combine the scaling with the first iteration of the k loop? But note: the scaling must be done exactly once per element. Alternatively, we can do the scaling in the same loop as the multiplication by initializing the accumulation with beta * C[i][j] and then adding the products. However, the original code does:
        C[i][j] = beta * C[i][j] + ... (over k)
   So we can rewrite the inner part as:
        C[i][j] = beta * C[i][j] + alpha * (sum over k of A[i][k] * B[k][j])

   But note: the original code does the scaling in a separate loop. We can fuse by doing:

        for (i ...) {
            for (j ...) {
                C[i][j] *= beta;   // This is the scaling
                for (k ...) {
                    C[i][j] += alpha * A[i][k] * B[k][j];
                }
            }
        }

   However, this changes the order: originally the scaling was done for all j for a fixed i, then the k loop. Now we are doing the scaling for each (i,j) and then immediately the k loop for that (i,j). This is acceptable because the scaling is independent of k.

   But note: the original code does the scaling for the entire matrix C first, then the multiplication. The fused version does the scaling and then the multiplication for each element. This is semantically equivalent.

   However, we must be cautious: the original code does the scaling in a separate loop that is outside the k loop. The fused version does the scaling inside the i and j loops and then the k loop. This is the same.

   But wait: the original code does:
        for i: 
            for j: C[i][j] *= beta
            for k: 
                for j: ... 

   The fused version would be:
        for i:
            for j:
                C[i][j] *= beta;
                for k:
                    ... 

   This is acceptable and actually improves locality because we are working on the same C[i][j] for the scaling and then the entire k loop.

   However, note that the original code does the scaling for all j for a fixed i, then the k loop for that i. The fused version does the same: for a fixed i, we do j and then k. So the order of memory accesses is similar.

   But the fused version might be better because it avoids writing and then reading C[i][j] again (the scaling writes, then the multiplication reads and writes again). In the fused version, we write C[i][j] (scaling) and then we read and write it in the k loop. However, the scaling is done once and then the k loop does multiple updates. This is the same as the original.

   However, we can do even better: we can combine the scaling and the multiplication into one expression. We can write:

        temp = beta * C[i][j];
        for (k=0; k<NK; k++)
            temp += alpha * A[i][k] * B[k][j];
        C[i][j] = temp;

   This avoids the multiple updates to C[i][j] in the inner loop. This is a common optimization for register reuse.

   So we can do:
        for (i ...) {
            for (j ...) {
                DATA_TYPE temp = beta * C[i][j];
                for (k ...) {
                    temp += alpha * A[i][k] * B[k][j];
                }
                C[i][j] = temp;
            }
        }

   This is a form of loop fusion and also introduces a temporary variable to hold the intermediate result.

   c) Loop Skewing: There are no carried dependencies that require skewing in this code. The inner loops are parallel.

 Step 2: OpenMP Directives

   We can parallelize the outer loop (i). We can also consider collapsing the i and j loops? But note: the inner loop is k, and we are using a temporary variable per (i,j). So we can parallelize the i loop and then the j loop can be done in parallel? Actually, the j loop is independent for each i. So we can do:

        #pragma omp parallel for private(j,k) collapse(2)  ... 

   But note: the collapse(2) would collapse the i and j loops. However, the inner loop is k and we have a temporary variable for each (i,j). So collapsing i and j is safe.

   Alternatively, we can do:

        #pragma omp parallel for private(j,k) 
        for (i=0; i<NI; i++) {
            for (j=0; j<NJ; j++) {
                ... 
            }
        }

   But without collapse, we are parallelizing only the i loop. The j loop is sequential. We can also try to parallelize the j loop? But note: the j loop is inside i, and we are using a temporary per (i,j). So we can collapse the two loops to have more parallelism.

   However, the problem sizes might be large, and collapsing two loops might create too many iterations. We can use a schedule and also consider the tile sizes.

   Also, we can use SIMD for the innermost loop (k). But note: the innermost loop is k, and we are accumulating into a scalar (temp). This is a reduction. We can use a SIMD reduction? But the inner loop is over k and we are doing a dot product. We can use SIMD for the k loop? However, the k loop is not easily vectorized because of the reduction. But modern compilers can vectorize reductions.

   Alternatively, we can interchange the j and k loops? Then the inner loop would be j, which is a vector operation. But note: the original order is i, then j (for scaling) and then k and j. In the fused version we have i, j, k. If we interchange j and k, we get i, k, j. Then the inner loop is j, which is contiguous in memory for B[k][j] and C[i][j]? Actually, C[i][j] is row-major, so j is contiguous. B[k][j] is also row-major. So we can do:

        for (i ...) {
            for (k ...) {
                for (j ...) {
                    temp = ... 
                }
            }
        }

   But note: we have a temporary per (i,j). If we interchange k and j, then we break the accumulation for a fixed (i,j). We would have to reorganize the temporary. We can do:

        for (i ...) {
            for (j ...) {
                temp = beta * C[i][j];
            }
            for (k ...) {
                for (j ...) {
                    temp += ...   // but now temp is for a fixed (i,j) and we are iterating k? 
                }
            }
        }

   This doesn't work because the temporary is per (i,j) and we are iterating k in the outer loop. We would have to store the entire matrix of temporaries? That would be too expensive.

   So we stick to the order: i, j, k.

   However, we can tile the i, j, k loops to improve cache locality. We can do:

        for (ii = 0; ii < NI; ii += TI) {
            for (jj = 0; jj < NJ; jj += TJ) {
                for (kk = 0; kk < NK; kk += TK) {
                    for (i = ii; i < ii+TI; i++) {
                        for (j = jj; j < jj+TJ; j++) {
                            // We cannot do the entire k loop in one go because we are tiling k as well.
                            // We break the k loop into tiles of size TK.
                            // We need to accumulate over k, so we do:
                            //   temp = beta * C[i][j]   (only for the first tile in k)
                            //   then for each tile in k: 
                            //        for k in [kk, kk+TK): temp += alpha * A[i][k] * B[k][j]
                            //   then at the end: C[i][j] = temp
                            // But note: the scaling by beta is only done once at the beginning.

                            // So we do:
                            DATA_TYPE temp = beta * C[i][j];
                            for (k = 0; k < NK; k++) {   // We are not tiling k in the inner loop? 
                                temp += alpha * A[i][k] * B[k][j];
                            }
                            C[i][j] = temp;
                        }
                    }
                }
            }
        }

   But wait, we are tiling i and j, but the k loop is not tiled in the above. We can tile k as well to improve cache locality for A and B. We can do:

        for (ii = 0; ii < NI; ii += TI) {
            for (jj = 0; jj < NJ; jj += TJ) {
                for (kk = 0; kk < NK; kk += TK) {
                    for (i = ii; i < ii+TI; i++) {
                        for (j = jj; j < jj+TJ; j++) {
                            // We are going to do a tile of k: from kk to kk+TK
                            // But note: we are accumulating over k, so we must do all k. We break the k loop into tiles.
                            // We can do:
                            //   if (kk == 0) {
                            //       temp = beta * C[i][j];
                            //   } else {
                            //       temp = C[i][j];   // because we are accumulating and we have already done the previous tiles
                            //   }
                            //   Then for k in [kk, kk+TK): temp += alpha * A[i][k] * B[k][j]
                            //   Then set C[i][j] = temp (if we are at the last tile) or just update for the next tile?
                            // But note: we cannot update C[i][j] until we have done all k. So we must do the entire k loop without updating C until the end.

                            // Therefore, we cannot break the k loop into tiles without storing the entire tile of C in a temporary array? 
                            // Alternatively, we can do the entire k loop without tiling? 

                            // However, we can restructure: we do the entire k loop for a tile of i and j. We can do:

                            //   for (i in tile_ii, j in tile_jj): 
                            //        temp = beta * C[i][j]
                            //        for (k=0; k<NK; k++) {
                            //            temp += alpha * A[i][k] * B[k][j]
                            //        }
                            //        C[i][j] = temp

                            // This is the same as the fused version without tiling k. But we are tiling i and j to have better spatial locality for C and also for A and B? 

                            // Actually, by tiling i and j, we are working on a block of C. We can also tile k to improve the reuse of A and B. We can do:

                            //   for (i in tile_ii) {
                            //        for (k in tile_kk) {
                            //            for (j in tile_jj) {
                            //                ... 
                            //            }
                            //        }
                            //   }

                            // But then we break the accumulation for a fixed (i,j). We would have to store the entire block of C in a temporary array? 

                            // Given the complexity, and since the problem does not specify the tile sizes, we might choose to tile only i and j and leave k as the innermost loop. This is a common tiling for matrix multiplication.

                            // We'll do:

                            DATA_TYPE temp = beta * C[i][j];
                            for (k = 0; k < _PB_NK; k++) {
                                temp += alpha * A[i][k] * B[k][j];
                            }
                            C[i][j] = temp;
                        }
                    }
                }
            }
        }

   But note: the above tiling of i and j does not change the inner loop order. It just breaks the i and j loops into tiles. This can help by keeping a tile of C in cache. Also, we are accessing A[i][k] for a fixed i and a range of k, and B[k][j] for a fixed j and a range of k. The tiling of j might help with B because we are accessing B[k][j] for consecutive j? But note: B is stored in row-major, so consecutive j is contiguous. So a tile in j will be contiguous. Similarly, A[i][k] for consecutive k is contiguous. So tiling i and j is beneficial.

   However, we must choose tile sizes. We can choose TI and TJ to be small enough to fit a tile of C (TI * TJ) in L1 cache, and also to have good reuse of A and B. But without specific cache sizes, we can choose arbitrary tile sizes, say 64x64.

   But note: the problem says "choose legal tile sizes". We can choose:

        #define TI 64
        #define TJ 64

   However, we must ensure that the tile sizes divide the problem sizes? Not necessarily, we can handle the boundaries.

   Alternatively, we can use a dynamic tile size? But the code must be fixed. We'll use fixed tile sizes.

   But note: the original code uses _PB_NI, _PB_NJ, _PB_NK. We must use these.

   We can do:

        int i, j, k, ii, jj;

        for (ii = 0; ii < _PB_NI; ii += TI) {
            for (jj = 0; jj < _PB_NJ; jj += TJ) {
                int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
                int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
                for (i = ii; i < i_end; i++) {
                    for (j = jj; j < j_end; j++) {
                        DATA_TYPE temp = beta * C[i][j];
                        for (k = 0; k < _PB_NK; k++) {
                            temp += alpha * A[i][k] * B[k][j];
                        }
                        C[i][j] = temp;
                    }
                }
            }
        }

   This is the tiled and fused version.

 Step 3: OpenMP

   We can parallelize the outer loop (ii). We can do:

        #pragma omp parallel for private(jj, i, j, k, i_end, j_end) 
        for (ii = 0; ii < _PB_NI; ii += TI) {
            for (jj = 0; jj < _PB_NJ; jj += TJ) {
                int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
                int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
                for (i = ii; i < i_end; i++) {
                    for (j = jj; j < j_end; j++) {
                        ... 
                    }
                }
            }
        }

   We can also consider collapsing the ii and jj loops? But note: the jj loop is inside ii. We can collapse two levels? 

        #pragma omp parallel for collapse(2) private(i, j, k, i_end, j_end)

   But note: the ii and jj loops are not perfectly nested because we have the private variables i_end and j_end that depend on ii and jj. We can restructure to have the ii and jj loops without the inner loops? Actually, we can do:

        #pragma omp parallel for collapse(2) private(i, j, k) 
        for (ii = 0; ii < _PB_NI; ii += TI) {
            for (jj = 0; jj < _PB_NJ; jj += TJ) {
                int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
                int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
                for (i = ii; i < i_end; i++) {
                    for (j = jj; j < j_end; j++) {
                        ... 
                    }
                }
            }
        }

   But note: the collapse(2) will collapse the two outer loops (ii and jj). This is safe because the inner loops (i and j) are independent for each (ii,jj). However, the number of iterations of the ii loop is _PB_NI / TI (approximately) and similarly for jj. The total number of iterations is (NI/TI) * (NJ/TJ). This might be too large if TI and TJ are small? But we are using fixed tile sizes, so it should be acceptable.

   Alternatively, we can parallelize only the ii loop and then the jj loop is sequential. This might be better if the number of tiles in ii is large enough to keep all threads busy.

   We'll choose to collapse the two outer loops to have more parallelism.

   Also, we can add a schedule. We can use static scheduling because the work per tile is uniform.

        schedule(static)

   We can also consider SIMD for the innermost loop (k). We can add:

        #pragma omp simd reduction(+:temp)

   But note: the inner loop is a reduction. We can do:

        DATA_TYPE temp = beta * C[i][j];
        #pragma omp simd reduction(+:temp)
        for (k = 0; k < _PB_NK; k++) {
            temp += alpha * A[i][k] * B[k][j];
        }

   This tells the compiler to vectorize the k loop and to use a reduction for temp.

   However, we must be cautious: the reduction requires that the variable 'temp' is private to the SIMD lane. The compiler will handle that.

   Also, we can align the arrays? But we don't have control over the allocation. We can add:

        #pragma omp simd aligned(C, A, B: 64) ...

   But the problem does not specify the alignment. We can skip the aligned clause.

   We can also use 'safelen' to specify the safe vector length? But we don't know the architecture. We can leave it to the compiler.

   So the inner loop becomes:

        DATA_TYPE temp = beta * C[i][j];
        #pragma omp simd reduction(+:temp)
        for (k = 0; k < _PB_NK; k++) {
            temp += alpha * A[i][k] * B[k][j];
        }
        C[i][j] = temp;

   We must declare 'temp' inside the j loop so that it is private to each (i,j).

   Also, we must declare the tile sizes. We'll define:

        #define TI 64
        #define TJ 64

   But note: the problem says "choose legal tile sizes". We assume 64 is a reasonable tile size for L1 cache.

   However, we must be cautious: the tile sizes must divide the problem sizes? We are handling the boundaries with i_end and j_end.

   We'll write the code accordingly.

   But note: the original code uses _PB_NI, _PB_NJ, _PB_NK. We must use these.

   Also, we must preserve the original variable names and types.

   We'll write the optimized code within the scop markers.

   Let's code accordingly.

   However, note: the original code uses:
        for (i = 0; i < _PB_NI; i++) {
            for (j = 0; j < _PB_NJ; j++)
                C[i][j] *= beta;
            for (k = 0; k < _PB_NK; k++) {
                for (j = 0; j < _PB_NJ; j++)
                    C[i][j] += alpha * A[i][k] * B[k][j];
            }
        }

   We are replacing it with the tiled and fused version.

   We must be cautious: the original code does the scaling for all j for a fixed i, then the k loop. Our version does the scaling and then the entire k loop for each (i,j). This is equivalent.

   We'll now write the code.

   Note: We must not change the function signature or the code outside the scop markers.

   We'll define the tile sizes as macros at the beginning of the scop region? But the problem says to preserve the marker lines. We can define them inside the scop region? But the problem says to preserve the markers. We can define them after the '#pragma scop' and before the code.

   Alternatively, we can hardcode the numbers. We'll use:

        const int TI = 64;
        const int TJ = 64;

   But note: we are inside a function. We can declare them as constants.

   However, the problem does not allow changing the code outside the scop markers. We are only allowed to change the code between the markers.

   We'll write:

#pragma scop
  int i, j, k, ii, jj;
  const int TI = 64;
  const int TJ = 64;

  #pragma omp parallel for collapse(2) private(i, j, k, ii, jj) schedule(static)
  for (ii = 0; ii < _PB_NI; ii += TI) {
    for (jj = 0; jj < _PB_NJ; jj += TJ) {
      int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
      int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
      for (i = ii; i < i_end; i++) {
        for (j = jj; j < j_end; j++) {
          DATA_TYPE temp = beta * C[i][j];
          #pragma omp simd reduction(+:temp)
          for (k = 0; k < _PB_NK; k++) {
            temp += alpha * A[i][k] * B[k][j];
          }
          C[i][j] = temp;
        }
      }
    }
  }
#pragma endscop

   But note: the variables ii and jj are declared in the outer scope. We are using them in the parallel region. We must make them private? Actually, the OpenMP directive has 'private(i, j, k, ii, jj)'. But note: ii and jj are loop variables of the outer loops. They are private by default in a parallel for? Actually, in a 'parallel for', the loop variable of the for-loop is private. But we have two nested loops (ii and jj) and we are collapsing them. The OpenMP standard says that the loop variables of the collapsed loops are private. So we don't need to list them? But we are also using ii and jj inside the inner loops. We must make them private? Actually, the OpenMP specification says that the loop variables of the associated loops are private. So we are safe.

   However, we are also using i_end and j_end. We must make them private? Yes, because they are declared inside the parallel region and are local to each thread. We have declared them inside the parallel region, so each thread will have its own. But note: we are using 'private' only for i, j, k, ii, jj. We did not list i_end and j_end. We must add them? Actually, they are declared inside the parallel region, so they are private by default? No, in C, variables declared inside a block are local to that block. Since the parallel region creates multiple threads, each thread will execute the entire block. So each thread will have its own i_end and j_end. So we don't need to list them.

   But note: the OpenMP directive is on the outer loop. The variables i_end and j_end are declared inside the loop body. They are local to each iteration of the outer loop. So they are private.

   We'll leave as is.

   However, we must check: the variables i, j, k are used in the inner loops. We have listed them as private. That is correct.

   Also, note: the variables alpha, beta, C, A, B are shared. We don't need to specify because by default they are shared. But to be safe, we can add 'default(none)' and then list all shared variables? The problem says we can use 'default(none)'. We'll do:

        #pragma omp parallel for collapse(2) private(i, j, k) default(none) shared(_PB_NI, _PB_NJ, _PB_NK, TI, TJ, C, A, B, alpha, beta) schedule(static)

   But note: TI and TJ are constants declared inside the scop. They are shared. Also, the problem sizes _PB_NI, _PB_NJ, _PB_NK are shared.

   We'll add 'default(none)' and list all shared variables.

   However, note: the variables ii and jj are loop variables and are private by the collapse(2). We don't need to list them as private? Actually, the OpenMP specification says that the loop variables of the collapsed loops are private. So we don't need to list them. But we are using them in the inner block. We must ensure they are private. The collapse(2) makes them private.

   We'll write:

        #pragma omp parallel for collapse(2) default(none) \
            shared(_PB_NI, _PB_NJ, _PB_NK, TI, TJ, C, A, B, alpha, beta) \
            private(i, j, k, i_end, j_end) \
            schedule(static)

   But note: we are using i_end and j_end. We must declare them as private? Yes, because they are used in the inner loops and are local to the thread. We have declared them inside the loop, but we must list them as private? Actually, the OpenMP specification says that variables declared inside the loop body are private. So we don't need to list them. But to be safe, we list them.

   Alternatively, we can avoid listing them by not declaring them until inside the loop? We are already doing that.

   We'll list only the variables that are declared outside the parallel region and are used inside: i, j, k. And we'll list the shared variables.

   But note: the variables ii and jj are loop variables and are private by the collapse. We don't list them.

   We'll write:

        #pragma omp parallel for collapse(2) default(none) \
            shared(_PB_NI, _PB_NJ, _PB_NK, TI, TJ, C, A, B, alpha, beta) \
            private(i, j, k) \
            schedule(static)

   And then inside the loop, we declare i_end and j_end. They are private by the fact that they are declared inside the loop.

   This should be safe.

   Let's write the final code.

   Note: We must not change the original variable declarations. The original code declares i, j, k at the beginning of the function. We are reusing them. We are also declaring ii, jj, i_end, j_end inside the scop. This is acceptable.

   We'll now output the entire code with the optimized scop region.

   IMPORTANT: We must preserve the rest of the code.

   We'll replace the scop region with the above.

   Let's write the code.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gemm.h"


/* Array initialization. */
static
void init_array(int ni, int nj, int nk,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      C[i][j] = (DATA_TYPE) ((i*j+1) % ni) / ni;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) (i*(j+1) % nk) / nk;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+2) % nj) / nj;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nj,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++) {
	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
    }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemm(int ni, int nj, int nk,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		 DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		 DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j, k;

//BLAS PARAMS
//TRANSA = 'N'
//TRANSB = 'N'
// => Form C := alpha*A*B + beta*C,
//A is NIxNK
//B is NKxNJ
//C is NIxNJ
#pragma scop
  const int TI = 64;
  const int TJ = 64;
  #pragma omp parallel for collapse(2) default(none) \
      shared(_PB_NI, _PB_NJ, _PB_NK, TI, TJ, C, A, B, alpha, beta) \
      private(i, j, k) \
      schedule(static)
  for (int ii = 0; ii < _PB_NI; ii += TI) {
    for (int jj = 0; jj < _PB_NJ; jj += TJ) {
      int i_end = (ii + TI) < _PB_NI ? (ii + TI) : _PB_NI;
      int j_end = (jj + TJ) < _PB_NJ ? (jj + TJ) : _PB_NJ;
      for (i = ii; i < i_end; i++) {
        for (j = jj; j < j_end; j++) {
          DATA_TYPE temp = beta * C[i][j];
          #pragma omp simd reduction(+:temp)
          for (k = 0; k < _PB_NK; k++) {
            temp += alpha * A[i][k] * B[k][j];
          }
          C[i][j] = temp;
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

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);

  /* Initialize array(s). */
  init_array (ni, nj, nk, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemm (ni, nj, nk,
	       alpha, beta,
	       POLYBENCH_ARRAY(C),
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nj,  POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}