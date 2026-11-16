We are given a code with two separate nested loops. The first loop computes x1[i] += A[i][j] * y_1[j] for all i and j.
The second loop computes x2[i] += A[j][i] * y_2[j] for all i and j.

We are to optimize the region between #pragma scop and #pragma endscop.

First, we consider the highest priority loop transformations:

1. Loop tiling: We have two separate double-nested loops. We can consider tiling each of them. However, note that the two loops are independent and operate on different arrays (x1 and x2). Also, the access patterns are different: the first loop accesses A[i][j] and the second accesses A[j][i]. 

   For the first loop: 
      for i in [0, N-1], for j in [0, N-1]: x1[i] += A[i][j] * y_1[j]
   This is a matrix-vector product. The inner loop (j) has a reduction on x1[i] and accesses A by rows and y_1 by elements.

   For the second loop:
      for i in [0, N-1], for j in [0, N-1]: x2[i] += A[j][i] * y_2[j]
   This is also a matrix-vector product but the matrix A is accessed by columns (because A[j][i] is the same as A^T[i][j]).

   Tiling the first loop: We can tile the i-loop and j-loop to improve cache locality. Since the inner loop is over j and we are accessing A[i][j] (row-wise) and y_1[j] (which is contiguous), we can tile the j-loop to have better cache reuse for A and y_1. However, note that the reduction on x1[i] is over the entire j-loop. We can tile the j-loop and then accumulate the reduction in a temporary variable for the tile and then add to x1[i] at the end of the tile.

   Similarly, for the second loop: we are accessing A[j][i] (column-wise) and y_2[j] (contiguous). Tiling the j-loop would help in reusing the same column of A multiple times? Actually, in the inner loop j, we are accessing A[j][i] for a fixed i and varying j. This is a column of A. But note that the inner loop is over j and we are accessing A[j][i] which is a column. However, in C, matrices are stored row-major, so accessing a column is not contiguous. Therefore, tiling the j-loop might not be as effective for the second loop because of the non-contiguous access. Alternatively, we can interchange the loops? But note: the reduction is on x2[i] and the inner loop is j. If we interchange, we would have to change the reduction to be over i? That would break the reduction because then we would be accumulating for different i's in the inner loop. So we cannot interchange without changing the reduction pattern.

   However, note that the two loops are independent. We can also consider fusing them? But they are two separate reductions and they use different arrays. So fusion is not straightforward.

   Let's consider tiling the first loop:

   We can do:
      for (ii = 0; ii < N; ii += TI)
         for (i = ii; i < ii+TI; i++)
            for (jj = 0; jj < N; jj += TJ)
               for (j = jj; j < jj+TJ; j++)
                  x1[i] += A[i][j] * y_1[j];

   But note: the inner two loops (j and then the tile of j) are over the same j. We can also tile the j-loop and then use a temporary variable for the reduction in the tile? Actually, we can do:

      for (i = 0; i < N; i++) {
         DATA_TYPE tmp = 0.0;
         for (j = 0; j < N; j++) {
            tmp += A[i][j] * y_1[j];
         }
         x1[i] += tmp;
      }

   But that is just a rewrite and doesn't change the memory access pattern. Alternatively, we can tile the j-loop and use a temporary for the tile:

      for (i = 0; i < N; i++) {
         DATA_TYPE tmp = 0.0;
         for (jj = 0; jj < N; jj += TJ) {
            for (j = jj; j < jj+TJ && j < N; j++) {
               tmp += A[i][j] * y_1[j];
            }
         }
         x1[i] += tmp;
      }

   This doesn't change the fact that we are accessing A row-wise and y_1 contiguously. The tiling in j might help if we choose TJ to be a multiple of the cache line size? But the main benefit of tiling is for 2D locality. Here, we are only accessing one row of A and one element of y_1 per j. The row of A is contiguous, so we are already getting good spatial locality. The problem is that we are accessing the same row of A and the same y_1 array for the entire j-loop. So tiling the j-loop doesn't add much because we are already accessing the row of A in order.

   However, note that the entire row of A might not fit in cache. If we break the j-loop into tiles, then we can reuse the same cache lines for the tile? Actually, we are accessing the row of A in order, so without tiling we are already getting good spatial locality. The issue is that the row might be too long and cause cache misses. Tiling the j-loop would break the row into chunks that fit in cache? But note: we are accessing the row in order, so the hardware prefetcher might do a good job. Therefore, tiling the j-loop might not be necessary for the first loop.

   For the second loop, we are accessing A by columns. This is non-contiguous and will cause cache misses. We can tile the j-loop? But note: the inner loop is j and we are accessing A[j][i] for a fixed i. This is a column. We can tile the j-loop to break the column into tiles? But the column is non-contiguous, so we are accessing one element per cache line. Tiling the j-loop would not change the fact that we are accessing non-contiguously. Alternatively, we can interchange the loops? But then we would have:

        for (j = 0; j < N; j++)
          for (i = 0; i < N; i++)
             x2[i] += A[j][i] * y_2[j];

   Now, the inner loop is i. For a fixed j, we are accessing A[j][i] for all i. This is a row of A, which is contiguous. And we are also accessing x2[i] for all i (contiguous) and y_2[j] (which is fixed for the entire inner loop). This is much better! And note that the reduction on x2[i] is now over j, but we are accumulating for each i in the inner loop. This is safe because the reduction is associative. So we can interchange the loops for the second loop.

   Therefore, we can do:

        for (j = 0; j < N; j++)
          for (i = 0; i < N; i++)
             x2[i] += A[j][i] * y_2[j];

   This transformation is legal because the order of the reductions does not matter (addition is commutative). This improves the memory access pattern for the second loop.

   Now, after interchanging the second loop, we have:

        First loop: 
          for i in [0, N-1]: 
             for j in [0, N-1]: x1[i] += A[i][j] * y_1[j]

        Second loop (interchanged):
          for j in [0, N-1]:
             for i in [0, N-1]: x2[i] += A[j][i] * y_2[j]

   Now, both loops have an inner loop that accesses contiguous memory: the first loop accesses A[i][j] (row) and y_1[j] (contiguous) and the second loop accesses A[j][i] (which is the same as A[j][i] and in the inner loop i, so for fixed j, we are accessing a row of A) and y_2[j] (fixed for the inner loop) and x2[i] (contiguous).

   Now, we can consider tiling both loops? But note that the two loops are independent and we can also consider fusing them? However, they are two separate reductions and they use different arrays. So fusion is not straightforward.

   Alternatively, we can tile each loop independently. For the first loop, we can tile the j-loop? But as discussed, it might not be necessary. For the second loop, we can tile the i-loop? Because the inner loop is over i and we are accessing a row of A (which is contiguous) and x2[i] (contiguous). Tiling the i-loop would break the row into tiles? But the row is contiguous and we are accessing it in order. So again, the hardware prefetcher might do well. However, if the row is very long, we might get cache misses. We can tile the i-loop to break it into chunks that fit in cache? But note: we are accessing the same row of A for the entire inner loop. So without tiling, we are already getting good spatial locality. Therefore, tiling might not be necessary.

   However, we are also going to parallelize. We want to have coarse-grained parallelism. We can use tiling to create tiles that can be assigned to threads? But note: the reductions are per i (for the first loop) and per i (for the second loop). So we can parallelize the outer loop of the first loop (over i) and the outer loop of the second loop (over j) because the reductions are independent for each i (in the first loop) and for each i (in the second loop) but note: in the second loop, the reduction for a fixed i is over j, and we are doing the inner loop over i. So if we parallelize the outer j-loop, then each thread would do the entire inner i-loop for a fixed j. But then the reduction for x2[i] would be done by multiple threads (each thread for a different j) and we would need to protect the reduction. Alternatively, we can change the reduction to be private and then do a final reduction? But that would require extra storage.

   Alternatively, we can change the second loop to have the reduction on the outer loop? Actually, we can restructure the second loop to have the reduction on the outer loop? But note: the reduction is over j for each i. We can do:

        for (i = 0; i < N; i++) {
           DATA_TYPE tmp = 0.0;
           for (j = 0; j < N; j++) {
              tmp += A[j][i] * y_2[j];
           }
           x2[i] += tmp;
        }

   This is the same as the original second loop but without the interchange. But we know that the original second loop had a bad access pattern. However, if we do this, then we can parallelize the outer i-loop. And the inner j-loop would access A[j][i] (column) which is bad. So we don't want that.

   Therefore, we stick with the interchanged version for the second loop.

   Now, for parallelization:

   For the first loop: we can parallelize the outer i-loop. Each thread would do the entire inner j-loop for a set of i's. The reduction for x1[i] is per i, so we can use a reduction clause? But note: the reduction is on x1[i] and we are accumulating for each i independently. So we can do:

        #pragma omp parallel for
        for (i = 0; i < N; i++) {
           for (j = 0; j < N; j++) {
              x1[i] += A[i][j] * y_1[j];
           }
        }

   But note: the inner j-loop is a reduction on x1[i]. We can use a temporary variable to hold the reduction for the inner j-loop and then assign to x1[i] at the end. This is more efficient because it avoids false sharing and reduces the number of writes to x1[i]. So:

        #pragma omp parallel for
        for (i = 0; i < N; i++) {
           DATA_TYPE tmp = 0.0;
           for (j = 0; j < N; j++) {
              tmp += A[i][j] * y_1[j];
           }
           x1[i] += tmp;
        }

   Similarly, for the second loop (interchanged):

        #pragma omp parallel for
        for (j = 0; j < N; j++) {
           for (i = 0; i < N; i++) {
              x2[i] += A[j][i] * y_2[j];
           }
        }

   But note: in the second loop, the reduction for x2[i] is over j. However, we are parallelizing over j. That means that multiple threads (each with a different j) will update the same x2[i]. This is a race condition. We cannot do that.

   We have two options for the second loop:

   Option 1: Change the loop order so that the reduction is on the outer loop. We can do:

        for (i = 0; i < N; i++) {
           DATA_TYPE tmp = 0.0;
           for (j = 0; j < N; j++) {
              tmp += A[j][i] * y_2[j];
           }
           x2[i] += tmp;
        }

   Then we can parallelize the outer i-loop. But the inner j-loop accesses A[j][i] (column) which is bad for cache.

   Option 2: Keep the interchanged version and use a reduction per i. But we are parallelizing over j, so we cannot update x2[i] from multiple threads without synchronization. We can use an array of temporaries? We can do:

        #pragma omp parallel
        {
           DATA_TYPE *tmp = (DATA_TYPE*) calloc(N, sizeof(DATA_TYPE));
           #pragma omp for
           for (j = 0; j < N; j++) {
              for (i = 0; i < N; i++) {
                 tmp[i] += A[j][i] * y_2[j];
              }
           }
           #pragma omp critical
           {
              for (i = 0; i < N; i++) {
                 x2[i] += tmp[i];
              }
           }
           free(tmp);
        }

   But this requires extra memory and a critical section. Also, the critical section would serialize the updates.

   Option 3: Use an array of reductions? But OpenMP reduction clause for arrays is not straightforward. We can use:

        #pragma omp parallel for reduction(+:x2[:N])

   But this is only available in OpenMP 4.5 and later. And it requires that the array is contiguous and we are reducing the entire array. However, note that we are only updating x2[i] for i in [0, N-1]. This is the entire array. So we can do:

        #pragma omp parallel for reduction(+:x2[:N])
        for (j = 0; j < N; j++) {
           for (i = 0; i < N; i++) {
              x2[i] += A[j][i] * y_2[j];
           }
        }

   But note: the reduction clause for arrays in OpenMP 4.5 and later is supported. However, we must be cautious because the original code uses a macro for the array (POLYBENCH_1D). We assume that the array is contiguous and of length N.

   Alternatively, we can avoid the array reduction by restructuring the second loop to have the reduction on the outer loop? But that would require the bad memory access pattern.

   Given the trade-off, we choose to use the array reduction if we are targeting OpenMP 4.5 or later. But the problem does not specify the OpenMP version. However, the original code uses a standard that we don't know. We must be safe.

   Alternatively, we can change the second loop to have the reduction on the outer loop and then use tiling to improve the cache behavior? We can tile the j-loop? But the inner loop would be over j and we are accessing A[j][i] (column) which is bad. We can tile the i-loop? But the outer loop is i and the inner loop is j. We are accessing A[j][i] for a fixed i and varying j. This is a column. We can tile the j-loop to break the column into tiles? But the column is non-contiguous. We can do:

        for (i = 0; i < N; i++) {
           DATA_TYPE tmp = 0.0;
           for (jj = 0; jj < N; jj += TJ) {
              for (j = jj; j < jj+TJ && j < N; j++) {
                 tmp += A[j][i] * y_2[j];
              }
           }
           x2[i] += tmp;
        }

   This doesn't change the access pattern: we are still accessing the column of A non-contiguously. But we are reusing the same cache line for the tile? Actually, we are accessing one element per cache line, so tiling doesn't help much.

   Therefore, we decide to use the interchanged version and use an array reduction. We assume OpenMP 4.5 or later.

   But note: the problem says "Do not add dependencies, data races, invalid OpenMP directives, or unsupported transformations." and "If no safe optimization is available, return the original code unchanged.".

   We must check if the array reduction is supported. The problem does not specify the OpenMP version. However, the PolyBench benchmarks are often used with recent compilers. We will use the array reduction.

   Alternatively, we can avoid the array reduction by using a different approach: we can parallelize the outer i-loop of the second loop without interchanging? But then we have the bad access pattern. We can also consider using SIMD? But that is a lower priority.

   Given the above, we decide:

   - For the first loop: we restructure to use a temporary for the inner j-loop and then parallelize the outer i-loop.
   - For the second loop: we interchange the loops to get contiguous access and then use an array reduction on x2.

   Now, we write the code for the scop region:

        #pragma scop
        // First loop: x1 = x1 + A * y1
        #pragma omp parallel for
        for (i = 0; i < _PB_N; i++) {
           DATA_TYPE tmp = 0.0;
           for (j = 0; j < _PB_N; j++) {
              tmp += A[i][j] * y_1[j];
           }
           x1[i] += tmp;
        }

        // Second loop: x2 = x2 + A^T * y2
        #pragma omp parallel for reduction(+:x2[:_PB_N])
        for (j = 0; j < _PB_N; j++) {
           for (i = 0; i < _PB_N; i++) {
              x2[i] += A[j][i] * y_2[j];
           }
        }
        #pragma endscop

   But note: the array reduction in OpenMP requires that the array is contiguous and we are reducing the entire array. We are reducing x2[0:_PB_N-1]. This is the entire array. So it is safe.

   However, the problem says: "Do not add dependencies, data races, invalid OpenMP directives, or unsupported transformations." and we are using an array reduction which is OpenMP 4.5. We must be cautious. If we are not sure, we can avoid the array reduction by using a different method.

   Alternative for the second loop without array reduction:

        We can use a private array for each thread and then do a final reduction. We can do:

        #pragma omp parallel
        {
           DATA_TYPE *local_x2 = (DATA_TYPE*) calloc(_PB_N, sizeof(DATA_TYPE));
           #pragma omp for
           for (j = 0; j < _PB_N; j++) {
              for (i = 0; i < _PB_N; i++) {
                 local_x2[i] += A[j][i] * y_2[j];
              }
           }
           #pragma omp critical
           {
              for (i = 0; i < _PB_N; i++) {
                 x2[i] += local_x2[i];
              }
           }
           free(local_x2);
        }

   This avoids the array reduction but introduces extra memory and a critical section. The critical section might become a bottleneck. We can use a reduction with a lock per element? That would be too heavy.

   We decide to use the array reduction because it is more efficient and we assume OpenMP 4.5 is available.

   Now, we also consider other OpenMP clauses:

   - We can use 'schedule(static)' for the parallel for loops? The default is static, so we don't need to specify.

   - We can use SIMD? The inner loops are contiguous and we can use SIMD. We can add 'simd' to the inner loops? But note: the outer loop is parallelized and the inner loop is a reduction. We can do:

        For the first loop: 
          #pragma omp parallel for
          for (i=...) {
             tmp = 0.0;
             #pragma omp simd reduction(+:tmp)
             for (j=...) {
                 tmp += ...;
             }
             x1[i] += tmp;
          }

        Similarly, for the second loop, we can add SIMD to the inner i-loop? But note: the inner i-loop is a reduction on an array (x2) and we are using an array reduction at the outer level. We cannot use a reduction on the inner loop for the entire array. We can do:

          #pragma omp parallel for reduction(+:x2[:_PB_N])
          for (j=...) {
             #pragma omp simd
             for (i=...) {
                 x2[i] += ...;
             }
          }

        But note: the inner loop is a reduction on x2[i] for each i. We can use a SIMD reduction? However, the reduction is over j for each i, and we are in the inner loop over i. We cannot use a reduction on the inner loop because the reduction is not for a scalar but for an array. We can do:

          #pragma omp parallel for reduction(+:x2[:_PB_N])
          for (j=...) {
             #pragma omp simd
             for (i=...) {
                 x2[i] += A[j][i] * y_2[j];
             }
          }

        This is safe because the inner loop is just a sequence of independent updates to x2[i] (but note: the same x2[i] is updated for every j, but in the inner loop we are updating for a fixed j and all i). However, without a reduction on the inner loop, the updates to x2[i] are not atomic. But in the context of the outer parallel loop, we are using an array reduction that will handle the reduction across threads. The inner loop is just a serial loop for a fixed j. So we can use SIMD for the inner loop without a reduction because the inner loop is not parallel. The SIMD will vectorize the independent updates to x2[i] for different i. This is safe.

        However, note that the inner loop updates x2[i] for all i. The SIMD will do multiple i's at the same time. This is safe because the updates are independent.

   Therefore, we add SIMD to the inner loops.

   Also, we can use 'aligned' if we know the arrays are aligned? The problem does not specify. We skip.

   We can also use 'collapse'? The loops are not perfectly nested in a way that we can collapse because the outer loop is parallel and the inner loop is SIMD. We don't collapse.

   We can use 'schedule'? We leave as static.

   We can use 'num_threads'? We don't specify.

   We can use 'proc_bind'? We don't specify.

   So the final code for the scop region:

        #pragma scop
        #pragma omp parallel for
        for (i = 0; i < _PB_N; i++) {
           DATA_TYPE tmp = 0.0;
           #pragma omp simd reduction(+:tmp)
           for (j = 0; j < _PB_N; j++) {
              tmp += A[i][j] * y_1[j];
           }
           x1[i] += tmp;
        }

        #pragma omp parallel for reduction(+:x2[:_PB_N])
        for (j = 0; j < _PB_N; j++) {
           #pragma omp simd
           for (i = 0; i < _PB_N; i++) {
              x2[i] += A[j][i] * y_2[j];
           }
        }
        #pragma endscop

   But note: the array reduction for x2 might be implemented by the compiler by creating a private array for each thread and then combining. This is efficient.

   However, we must be cautious: the array reduction in OpenMP 4.5 requires that the array is contiguous and we are reducing the entire array. We are reducing x2[0:_PB_N-1] which is the entire array. So it is safe.

   Also, note that the original code uses _PB_N which is a macro for the problem size.

   We must preserve the original semantics. The transformations we did are:

     - For the first loop: we introduced a temporary variable to hold the reduction for the inner j-loop. This is safe and does not change the result.

     - For the second loop: we interchanged the loops and then used an array reduction. The interchange is safe because the reduction is associative. The array reduction is safe because it is equivalent to the original reduction.

   Therefore, we output the code with these changes.

   However, note: the problem says "Do not modify code outside the marked regions". We are only modifying the code between #pragma scop and #pragma endscop.

   We write the entire function kernel_mvt with the scop region replaced.

   But note: the original code uses variables i and j. We are reusing them. We must declare them? The original code declares them at the beginning of the function. We are using the same.

   We output the entire kernel_mvt function with the scop region replaced.

   However, the problem says: "Return ONLY the complete optimized C code." and "Preserve the '#pragma scope' and '#pragma endscop' markers exactly as supplied.".

   We note that the original code has:

        static
        void kernel_mvt(int n,
                DATA_TYPE POLYBENCH_1D(x1,N,n),
                DATA_TYPE POLYBENCH_1D(x2,N,n),
                DATA_TYPE POLYBENCH_1D(y_1,N,n),
                DATA_TYPE POLYBENCH_1D(y_2,N,n),
                DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
        {
          int i, j;

        #pragma scop
          ... original loops ...
        #pragma endscop

        }

   We replace the scop region with our optimized code.

   We also note that the original code uses _PB_N. We keep that.

   Final code for the scop region:

        #pragma scop
        #pragma omp parallel for
        for (i = 0; i < _PB_N; i++) {
           DATA_TYPE tmp = 0.0;
           #pragma omp simd reduction(+:tmp)
           for (j = 0; j < _PB_N; j++) {
              tmp += A[i][j] * y_1[j];
           }
           x1[i] += tmp;
        }

        #pragma omp parallel for reduction(+:x2[:_PB_N])
        for (j = 0; j < _PB_N; j++) {
           #pragma omp simd
           for (i = 0; i < _PB_N; i++) {
              x2[i] += A[j][i] * y_2[j];
           }
        }
        #pragma endscop

   We output the entire file? The problem says: "Return ONLY the complete optimized C code." and we are to preserve the rest of the code. But note: the problem says "Optimize only code in regions delimited by ...". So we only change the scop region.

   However, the problem says: "Return ONLY the complete optimized C code." meaning the entire file? But we are only allowed to change the scop region. So we output the entire file with the scop region replaced.

   But note: the problem says: "Do not modify code outside the marked regions.".

   Therefore, we output the entire file with the scop region in kernel_mvt replaced by the above code.

   However, the problem says: "Return ONLY the complete optimized C code." and we are to return the entire code? Yes.

   But note: the original code has a main function and other functions. We must return the entire file.

   We write the entire file with the scop region in kernel_mvt replaced.

   However, the problem says: "Preserve the '#pragma scope' and '#pragma endscop' markers exactly as supplied." so we keep the exact markers.

   We note that the original code uses the macro _PB_N. We keep that.

   We also note that the original code uses the type DATA_TYPE and the macros for array declarations. We keep them.

   We output the entire code.

   But note: the problem says "Output requirements: Return ONLY the complete optimized C code." and we are to return the entire file? Yes.

   However, the problem says: "No explanations, Markdown code fences, or surrounding text.".

   So we output the entire C code as a single code block.

   We write:
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "mvt.h"


/* Array initialization. */
static
void init_array(int n,
		DATA_TYPE POLYBENCH_1D(x1,N,n),
		DATA_TYPE POLYBENCH_1D(x2,N,n),
		DATA_TYPE POLYBENCH_1D(y_1,N,n),
		DATA_TYPE POLYBENCH_1D(y_2,N,n),
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      x1[i] = (DATA_TYPE) (i % n) / n;
      x2[i] = (DATA_TYPE) ((i + 1) % n) / n;
      y_1[i] = (DATA_TYPE) ((i + 3) % n) / n;
      y_2[i] = (DATA_TYPE) ((i + 4) % n) / n;
      for (j = 0; j < n; j++)
	A[i][j] = (DATA_TYPE) (i*j % n) / n;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(x1,N,n),
		 DATA_TYPE POLYBENCH_1D(x2,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("x1");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, x1[i]);
  }
  POLYBENCH_DUMP_END("x1");

  POLYBENCH_DUMP_BEGIN("x2");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, x2[i]);
  }
  POLYBENCH_DUMP_END("x2");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_mvt(int n,
		DATA_TYPE POLYBENCH_1D(x1,N,n),
		DATA_TYPE POLYBENCH_1D(x2,N,n),
		DATA_TYPE POLYBENCH_1D(y_1,N,n),
		DATA_TYPE POLYBENCH_1D(y_2,N,n),
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

#pragma scop
  #pragma omp parallel for
  for (i = 0; i < _PB_N; i++) {
     DATA_TYPE tmp = 0.0;
     #pragma omp simd reduction(+:tmp)
     for (j = 0; j < _PB_N; j++) {
        tmp += A[i][j] * y_1[j];
     }
     x1[i] += tmp;
  }
  
  #pragma omp parallel for reduction(+:x2[:_PB_N])
  for (j = 0; j < _PB_N; j++) {
     #pragma omp simd
     for (i = 0; i < _PB_N; i++) {
        x2[i] += A[j][i] * y_2[j];
     }
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(x1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x2, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y_1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y_2, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n,
	      POLYBENCH_ARRAY(x1),
	      POLYBENCH_ARRAY(x2),
	      POLYBENCH_ARRAY(y_1),
	      POLYBENCH_ARRAY(y_2),
	      POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_mvt (n,
	      POLYBENCH_ARRAY(x1),
	      POLYBENCH_ARRAY(x2),
	      POLYBENCH_ARRAY(y_1),
	      POLYBENCH_ARRAY(y_2),
	      POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(x1), POLYBENCH_ARRAY(x2)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(x1);
  POLYBENCH_FREE_ARRAY(x2);
  POLYBENCH_FREE_ARRAY(y_1);
  POLYBENCH_FREE_ARRAY(y_2);

  return 0;
}