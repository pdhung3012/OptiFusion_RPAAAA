We are given the covariance.c code and are to optimize the region between #pragma scop and #pragma endscop.
 The code has three main parts:
  1. Compute mean for each column j: mean[j] = (sum_{i} data[i][j]) / float_n
  2. Subtract the mean from each element: data[i][j] -= mean[j]
  3. Compute the covariance matrix: 
        cov[i][j] = (sum_{k} (data[k][i] * data[k][j])) / (float_n - 1.0)
     and then set cov[j][i] = cov[i][j] (since covariance matrix is symmetric)

 We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

 Step 1: Loop Transformations
   - Tiling: We note that the third part (covariance computation) is a symmetric matrix update and involves a reduction over k. 
     This is similar to a matrix multiplication (data^T * data) scaled by 1/(n-1). Tiling can be applied to the i,j loops to improve cache locality.
     However, note that the first two parts are relatively simple and might not require tiling.

   - Fusion: The first two parts are adjacent and have compatible iteration spaces? 
        Part1: j from 0 to M-1, then i from 0 to N-1 -> mean[j] = sum(data[i][j]) / n
        Part2: i from 0 to N-1, j from 0 to M-1 -> data[i][j] -= mean[j]
     We cannot fuse these because the first part must complete for all j before the second part can start (since the second part uses mean[j] for every j).

   - Skewing: There are no carried dependencies that require skewing in the first two parts. The third part has a reduction over k and then a symmetric assignment. 
     The reduction is independent for each (i,j) pair, so no skewing is needed.

   Therefore, we focus on tiling the third part (the covariance computation) because it is the most compute-intensive and has a matrix-multiplication-like pattern.

   Proposed tiling for the third part:
     We tile the i and j loops (which are the outer loops for the covariance computation) with a tile size of, say, 64 (a common cache-friendly block size). 
     We note that the inner loop is over k (from 0 to N-1). We can also tile the k loop? But note that the inner loop is a reduction and we are accumulating into cov[i][j]. 
     However, the reduction is over k and we are not tiling the reduction dimension in the same way as in matrix multiplication. Instead, we are tiling the two outer dimensions.

   How to tile:
     We break the i and j loops into blocks of size TI and TJ. We then iterate over the blocks. For each block, we load a tile of data (for the current block of i and j) and compute the covariance for that tile.

   But note: the inner loop is over k, and we are accessing data[k][i] and data[k][j]. We want to reuse the data for a block of i and j for multiple k? Actually, we are iterating k from 0 to N-1. 
   We can also tile the k loop? However, the reduction is over the entire k, so we cannot break the reduction. Instead, we can tile the i and j and then for each tile, we do the entire k loop.

   Alternatively, we can do:
        for (ii = 0; ii < M; ii += TI)
          for (jj = ii; jj < M; jj += TJ)   // note: we are only doing the lower triangle including the diagonal
            for (i = ii; i < min(ii+TI, M); i++)
              for (j = jj; j < min(jj+TJ, M); j++)
                if (j >= i) {   // because we are only computing the lower triangle? Actually, the original code does j from i to M-1.
                  ... the inner k loop ...

   However, note the original code: 
        for (i=0; i<M; i++)
          for (j=i; j<M; j++)

   We are only computing the lower triangle (including the diagonal). Then we set the symmetric element.

   We can tile the i and j loops. We choose a tile size that fits in L1 cache. Since we are accessing two columns of data (for a fixed i and j) over k, and we are accumulating a tile of cov, we want the tile of cov to fit in cache and also the two columns of data (which are of length N) to be in cache? Actually, we are iterating k, so we are streaming through the columns. 

   Alternatively, we can interchange the loops to have k in the outer and then tile i and j? But note: the reduction is over k, so we cannot change the order arbitrarily.

   We decide to tile the i and j loops. We choose a tile size of 64 for both i and j (if M is large enough). We must be cautious because the inner loop is over k and we are doing a reduction. We can do:

        #pragma scop
        ... first two parts unchanged ...

        // Tiled version of the third part
        const int TI = 64;
        const int TJ = 64;
        for (int ii = 0; ii < _PB_M; ii += TI) {
          int i_end = min(ii+TI, _PB_M);
          for (int jj = ii; jj < _PB_M; jj += TJ) {   // note: we start jj from ii because of the lower triangle
            int j_end = min(jj+TJ, _PB_M);
            for (int i = ii; i < i_end; i++) {
              for (int j = jj; j < j_end; j++) {
                if (j < i) continue;   // because we are only computing j>=i? Actually, we are iterating j from jj to j_end, but we must ensure j>=i? 
                // But note: our jj starts at ii and we are iterating i from ii to i_end and j from jj to j_end. 
                // However, if jj is less than ii, then we might have j < i? Actually, we set jj starting at ii, so j>=ii and i>=ii -> j>=i? 
                // But wait: we are iterating i from ii to i_end and j from jj to j_end. Since jj>=ii, then j>=jj>=ii<=i -> so j>=i? 
                // Actually, no: because i might be greater than j? For example, if ii=0, jj=0, then i and j both start at 0. 
                // But we are only computing the lower triangle? The original code does j from i to M-1. 
                // We are iterating i from ii to i_end and j from jj to j_end. We must ensure that we only compute when j>=i? 
                // However, note that we are starting j at jj which is at least ii, and i is at least ii. But if j < i, then we skip? 
                // Actually, we are iterating j from jj to j_end, and i from ii to i_end. We want to compute only when j>=i? 
                // But note: the original code does for i, then j from i to M-1. So we are only computing the lower triangle. 
                // In our tiled version, we are iterating over a block of i and j. We must compute only the part of the block that is in the lower triangle? 
                // We can do: if (j < i) then skip? But note: we are starting j at jj which is at least ii, and i is at least ii. 
                // However, if we are in a block that is above the diagonal? Actually, we are starting jj at ii, so the block is from (ii, ii) to (i_end, j_end). 
                // This block is entirely in the lower triangle? Not exactly: because we are iterating i from ii to i_end and j from jj to j_end, and we require j>=i? 
                // Actually, we are only computing the lower triangle, so we should only compute when j>=i. 
                // But note: we are iterating j from jj to j_end, and i from ii to i_end. We must skip if j < i? 
                // However, we are starting j at jj which is at least ii, and i is at least ii. But if j < i, then we skip? 
                // Example: ii=0, jj=0, then we have i=0, j=0 -> j>=i -> compute. Then i=0, j=1 -> j>=i -> compute. Then i=1, j=0 -> j<i -> skip? 
                // But wait, we are iterating j from jj (which is 0) to j_end, and i from ii (0) to i_end. So we are iterating over the entire block? 
                // We must only compute the part of the block that is in the lower triangle? 
                // We can do: if (j < i) continue; 
                // But note: the original code does j from i to M-1, so we are only computing the lower triangle. 
                // Therefore, we must skip when j < i.

                // However, we are starting j at jj which is at least ii, and i is at least ii. But if we are in a block that is above the diagonal? 
                // Actually, we are iterating i from ii to i_end and j from jj to j_end. We are only interested in j>=i. 
                // So we do: if (j < i) continue;

                // But note: we are also going to set cov[j][i] = cov[i][j] at the end. However, we are only computing the lower triangle? 
                // Actually, the original code does:
                //   for (i=0; i<M; i++)
                //     for (j=i; j<M; j++)
                //        ... compute cov[i][j] ...
                //        cov[j][i] = cov[i][j];
                // So we are computing the lower triangle and then setting the symmetric element.

                // Therefore, we must only compute when j>=i. So we skip if j < i.

                // However, we are iterating j from jj to j_end and i from ii to i_end. We can avoid the condition by adjusting the loop bounds? 
                // We can set the inner j loop to start at max(jj, i) and then we don't need the condition? 
                // But note: we are tiling, so we want to keep the inner loops simple. We'll use the condition.

                // Alternatively, we can restructure the j loop to start at i? But then we break the tiling? 
                // We decide to keep the condition.

                // But note: we are also going to set the symmetric element. We must be cautious: we are only computing the lower triangle? 
                // Actually, we are computing cov[i][j] for j>=i, and then we set cov[j][i] = cov[i][j]. 
                // So we are computing the entire matrix? But we are only doing the lower triangle and then setting the upper by symmetry.

                // Therefore, we do:
                if (j < i) continue;

                cov[i][j] = SCALAR_VAL(0.0);
                for (k = 0; k < _PB_N; k++) {
                  cov[i][j] += data[k][i] * data[k][j];
                }
                cov[i][j] /= (float_n - SCALAR_VAL(1.0));
                if (i != j) {
                  cov[j][i] = cov[i][j];
                }
              }
            }
          }
        }

   However, note: the above tiling does not change the order of the k loop. We are still doing the entire k loop for each (i,j). 
   This tiling improves locality for the cov array: we are writing a tile of cov and then we set the symmetric part. 
   But note: the data array is accessed by columns: data[k][i] and data[k][j]. We are streaming through k, so we are accessing two columns (i and j) for the entire k. 
   The tiling of i and j does not help with the data array because we are still streaming through k. 

   Alternatively, we can interchange the loops to have k in the outer and then tile i and j? But note: the reduction is over k, so we cannot break the reduction. 
   We can do:

        for (k = 0; k < _PB_N; k++) {
          for (i = 0; i < _PB_M; i++) {
            for (j = i; j < _PB_M; j++) {
              cov[i][j] += data[k][i] * data[k][j];
            }
          }
        }

   But then we are not doing the division and the symmetric assignment until the end? And we are accumulating. This is a reduction. 
   We can tile the i and j loops inside the k loop? 

   However, the original code does the entire reduction for one (i,j) and then moves to the next. We are changing the order. 
   This is a reduction, so it is safe to change the order? Yes, because addition is associative. 

   But note: floating point addition is not associative, but the problem does not specify strict floating point semantics? 
   We are allowed to reassociate? The problem says: "Do not change observable behavior." and "preserve original semantics". 
   We must be cautious. However, the original code does not specify the order of summation? It is a simple loop. 
   We are allowed to change the order? The problem says: "Do not change observable behavior." and we are using the same operations. 
   But floating point reassociation can change the result. 

   Therefore, we avoid changing the order of the reduction. We stick to the original order: for each (i,j), we do the entire k loop.

   So we keep the tiling of i and j as above.

   But note: the tiling of i and j might not be beneficial for the data array because we are still streaming through k. 
   However, it does help with the cov array: we are writing a tile of cov and then we set the symmetric part. 
   Also, we are reusing the same tile of cov for multiple k? No, because we are accumulating over k. 

   Actually, the tiling of i and j does not help with the data array because we are accessing two entire columns (i and j) for each (i,j) in the tile. 
   And we are doing that for every k. 

   We might consider tiling the k loop? But note: the reduction is over the entire k, so we cannot break the reduction. 

   Alternatively, we can use a different tiling strategy: tile the k loop and then for each tile of k, we update the entire cov matrix? 
   But that would require storing partial sums for the entire cov matrix, which is O(M^2) and might be too big. 

   Given the constraints, we decide to tile the i and j loops as above. We choose a tile size of 64. 

   However, note: the first two parts are not tiled. They are relatively simple and we leave them as is.

 Step 2: OpenMP directives
   We now add OpenMP directives to the three parts.

   Part1: 
        for (j = 0; j < _PB_M; j++) {
          mean[j] = SCALAR_VAL(0.0);
          for (i = 0; i < _PB_N; i++)
            mean[j] += data[i][j];
          mean[j] /= float_n;
        }

     We can parallelize the outer j loop? 
        Each j is independent -> yes.
        We can do: #pragma omp parallel for

     But note: the inner loop is a reduction? Actually, we are accumulating into mean[j]. 
        We can use a reduction clause? But note: the reduction is over i for each j. 
        We can do: 
            #pragma omp parallel for private(i) reduction(+:mean) 
        But wait: mean is an array. We cannot use reduction on an array? 
        Alternatively, we can make mean private and then do a critical section? That would be inefficient.

     Instead, we note that each j iteration is independent and writes to a separate element of mean. 
        So we can do:
            #pragma omp parallel for
            for (j=0; j<_PB_M; j++) {
                ... 
            }

        But note: the inner loop is a reduction for a single element. We don't need a reduction clause because each thread writes to its own j.

     However, we must declare i as private? Yes, because it is a loop index.

     We can do:
        #pragma omp parallel for private(i)
        for (j=0; j<_PB_M; j++) {
          mean[j] = SCALAR_VAL(0.0);
          for (i=0; i<_PB_N; i++) {
            mean[j] += data[i][j];
          }
          mean[j] /= float_n;
        }

   Part2:
        for (i=0; i<_PB_N; i++)
          for (j=0; j<_PB_M; j++)
            data[i][j] -= mean[j];

     We can parallelize the outer i loop? 
        Each i is independent -> yes.
        We can do: #pragma omp parallel for

        But note: the inner loop is over j, and we are reading mean[j] (which is shared) and writing data[i][j] (which is independent per i).
        We can do:
            #pragma omp parallel for private(j)
            for (i=0; i<_PB_N; i++) {
              for (j=0; j<_PB_M; j++) {
                data[i][j] -= mean[j];
              }
            }

   Part3: the tiled version we created above.

        We can parallelize the outer tile loops? 
          We have two outer loops: ii and jj. 
          We note that the tiles are independent? 
          But note: we are writing to cov[i][j] and then setting cov[j][i] = cov[i][j]. 
          However, we are only writing to the lower triangle and then setting the symmetric element. 
          But note: when we set cov[j][i] = cov[i][j], we are writing to a location that might be in a different tile? 
          However, we are iterating over tiles in a way that the lower triangle is computed. 
          And we are setting the symmetric element only when i != j. 
          But note: if we are in a tile (ii, jj) and we compute cov[i][j] for i in [ii, ii+TI) and j in [jj, jj+TJ) and j>=i, then we set cov[j][i] = cov[i][j]. 
          The location (j,i) might be in a different tile? 
          Actually, we are iterating over the entire lower triangle. The symmetric element (j,i) is in the upper triangle and we are not computing it in the lower triangle. 
          But note: we are setting it from the computed value. 
          However, we are not computing the upper triangle. 

          The problem: if two tiles are being computed in parallel and they both try to write to the same symmetric element? 
          For example, tile (ii, jj) might compute cov[i][j] and then set cov[j][i]. Meanwhile, tile (jj, ii) might compute cov[j][i] and then set cov[i][j]? 
          But note: we are only iterating jj from ii to M-1. So we are only computing the lower triangle. The tile (jj, ii) would not be computed because jj>=ii and ii is the starting row? 
          Actually, we are iterating ii from 0 to M-1 and for each ii, we iterate jj from ii to M-1. So we are only covering the lower triangle. 
          Therefore, the symmetric element (j,i) is in the upper triangle and we are not computing it in the lower triangle. 
          But note: we are setting it from the computed value. 

          However, the tile that computes (i,j) (with j>=i) will set cov[j][i] = cov[i][j]. 
          And the tile that would compute (j,i) (if we were computing the entire matrix) is not computed because we skip when j<i? 
          Actually, we are not computing the entire matrix. We are only computing the lower triangle. 

          Therefore, the write to cov[j][i] is safe because no other tile will write to that location? 
          But note: the tile that computes (i,j) writes to cov[i][j] and then to cov[j][i]. 
          And the tile that computes (j,i) would be skipped because j>i? So no other tile writes to cov[j][i]. 

          Therefore, the tiles are independent? 

          However, note: the tile (ii, jj) and the tile (ii, jj') (with jj' != jj) might write to the same row? 
          But they are writing to different j in the same row? 
          Actually, they are writing to different (i,j) pairs. 

          So we can parallelize the outer tile loops? 

          We can do:
            #pragma omp parallel for private(jj, j_end, i, j, k) 
            for (ii=0; ii<_PB_M; ii+=TI) {
              ...
            }

          But note: the inner loops are over jj, then i, then j, then k. 
          We must be cautious: the inner loops are nested and we are using min functions. 

          Alternatively, we can collapse the two outer tile loops? 
            #pragma omp parallel for collapse(2) 
          But note: the inner loop for jj starts at ii, so the iteration space is triangular. 
          We cannot collapse because the inner loop bound depends on the outer loop index. 

          Therefore, we parallelize only the outer ii loop? 

          We do:
            #pragma omp parallel for private(jj, j_end, i, j, k) 
            for (ii=0; ii<_PB_M; ii+=TI) {
              int i_end = min(ii+TI, _PB_M);
              for (jj=ii; jj<_PB_M; jj+=TJ) {
                ...
              }
            }

          But note: the inner jj loop is not parallelized. We are only parallelizing the ii loop. 
          This might not be enough parallelism if TI is large and the number of tiles in ii is small. 

          Alternatively, we can parallelize the jj loop? But note: the jj loop is inside the ii loop and the jj loop bound depends on ii. 
          We can use the OpenMP 'for' directive on the jj loop? But then we would have nested parallelism? 
          We avoid nested parallelism because it can lead to too many threads. 

          We decide to parallelize the outer ii loop. 

          However, we note that the work per ii tile might be uneven because the jj loop goes from ii to M-1. 
          The first tiles (small ii) have more jj tiles. 

          We can use dynamic scheduling? But static is usually better for regular work. 

          We'll use static scheduling. 

   We also note that the inner k loop is a reduction? We are accumulating cov[i][j] over k. 
        We are not using a reduction clause because we are not parallelizing the k loop. 

   We must declare the inner loop variables as private? 
        We have: jj, j_end, i, j, k. 
        We are using them in the inner loops. We must make them private to the thread. 

   We can do:
        #pragma omp parallel for private(jj, j_end, i, j, k) 
        for (ii=0; ii<_PB_M; ii+=TI) {
          int i_end = min(ii+TI, _PB_M);
          for (jj=ii; jj<_PB_M; jj+=TJ) {
            int j_end = min(jj+TJ, _PB_M);
            for (i=ii; i<i_end; i++) {
              for (j=jj; j<j_end; j++) {
                if (j < i) continue;
                ... 
              }
            }
          }
        }

   But note: the variables jj, j_end, i, j, k are declared in the outer scope? We are using them in the inner loops. 
        We must declare them as private in the OpenMP directive. 

   Alternatively, we can declare them inside the parallel region? 
        We can do:
            #pragma omp parallel
            {
              for (ii=0; ii<_PB_M; ii+=TI) {
                int i_end = min(ii+TI, _PB_M);
                for (jj=ii; jj<_PB_M; jj+=TJ) {
                  ... 
                }
              }
            }

        But then we are not worksharing. 

   We want worksharing on the ii loop. 

   We do:
        #pragma omp parallel for private(jj, j_end, i, j, k) 
        for (ii=0; ii<_PB_M; ii+=TI) {
          int i_end = min(ii+TI, _PB_M);
          for (jj=ii; jj<_PB_M; jj+=TJ) {
            int j_end = min(jj+TJ, _PB_M);
            for (i=ii; i<i_end; i++) {
              for (j=jj; j<j_end; j++) {
                if (j < i) continue;
                cov[i][j] = SCALAR_VAL(0.0);
                for (k=0; k<_PB_N; k++) {
                  cov[i][j] += data[k][i] * data[k][j];
                }
                cov[i][j] /= (float_n - SCALAR_VAL(1.0));
                if (i != j) {
                  cov[j][i] = cov[i][j];
                }
              }
            }
          }
        }

   We must be cautious: the variables jj, j_end, i, j, k are used in the inner loops and must be private to each thread. 
        We list them in the private clause. 

   However, note: the variable k is declared in the inner loop? Actually, we are declaring it in the for loop? 
        In C99, we can declare it in the for loop. But we are using it in the inner loop. 
        We must make sure that each thread has its own k. 

   We are using private(k) so that each thread has a private copy. 

   But note: the inner loop for k is inside the j loop, which is inside the i loop, which is inside the jj loop, which is inside the ii loop. 
        And we are parallelizing the ii loop. So each thread will have its own k. 

   We also note that the division and the symmetric assignment are done after the k loop. 

   We must be cautious: the symmetric assignment cov[j][i] = cov[i][j] might be writing to a location that is in a different tile? 
        But as argued, it is safe because no other tile will write to that location. 

   We also note that the first two parts are parallelized independently. 

   We must add the necessary OpenMP directives to the first two parts as described.

   We also note: we are using min, so we must include <stdio.h> or <stdlib.h>? The original code includes <stdio.h> and <math.h> (which might include stdlib.h). 
        We assume min is available. 

   We define TI and TJ as constants. We choose 64. 

   We must be cautious: if M is less than 64, then the tile will be smaller. 

   We write the code accordingly.

   We also note: the original code uses _PB_M and _PB_N. We use the same.

   We now write the entire scop region.

   Important: We must preserve the #pragma scop and #pragma endscop markers.

   We also note: the original code uses SCALAR_VAL(0.0) and SCALAR_VAL(1.0). We use the same.

   We write the code accordingly.

   Let's code:

        #pragma scop
        // Part1: compute mean
        #pragma omp parallel for private(i)
        for (j = 0; j < _PB_M; j++) {
          mean[j] = SCALAR_VAL(0.0);
          for (i = 0; i < _PB_N; i++) {
            mean[j] += data[i][j];
          }
          mean[j] /= float_n;
        }

        // Part2: subtract mean
        #pragma omp parallel for private(j)
        for (i = 0; i < _PB_N; i++) {
          for (j = 0; j < _PB_M; j++) {
            data[i][j] -= mean[j];
          }
        }

        // Part3: compute covariance (tiled and parallel)
        const int TI = 64;
        const int TJ = 64;
        #pragma omp parallel for private(jj, j_end, i, j, k)
        for (int ii = 0; ii < _PB_M; ii += TI) {
          int i_end = (ii + TI > _PB_M) ? _PB_M : ii + TI;
          for (jj = ii; jj < _PB_M; jj += TJ) {
            int j_end = (jj + TJ > _PB_M) ? _PB_M : jj + TJ;
            for (i = ii; i < i_end; i++) {
              for (j = jj; j < j_end; j++) {
                if (j < i) continue;
                cov[i][j] = SCALAR_VAL(0.0);
                for (k = 0; k < _PB_N; k++) {
                  cov[i][j] += data[k][i] * data[k][j];
                }
                cov[i][j] /= (float_n - SCALAR_VAL(1.0));
                if (i != j) {
                  cov[j][i] = cov[i][j];
                }
              }
            }
          }
        }
        #pragma endscop

   We use the ternary operator for min to avoid including <algorithm> or <stdlib.h> if not already included.

   We note: the original code uses _PB_M and _PB_N, which are defined in the polybench header.

   We assume that the code is compiled with C99 or later so that we can declare loop variables in the for loop.

   We also note: the variable k is declared in the inner loop? We are declaring it in the for loop: for (k=0; ...). 
        This is allowed in C99.

   We also note: the variables ii, jj, i, j, k are declared in the scop region. We must not have any conflicts.

   We declare TI and TJ as constants inside the scop region.

   We also note: the original code uses int for i, j, k. We use the same.

   We write the code accordingly.

   However, note: the original code uses:
        for (j = 0; j < _PB_M; j++)
        for (i = 0; i < _PB_N; i++)
        for (i = 0; i < _PB_N; i++)
        for (j = 0; j < _PB_M; j++)
        for (i = 0; i < _PB_M; i++)
        for (j = i; j < _PB_M; j++)

   We are changing the order of the loops in the third part? We are tiling the i and j loops. 

   We are also changing the inner loop order? We are keeping the k loop innermost. 

   We are not changing the order of the k loop relative to the accumulation. 

   We are not changing the order of the reduction. 

   We are not changing the order of the first two parts. 

   We are not changing the order of the symmetric assignment. 

   Therefore, we believe the transformation is safe.

   We output the entire code with the scop region replaced.

   Note: We must preserve the function signature and the rest of the code.

   We only change the code between #pragma scop and #pragma endscop.

   We write the entire kernel_covariance function with the optimized scop region.

   We also note: the original code uses POLYBENCH_2D and POLYBENCH_1D macros. We do not change them.

   We write the code accordingly.

   Let's write the entire function.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "covariance.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)n;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = ((DATA_TYPE) i*j) / M;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(cov,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("cov");
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, cov[i][j]);
    }
  POLYBENCH_DUMP_END("cov");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_covariance(int m, int n,
		       DATA_TYPE float_n,
		       DATA_TYPE POLYBENCH_2D(data,N,M,n,m),
		       DATA_TYPE POLYBENCH_2D(cov,M,M,m,m),
		       DATA_TYPE POLYBENCH_1D(mean,M,m))
{
  int i, j, k;

#pragma scop
  #pragma omp parallel for private(i)
  for (j = 0; j < _PB_M; j++) {
    mean[j] = SCALAR_VAL(0.0);
    for (i = 0; i < _PB_N; i++) {
      mean[j] += data[i][j];
    }
    mean[j] /= float_n;
  }
  
  #pragma omp parallel for private(j)
  for (i = 0; i < _PB_N; i++) {
    for (j = 0; j < _PB_M; j++) {
      data[i][j] -= mean[j];
    }
  }
  
  const int TI = 64;
  const int TJ = 64;
  #pragma omp parallel for private(jj, j_end, i, j, k)
  for (int ii = 0; ii < _PB_M; ii += TI) {
    int i_end = (ii + TI > _PB_M) ? _PB_M : ii + TI;
    for (jj = ii; jj < _PB_M; jj += TJ) {
      int j_end = (jj + TJ > _PB_M) ? _PB_M : jj + TJ;
      for (i = ii; i < i_end; i++) {
        for (j = jj; j < j_end; j++) {
          if (j < i) continue;
          cov[i][j] = SCALAR_VAL(0.0);
          for (k = 0; k < _PB_N; k++) {
            cov[i][j] += data[k][i] * data[k][j];
          }
          cov[i][j] /= (float_n - SCALAR_VAL(1.0));
          if (i != j) {
            cov[j][i] = cov[i][j];
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
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE float_n;
  POLYBENCH_2D_ARRAY_DECL(data,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(cov,DATA_TYPE,M,M,m,m);
  POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);


  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_covariance (m, n, float_n,
		     POLYBENCH_ARRAY(data),
		     POLYBENCH_ARRAY(cov),
		     POLYBENCH_ARRAY(mean));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(cov)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(cov);
  POLYBENCH_FREE_ARRAY(mean);

  return 0;
}