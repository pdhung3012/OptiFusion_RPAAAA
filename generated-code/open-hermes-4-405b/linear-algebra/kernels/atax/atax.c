We are given the task to optimize the code within the '#pragma scop' and '#pragma endscop' markers in the kernel_atax function.

The original code in the scop region:

  for (i = 0; i < _PB_N; i++)
    y[i] = 0;
  for (i = 0; i < _PB_M; i++)
    {
      tmp[i] = SCALAR_VAL(0.0);
      for (j = 0; j < _PB_N; j++)
	tmp[i] = tmp[i] + A[i][j] * x[j];
      for (j = 0; j < _PB_N; j++)
	y[j] = y[j] + A[i][j] * tmp[i];
    }

We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

Step 1: Loop Transformations

1. Loop Tiling: We have nested loops. The inner loops are over j. However, note that the outer loop for the first part (initializing y) is independent. The main computation has two inner loops: one for computing tmp[i] and then one for updating y. We might consider tiling the j loops? But note that the two inner loops are adjacent and both over j. Also, the outer loop for i (from 0 to _PB_M) is the main loop.

   However, the two inner loops (the j loops) are not nested but sequential. We cannot tile across them without fusion? Also, the first loop (initializing y) is a separate loop.

   Let's consider the two inner loops in the i-loop: they are both over j and have the same bounds. We can try to fuse them? But note: the first inner loop (for tmp) uses A[i][j] and x[j] to accumulate into tmp[i]. The second inner loop uses A[i][j] and tmp[i] to update y[j]. There is a dependence: the second loop uses tmp[i] which is computed in the first loop. So we cannot fuse the two inner loops arbitrarily.

   However, we can try to interchange the order? Actually, the two inner loops are independent in terms of j? But note: the first inner loop is a reduction on tmp[i] and the second is a reduction on y[j]. We cannot change the order of the two inner loops because the second loop uses the result of the first (tmp[i] is computed in the first and used in the second). So we cannot fuse them without changing the order? Actually, we can fuse them by doing both operations in the same j loop? But note: the first loop must complete for a fixed i before the second loop starts because the second loop uses the entire tmp[i] (which is the result of the first loop). Therefore, we cannot fuse the two inner loops.

   Alternatively, we can consider tiling the outer loop? The outer loop is over i (from 0 to _PB_M). We can tile the i-loop? But note that the two inner loops are over j and are independent of the outer loop except for the index i. Also, the tmp array is indexed by i and the y array by j. There is no immediate benefit in tiling the i-loop because the inner loops are over j and the data accessed in the inner loops (A[i][j] and x[j]) are contiguous in j. The tmp array is only of size M (which is the outer loop) and is private to each i. The y array is of size N and is updated in the second inner loop.

   However, note that the second inner loop updates y[j] for each j and for each i. This is a reduction over i for each j. The access pattern for y is: for each i, we update every j. This is a reduction and the order of updates does not matter. But the first inner loop (for tmp) is a dot product for each i.

   We might consider tiling the j loops? But the j loops are inner and already have good locality? The problem is that the two inner loops are separate and we are traversing the same A[i][j] twice. We can try to fuse the two inner loops? But as argued, we cannot because the second loop depends on the entire first loop.

   Alternatively, we can try to change the order of the loops? The outer loop is i, then two inner loops over j. We can try to interchange the i and j loops? But note: the first inner loop (for tmp) is a reduction over j for each i. The second inner loop is a reduction over i for each j? Actually, the second inner loop is: for a fixed i, we update y[j] for each j. So if we interchange the i and j loops, we would have:

     for j=0 to N-1
        for i=0 to M-1
          ... 

   But then the first inner loop (which is now the j loop) would be: for a fixed j, we do the first part? Actually, the first inner loop (the one that computes tmp) is:

        tmp[i] = 0;
        for j=0 to N-1
            tmp[i] += A[i][j] * x[j];

   If we interchange the i and j loops, we would have:

        for j=0 to N-1
            for i=0 to M-1
                tmp[i] += A[i][j] * x[j];

   But note: the initialization of tmp[i] to 0 must be done for each i. We cannot do that inside the j loop because then we would reset tmp[i] for every j. So we would have to move the initialization outside the j loop? Actually, we can do:

        for i=0 to M-1
            tmp[i] = 0;

        for j=0 to N-1
            for i=0 to M-1
                tmp[i] += A[i][j] * x[j];

   Then the second part (updating y) is:

        for j=0 to N-1
            for i=0 to M-1
                y[j] += A[i][j] * tmp[i];

   But note: the second part uses tmp[i] which is computed in the first part. And the first part now is a nested loop: j then i. The second part is also j then i. And the two parts are independent? Actually, the entire tmp array is computed in the first part and then used in the second part. So we can do:

        // Initialize tmp
        for (i=0; i<_PB_M; i++)
            tmp[i] = 0;

        // Compute tmp: now j outer, i inner
        for (j=0; j<_PB_N; j++)
            for (i=0; i<_PB_M; i++)
                tmp[i] += A[i][j] * x[j];

        // Update y: j outer, i inner
        for (j=0; j<_PB_N; j++)
            for (i=0; i<_PB_M; i++)
                y[j] += A[i][j] * tmp[i];

   This transformation is legal? The original order for the tmp computation was: for each i, for each j: tmp[i] += A[i][j]*x[j]. The new order: for each j, for each i: same operation. The order of additions is changed, but since addition is commutative, it is safe. Similarly, the update of y: originally for each i, for each j: y[j] += ... . Now for each j, for each i: same. Again, the order of additions is changed but safe.

   This transformation is called loop interchange. It is not one of the three (tiling, fusion, skewing) but note that the problem says "loop transformations" and we are allowed to restructure loops. Also, note that the problem says "examine tiling, fusion, skewing" but we are allowed to do other transformations if they are beneficial? However, the problem says "HIGHEST PRIORITY — LOOP TRANSFORMATIONS (examine FIRST, in this order): 1. Loop tiling / blocking, 2. Loop fusion, 3. Loop skewing". So we must try these first.

   But note: after the interchange, we have two separate loops: one for computing tmp and one for updating y. And both are nested with j outer and i inner. Now, we can try to fuse these two loops? Because they have the same outer loop (j) and the same inner loop (i). And the two inner loops (over i) are independent? Actually, the first inner loop (for tmp) updates tmp[i] and the second inner loop (for y) updates y[j] and uses tmp[i]. But note: the entire tmp array is computed in the first part and then used in the second part. We cannot fuse the two j loops because the second j loop uses the entire tmp array that is computed in the first j loop. So we cannot fuse the two outer j loops.

   However, we can try to tile the j loop? The j loop is now the outer loop for both parts. We can tile the j loop to improve locality? The data accessed: in the first part, for a fixed j, we access A[i][j] for all i and x[j] (which is the same for all i). Then in the second part, for the same j, we access A[i][j] again and tmp[i] (which is the same for all j). So if we tile the j loop, we can have a tile of j values and then for each tile, we do:

        for j in tile
            for i
                tmp[i] += A[i][j] * x[j];

        for j in tile
            for i
                y[j] += A[i][j] * tmp[i];

   But note: the second part uses tmp[i] which is updated in the first part. And the first part updates tmp[i] for the entire tile of j. Then the second part uses the updated tmp[i] for the same tile of j. This is correct? Actually, the original code computed tmp[i] for all j (the entire range) and then used it for all j. Now, we are breaking the j loop into tiles. In the first part of the tile, we update tmp[i] for a set of j's. Then in the second part, we use that partially updated tmp[i] for the same set of j's. This is not the same as the original! The original used the entire tmp[i] (which is the sum over all j) for every j in the update of y. But now, in the second part of the tile, we are using a tmp[i] that only has the contributions from the current tile of j. This is incorrect.

   Therefore, we cannot tile the j loop in this way.

   Alternatively, we can try to tile the i loop? The i loop is the inner loop. We can tile the i loop? But note: the inner loop is over i and the outer over j. The data accessed: for a fixed j, we access A[i][j] for a range of i. The array A is stored in row-major order? The problem does not specify, but typically in C, 2D arrays are row-major. So A[i][j] for fixed j and varying i is a column, which is not contiguous. Therefore, tiling the i loop might not help because the access to A is strided.

   Given the complexity and the fact that we are not seeing an obvious tiling opportunity, we might skip tiling for now.

2. Loop Fusion: We have two adjacent loops: the first initializes y to 0, and the second is the main computation. We can try to fuse the initialization of y with the main computation? But note: the main computation updates y[j] for every i. The initialization sets y[j] to 0. We can move the initialization inside the main loop? Actually, we can set y[j] to 0 at the beginning of the j loop? But the main computation has two parts: the tmp computation and the y update. And the y update is done for every i. We cannot set y[j] to 0 in the middle of the i loop because we are accumulating.

   Alternatively, we can do:

        for (j=0; j<_PB_N; j++)
            y[j] = 0;

        for (i=0; i<_PB_M; i++)
        {
            tmp[i] = 0;
            for (j=0; j<_PB_N; j++)
                tmp[i] += A[i][j] * x[j];
            for (j=0; j<_PB_N; j++)
                y[j] += A[i][j] * tmp[i];
        }

   This is the same as the original. We cannot fuse the initialization of y with the main loop because the main loop updates y for every i.

   What about fusing the two inner loops? We already discussed that we cannot because of the dependence.

   So loop fusion does not seem applicable.

3. Loop Skewing: This is typically for loops with carried dependencies to expose parallelism. In our code, the main loop (i) has two inner loops. The first inner loop (j) is a reduction and the second is also a reduction. There are no carried dependencies that prevent parallelism? Actually, the outer i loop is parallel? Let's check:

        for (i=0; i<_PB_M; i++)
        {
            tmp[i] = 0;
            for (j=0; j<_PB_N; j++)
                tmp[i] += A[i][j] * x[j];   // This is a reduction on tmp[i] -> private to i
            for (j=0; j<_PB_N; j++)
                y[j] += A[i][j] * tmp[i];   // This is a reduction on y[j] -> shared, but reduction
        }

   The outer i loop is parallel because:
        - tmp[i] is written only by iteration i -> private.
        - y[j] is updated by multiple i's, but it is a reduction (commutative and associative).

   So we don't need skewing to expose parallelism.

   Therefore, we might not need any of the three transformations? But note: the problem says to examine them in order. We did and found none applicable? Then we move to OpenMP.

Step 2: OpenMP Directives

We have two parts: the initialization of y and the main i loop.

The initialization of y:

        for (i=0; i<_PB_N; i++)
            y[i] = 0;

   This is a parallel loop? Yes, because each iteration is independent.

   The main i loop:

        for (i=0; i<_PB_M; i++)
        {
            tmp[i] = 0;
            for (j=0; j<_PB_N; j++)
                tmp[i] += A[i][j] * x[j];
            for (j=0; j<_PB_N; j++)
                y[j] += A[i][j] * tmp[i];
        }

   We can parallelize the outer i loop. We note:
        - tmp[i] is private to each i.
        - y[j] is shared and requires a reduction.

   However, note that the inner j loops are reductions. We can use a reduction clause for y? But the reduction is over the entire array y? OpenMP reduction clause for arrays is supported in OpenMP 4.5 and later? But the problem does not specify the OpenMP version. Alternatively, we can use a critical section or atomic operations? But that would be inefficient.

   Alternatively, we can use an array reduction? But the problem says we can use 'reduction(+:x)' etc. and consider FP reassociation. We can do:

        #pragma omp parallel for private(j) reduction(+: y[:_PB_N])

   But note: the reduction for an array is supported in OpenMP 4.5 and later. However, the problem does not specify the version. We must be cautious.

   Alternatively, we can restructure the code to avoid the reduction on the entire array? We can use a private array for y and then do a reduction at the end? But that would require extra memory and a final reduction step.

   However, the problem says we can use reduction for arrays? The example in the problem says 'reduction(+:x)' and x is a scalar. But the standard allows array reductions? Actually, in OpenMP 3.1 and later, array reduction is allowed for C/C++.

   But note: the original code uses y as a shared array. We can do:

        #pragma omp parallel for private(j) reduction(+: y[0:_PB_N])

   However, the syntax for array sections in reduction is available in OpenMP 4.0 and later? We must check the standard. Alternatively, we can use:

        #pragma omp parallel for private(j) reduction(+: y)

   But that would reduce the entire array? Actually, the reduction clause for an array variable in C is allowed and it will create a private copy for each thread and then combine them by element-wise addition.

   So we can do:

        #pragma omp parallel for private(j) reduction(+: y)

   But note: the array y is of size _PB_N. This is acceptable.

   However, the inner j loops: the first inner loop (for tmp) is a reduction on tmp[i] and is private to the i loop. The second inner loop updates y[j] and we are using a reduction on the entire y array.

   Also, we must consider the initialization of tmp[i] to 0. That is done inside the parallel region and is private to each thread? Actually, tmp[i] is indexed by i, and i is the loop variable of the parallel loop. So each thread will have its own set of i's. We can make tmp private? But note: tmp is an array of size M. We cannot make the entire array private because that would be too big. Instead, we note that each iteration i writes only to tmp[i]. So we can use:

        #pragma omp parallel for private(j) reduction(+: y)
        for (i=0; i<_PB_M; i++)
        {
            tmp[i] = 0;   // This is safe because each i is done by one thread.
            for (j=0; j<_PB_N; j++)
                tmp[i] += A[i][j] * x[j];
            for (j=0; j<_PB_N; j++)
                y[j] += A[i][j] * tmp[i];
        }

   But note: the reduction on y is done by element-wise addition. This is safe.

   However, the inner j loops: we can also try to parallelize the j loops? But note: the first inner loop is a reduction on tmp[i] and is not parallelizable without an inner reduction. The second inner loop is a reduction on y[j] and we are already using an outer reduction. We can try to collapse the i and j loops? But the two inner loops are not the same. We have two separate inner loops.

   Alternatively, we can try to use SIMD for the inner j loops? We can add 'simd' to the inner loops? But note: the inner loops are reductions. We can use:

        for (j=0; j<_PB_N; j++)
            tmp[i] += A[i][j] * x[j];

        can be:

        #pragma omp simd reduction(+: tmp[i])
        for (j=0; j<_PB_N; j++)
            tmp[i] += A[i][j] * x[j];

   Similarly, the second inner loop:

        #pragma omp simd
        for (j=0; j<_PB_N; j++)
            y[j] += A[i][j] * tmp[i];

   But note: the second inner loop updates y[j] which is shared. We are already in a parallel region for i, and we are using a reduction on y for the outer parallel loop. The inner SIMD loop would be done by one thread (the one that owns this i) and it is safe to use SIMD for the reduction? Actually, the reduction in the outer parallel loop will take care of combining the private copies of y. The inner SIMD loop is just a sequential reduction within the thread.

   However, the second inner loop is a reduction on y[j] and we are using SIMD. We can use a reduction in the SIMD? But note: the reduction in the outer parallel loop already creates a private copy of y for the thread. The inner SIMD loop is just updating that private copy. So we don't need an inner reduction for the SIMD? We can just use:

        #pragma omp simd
        for (j=0; j<_PB_N; j++)
            y_private[j] += A[i][j] * tmp[i];

   But wait, in the outer parallel loop, the reduction on y will create a private copy for the thread. We don't have to do anything special for the inner SIMD.

   However, the problem says we can use 'simd' and 'reduction' together? Actually, the reduction in the outer parallel loop is for the entire array y. The inner SIMD loop is just a sequential loop that updates the private copy. So we can use SIMD without an inner reduction.

   But note: the inner loop for the second part is a reduction on y[j] and we are using SIMD. We can use:

        #pragma omp simd
        for (j=0; j<_PB_N; j++)
            y[j] += A[i][j] * tmp[i];

   This is safe because the thread has its own private copy of y.

   Similarly, the first inner loop is a reduction on tmp[i] and we can use:

        #pragma omp simd reduction(+: tmp[i])
        for (j=0; j<_PB_N; j++)
            tmp[i] += A[i][j] * x[j];

   This is safe because tmp[i] is private to the thread (since i is the outer loop index and each thread has its own i).

   Now, we can also consider collapsing the outer i loop and the inner j loops? But we have two inner j loops. We cannot collapse the outer i loop with the inner j loops because the inner j loops are two separate loops. We could try to fuse the two inner j loops? But we cannot because of the dependence.

   Alternatively, we can try to collapse the outer i loop and the first inner j loop? But the first inner j loop is a reduction and the outer i loop is parallel. We can do:

        #pragma omp parallel for collapse(2) private(i,j) reduction(+: y) reduction(+: tmp[i])
        for (i=0; i<_PB_M; i++)
            for (j=0; j<_PB_N; j++)
                tmp[i] += A[i][j] * x[j];

   But note: the reduction on tmp[i] is per i. We cannot collapse because the reduction on tmp[i] is not a scalar but an array element. The reduction clause in OpenMP for an array element is not standard. We can only reduce the entire array? But we don't want to reduce the entire tmp array because each i is independent.

   Alternatively, we can remove the collapse and just use the outer parallel for and inner SIMD.

   We can also consider the scheduling? The outer i loop: we can use static scheduling because the work per iteration is the same.

   Also, we can use 'default(none)' to be safe.

   And we can use 'aligned' if we know the alignment? But the problem does not specify.

   Let's write the OpenMP version:

        // Initialize y
        #pragma omp parallel for
        for (int i=0; i<_PB_N; i++)
            y[i] = 0;

        // Main computation
        #pragma omp parallel for private(j) reduction(+: y) default(none) shared(A, x, tmp, _PB_M, _PB_N)
        for (int i=0; i<_PB_M; i++)
        {
            tmp[i] = 0;
            #pragma omp simd reduction(+: tmp[i])
            for (j=0; j<_PB_N; j++)
                tmp[i] += A[i][j] * x[j];

            #pragma omp simd
            for (j=0; j<_PB_N; j++)
                y[j] += A[i][j] * tmp[i];
        }

   But note: the variable j is declared in the outer scope? We have to declare it inside the parallel region? Actually, we can declare it in the for loop? We are using C99? The original code uses:

        int i, j;

   We can change to:

        for (int i=0; ...   and for (int j=0; ...

   But the problem says we can restructure the code. We can change the declaration of i and j to be inside the loops? That is allowed in C99.

   Alternatively, we can leave j as declared in the outer scope and use 'private(j)' in the OpenMP directive.

   We'll do:

        #pragma omp parallel for private(j) ... 
        for (i=0; i<_PB_M; i++)
        {
            ...
            for (j=0; ...)   // j is private by the clause
        }

   Also, note: the reduction on y is for the entire array. We are using 'reduction(+: y)'. This is acceptable.

   But note: the array y is shared and we are reducing it. The OpenMP runtime will create a private copy for each thread and then combine them by element-wise addition.

   We must be cautious: the array y is of size _PB_N. This is acceptable.

   Also, we can add 'schedule(static)' for the outer loop? The default is static, so we can omit.

   We can also consider using 'num_threads' and 'proc_bind'? But without specific requirements, we can leave it to the runtime.

   Now, the initialization of y: we can also use SIMD? We can do:

        #pragma omp parallel for simd
        for (int i=0; i<_PB_N; i++)
            y[i] = 0;

   This will use multiple threads and each thread will use SIMD for the stores.

   But note: the stores are to contiguous memory? Yes, because y is a 1D array. So SIMD is beneficial.

   So we can change the initialization to:

        #pragma omp parallel for simd
        for (int i=0; i<_PB_N; i++)
            y[i] = 0;

   Now, we must check the data sharing: in the main computation, we have shared(A, x, tmp, _PB_M, _PB_N). But note: tmp is written by each thread and is indexed by i. We don't need to share tmp? Actually, tmp is an output array. We are writing to tmp[i] and then reading it. But each thread writes to its own set of i's. So it is safe. We can mark it as shared? Or we can mark it as shared because the entire array is shared? But each thread only writes to its own part. So shared is correct.

   However, note: the reduction on y is done by the OpenMP runtime, so we don't have to worry about the sharing of y.

   We can also use 'aligned' for the arrays? But without knowing the alignment, we can skip.

   Finally, we must preserve the original semantics. We are using reductions and the order of floating point operations might change? The problem says "consider FP reassociation". We are allowed to change the order? The original code did the reductions in a specific order. Our transformation changes the order of additions. This is acceptable because the problem says "consider FP reassociation", meaning we are allowed to change the order.

   Therefore, we can write the optimized code.

   But note: the original code uses _PB_N and _PB_M. We must use the same.

   Let's write the code within the scop region:

        #pragma scop
        #pragma omp parallel for simd
        for (int i = 0; i < _PB_N; i++)
            y[i] = 0;

        #pragma omp parallel for private(j) reduction(+: y) default(none) shared(A, x, tmp, _PB_M, _PB_N)
        for (int i = 0; i < _PB_M; i++)
        {
            tmp[i] = SCALAR_VAL(0.0);
            #pragma omp simd reduction(+: tmp[i])
            for (j = 0; j < _PB_N; j++)
                tmp[i] = tmp[i] + A[i][j] * x[j];

            #pragma omp simd
            for (j = 0; j < _PB_N; j++)
                y[j] = y[j] + A[i][j] * tmp[i];
        }
        #pragma endscop

   We changed the inner loop variables to declare j in the outer scope? Actually, we are using the same j for both inner loops. We can declare j inside the parallel region? We can do:

        for (int i=0; i<_PB_M; i++)
        {
            tmp[i] = 0;
            for (int j=0; j<_PB_N; j++)
                ...
            for (int j=0; j<_PB_N; j++)
                ...
        }

   This is better because j is local to the block. And we don't need the 'private(j)' clause because j is declared inside the loop. So we can remove the 'private(j)'.

   Revised:

        #pragma omp parallel for reduction(+: y) default(none) shared(A, x, tmp, _PB_M, _PB_N)
        for (int i = 0; i < _PB_M; i++)
        {
            tmp[i] = SCALAR_VAL(0.0);
            #pragma omp simd reduction(+: tmp[i])
            for (int j = 0; j < _PB_N; j++)
                tmp[i] = tmp[i] + A[i][j] * x[j];

            #pragma omp simd
            for (int j = 0; j < _PB_N; j++)
                y[j] = y[j] + A[i][j] * tmp[i];
        }

   This is cleaner.

   Also, note: the initialization of tmp[i] to 0 is done by the assignment. We can remove the SCALAR_VAL(0.0) and just use 0? The original uses SCALAR_VAL(0.0) which is a macro for 0.0. We can keep it.

   We are done.

   But note: the problem says to preserve the marker lines. We must not change the '#pragma scop' and '#pragma endscop'.

   Therefore, the final code for the scop region is:

        #pragma scop
        #pragma omp parallel for simd
        for (int i = 0; i < _PB_N; i++)
            y[i] = 0;

        #pragma omp parallel for reduction(+: y) default(none) shared(A, x, tmp, _PB_M, _PB_N)
        for (int i = 0; i < _PB_M; i++)
        {
            tmp[i] = SCALAR_VAL(0.0);
            #pragma omp simd reduction(+: tmp[i])
            for (int j = 0; j < _PB_N; j++)
                tmp[i] = tmp[i] + A[i][j] * x[j];

            #pragma omp simd
            for (int j = 0; j < _PB_N; j++)
                y[j] = y[j] + A[i][j] * tmp[i];
        }
        #pragma endscop

   We removed the 'private(j)' because j is now declared inside the loop.

   We also changed the loop variables to be declared in the for loops.

   This should be safe and correct.

   Let's check the data sharing: 
        - A, x, tmp, _PB_M, _PB_N are shared.
        - y is shared and reduced.
        - i and j are private by the loop constructs.

   We use 'default(none)' to force explicit sharing.

   We use SIMD for the inner loops and for the initialization.

   We use parallel for the outer loops.

   We use reduction for y in the outer parallel loop and for tmp[i] in the inner SIMD loop.

   This meets the requirements.

   We output the entire code? The problem says: "Return ONLY the complete optimized C code." and "Preserve the '#pragma scope' and '#pragma endscop' markers exactly as supplied."

   We only change the code within the scop region. The rest of the code remains the same.

   Therefore, we output the entire code with the scop region replaced.

   Note: the original code uses:

        for (i = 0; i < _PB_N; i++)
            y[i] = 0;

        for (i = 0; i < _PB_M; i++)
        {
            tmp[i] = SCALAR_VAL(0.0);
            for (j = 0; j < _PB_N; j++)
                tmp[i] = tmp[i] + A[i][j] * x[j];
            for (j = 0; j < _PB_N; j++)
                y[j] = y[j] + A[i][j] * tmp[i];
        }

   We changed the loop variables to be declared in the for loops. This is allowed.

   We also added OpenMP directives.

   We are done.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "atax.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n))
{
  int i, j;
  DATA_TYPE fn;
  fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
      x[i] = 1 + (i / fn);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A[i][j] = (DATA_TYPE) ((i+j) % n) / (5*m);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("y");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, y[i]);
  }
  POLYBENCH_DUMP_END("y");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_atax(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n),
		 DATA_TYPE POLYBENCH_1D(tmp,M,m))
{
  int i, j;

#pragma scop
  #pragma omp parallel for simd
  for (int i = 0; i < _PB_N; i++)
    y[i] = 0;
  
  #pragma omp parallel for reduction(+: y) default(none) shared(A, x, tmp, _PB_M, _PB_N)
  for (int i = 0; i < _PB_M; i++)
    {
      tmp[i] = SCALAR_VAL(0.0);
      #pragma omp simd reduction(+: tmp[i])
      for (int j = 0; j < _PB_N; j++)
	tmp[i] = tmp[i] + A[i][j] * x[j];
      #pragma omp simd
      for (int j = 0; j < _PB_N; j++)
	y[j] = y[j] + A[i][j] * tmp[i];
    }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, M, N, m, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, M, m);

  /* Initialize array(s). */
  init_array (m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(x));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_atax (m, n,
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(x),
	       POLYBENCH_ARRAY(y),
	       POLYBENCH_ARRAY(tmp));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(tmp);

  return 0;
}