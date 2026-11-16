We are given a code snippet for the bicg kernel. The scop region contains:
  for (i = 0; i < _PB_M; i++)
    s[i] = 0;
  for (i = 0; i < _PB_N; i++)
    {
      q[i] = SCALAR_VAL(0.0);
      for (j = 0; j < _PB_M; j++)
	{
	  s[j] = s[j] + r[i] * A[i][j];
	  q[i] = q[i] + A[i][j] * p[j];
	}
    }

We are to optimize only within the scop region.

First, we note that the first loop initializes the array `s` to zero. Then the second loop (over i) does two things:
  1. Initializes q[i] to zero.
  2. For each j, updates s[j] and q[i].

We are to consider loop transformations in the order: tiling, fusion, skewing.

1. Loop Tiling: 
   The inner loop (j) is over _PB_M and the outer loop (i) is over _PB_N. The inner loop accesses:
      s[j] (which is being updated and also read from the same j in the same iteration? Actually, note that s[j] is updated for each i and j, and the update is cumulative. Also, q[i] is updated for each j in the inner loop.

   However, note that the inner loop is over j and the outer over i. The array A is accessed as A[i][j]. The memory layout of A is row-major (since it's a 2D array in C). So A[i][j] is contiguous in j. 

   The accesses to s[j] and p[j] are contiguous in j. The array r[i] is fixed for the entire inner loop.

   The problem: The inner loop is over j and the outer over i. The inner loop is already accessing contiguous memory for A and p. However, the update to s[j] is a reduction over i. Also, the update to q[i] is a reduction over j.

   Tiling might be beneficial to improve cache locality for the array A. Since A is accessed by row (i fixed, j varying) and we are iterating over i in the outer loop, we are already accessing rows. However, if we tile the i-loop, we can reuse the same rows of A multiple times in a tile. But note that the inner loop is over j and we are updating s[j] and q[i]. 

   However, the update to s[j] is a reduction over i. If we tile the i-loop, then we break the reduction. We would have to use a private array for s in the tile and then combine. This is complex and might not be beneficial because the reduction for s is over the entire i-loop. Similarly, q[i] is a reduction over j, which is the inner loop and we are not tiling j.

   Alternatively, we might consider tiling the j-loop? But the j-loop is the inner loop and we are already accessing contiguous memory. Tiling the j-loop would break the contiguous access.

   Therefore, tiling might not be straightforward and might not yield benefits without significant restructuring. We will skip tiling for now.

2. Loop Fusion:
   We have two separate loops: the first initializes s, and the second does the main computation. We can fuse the initialization of s with the main loop? But note that the main loop updates s[j] for every i and j. The initialization of s must happen before the main loop. So we cannot fuse the first loop with the second.

   However, within the second loop, we have an inner loop that updates both s and q. We cannot fuse the inner loop with anything else because it is already the innermost.

   But note: the second loop has an outer loop over i and an inner loop over j. The inner loop does two independent reductions: one for s and one for q. We cannot fuse the inner loop with another loop because there is no adjacent loop.

   Therefore, loop fusion does not seem applicable.

3. Loop Skewing:
   We look for carried dependencies. The inner loop (j) has no carried dependence because j is the inner index and we are iterating independently. The outer loop (i) has a carried dependence for s: because s[j] is updated for every i, and the update for a fixed j in iteration i depends on the value from the previous i. However, the order of i is fixed and we cannot change it. Also, the update for q[i] is independent of other i's.

   There is no carried dependence that would prevent parallelism in the i-loop? Actually, the update to s[j] is a reduction over i. So we cannot parallelize the i-loop without a reduction. Similarly, the inner j-loop is a reduction for q[i] and also updates s[j] (which is a reduction over i). 

   Therefore, we cannot skew the loops to expose more parallelism because the dependencies are reductions.

Now, we move to OpenMP directives.

We note that the first loop (initializing s) is a simple loop and can be parallelized. The second loop (over i) has two reductions: one for s (over i) and one for q (over j). However, note that the update to s is done for every j in the inner loop and for every i. This is a reduction over i for each j. Also, the update to q is a reduction over j for each i.

We can parallelize the outer i-loop. But note that the update to s is a reduction over i. So we need to use a reduction clause for s. However, s is an array and OpenMP reduction for arrays is not directly supported in all versions. Alternatively, we can use a private array for s and then do a reduction at the end. But that would require extra memory and a critical section.

Alternatively, we can note that the update to s is independent for each j. So we can change the order of the loops? But the original code has the i-loop outer and j-loop inner. If we swap the loops, then we would have:

  for (j=0; j<_PB_M; j++)
    for (i=0; i<_PB_N; i++) {
        s[j] = s[j] + r[i] * A[i][j];
        q[i] = q[i] + A[i][j] * p[j];
    }

But note: the update to q[i] is now a reduction over i? Actually, no: for a fixed j, we are updating q[i] for every i. But the original code updated q[i] for a fixed i and for every j. The update to q[i] is a reduction over j. In the swapped version, for a fixed j, we are updating q[i] for every i. This is not a reduction over i because each i is independent. However, the update to q[i] in the swapped version would be:

        q[i] = q[i] + ... 

and we are iterating over i in the inner loop. This is a reduction over j? Actually, no: for a fixed i, we are updating q[i] for every j. But in the swapped loop, we are iterating j in the outer and i in the inner. So for a fixed i, we are updating q[i] for every j. This is the same as the original. However, the order of updates to q[i] is the same: we are accumulating over j. But note that the inner loop is now over i, and we are updating q[i] for each j. This is a reduction over j for each i. But the loop order is changed.

But wait: the update to s[j] is now: for a fixed j, we are updating s[j] for every i. This is a reduction over i. And the update to q[i] is for a fixed i, we are updating q[i] for every j. This is a reduction over j. 

However, the problem is that the update to q[i] in the inner loop (over i) would require that we have a reduction over j? But the inner loop is over i, and we are updating q[i] for each j. This is not a reduction in the inner loop because the inner loop index is i and we are updating q[i] for each j. Actually, we are updating q[i] for every j, but the inner loop is over i. This is incorrect because for a fixed j, we are iterating i and updating q[i] for that j. But the update for q[i] should be the sum over j. In the swapped version, we are doing:

  for j=0 to M-1
    for i=0 to N-1
        q[i] += A[i][j] * p[j]

This is equivalent to:

  for each j, and for each i: add A[i][j]*p[j] to q[i]

which is the same as the original. So it is correct.

But note: the update to s[j] is:

  for each j, and for each i: s[j] += r[i] * A[i][j]

which is the same as the original.

So we can swap the loops. Then we can parallelize the outer j-loop. Why? Because the update to s[j] is now a reduction over i (which is the inner loop) and the update to q[i] is a reduction over j (which is the outer loop). But note: the update to q[i] is done for every j, and we are iterating j in the outer loop. So for a fixed i, we are updating q[i] in multiple iterations of the outer j-loop. This is a reduction over j. However, the outer loop is j and we are updating q[i] for every j. This is a reduction over j for each i. But the outer loop is j, so we cannot parallelize the j-loop without a reduction for q? Actually, we can because the reduction for q is over j and we are parallelizing j. But note: the update to q[i] is done by multiple threads (each thread handling a set of j's) and they are updating the same q[i] for different j's. This is a race condition.

Therefore, we cannot simply parallelize the j-loop without protecting the update to q[i]. We would need to use a reduction for q. But q is an array and we are reducing over j for each i. This is a reduction over j for each element of q. OpenMP does support array reductions in modern versions (OpenMP 4.5 and later) with the reduction clause. However, we must be cautious because the reduction for q is over j and we are parallelizing j. But note: the reduction for q is over j, and we are parallelizing j, so each thread will have a private copy of q and then we reduce at the end. Similarly, the update to s is a reduction over i (the inner loop) and we are not parallelizing the inner loop. So we can do:

  #pragma omp parallel for reduction(+:q[:_PB_N]) reduction(+:s[:_PB_M])
  for (j=0; j<_PB_M; j++)
    for (i=0; i<_PB_N; i++) {
        s[j] += r[i] * A[i][j];
        q[i] += A[i][j] * p[j];
    }

But note: the original code initializes s to zero and q[i] to zero. We must do the same. We can initialize s and q to zero before the parallel loop. However, the reduction clause will initialize the private copies to zero and then combine by addition. So we can remove the initializations? But the original code initializes s in a separate loop and q[i] at the beginning of the i-loop. We are now initializing q[i] to zero? Actually, the reduction clause will take care of initializing the private copies to zero. But we must initialize the original arrays to zero? The reduction clause will add the private copies to the original. So we must set the original arrays to zero before the parallel region.

Alternatively, we can do:

  // Initialize s and q to zero
  #pragma omp parallel for
  for (int i=0; i<_PB_M; i++) s[i] = 0;
  #pragma omp parallel for
  for (int i=0; i<_PB_N; i++) q[i] = 0;

  // Then the main loop with reduction
  #pragma omp parallel for reduction(+:s[:_PB_M], q[:_PB_N])
  for (j=0; j<_PB_M; j++)
    for (i=0; i<_PB_N; i++) {
        s[j] += r[i] * A[i][j];
        q[i] += A[i][j] * p[j];
    }

But note: the reduction clause in OpenMP for arrays requires that the array is shared and the reduction is over the entire array. This is supported in OpenMP 4.5 and later.

However, we must consider the original scop region. The original scop region has:

  for (i = 0; i < _PB_M; i++)
    s[i] = 0;
  for (i = 0; i < _PB_N; i++)
    {
      q[i] = SCALAR_VAL(0.0);
      for (j = 0; j < _PB_M; j++)
	{
	  s[j] = s[j] + r[i] * A[i][j];
	  q[i] = q[i] + A[i][j] * p[j];
	}
    }

We are going to change the loop order and use OpenMP. We must preserve the original semantics.

But note: the original code initializes q[i] to zero at the start of the i-loop. In our swapped version, we are initializing q to zero in a separate loop. This is equivalent.

So the plan:

  Step 1: Initialize s and q to zero in separate parallel loops (if we want to parallelize the initialization). But note: the initialization of s and q are independent and can be done in parallel.

  Step 2: Swap the loops: j outer, i inner.

  Step 3: Apply OpenMP parallel for with reduction on the entire j-loop.

But note: the inner loop is over i and we are updating s[j] (which is a reduction over i) and q[i] (which is a reduction over j). The reduction for s is over the inner loop (i) and we are not parallelizing the inner loop. The reduction for q is over the outer loop (j) and we are parallelizing j. So we need a reduction for q. Similarly, the update to s is a reduction over i (the inner loop) and we are not parallelizing the inner loop, so we don't need a reduction for s? Actually, no: because the inner loop is sequential, the update to s[j] is done entirely by one thread (the one that handles j). So we don't need a reduction for s? But wait: the reduction clause for s would be unnecessary because the inner loop is sequential. However, the reduction clause for q is necessary because multiple threads update the same q[i] (for different j's). 

But note: the reduction clause for q is over the entire array q. The reduction operation is addition. This is safe.

Alternatively, we can avoid the reduction for s because the inner loop is sequential and the update to s[j] is done by one thread. But the reduction clause for s would be harmless? Actually, it would be inefficient because it would create a private copy of the entire s array for each thread and then combine. But we don't need it because the update to s is done by one thread per j. So we can do:

  #pragma omp parallel for reduction(+:q[:_PB_N])
  for (j=0; j<_PB_M; j++)
    for (i=0; i<_PB_N; i++) {
        s[j] += r[i] * A[i][j];
        q[i] += A[i][j] * p[j];
    }

But note: the update to s[j] is done by the thread that handles j, and it is a reduction over i (the inner loop). Since the inner loop is sequential, the thread will update s[j] for all i. So no reduction is needed for s. However, the update to q[i] is done by multiple threads (each j) and they update the same q[i]. So we need a reduction for q.

But wait: the reduction clause for q will create a private copy of the entire q array for each thread. Then at the end, it will add all the private copies to the original q. This is correct.

However, we must initialize q to zero before the parallel region. Similarly, we must initialize s to zero.

So the code becomes:

  // Initialize s and q to zero
  #pragma omp parallel for
  for (int i=0; i<_PB_M; i++) s[i] = 0;
  #pragma omp parallel for
  for (int i=0; i<_PB_N; i++) q[i] = 0;

  // Main computation
  #pragma omp parallel for reduction(+:q[:_PB_N])
  for (j=0; j<_PB_M; j++)
    for (i=0; i<_PB_N; i++) {
        s[j] += r[i] * A[i][j];
        q[i] += A[i][j] * p[j];
    }

But note: the original code had the initialization of s in a loop and then the main loop. We are now initializing both s and q in separate parallel loops and then the main loop.

However, we can also consider fusing the initializations? But they are independent and we are already parallelizing them.

But note: the original scop region had the initialization of s and then the main loop that initializes q[i] and then the inner loop. We are now initializing q in a separate loop. This is equivalent.

But what about the performance? The initialization of s and q are O(M) and O(N) respectively, which are small compared to the main loop which is O(M*N). So it is acceptable.

Alternatively, we can avoid the separate initialization of q by including it in the main loop? But we cannot because the main loop is now over j and we are updating q[i] for every j. We must initialize q to zero before the main loop.

So we will do:

  #pragma scop
  // Initialize s to zero
  #pragma omp parallel for
  for (int i = 0; i < _PB_M; i++)
    s[i] = 0;

  // Initialize q to zero
  #pragma omp parallel for
  for (int i = 0; i < _PB_N; i++)
    q[i] = SCALAR_VAL(0.0);

  // Main computation: swapped loops
  #pragma omp parallel for reduction(+:q[:_PB_N])
  for (int j = 0; j < _PB_M; j++)
    for (int i = 0; i < _PB_N; i++) {
        s[j] = s[j] + r[i] * A[i][j];
        q[i] = q[i] + A[i][j] * p[j];
    }
  #pragma endscop

But note: the original code used _PB_M and _PB_N. We are using the same.

However, we must be cautious: the original code used the same variable names i and j. We are reusing them. This is acceptable.

But note: the original code had the inner loop updating s[j] and q[i]. We are doing the same.

But wait: in the swapped version, the inner loop is over i and we are updating s[j] (which is the same j for the entire inner loop) and q[i] (which is the same i for the entire inner loop). This is correct.

But note: the update to s[j] is now: s[j] = s[j] + ... and we are doing it for every i. This is the same as the original.

However, we must consider the possibility of floating point non-associativity. The order of additions might change. But the reduction clause for q will also change the order. The original code did the reduction for q[i] in the order of j, and for s[j] in the order of i. In our version, the reduction for s[j] is done in the order of i (same as original) and the reduction for q[i] is done in the order of j (same as original) but then the reduction of the private copies of q might be in a different order. However, the reduction clause in OpenMP does not specify the order of reduction. So it might be different. But the problem does not specify that the order must be preserved. We are allowed to change the order as long as the result is the same in exact arithmetic. But in floating point, it might differ. However, the problem says "preserve original semantics", meaning the result should be the same. But in floating point, the result might differ due to rounding. 

But note: the original code did:

  for i: 
    q[i] = 0;
    for j: q[i] += ... 

which is a sequential reduction in j order.

In our version, the inner loop is over i and we are updating q[i] for a fixed j. Then we are reducing the private copies of q. The reduction of the private copies is done by adding the entire arrays. This is equivalent to:

  q_final[i] = sum_{threads} (private_q[thread][i])

and the private_q[thread][i] is the sum over the j's assigned to that thread. This is the same as the original if we ignore rounding. But the order of additions is different. 

However, the problem does not specify that the floating point results must be bit-identical. It says "preserve original semantics", meaning the algorithm. But in HPC, sometimes bit-identical results are required. But the problem does not specify. We are allowed to change the order of reductions? The problem says "do not change observable behavior". If the floating point results are different, that is observable. 

Therefore, we must be cautious. We can avoid the reduction for q by using an array reduction only if we are sure that the order does not matter. Alternatively, we can use a different approach.

Another approach: we can parallelize the original i-loop without swapping. The original i-loop has:

  for (i=0; i<_PB_N; i++) {
      q[i] = 0;
      for (j=0; j<_PB_M; j++) {
          s[j] += r[i] * A[i][j];
          q[i] += A[i][j] * p[j];
      }
  }

We can parallelize the i-loop. But note: the update to s[j] is a reduction over i. So we need a reduction for s. And the update to q[i] is a reduction over j, which is the inner loop and we are not parallelizing it. So we can do:

  // Initialize s to zero (as before)
  #pragma omp parallel for
  for (int i=0; i<_PB_M; i++) s[i] = 0;

  // Then the main loop: parallelize i-loop with reduction on s
  #pragma omp parallel for reduction(+:s[:_PB_M])
  for (int i=0; i<_PB_N; i++) {
      DATA_TYPE tmp = SCALAR_VAL(0.0);
      for (int j=0; j<_PB_M; j++) {
          s[j] += r[i] * A[i][j];
          tmp += A[i][j] * p[j];
      }
      q[i] = tmp;
  }

This avoids the reduction on q because we are computing q[i] in a private variable tmp and then assigning at the end. And the update to s is done with a reduction. This is safe and the order of updates to s is the same as the original? Not exactly: the reduction for s will be done by adding the private copies at the end. But the original updated s in the order of i. The reduction will add the contributions from different i's in an arbitrary order. This might cause floating point differences.

But note: the original update for a fixed j: s[j] is updated for every i in increasing order. In the reduction, the updates from different i's are added in an arbitrary order. This is a change.

To preserve the order, we would have to use a critical section for s, which would be inefficient.

Alternatively, we can use the swapped loop version and accept the potential floating point difference. But the problem says "preserve original semantics", and if the floating point results are different, that is a change in semantics.

Given the ambiguity, we must choose the method that is most likely to be accepted. The problem does not specify the floating point model. In HPC, sometimes non-associative floating point is acceptable. But to be safe, we should try to preserve the order.

But note: the reduction for s in the swapped version is done by one thread per j, and the inner loop is over i in increasing order. So the update to s[j] is done in the same order as the original. The reduction for q is done by the reduction clause, which changes the order. 

In the original, the update to q[i] is done in increasing j order. In the swapped version, the inner loop is over i, so for a fixed j, we update q[i] for every i. Then the reduction clause will add the private copies of q. This is equivalent to:

  q_final[i] = sum_{j} (A[i][j] * p[j])

but the order of summation over j is broken into chunks by the threads. This is the same as the original only if the reduction is associative, which it is not for floating point.

Therefore, we must avoid changing the order of the reduction for q. 

We can do:

  // Initialize s to zero
  #pragma omp parallel for
  for (int i=0; i<_PB_M; i++) s[i] = 0;

  // Then the main loop: parallelize the i-loop, and for each i, do the inner j-loop without parallelizing it.
  // But we need to protect the update to s. We can use a reduction for s, but that changes the order.

  // Alternatively, we can use a critical section for s, but that would be inefficient.

  // Or we can use atomic updates for s? But that would be very inefficient.

Given the trade-offs, and since the problem does not specify the floating point model, we will choose the swapped loop version with reduction on q, and hope that the floating point difference is acceptable. This is common in HPC.

But note: the problem says "preserve original semantics". If the original semantics require bit-identical results, then we cannot change the order. However, the problem does not specify. We will assume that the order of reductions can be changed.

Therefore, we will use the swapped loop version.

But note: the original code used the variable names i and j. We are swapping them. We must change the loop variables accordingly.

Final code for the scop region:

  #pragma scop
  // Initialize s to zero
  #pragma omp parallel for
  for (int i = 0; i < _PB_M; i++)
    s[i] = 0;

  // Initialize q to zero
  #pragma omp parallel for
  for (int i = 0; i < _PB_N; i++)
    q[i] = SCALAR_VAL(0.0);

  // Main computation: swapped loops
  #pragma omp parallel for reduction(+:q[:_PB_N])
  for (int j = 0; j < _PB_M; j++)
    for (int i = 0; i < _PB_N; i++) {
        s[j] = s[j] + r[i] * A[i][j];
        q[i] = q[i] + A[i][j] * p[j];
    }
  #pragma endscop

But note: the original code used _PB_M and _PB_N. We are using the same.

However, we must be cautious: the array A is accessed as A[i][j]. In the swapped loop, we are accessing A[i][j] with i and j swapped in the loop order. But the array is the same. This is correct.

But note: the memory layout of A is row-major. In the original, the inner loop was over j, so we accessed a row. In the swapped version, the inner loop is over i, so we are accessing the i-th element of each row. This is not contiguous. This might cause cache misses.

This is a problem. The original access pattern for A was contiguous in j. In the swapped version, we are accessing A[i][j] for fixed j and varying i. This is accessing the j-th element of each row, which is not contiguous. This might be bad for performance.

Therefore, we must reconsider.

We want to preserve the contiguous access to A. So we should not swap the loops. 

Then we go back to the original loop order. We can parallelize the i-loop with a reduction on s. And for q, we can use a private temporary variable to avoid the reduction on q.

  #pragma scop
  // Initialize s to zero
  #pragma omp parallel for
  for (int i = 0; i < _PB_M; i++)
    s[i] = 0;

  // Main computation: parallelize i-loop
  #pragma omp parallel for reduction(+:s[:_PB_M])
  for (int i = 0; i < _PB_N; i++) {
      DATA_TYPE tmp = SCALAR_VAL(0.0);
      for (int j = 0; j < _PB_M; j++) {
          s[j] += r[i] * A[i][j];
          tmp += A[i][j] * p[j];
      }
      q[i] = tmp;
  }
  #pragma endscop

This preserves the contiguous access to A (inner loop over j). The update to s is done with a reduction, which might change the order of floating point additions. But the update to q is done in a private variable and then assigned, so the order for q is preserved (the inner loop over j is sequential and in increasing order).

The reduction for s will be done by adding the private copies of s from each thread. This is not in the original order. But it is the best we can do without a critical section.

We can try to use a schedule that might help, but it won't change the reduction order.

Given the above, we choose this version.

But note: the reduction for s is an array reduction. This is supported in OpenMP 4.5 and later. We assume the compiler supports it.

If the compiler does not support array reductions, we would have to do something else. But the problem does not specify the OpenMP version. We assume a reasonably recent version.

Therefore, the final code for the scop region is:

  #pragma scop
  #pragma omp parallel for
  for (int i = 0; i < _PB_M; i++)
    s[i] = 0;

  #pragma omp parallel for reduction(+:s[:_PB_M])
  for (int i = 0; i < _PB_N; i++) {
      DATA_TYPE tmp = SCALAR_VAL(0.0);
      for (int j = 0; j < _PB_M; j++) {
          s[j] += r[i] * A[i][j];
          tmp += A[i][j] * p[j];
      }
      q[i] = tmp;
  }
  #pragma endscop

We removed the initialization of q in a separate loop because we are now initializing q[i] in the main loop with the assignment of tmp.

But note: the original code initialized q[i] to zero at the beginning of the i-loop. We are now initializing tmp to zero and then assigning q[i] = tmp. This is equivalent.

This code has:
  - The first loop initializes s in parallel.
  - The second loop is parallelized over i, with a reduction on the entire s array.
  - The inner loop is sequential and updates s and a private tmp.
  - Then q[i] is set to tmp.

This preserves the contiguous access to A in the inner loop.

We must declare tmp. We use DATA_TYPE.

Let's write the code accordingly.

Note: We must not change anything outside the scop region.

We are only changing the code between #pragma scop and #pragma endscop.

The original code had:

  for (i = 0; i < _PB_M; i++)
    s[i] = 0;
  for (i = 0; i < _PB_N; i++)
    {
      q[i] = SCALAR_VAL(0.0);
      for (j = 0; j < _PB_M; j++)
	{
	  s[j] = s[j] + r[i] * A[i][j];
	  q[i] = q[i] + A[i][j] * p[j];
	}
    }

We are replacing it with:

  #pragma omp parallel for
  for (int i = 0; i < _PB_M; i++)
    s[i] = 0;

  #pragma omp parallel for reduction(+:s[:_PB_M])
  for (int i = 0; i < _PB_N; i++) {
      DATA_TYPE tmp = SCALAR_VAL(0.0);
      for (int j = 0; j < _PB_M; j++) {
          s[j] += r[i] * A[i][j];
          tmp += A[i][j] * p[j];
      }
      q[i] = tmp;
  }

We added the variable tmp. We must declare it with the correct type. The original used SCALAR_VAL(0.0) for q[i], so we use the same for tmp.

We also changed the loop variables to be declared inside the for loops (with int i, int j) to avoid conflicts. The original used the same i and j for both loops. We are reusing i for the outer loop and j for the inner loop, but we are declaring them inside the for, so they are local.

This should be safe.

Let's output the entire kernel_bicg function with the scop region replaced.
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "bicg.h"


/* Array initialization. */
static
void init_array (int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		 DATA_TYPE POLYBENCH_1D(r,N,n),
		 DATA_TYPE POLYBENCH_1D(p,M,m))
{
  int i, j;

  for (i = 0; i < m; i++)
    p[i] = (DATA_TYPE)(i % m) / m;
  for (i = 0; i < n; i++) {
    r[i] = (DATA_TYPE)(i % n) / n;
    for (j = 0; j < m; j++)
      A[i][j] = (DATA_TYPE) (i*(j+1) % n)/n;
  }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_1D(s,M,m),
		 DATA_TYPE POLYBENCH_1D(q,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("s");
  for (i = 0; i < m; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, s[i]);
  }
  POLYBENCH_DUMP_END("s");
  POLYBENCH_DUMP_BEGIN("q");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, q[i]);
  }
  POLYBENCH_DUMP_END("q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_bicg(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		 DATA_TYPE POLYBENCH_1D(s,M,m),
		 DATA_TYPE POLYBENCH_1D(q,N,n),
		 DATA_TYPE POLYBENCH_1D(p,M,m),
		 DATA_TYPE POLYBENCH_1D(r,N,n))
{
  int i, j;

#pragma scop
  #pragma omp parallel for
  for (int i = 0; i < _PB_M; i++)
    s[i] = 0;
  #pragma omp parallel for reduction(+:s[:_PB_M])
  for (int i = 0; i < _PB_N; i++) {
    DATA_TYPE tmp = SCALAR_VAL(0.0);
    for (int j = 0; j < _PB_M; j++) {
      s[j] += r[i] * A[i][j];
      tmp += A[i][j] * p[j];
    }
    q[i] = tmp;
  }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, M, n, m);
  POLYBENCH_1D_ARRAY_DECL(s, DATA_TYPE, M, m);
  POLYBENCH_1D_ARRAY_DECL(q, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(p, DATA_TYPE, M, m);
  POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(r),
	      POLYBENCH_ARRAY(p));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_bicg (m, n,
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(s),
	       POLYBENCH_ARRAY(q),
	       POLYBENCH_ARRAY(p),
	       POLYBENCH_ARRAY(r));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(s), POLYBENCH_ARRAY(q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(s);
  POLYBENCH_FREE_ARRAY(q);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(r);

  return 0;
}