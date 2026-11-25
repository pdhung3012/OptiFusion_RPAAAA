/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* durbin.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "durbin.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_1D(r,N,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      r[i] = (n+1-i);
    }
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
void kernel_durbin(int n,
		   DATA_TYPE POLYBENCH_1D(r,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n))
{
 DATA_TYPE z[N];
 DATA_TYPE alpha;
 DATA_TYPE beta;
 DATA_TYPE sum;

 int i,k;

#pragma scop
We are given the Durbin kernel code. The task is to optimize the code within the '#pragma scop' and '#pragma endscop' markers.

Original code in the scop region:

  y[0] = -r[0];
  beta = SCALAR_VAL(1.0);
  alpha = -r[0];

  for (k = 1; k < _PB_N; k++) {
    beta = (1-alpha*alpha)*beta;
    sum = SCALAR_VAL(0.0);
    for (i=0; i<k; i++) {
       sum += r[k-i-1]*y[i];
    }
    alpha = - (r[k] + sum)/beta;

    for (i=0; i<k; i++) {
       z[i] = y[i] + alpha*y[k-i-1];
    }
    for (i=0; i<k; i++) {
      y[i] = z[i];
    }
    y[k] = alpha;
  }

We are to apply loop transformations (tiling, fusion, skewing) first, then OpenMP directives.

Analysis:

1. Loop Tiling: The inner loops are over 'i' and are of length k (which grows). The outer loop is over k. The inner loops are not easily tileable because the inner loop bounds depend on the outer loop index. Also, the inner loops are short for small k. Tiling might not be beneficial and could complicate the code.

2. Loop Fusion: There are two adjacent loops over the same index range (i from 0 to k-1) that write to z and then copy from z to y. We can fuse these two loops? Actually, the first loop writes to z, and the second loop reads from z and writes to y. We cannot fuse them arbitrarily because the second loop uses the entire z array. However, note that the second loop is just a copy. We can eliminate the copy by writing directly to y? But wait, the first loop uses y in the computation (both y[i] and y[k-i-1]). We cannot overwrite y[i] until we have computed all the new values. Therefore, we must use a temporary array (z) to hold the new values and then copy. Alternatively, we can use a temporary array for the entire y? But the problem is that the update for y[i] depends on the old y[i] and y[k-i-1]. So we cannot update in-place without a temporary.

However, note that the two loops are:

  for (i=0; i<k; i++) {
      z[i] = y[i] + alpha*y[k-i-1];
  }
  for (i=0; i<k; i++) {
      y[i] = z[i];
  }

We can fuse these two loops by writing directly to y? But that would be incorrect because the second loop uses the old y for the entire range. However, observe that the first loop writes to z and the second loop reads from z and writes to y. We can avoid the copy by using a temporary array for the entire y? Actually, we are already using z as a temporary. The copy is necessary because we are updating the entire y array for indices 0 to k-1 and we need the old values for the entire range.

But note: the update for y[i] uses y[k-i-1]. Since k-i-1 might be less than i, we cannot update in-place without a temporary. So the two loops are necessary.

However, we can fuse the two loops by writing the result of the first loop directly to y? That would be incorrect because the value of y[k-i-1] might be overwritten by an earlier iteration. Therefore, we must use a temporary array.

But we can avoid the copy by using a double buffer? The algorithm already uses two arrays: y and z. We are using z as a temporary for the update. Then we copy z to y. We can avoid the copy by swapping the roles of y and z? However, the next iteration uses y again. So we cannot avoid the copy without changing the entire algorithm.

Alternatively, we can use a temporary array for the entire y? But that would require O(n) extra space and we are already using z (which is of size n). So we are already using O(n) extra space.

Therefore, we cannot fuse the two loops without changing the algorithm. We leave them as is.

3. Loop Skewing: The inner loops are simple and do not have carried dependencies that would benefit from skewing. The outer loop has a carried dependence on beta and alpha, and the inner loops are reductions and then a simple update. Skewing does not seem applicable.

So we move to OpenMP.

OpenMP opportunities:

- The outer loop (k) is not parallel because each iteration depends on the previous (beta and alpha are carried, and y is updated and used in the next iteration).

- The inner loops (over i) are parallel? Let's check:

  a. The first inner loop (sum reduction) is a reduction. We can parallelize it with a reduction clause.

  b. The second inner loop (computing z) is a parallel loop because each iteration writes to a different z[i] and reads from y (which is not written in this loop). However, note that the same y array is read by multiple iterations. But since we are only reading, it is safe.

  c. The third inner loop (copying z to y) is also parallel because each iteration writes to a different y[i] and reads from z[i].

But note: the second and third loops are over the same range and we can parallelize them together? Actually, we can parallelize the second and third loops separately.

However, the first inner loop (the reduction) is the most expensive and a good candidate for parallelization.

Plan:

- Parallelize the reduction loop (the first inner loop) with OpenMP and a reduction clause.

- The other two inner loops (the update and the copy) are also parallel and we can parallelize them as well.

But note: the reduction loop is the only one that is computationally heavy? The other two are memory-bound and might not benefit as much. However, we can try to parallelize all three inner loops.

But we must be cautious: the reduction loop is a sum reduction. We can use:

  #pragma omp parallel for reduction(+:sum)

For the other two loops, we can use:

  #pragma omp parallel for

But note: the two loops (the update and the copy) are independent and we can run them in parallel? Actually, they are sequential: the update must complete before the copy. So we cannot run them at the same time.

Alternatively, we can fuse the update and the copy? But we cannot because the update writes to z and the copy reads from z. We can avoid the copy by writing directly to y? As discussed, that is not safe.

But note: the update loop and the copy loop are both over the same range and we can combine them? Actually, we can write:

  for (i=0; i<k; i++) {
      y[i] = y[i] + alpha*y[k-i-1];
  }

But that is incorrect because we are overwriting y[i] and then using the old y[k-i-1] for other i. However, if we do:

  for (i=0; i<k; i++) {
      double temp = y[i] + alpha*y[k-i-1];
      y[i] = temp;
  }

This is still incorrect because when i is small, k-i-1 is large and might be greater than i. Then when we update y[i], we are using the old y[k-i-1] which is correct. But when we update y[k-i-1] (if k-i-1 < k) then we are using the old y[i] which is also correct. However, if k-i-1 < i, then when we update y[i] we are using the old y[k-i-1] (which is correct) and when we update y[k-i-1] we are using the old y[i] (which is also correct). So it is safe? 

But note: the update for y[i] uses y[k-i-1]. The update for y[k-i-1] uses y[i]. So if we update y[i] first, then when we update y[k-i-1] we are using the updated y[i]? That would be incorrect because the original algorithm uses the old y[i] for both.

Therefore, we must use a temporary array. So we cannot avoid the two loops.

But we can avoid the copy by using a temporary array for the entire y? We are already using z. So we are stuck with the two loops.

So we will parallelize the three inner loops separately.

However, the reduction loop is the most important. The other two loops are memory-bound and might not scale well. But we can try.

But note: the reduction loop is a sum over i. We can parallelize it with a reduction.

The update loop (writing to z) is parallel because each iteration writes to a different z[i] and reads shared y (which is not written in this loop). So we can use:

  #pragma omp parallel for

Similarly, the copy loop is parallel.

But we must be cautious: the reduction loop is the only one that requires a reduction. The others are straightforward.

Also, we must consider the overhead of multiple parallel regions. We are inside a k loop. For each k, we would have three parallel regions. That might be expensive.

We can use a single parallel region for the entire k iteration? But the k loop is sequential. We can use:

  #pragma omp parallel
  {
    for (k=1; k<_PB_N; k++) {
      ... // but then we have to manage the inner loops with OpenMP worksharing inside
    }
  }

But the k loop is sequential, so we cannot parallelize it. We can only parallelize the inner loops.

Alternatively, we can use a single parallel region and then use worksharing for the inner loops. This avoids the overhead of creating a new team for each inner loop.

So:

  #pragma omp parallel
  {
    for (k=1; k<_PB_N; k++) {
      beta = (1-alpha*alpha)*beta;
      sum = 0.0;
      #pragma omp for reduction(+:sum)
      for (i=0; i<k; i++) {
        sum += r[k-i-1]*y[i];
      }
      alpha = - (r[k] + sum)/beta;

      #pragma omp for
      for (i=0; i<k; i++) {
        z[i] = y[i] + alpha*y[k-i-1];
      }
      #pragma omp for
      for (i=0; i<k; i++) {
        y[i] = z[i];
      }
      y[k] = alpha;
    }
  }

But note: the variables beta, alpha, and sum are shared? We must be cautious.

- beta: updated in each k iteration and used in the next. It is carried. We must make sure that the update to beta is done by one thread and then broadcast? Actually, the update to beta is done at the beginning of the k iteration and then used in the computation of alpha. The update to beta is:

      beta = (1-alpha*alpha)*beta;

This uses the old beta and the old alpha (from the previous k). Then we compute sum and then alpha. Then we update y and set y[k]=alpha.

The variables beta and alpha are shared and updated in the sequential part. We must make sure that the sequential part is executed by only one thread. We can use a master or single directive? But note: the entire k loop is inside the parallel region and we want the sequential parts to be executed by the master thread.

We can restructure:

  #pragma omp parallel
  {
    for (k=1; k<_PB_N; k++) {
      // Sequential part: update beta and compute alpha
      // We do this in a single thread
      #pragma omp single
      {
        beta = (1-alpha*alpha)*beta;
        sum = 0.0;
        // But note: the reduction loop is parallel and we are about to do it. We cannot do the reduction in the single.
        // So we break the sequential part: we update beta, then we do the reduction in parallel, then we compute alpha in the single.
        // Then we do the two parallel loops, and then set y[k] in the single.
      }

      // But the reduction loop must be done by all threads? Actually, we can do:
      //   update beta (single)
      //   then reduction loop (parallel for)
      //   then compute alpha (single)
      //   then update z (parallel for)
      //   then copy to y (parallel for)
      //   then set y[k] (single)

      // However, the reduction loop uses y which is being updated in the previous k? But note: the reduction loop uses the old y (from the previous k) because we haven't updated y for this k yet. So it is safe.

      // But the reduction loop uses r and y. They are shared and read-only in the reduction.

      // So:

      #pragma omp single
      {
        beta = (1-alpha*alpha)*beta;
      }
      sum = 0.0;
      #pragma omp for reduction(+:sum)
      for (i=0; i<k; i++) {
        sum += r[k-i-1]*y[i];
      }
      #pragma omp single
      {
        alpha = - (r[k] + sum)/beta;
      }

      #pragma omp for
      for (i=0; i<k; i++) {
        z[i] = y[i] + alpha*y[k-i-1];
      }
      #pragma omp for
      for (i=0; i<k; i++) {
        y[i] = z[i];
      }
      #pragma omp single
      {
        y[k] = alpha;
      }
    }
  }

But note: the variable 'sum' is now shared? We are using a reduction, so that is safe.

However, the variable 'alpha' is computed in the single and then used in the parallel loops. We must make sure that the parallel loops see the updated alpha. Since the single is executed by one thread and then the parallel loops are started, and the alpha is written before the parallel loops, it should be visible.

But we must declare the variables appropriately. We are using:

  DATA_TYPE z[N];
  DATA_TYPE alpha;
  DATA_TYPE beta;
  DATA_TYPE sum;

These are declared outside the parallel region. We must use shared for these? But note: the parallel region is inside the function and these variables are shared by default? In OpenMP, variables in the outer scope are shared by default.

But we must be cautious: the reduction variable 'sum' is private in the reduction, but we are using it in the single after the reduction. The reduction clause will make the reduction variable shared and then combine the private copies. But note: the reduction clause on a for loop creates a private copy for each thread and then combines them. Then the shared variable is updated. So after the parallel for, the shared 'sum' will have the total.

But we are using the same 'sum' variable. So it is safe.

However, the variables alpha and beta are updated in the single and then used in the parallel loops. They are shared and the updates in the single are visible to all threads.

But note: the initial values of alpha and beta are set outside the parallel region. We must make sure they are shared.

Also, the arrays r, y, z are shared.

But note: the array y is updated in the parallel loop (the copy loop). We must be cautious about false sharing? But the copy loop is parallel and each thread writes to a different y[i]. So it is safe.

However, the update loop (writing to z) is also parallel and safe.

But note: the reduction loop reads y and r. They are not written in that loop, so it is safe.

But there is a problem: the reduction loop uses y, and then the update loop uses y. But the update loop uses the same y that was used in the reduction? Yes, because we haven't updated y for this k yet. So it is safe.

But note: the reduction loop uses y[i] for i in [0, k-1]. The update loop uses y[i] and y[k-i-1] for i in [0, k-1]. The copy loop then updates y[i] for i in [0, k-1]. So the reduction and update loops use the old y, and the copy loop writes the new y.

So the code is correct.

But we must consider the overhead of the multiple worksharing directives inside the k loop. For each k, we have:

  single
  for (reduction)
  single
  for (update)
  for (copy)
  single

This might be heavy for small k. But for large k, the inner loops are long and the overhead is amortized.

Alternatively, we can use a single parallel region and then use worksharing for the inner loops without the single for the sequential parts? But the sequential parts (beta update, alpha computation, and y[k] assignment) must be done by one thread.

We can do:

  #pragma omp parallel
  {
    for (k=1; k<_PB_N; k++) {
      // Update beta: must be done by one thread
      #pragma omp single
      beta = (1-alpha*alpha)*beta;

      // Reduction loop
      sum = 0.0;
      #pragma omp for reduction(+:sum)
      for (i=0; i<k; i++) {
        sum += r[k-i-1]*y[i];
      }

      // Compute alpha: must be done by one thread
      #pragma omp single
      alpha = - (r[k] + sum)/beta;

      // Update z: parallel
      #pragma omp for
      for (i=0; i<k; i++) {
        z[i] = y[i] + alpha*y[k-i-1];
      }

      // Copy to y: parallel
      #pragma omp for
      for (i=0; i<k; i++) {
        y[i] = z[i];
      }

      // Set y[k]: must be done by one thread
      #pragma omp single
      y[k] = alpha;
    }
  }

But note: the variable 'sum' is shared and we are using a reduction. The reduction clause will take care of the private copies and the final shared value.

However, the variable 'alpha' is updated in the single and then used in the update loop. The update loop is a worksharing loop that uses alpha. Since the single is executed by one thread and then the worksharing loop is started, the value of alpha should be visible to all threads.

But we must declare the variables appropriately. We are using default shared? We should use default(none) to be safe? But the problem says we can use data sharing clauses.

We can add:

  #pragma omp parallel default(none) shared(r, y, z, beta, alpha, sum, _PB_N)

But note: the variables k and i are declared in the outer scope. We must make them private? Actually, k is the loop variable of the outer loop and i is the loop variable of the inner loops. We can make them private by declaring them inside the parallel region? But they are declared outside.

We can do:

  int i, k;

  #pragma omp parallel default(none) shared(r, y, z, beta, alpha, sum, _PB_N) private(i, k)
  {
    for (k=1; k<_PB_N; k++) {
      ...
    }
  }

But note: the outer loop variable k is private by default in a for loop? Actually, in a for loop, the loop variable is private. But we are using a for loop inside the parallel region without a worksharing directive? The outer loop is not a worksharing loop. We are using a for loop that is executed by all threads? That would be a problem.

We must make the outer loop (k) executed by only one thread? Or we must use a worksharing loop for k? But k is sequential.

We can use:

  #pragma omp parallel
  {
    #pragma omp single
    {
      for (k=1; k<_PB_N; k++) {
        ... // the entire body, including the worksharing for the inner loops
      }
    }
  }

But then the inner worksharing loops would be inside the single, and only one thread would execute them? That is not what we want.

Alternatively, we can use a worksharing loop for the outer k? But k is sequential. We cannot.

So we must have the outer k loop executed by all threads? But that would be a race because each thread would try to update beta, alpha, and y[k]. We don't want that.

Therefore, we must have the outer k loop executed by one thread, and then the inner loops are workshared. We can do:

  #pragma omp parallel
  {
    #pragma omp single
    {
      for (k=1; k<_PB_N; k++) {
        // Update beta: done by the single thread
        beta = (1-alpha*alpha)*beta;

        // Reduction loop: we want to parallelize this
        // But we are in a single, so we cannot use a worksharing directive here? 
        // We can use a parallel for inside the single? That would create a new team? We don't want that.

        // Instead, we can use a worksharing directive that is executed by the team of threads that are in the parallel region? 
        // But we are in a single, so only one thread is active. We need to use a worksharing directive that is outside the single.

        // We can restructure: do the sequential parts in the single, and then use a worksharing directive for the inner loops that are outside the single? But the inner loops are inside the k loop.

        // We can do:

        //   update beta (in the single)
        //   then do the reduction loop with a worksharing for (but we are in a single, so we cannot)

        // This is a problem.

        // Alternative: we can use a master directive for the sequential parts and then worksharing for the inner loops without a single for the entire k loop.

        // We can do:

        //   #pragma omp master
        //   {
        //      beta = ...;
        //   }
        //   #pragma omp barrier   // to make sure beta is updated
        //   ... reduction loop ...

        // But that is messy.

        // Another idea: we can use a single for the entire k loop body? But then the inner worksharing loops would be executed by the single thread? That is not what we want.

        // We must break the k loop body into sequential and parallel parts and use worksharing directives that are not nested in a single.

        // We can do:

        //   for (k=1; k<_PB_N; k++) {
        //      #pragma omp master
        //      {
        //          beta = ...;
        //      }
        //      #pragma omp barrier
        //      ... reduction loop with worksharing ...
        //      #pragma omp master
        //      {
        //          alpha = ...;
        //      }
        //      #pragma omp barrier
        //      ... update loop with worksharing ...
        //      ... copy loop with worksharing ...
        //      #pragma omp master
        //      {
        //          y[k] = alpha;
        //      }
        //   }

        // But the barriers are expensive.

        // Given the complexity, and since the k loop is sequential, it might be simpler to parallelize only the reduction loop and leave the other two loops sequential? Or parallelize the reduction loop and the update loop, but not the copy loop? The copy loop is O(k) and might be cheap.

        // But the problem says to prioritize parallelism.

        // We decide to parallelize only the reduction loop and the update loop, and leave the copy loop sequential? Or we can parallelize the reduction loop and then do the update and copy in one parallel loop? But the update and copy are two separate loops.

        // Alternatively, we can use a single parallel region for the entire k loop and use worksharing for the inner loops without the single for the sequential parts, but then we must protect the sequential parts with a master or single and use barriers.

        // We choose to use:

        #pragma omp parallel
        {
          for (k=1; k<_PB_N; k++) {
            #pragma omp master
            {
              beta = (1-alpha*alpha)*beta;
            }
            #pragma omp barrier

            sum = 0.0;
            #pragma omp for reduction(+:sum)
            for (i=0; i<k; i++) {
              sum += r[k-i-1]*y[i];
            }
            #pragma omp barrier

            #pragma omp master
            {
              alpha = - (r[k] + sum)/beta;
            }
            #pragma omp barrier

            #pragma omp for
            for (i=0; i<k; i++) {
              z[i] = y[i] + alpha*y[k-i-1];
            }
            #pragma omp barrier

            // We can do the copy in parallel too
            #pragma omp for
            for (i=0; i<k; i++) {
              y[i] = z[i];
            }
            #pragma omp barrier

            #pragma omp master
            {
              y[k] = alpha;
            }
            #pragma omp barrier
          }
        }

        // But the barriers are very expensive, especially for small k.

        // Given the overhead, and since the reduction loop is the most expensive, we might only parallelize the reduction loop and leave the rest sequential.

        // We decide to parallelize only the reduction loop.

        // So:

        for (k=1; k<_PB_N; k++) {
          beta = (1-alpha*alpha)*beta;
          sum = 0.0;
          #pragma omp parallel for reduction(+:sum)
          for (i=0; i<k; i++) {
            sum += r[k-i-1]*y[i];
          }
          alpha = - (r[k] + sum)/beta;

          for (i=0; i<k; i++) {
            z[i] = y[i] + alpha*y[k-i-1];
          }
          for (i=0; i<k; i++) {
            y[i] = z[i];
          }
          y[k] = alpha;
        }

        // This creates a new parallel region for each k. For large k, the reduction loop is long and the overhead is amortized. For small k, the overhead might be high, but the reduction loop is short anyway.

        // We can use a schedule for the reduction loop? The default is static, which is fine.

        // Also, we can use SIMD for the inner loops? The compiler might vectorize the reduction loop and the update and copy loops.

        // We can add:

        #pragma omp parallel for reduction(+:sum) simd
        for (i=0; i<k; i++) {
          sum += r[k-i-1]*y[i];
        }

        // And for the update and copy loops, we can use:

        #pragma omp simd
        for (i=0; i<k; i++) {
          z[i] = y[i] + alpha*y[k-i-1];
        }
        #pragma omp simd
        for (i=0; i<k; i++) {
          y[i] = z[i];
        }

        // But note: the reduction loop is already parallelized with OpenMP, and we are adding simd. This is a parallel for simd.

        // However, the update and copy loops are not parallelized with OpenMP, only with SIMD.

        // This might be a good compromise.

        // But the problem says to prioritize parallelism. The reduction loop is the most expensive, so we focus on that.

        // We decide to use:

        for (k=1; k<_PB_N; k++) {
          beta = (1-alpha*alpha)*beta;
          sum = 0.0;
          #pragma omp parallel for reduction(+:sum) simd
          for (i=0; i<k; i++) {
            sum += r[k-i-1]*y[i];
          }
          alpha = - (r[k] + sum)/beta;

          #pragma omp simd
          for (i=0; i<k; i++) {
            z[i] = y[i] + alpha*y[k-i-1];
          }
          #pragma omp simd
          for (i=0; i<k; i++) {
            y[i] = z[i];
          }
          y[k] = alpha;
        }

        // We also add 'simd' to the reduction loop to allow vectorization.

        // We must ensure that the arrays are aligned? We can use 'aligned' if we know the alignment, but we don't. We can omit it.

        // Also, we can use 'safelen' if needed, but we omit.

        // This version creates a new parallel region for each k. We can use a single parallel region for the entire k loop? But then we would have to use worksharing for the reduction loop and then the sequential parts would require synchronization.

        // Given the trade-offs, we choose the version with a parallel region for each reduction loop.

        // But note: the reduction loop is the only one we are parallelizing with OpenMP. The other two loops are only SIMD.

        // We also note that the reduction loop might be the bottleneck.

        // We also add 'default(none)' to the parallel for to be safe? But the variables are in the outer scope. We can specify:

          #pragma omp parallel for reduction(+:sum) default(none) shared(r,y,k) private(i)

        // But k is the outer loop index, and it is shared? Actually, k is the same for all threads in this parallel region. And i is private.

        // We can do:

          #pragma omp parallel for reduction(+:sum) default(none) shared(k, r, y) private(i)

        // But note: the variables r and y are shared, and k is shared.

        // However, the variable 'sum' is the reduction variable.

        // But the initial value of sum is set to 0.0 before the parallel region. The reduction will add to that.

        // This should be safe.

        // But note: the variable 'k' is used in the loop bound and in the index 'k-i-1'. It is shared and read-only.

        // So we write:

          #pragma omp parallel for reduction(+:sum) default(none) shared(k, r, y) private(i)

        // But the variable 'k' is declared in the outer scope. We are capturing it by sharing.

        // Alternatively, we can use firstprivate for k? But it is not necessary because it is not modified.

        // We decide to use shared.

        // However, the OpenMP standard says that variables in the outer scope are shared by default. We are using default(none) so we must list them.

        // We also note that the variable 'i' is private.

        // So the final code for the scop region is:

          y[0] = -r[0];
          beta = SCALAR_VAL(1.0);
          alpha = -r[0];

          for (k = 1; k < _PB_N; k++) {
            beta = (1-alpha*alpha)*beta;
            sum = SCALAR_VAL(0.0);
            #pragma omp parallel for reduction(+:sum) default(none) shared(k, r, y) private(i)
            for (i=0; i<k; i++) {
              sum += r[k-i-1]*y[i];
            }
            alpha = - (r[k] + sum)/beta;

            #pragma omp simd
            for (i=0; i<k; i++) {
              z[i] = y[i] + alpha*y[k-i-1];
            }
            #pragma omp simd
            for (i=0; i<k; i++) {
              y[i] = z[i];
            }
            y[k] = alpha;
          }

        // We also add SIMD to the reduction loop? We already have 'simd' in the pragma? Actually, we have 'parallel for simd'? We can combine:

          #pragma omp parallel for reduction(+:sum) simd

        // But the 'simd' clause is for the worksharing loop. We can do:

          #pragma omp parallel for reduction(+:sum) simd

        // This allows the loop to be both parallelized and vectorized.

        // But note: the loop is parallelized across threads and within each thread it is vectorized.

        // We can also add 'aligned' if we know the alignment, but we don't.

        // We decide to use:

          #pragma omp parallel for reduction(+:sum) simd

        // without default(none) for simplicity? But the problem says to use data sharing clauses. We can use default(none) to be safe.

        // We'll use:

          #pragma omp parallel for reduction(+:sum) default(none) shared(k, r, y) private(i) simd

        // But the 'simd' clause is not needed if the compiler can vectorize automatically? We add it to give a hint.

        // Alternatively, we can rely on the compiler to vectorize the inner loops without OpenMP SIMD? But the problem says to use SIMD directives.

        // We'll keep the SIMD directives.

        // However, the reduction loop is now both parallel and vectorized.

        // We also note that the two other loops are only vectorized.

        // This should be a good optimization.

        // Let's write the final code.

        Note: We must preserve the '#pragma scop' and '#pragma endscop' markers.

        We only change the code between them.

        We also note that the variable 'i' is declared in the outer scope. We are using it as private in the OpenMP loop, which is safe.

        We also note that the variable 'k' is the loop index of the outer loop, and we are using it in the inner loop. We are sharing it in the OpenMP region, which is safe because it is not modified.

        We also note that the arrays r and y are shared.

        We also note that the variable 'sum' is the reduction variable.

        We also note that the variable 'alpha' is used in the inner loops. It is computed before the inner loops and is shared. In the reduction loop, it is not used. In the update loop, it is used. But the update loop is after the reduction loop and is not in the parallel region. So it is safe.

        We also note that the variable 'beta' is updated in the outer loop and used in the reduction loop? No, it is not used in the reduction loop. It is used only in the sequential part.

        So the code is safe.

        We write the optimized code.
#pragma scop
 y[0] = -r[0];
 beta = SCALAR_VAL(1.0);
 alpha = -r[0];

 for (k = 1; k < _PB_N; k++) {
   beta = (1-alpha*alpha)*beta;
   sum = SCALAR_VAL(0.0);
   #pragma omp parallel for reduction(+:sum) default(none) shared(k, r, y) private(i) simd
   for (i=0; i<k; i++) {
      sum += r[k-i-1]*y[i];
   }
   alpha = - (r[k] + sum)/beta;

   #pragma omp simd
   for (i=0; i<k; i++) {
      z[i] = y[i] + alpha*y[k-i-1];
   }
   #pragma omp simd
   for (i=0; i<k; i++) {
     y[i] = z[i];
   }
   y[k] = alpha;
 }
#pragma endscop
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(r));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_durbin (n,
		 POLYBENCH_ARRAY(r),
		 POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(r);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
