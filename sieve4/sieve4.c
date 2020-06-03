/* 优化1：去掉偶数
   优化2：消除广播 
   优化3：cache优化 目标机E5-2660v4 L3缓存35MB  本机G4600 L3缓存3MB*/

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define CACHE_SIZE 36700160       /* 35MB */


int main(int argc, char *argv[])
{
    int   count;              /* Local prime count */
    double elapsed_time;      /* Parallel execution time */
    /*uint64_t    first;        /* Index of first multiple */
    int    global_count;      /* Global prime count */
    uint64_t    high_value;   /* Highest value on this proc */
    uint64_t    i;
    int    id;                /* Process ID number */
    uint64_t    index;        /* Index of current prime */
    uint64_t    low_value;    /* Lowest value on this proc */
    char  *marked;            /* Portion of 2,...,'n' */
    uint64_t    n;            /* Sieving from 2, ..., 'n' */
    int    p;                 /* Number of processes */
    uint64_t    proc0_size;   /* Size of proc 0's subarray */
    uint64_t    prime;        /* Current prime */
    uint64_t    size;         /* Elements in 'marked' */
    uint64_t    local_index;  /* Index of current local prime ofr each proc */
    uint64_t    local_prime;  /* Current local prime for each proc */
    char  *local_primes;      /* Prime numbers for each proc of 2,...,'sqrt(high_value)' */
    uint64_t    local_size;   /* Elements in 'local_primes' */

    int block_step;                /* Current block number */
    int block_step_num;            /* Number of blocks could 'n' be divided */
    uint64_t    block_first;       /* Index of first multiple in this block */
    uint64_t    block_size;        /* Size of a block's subarray */
    uint64_t    block_low_value;   /* Lowest value of this block */
    uint64_t    block_high_value;  /* Highest value of this block */


    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoi(argv[1]);

    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */

    low_value = 2 + id * (n - 1) / p;
    high_value = 1 + (id + 1)*(n - 1) / p;
    size = (high_value - low_value + 1) / 2;
    local_size = ((int)sqrt((double)high_value) + 1) / 2 - 1;  /* just like the sieve3.c */

    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = (n - 1) / p;

    if ((2 + proc0_size) < (int)sqrt((double)n)) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate this process's share of the array. */

    marked = (char *)malloc(size);
    local_primes = (char *)malloc(local_size);

    if (marked == NULL || local_primes == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    for (i = 0; i < size; i++) marked[i] = 0;
    for (i = 0; i < local_size; i++) local_primes[i] = 0;

    local_index = 0;
    local_prime = 3;
    do {
        for (i = (local_prime * local_prime - 3) / 2; i < local_size; i += local_prime)
            local_primes[i] = 1;
        while (local_primes[++local_index]);
        local_prime = 2 * local_index + 3;
    } while (local_prime * local_prime <= (int)sqrt((double)high_value));

    block_size = CACHE_SIZE / sizeof(uint64_t);
    block_step_num = size / block_size;
    if (size > block_step_num * block_size) {
        block_step_num += 1;
    }
    /*printf("size=%lld, block_size=%lld, block_step_num=%d\n", size, block_size, block_step_num);*/
    for (block_step = 0; block_step < block_step_num; block_step++) {
        block_low_value = low_value + block_size * block_step * 2;
        block_high_value = MIN(high_value, block_low_value + block_size * 2);
        /*printf("step %d, block_low_value=%lld, block_high_value=%lld\n", block_step, block_low_value, block_high_value);*/
        index = 0;
        prime = 3;
        do {
            if (prime * prime > block_low_value)
                block_first = (prime * prime - block_low_value) / 2;
            else {
                if (!(block_low_value % prime)) block_first = block_low_value % 2 ? 0 : prime / 2;
                else {
                    block_first = prime - block_low_value % prime;
                    if (!((block_low_value + block_first) % 2)) block_first = (prime + block_first) / 2;
                    else block_first /= 2;
                }
            }
            if (block_step < block_step_num - 1) {
                for (i = block_first; i < block_size; i += prime) {
                    marked[i + block_size * block_step] = 1;
                }
            }
            else {  /* the last block may be not full */
                for (i = block_first; i < block_size; i += prime) {
                    if (i + block_size * block_step >= size) {
                        break;
                    }
                    marked[i + block_size * block_step] = 1;
                }
            }
            while (local_primes[++index]);
            prime = index * 2 + 3;
        } while (prime * prime <= block_high_value);
    }

    count = 0;
    for (i = 0; i < size; i++) {
        if (!marked[i]) count++;
    }
    /* printf("process %d:\n", id);
     for (int i = 0; i < size; i++) printf("%d ", marked[i]);
     printf("\n");
     printf("count of process %d: %d    range:from %d to %d \n", id, count, low_value, high_value); */
    if (p > 1) MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
        0, MPI_COMM_WORLD);
    else global_count = count;
    global_count++;

    /* Stop the timer */

    elapsed_time += MPI_Wtime();


    /* Print the results */

    if (!id) {
        printf("There are %d primes less than or equal to %lld\n",
            global_count, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}
