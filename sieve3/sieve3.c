/* 优化1：去掉偶数
   优化2：消除广播 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main(int argc, char *argv[])
{
    int    count;        /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    int    first;        /* Index of first multiple */
    int    global_count; /* Global prime count */
    int    high_value;   /* Highest value on this proc */
    int    i;
    int    id;           /* Process ID number */
    int    index;        /* Index of current prime */
    int    low_value;    /* Lowest value on this proc */
    char  *marked;       /* Portion of 2,...,'n' */
    int    n;            /* Sieving from 2, ..., 'n' */
    int    p;            /* Number of processes */
    int    proc0_size;   /* Size of proc 0's subarray */
    int    prime;        /* Current prime */
    int    size;         /* Elements in 'marked' */
    int    local_index;  /* Index of current local prime ofr each proc */
    int    local_prime;  /* Current local prime for each proc */
    char  *local_primes; /* Prime numbers for each proc of 2,...,'sqrt(high_value)' */
    int    local_size;   /* Elements in 'local_primes' */

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
    local_size = ((int)sqrt((double)high_value) + 1) / 2 - 1;  /* 保留奇数，同优化2方法 */

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


    index = 0;
    prime = 3;
    do {
        if (prime * prime > low_value)
            first = (prime * prime - low_value) / 2;
        else {
            if (!(low_value % prime)) first = low_value % 2 ? 0 : prime / 2;
            else {
                first = prime - low_value % prime;
                if (!((low_value + first) % 2)) first = (prime + first) / 2;
                else first /= 2;
            }
        }
        for (i = first; i < size; i += prime) {
            marked[i] = 1;
        }
        /*if (!id) {
            while (marked[++index]);
            prime = index * 2 + 3;
        }*/
        while (local_primes[++index]);
        prime = index * 2 + 3;
        /* if (p > 1) MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);*/
    } while (prime * prime <= high_value);
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
        printf("There are %d primes less than or equal to %d\n",
            global_count, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}
