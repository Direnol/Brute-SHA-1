#include "sha-1.h"

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr, "Use %s str max_lenght alpha\n", argv[0]);
        return EXIT_FAILURE;
    }
    int rank, commsize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    double_t ttotal = -MPI_Wtime();

    // init parametrs
    size_t L; // max lenght
    sha_t sha;
    ull_t N;
    if (init_alpha(argv, &L, &N, &sha)) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // count of combinations

    char *maybe = calloc(L + 1, sizeof(char));
    if (!maybe) {
        perror("Alloc maybe");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    ull_t start_comb, comb_per_proc;
    comb_per_proc = get_block_size(N, rank, commsize);
    start_comb = get_start_block(N, (ull_t) rank, commsize);
    size_t n = comb(maybe, start_comb, L, &sha);

#if DEBUG
    {
        debug(printf("Start [%s] index %llu count %llu rank %d commsize %d\n",
                     maybe,
                     start_comb,
                     comb_per_proc,
                     rank,
                     commsize))
    }
#endif

    backtrack(maybe, &sha, n, 0, comb_per_proc);

#if DEBUG
    {
        debug(printf("prev [%d] n %lu %llu\n", rank, n, sha.c))
    }
#endif

    for (size_t i = n + 1; i <= L && !sha.found && sha.c < comb_per_proc; ++i) {
        InitString(maybe, i, sha.alphabet[0]);
        backtrack(maybe, &sha, i, 0, comb_per_proc);
    }


    ttotal += MPI_Wtime();

    if (sha.found) {
        printf("[%d|%d] %lf\n", rank, commsize, ttotal);
        freeAll(maybe, sha.alpha, sha.mem_cur_sha, sha.mem_sha, sha.alphabet, NULL);
        MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
    }

    double tmax;
    MPI_Reduce(&ttotal, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("%d %d %.6f\n", rank, commsize, ttotal);
    }
    freeAll(maybe, sha.alpha, sha.mem_cur_sha, sha.mem_sha, sha.alphabet, NULL);
    MPI_Finalize();
    return EXIT_SUCCESS;
}