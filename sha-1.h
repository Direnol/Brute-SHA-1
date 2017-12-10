#ifndef BRUTESHA_1_SHA_1_H
#define BRUTESHA_1_SHA_1_H

#include <math.h>
#include <openssl/sha.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <stdarg.h>

#define DEBUG 0
#define debug(func) int rrank, ccommsize;                                                    \
                    MPI_Comm_size(MPI_COMM_WORLD, &ccommsize);                               \
                    MPI_Comm_rank(MPI_COMM_WORLD, &rrank);                                  \
                    int prev = rrank ? rrank - 1 : MPI_PROC_NULL;                             \
                    int next = rrank != ccommsize - 1 ? rrank + 1 : MPI_PROC_NULL;             \
                    int tmp = 1;                                                            \
                    MPI_Recv(&tmp, 0, MPI_INT, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
                    func;                                                                    \
                    MPI_Send(&tmp, 0, MPI_INT, next, 0, MPI_COMM_WORLD);                    \
                    MPI_Barrier(MPI_COMM_WORLD);

typedef unsigned long long ull_t;

typedef struct {
    unsigned char *mem_sha;
    unsigned char *mem_cur_sha;
    size_t n;
    int *alpha; //alpha[c] == index -> alphabet[index] == c
    char *alphabet;
    char found; // 0 -  not found | 1 - found
    ull_t c;// local count of combination
} sha_t;

void freeAll(void *arg, ...);

void uniq(char *arr, size_t *n);

int init_alpha(char **argv, size_t *L, ull_t *N, sha_t *sha);

void InitString(char *str, size_t len, char c);

size_t comb(char *str, ull_t combination, size_t len, sha_t *sha);

ull_t get_block_size(ull_t n, int rank, int nprocs);

ull_t get_start_block(ull_t n, ull_t rank, int nprocs);

void backtrack(char *s, sha_t *sha, size_t len, int cur, long long int comb_per_proc);
#endif //BRUTESHA_1_SHA_1_H
