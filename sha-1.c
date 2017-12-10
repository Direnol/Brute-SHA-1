#include "sha-1.h"


void freeAll(void *arg, ...)
{
    va_list ap;
    va_start(ap, arg);
    while (arg != NULL) {
        free(arg);
        arg = va_arg(ap, void *);

    }
    va_end(ap);

}

void uniq(char *arr, size_t *n)
{
    for (size_t i = 0; i < (*n - 1); ++i) {
        if (arr[i] == arr[i + 1]) {
            int k = 0;
            while (arr[i] == arr[i + k + 1]) ++k;
            for (size_t j = i; j < *n - k + 1; ++j) {
                arr[j] = arr[j + k];
            }
            (*n) -= k;
        }
    }
}

int init_alpha(char **argv, size_t *L, ull_t *N, sha_t *sha)
{
    sha->alphabet = strdup(argv[3]);
    sha->n = strlen(sha->alphabet);
    qsort(sha->alphabet, sha->n, sizeof(char), (__compar_fn_t) strcmp);

    uniq(sha->alphabet, &sha->n);

    // atoi
    char *endptr;
    *L = (size_t) strtoll(argv[2], &endptr, 10);
    if (endptr == argv[2]) {
        perror("Strtoll:");
        return EXIT_FAILURE;
    }

    sha->alpha = calloc((size_t) (sha->alphabet[sha->n - 1] + 2), sizeof(char));

    if (!sha->alpha) {
        perror("Malloc alpha:");
        return EXIT_FAILURE;
    }
    for (size_t i = 0; i < sha->n; ++i) {
        sha->alpha[sha->alphabet[i]] = (int) i;
    }

    // get key
    sha->mem_sha = malloc(20 * sizeof(char));
    sha->mem_cur_sha = malloc(20 * sizeof(char));
    if (SHA1(argv[1], strlen(argv[1]), sha->mem_sha) == NULL) {
        free(sha->alpha);
        fprintf(stderr, "Hash is not created\n");
        return EXIT_FAILURE;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    *N = 0;
    for (int i = 1; i <= *L; ++i) *N += pow(sha->n, i);

    if (rank == 0) {
        printf("Alphabet [%s] | Max lenght [%lu] | N [%llu] | Input [%s]-->\n[",
               sha->alphabet, *L, *N, argv[1]);
        for (int i = 0; i < 20; ++i) printf("%x", sha->mem_sha[i]);
        printf("]\n");
    }

    sha->c = 0;
    sha->found = 0;
    return EXIT_SUCCESS;
}

void InitString(char *str, size_t len, char c)
{
    for (size_t i = 0; i < len; ++i) str[i] = c;
    str[len] = '\0';
}

size_t comb(char *str, ull_t combination, size_t len, sha_t *sha)
{
    ull_t sum = 0;
    ull_t cur_len = 0;
    for (ull_t i = 1; i <= len; ++i) {
        ull_t tmp = (ull_t) pow(sha->n, i);
        if (sum + tmp >= combination) {
            cur_len = i;
            break;
        }
        sum += tmp;
    }
    combination -= (sum);
    ull_t ind;
    ull_t p = (ull_t) pow(sha->n, cur_len - 1);
    for (ull_t i = 0; i < cur_len; ++i) {
        ind = (ull_t) floor(combination / p);
        combination -= (ind * p);
        str[i] = sha->alphabet[ind];
        p /= sha->n;
    }
    str[cur_len] = '\0';
    return cur_len;
}

ull_t get_block_size(ull_t n, int rank, int nprocs)
{
    ull_t s = n / nprocs;
    if (n % nprocs > rank)
        s++;
    return s;
}

ull_t get_start_block(ull_t n, ull_t rank, int nprocs)
{
    ull_t rem = n % nprocs;
    return n / nprocs * rank + ((rank >= rem) ? rem : rank);
}

void backtrack(char *s, sha_t *sha, size_t len, int cur, long long int comb_per_proc)
{
    if (cur == len || sha->found || sha->c >= comb_per_proc) return;

    int j;
    char check;
    if (cur < len - 1) {
        backtrack(s, sha, len, cur + 1, comb_per_proc);
        j = sha->alpha[s[cur]] + 1;
        check = 0;
    } else {
        j = sha->alpha[s[cur]];
        check = 1;
    }

    for (int k = j; k < sha->n && !sha->found && sha->c < comb_per_proc; ++k) {
        s[cur] = sha->alphabet[k];
        if (check) {
            ++sha->c;
//            printf("%s %llu\n", s, sha->c);
            if (SHA1(s, strlen(s), sha->mem_cur_sha) == NULL) {
                fprintf(stderr, "SHA in backtrack\n");
                exit(1);
            }
            if (strncmp(sha->mem_cur_sha, sha->mem_sha, 20) == 0) {
                sha->found = 1;
                printf("[%s] ", s);
                return;
            }
        }
        if (cur != len - 1) backtrack(s, sha, len, cur + 1, comb_per_proc);
    }
    s[cur] = sha->alphabet[0];
}
