#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define N 4

double A[N][N];
double A_flat[N*N];

void init_array() {
    int i, j;
    srand(123);
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            double ra = 1.0 + ((double)N * rand() / (RAND_MAX + 1.0));
            A[i][j] = ra;
        }
    }
    for (i = 0; i < N; i++) {
        A[i][i] = 0;
    }
}

void print_array() {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            fprintf(stderr, "%lf ", A[i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void flat_array() {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A_flat[i*N+j] = A[i][j];
        }
    }
}

void deflat_array() {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[i][j] = A_flat[i*N+j];
        }
    }

}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    init_array();

    flat_array();

    int rows_per_process = N / size;
    double *local_A = malloc(rows_per_process * N * sizeof(double));
    MPI_Scatter(A_flat, rows_per_process * N, MPI_DOUBLE, local_A, rows_per_process * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double *k_row = malloc(N * sizeof(double));
    for (int k = 0; k < N; k++) {
        int owner = k / rows_per_process;
        if (rank == owner) {
            int k_local = k % rows_per_process;
            for (int j = 0; j < N; j++) {
                k_row[j] = local_A[k_local * N + j];
            }
        }
        MPI_Bcast(k_row, N, MPI_DOUBLE, owner, MPI_COMM_WORLD);

        for (int i = 0; i < rows_per_process; i++) {
            for (int j = 0; j < N; j++) {
                local_A[i * N + j] = min(local_A[i * N + j], local_A[i * N + k] + k_row[j]);
            }
        }
    }

    MPI_Gather(local_A, rows_per_process * N, MPI_DOUBLE, A_flat, rows_per_process * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    deflat_array();

    if (rank == 0) {
        print_array();
    }

    free(local_A);
    free(k_row);
    MPI_Finalize();
    return 0;
}
