#include <stdio.h>
#include <stdlib.h>

#define N 1024
#define T 100

double A[2][N][N];

void init_array()
{
    int i, j;
    const int BASE = 1024;

    srand(42);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[0][i][j] = 1.0 * (rand() % BASE);
        }
    }
}

void print_array()
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            fprintf(stderr, "%0.2lf ", A[T%2][i][j]);
        }
        fprintf(stderr, "\n");
    }
}


int main()
{
    init_array();

#pragma scop
    for (int t = 0; t < T; t++) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                A[(t+1)%2][i][j] = 0.125 * (A[t%2][i+1][j] - 2.0 * A[t%2][i][j] + A[t%2][i-1][j])
                                 + 0.125 * (A[t%2][i][j+1] - 2.0 * A[t%2][i][j] + A[t%2][i][j-1])
                                 + A[t%2][i][j];
            }
        }
    }
#pragma endscop

    print_array();

    return 0;
}
