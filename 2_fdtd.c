#include <stdio.h>
#include <stdlib.h>

#define tmax 100
#define nx 1024
#define ny 1024

double ex[nx][ny+1];
double ey[nx+1][ny];
double hz[nx][ny];

void init_array()
{
    int i, j;

    for (i=0; i<nx+1; i++)  {
        for (j=0; j<ny; j++)  {
            ey[i][j] = 0;
        }
    }

    for (i=0; i<nx; i++)  {
        for (j=0; j<ny+1; j++)  {
            ex[i][j] = 0;
        }
    }

    for (j=0; j<ny; j++)  {
        ey[0][j] = ((double)j)/ny;
    }

    for (i=0; i<nx; i++)    {
        for (j=0; j<ny; j++)  {
            hz[i][j] = 0;
        }
    }
}

void print_array()
{
    int i, j;

    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            fprintf(stderr, "%lf ", hz[i][j]);
        }
        fprintf(stderr, "\n");
    }
}

int main()
{
	int t, i, j;

	init_array();

#pragma scop
    for(t = 0; t < tmax; t++) {
        for (j = 0; j < ny; j++)
            ey[0][j] = t;
        for (i = 1; i < nx; i++)
            for (j = 0; j < ny; j++)
                ey[i][j] = ey[i][j] - 0.5 * (hz[i][j] - hz[i-1][j]);
        for (i = 0; i < nx; i++)
            for (j = 1; j < ny; j++)
                ex[i][j] = ex[i][j] - 0.5 * (hz[i][j] - hz[i][j-1]);
        for (i = 0; i < nx; i++)
            for (j = 0; j < ny; j++)
                hz[i][j] = hz[i][j] - 0.7 * (ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
    }
#pragma endscop

    print_array();

    return 0;
}
