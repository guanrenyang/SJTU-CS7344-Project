#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 4

double A[N][N][N]; // k, j ,i

int start_k_0, end_k_0, start_i_0, end_i_0, start_j_0, end_j_0;
int start_k_1, end_k_1, start_i_1, end_i_1, start_j_1, end_j_1;

void init_array()
{
    int i, j, k;
    for (k=0; k<N; k++) {
        for (j=0; j<N; j++) {
            for (i=0; i<N; i++) {
                // A[k][j][i] = (1+(i*j+k)%1024)/3.0;
                A[k][j][i] = k*N*N + j*N + i;
            }
        }
    }
}

void print_array()
{
    int i, j, k;

    for (k=0; k<N; k++) {
        for (j=0; j<N; j++) {
            for (i=0; i<N; i++) {
                fprintf(stderr, "%lf ", A[k][j][i]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
    }
}

inline int index_plane(int i, int j){
    return i * N/2 + j;
}

void print_subarray(double A[N][N][N], int starts[3], int subsizes[3]) {
    printf("Subarray data:\n");
    for (int i = starts[0]; i < starts[0] + subsizes[0]; i++) {
        for (int j = starts[1]; j < starts[1] + subsizes[1]; j++) {
            for (int k = starts[2]; k < starts[2] + subsizes[2]; k++) {
                printf("%.2f ", A[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

void store_subarray_in_1D(double A[N][N][N], int starts[3], int subsizes[3], double* subarray_1D) {
    int idx = 0;
    for (int i = starts[0]; i < starts[0] + subsizes[0]; i++) {
        for (int j = starts[1]; j < starts[1] + subsizes[1]; j++) {
            for (int k = starts[2]; k < starts[2] + subsizes[2]; k++) {
                subarray_1D[idx++] = A[i][j][k];
            }
        }
    }
}

// Axis-dependent functions
void set_index_for_array(int *start_k_0, int *end_k_0, int *start_i_0, int *end_i_0, int *start_j_0, int *end_j_0, int *start_k_1, int *end_k_1, int *start_i_1, int *end_i_1, int *start_j_1, int *end_j_1, int rank, char dependent_axis) {
    if (dependent_axis=='i'){
        if (rank == 0) {
            *start_k_0 = 0;
            *end_k_0 = N/2;
            *start_i_0 = 0;
            *end_i_0 = N/2;
            *start_j_0 = N/2;
            *end_j_0 = N;
            *start_k_1 = N/2;
            *end_k_1 = N;
            *start_i_1 = N/2;
            *end_i_1 = N;
            *start_j_1 = 0;
            *end_j_1 = N/2;
        } else if (rank == 1){
            *start_k_0 = 0;
            *end_k_0 = N/2;
            *start_i_0 = 0;
            *end_i_0 = N/2;
            *start_j_0 = 0;
            *end_j_0 = N/2;
            *start_k_1 = N/2;
            *end_k_1 = N;
            *start_i_1 = N/2;
            *end_i_1 = N;
            *start_j_1 = N/2;
            *end_j_1 = N;
        } else if (rank == 2){
            *start_k_0 = N/2;
            *end_k_0 = N;
            *start_i_0 = 0;
            *end_i_0 = N/2;
            *start_j_0 = N/2;
            *end_j_0 = N;
            *start_k_1 = 0;
            *end_k_1 = N/2;
            *start_i_1 = N/2;
            *end_i_1 = N;
            *start_j_1 = 0;
            *end_j_1 = N/2;
        } else if (rank == 3){
            *start_k_0 = N/2;
            *end_k_0 = N;
            *start_i_0 = 0;
            *end_i_0 = N/2;
            *start_j_0 = 0;
            *end_j_0 = N/2;
            *start_k_1 = 0;
            *end_k_1 = N/2;
            *start_i_1 = N/2;
            *end_i_1 = N;
            *start_j_1 = N/2;
            *end_j_1 = N;
        }
    } else if (dependent_axis=='j') {
        if (rank == 0) {
            *start_k_0 = N/2;
            *end_k_0 = N;
            *start_i_0 = N/2;
            *end_i_0 = N;
            *start_j_0 = 0;
            *end_j_0 = N/2;
            *start_k_1 = 0;
            *end_k_1 = N/2;
            *start_i_1 = 0;
            *end_i_1 = N/2;
            *start_j_1 = N/2;
            *end_j_1 = N;
        } else if (rank == 1){
            *start_k_0 = 0;
            *end_k_0 = N/2;
            *start_i_0 = 0;
            *end_i_0 = N/2;
            *start_j_0 = 0;
            *end_j_0 = N/2;
            *start_k_1 = N/2;
            *end_k_1 = N;
            *start_i_1 = N/2;
            *end_i_1 = N;
            *start_j_1 = N/2;
            *end_j_1 = N;
        } else if (rank == 2){
            *start_k_0 = 0;
            *end_k_0 = N/2;
            *start_i_0 = N/2;
            *end_i_0 = N;
            *start_j_0 = 0;
            *end_j_0 = N/2;
            *start_k_1 = N/2;
            *end_k_1 = N;
            *start_i_1 = 0;
            *end_i_1 = N/2;
            *start_j_1 = N/2;
            *end_j_1 = N;
        } else if (rank == 3){
            *start_k_0 = N/2;
            *end_k_0 = N;
            *start_i_0 = 0;
            *end_i_0 = N/2;
            *start_j_0 = 0;
            *end_j_0 = N/2;
            *start_k_1 = 0;
            *end_k_1 = N/2;
            *start_i_1 = N/2;
            *end_i_1 = N;
            *start_j_1 = N/2;
            *end_j_1 = N;
        }
    } else if (dependent_axis=='k'){
        if (rank == 0) {
            *start_k_0 = 0;
            *end_k_0 = N/2;
            *start_i_0 = 0;
            *end_i_0 = N/2;
            *start_j_0 = N/2;
            *end_j_0 = N;
            *start_k_1 = N/2;
            *end_k_1 = N;
            *start_i_1 = N/2;
            *end_i_1 = N;
            *start_j_1 = 0;
            *end_j_1 = N/2;
        } else if (rank == 1){
            *start_k_0 = 0;
            *end_k_0 = N/2;
            *start_i_0 = 0;
            *end_i_0 = N/2;
            *start_j_0 = 0;
            *end_j_0 = N/2;
            *start_k_1 = N/2;
            *end_k_1 = N;
            *start_i_1 = N/2;
            *end_i_1 = N;
            *start_j_1 = N/2;
            *end_j_1 = N;
        } else if (rank == 2){
            *start_k_0 = 0;
            *end_k_0 = N/2;
            *start_i_0 = N/2;
            *end_i_0 = N;
            *start_j_0 = 0;
            *end_j_0 = N/2;
            *start_k_1 = N/2;
            *end_k_1 = N;
            *start_i_1 = 0;
            *end_i_1 = N/2;
            *start_j_1 = N/2;
            *end_j_1 = N;
        } else if (rank == 3){
            *start_k_0 = 0;
            *end_k_0 = N/2;
            *start_i_0 = N/2;
            *end_i_0 = N;
            *start_j_0 = N/2;
            *end_j_0 = N;
            *start_k_1 = N/2;
            *end_k_1 = N;
            *start_i_1 = 0;
            *end_i_1 = N/2;
            *start_j_1 = 0;
            *end_j_1 = N/2;
        }
    }
}

void set_partner_rank(int *partner_rank, int rank, char dependent_axis){
    if (dependent_axis=='i'){
        if (rank==0){
            *partner_rank = 3;
        } else if (rank==1){
            *partner_rank = 2;
        } else if (rank==2){
            *partner_rank = 1;
        } else if (rank==3){
            *partner_rank = 0;
        }
    } else if (dependent_axis=='j'){
        if (rank==0){
            *partner_rank = 1;
        } else if (rank==1){
            *partner_rank = 0;
        } else if (rank==2){
            *partner_rank = 3;
        } else if (rank==3){
            *partner_rank = 2;
        }
    } else if (dependent_axis=='k'){
        if (rank==0){
            *partner_rank = 2;
        } else if (rank==1){
            *partner_rank = 3;
        } else if (rank==2){
            *partner_rank = 0;
        } else if (rank==3){
            *partner_rank = 1;
        }
    }
}

void set_sizes_for_transmission(int *starts, int *subsizes, char dependent_axis, int rank){
    if (dependent_axis=='i'){
        starts[2] = N/2-1; // i = 511
        starts[1] = rank%2==0 ? N/2 : 0; // j = 512 for rank 0,2 and j = 0 for rank 1,3
        starts[0] = rank<2 ? 0 : N/2; // k = 0 for rank 0,1 and k = 512 for rank 2,3
        
        subsizes[2] = 1;
        subsizes[1] = N/2;
        subsizes[0] = N/2;
    } else if (dependent_axis=='j'){
        starts[0] = (rank==0||rank==3) ? N/2 : 0;
        starts[1] = N/2-1;
        starts[2] = rank%2==0 ? N/2 : 0;
        
        subsizes[0] = N/2;
        subsizes[1] = 1;
        subsizes[2] = N/2;
    } else if (dependent_axis=='k'){
        starts[0] = N/2-1;
        starts[1] = (rank==0 || rank==3) ? N/2 : 0;
        starts[2] = (rank<2) ? 0 : N/2;
        
        subsizes[0] = 1;
        subsizes[1] = N/2;
        subsizes[2] = N/2;
    }
}

void compute (int start_k, int end_k, int start_i, int end_i, int start_j, int end_j, double *boundary_buffer, char dependent_axis) {
    int i, j, k;
    if (dependent_axis=='i'){
        for (k = start_k; k < end_k; k++) {
            for (j = start_j; j < end_j; j++) {
                if (boundary_buffer)
                    A[k][j][start_i] = A[k][j][start_i] * 0.4 - boundary_buffer[index_plane(k-start_k, j-start_j)] * 0.6;
                for (i = start_i + 1; i < end_i; i++) {
                    A[k][j][i] = A[k][j][i] * 0.4 - A[k][j][i-1] * 0.6;
                }
            }
        }       
    } else if (dependent_axis=='j') {
        for (k = start_k; k < end_k; k++) {
            for (i = start_i; i < end_i; i++) {
                if (boundary_buffer)
                    A[k][start_j][i] = A[k][start_j][i] * 0.5 - boundary_buffer[index_plane(k-start_k, i-start_i)] * 0.5;
                for (j = start_j+1; j < end_j; j++) {
                    A[k][j][i] = A[k][j][i] * 0.5 - A[k][j-1][i] * 0.5;
                }
            }
        }   
    } else if (dependent_axis=='k') {
        for (j = start_j; j < end_j; j++) {
            for (i = start_i; i < end_i; i++) {
                if (boundary_buffer)
                    A[start_k][j][i] = A[start_k][j][i] * 0.6 - boundary_buffer[index_plane(j-start_j, i-start_i)] * 0.4;
                for (k = start_k+1; k < end_k; k++) {
                    A[k][j][i] = A[k][j][i] * 0.6 - A[k-1][j][i] * 0.4;
                }
            }
        }   
    }
}

void compute_axis_dependent_loop(char dependent_axis, int rank) {
    int i, j, k;
    int partner_rank;
    int starts[3], subsizes[3], bigsizes[3] = {N, N, N};

    double *send_buffer = malloc((N/2) * (N/2) * sizeof(double));
    double *recv_buffer = malloc((N/2) * (N/2) * sizeof(double));

    set_index_for_array(&start_k_0, &end_k_0, &start_i_0, &end_i_0, &start_j_0, &end_j_0, &start_k_1, &end_k_1, &start_i_1, &end_i_1, &start_j_1, &end_j_1, rank, dependent_axis);

    compute(start_k_0, end_k_0, start_i_0, end_i_0, start_j_0, end_j_0, NULL, dependent_axis);

    /* Transfer */
    set_partner_rank(&partner_rank, rank, dependent_axis);         
    set_sizes_for_transmission(starts, subsizes, dependent_axis, rank);
    store_subarray_in_1D(A, starts, subsizes, send_buffer);

    MPI_Request send_req, recv_req;
    MPI_Status recv_status;
    MPI_Isend(send_buffer, (N/2) * (N/2), MPI_DOUBLE, partner_rank, 0, MPI_COMM_WORLD, &send_req);
    MPI_Irecv(recv_buffer, (N/2) * (N/2), MPI_DOUBLE, partner_rank, 0, MPI_COMM_WORLD, &recv_req);
    MPI_Wait(&recv_req, &recv_status);
    /* Done Transfer*/

    compute(start_k_1, end_k_1, start_i_1, end_i_1, start_j_1, end_j_1, recv_buffer, dependent_axis);

    free(send_buffer);
    free(recv_buffer);
}


int main(int argc, char **argv)
{
    init_array();
    
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // check is size == 4
    if (size != 4) {
        fprintf(stderr, "size must be 4\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    compute_axis_dependent_loop('i', rank);
    compute_axis_dependent_loop('j', rank);
    compute_axis_dependent_loop('k', rank);
      
    MPI_Finalize();
}
