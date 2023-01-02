#include <mpi.h>
#include <stdio.h>
#include "stdlib.h"

#define ARRAY_SIZE 1000

int main(int argc, char** argv)
{
    double time;
    int size, rank;
    MPI_Status status;

    int* arr = new int[ARRAY_SIZE];

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank % 2 == 0) {
        if (rank != size - 1) {
            for (int i = 0; i < ARRAY_SIZE; i++) arr[i] = i;
            time = MPI_Wtime();
            MPI_Send(arr, ARRAY_SIZE, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            printf("-> Process %d sent message of %d elements to process %d \n", rank, ARRAY_SIZE, rank + 1);
            //printf("Sent time: %f", time);
        }
    }
    else {
        MPI_Recv(arr, ARRAY_SIZE, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        time = MPI_Wtime();
        printf("<- Process %d recieved message of %d elements from process %d \n", rank, ARRAY_SIZE, rank - 1);
        //printf("Delivered time: %f", time);
    }

    MPI_Finalize();
    return 0;
}
