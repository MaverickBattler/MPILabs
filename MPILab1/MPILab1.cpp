#include <iostream>
#include <stdio.h>
#include <mpi.h>

#define ARRAY_SIZE 100000

int main(int argc, char* argv[])
{
    double timeStart, timeFinish;
    int ProcNum, ProcRank, RecvRank;
    MPI_Status Status;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int* arr = new int[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; i++) arr[i] = i;

    if (ProcRank == 0)
    {
        timeStart = MPI_Wtime();
        for (int i = 1; i < ProcNum; i++)
        {
            MPI_Send(arr, ARRAY_SIZE, MPI_INT, i, 0, MPI_COMM_WORLD);
            std::cout << "-> Process 0 sent the message to process " << i << std::endl;
            MPI_Recv(&RecvRank, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            std::cout << "<- Process 0 recieved a confirmation message from process " << RecvRank << std::endl;
        }
        timeFinish = MPI_Wtime();
        printf("Message length: %d\n", ARRAY_SIZE);
        //printf("Start time: %f\n", timeStart);
        //printf("End time: %f\n", timeFinish);
        //printf("Total time: %f\n", timeFinish - timeStart);
    }
    else
    {
        MPI_Recv(arr, ARRAY_SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
        std::cout << "<- Process " << ProcRank << " recieved the message" << std::endl;
        MPI_Send(&ProcRank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        std::cout << "-> Process " << ProcRank << " sent the confirmation message" << std::endl;
    }

    MPI_Finalize();

    return 0;
}