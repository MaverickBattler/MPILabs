#include <iostream>
#include <stdio.h>
#include <mpi.h>

#define ARRAY_SIZE 100

void mergeArrays(int arr1[], int arr2[], int n1, int n2, int arr3[])
{
    int i = 0, j = 0, k = 0;
    while (i < n1) {
        arr3[k++] = arr1[i++];
    }
    while (j < n2) {
        arr3[k++] = arr2[j++];
    }
}

int main(int argc, char** argv)
{
    double timeStart, timeFinish;
    int procNum, procRank, recvRank;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    int arrSize = procNum - 1 + ARRAY_SIZE;
    int* arr = new int[arrSize];

    if (procRank == 0)
    {
        int* routePart = new int[procNum - 1];
        for (int i = 1; i < procNum; i++)
        {
            routePart[i - 1] = i;
        }
        int* informationPart = new int[ARRAY_SIZE];
        for (int i = 0; i < ARRAY_SIZE; i++)
        {
            informationPart[i] = i;
        }
        mergeArrays(routePart, informationPart, procNum - 1, ARRAY_SIZE, arr);

        int destProc = arr[0];
        timeStart = MPI_Wtime();
        MPI_Send(arr, arrSize, MPI_INT, destProc, 0, MPI_COMM_WORLD);
        std::cout << "-> Process 0 sent the message to process " << destProc << std::endl;
        MPI_Recv(&recvRank, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
        timeFinish = MPI_Wtime();
        printf("Message length: %d\n", ARRAY_SIZE);
        //printf("Start time: %f\n", timeStart);
        //printf("End time: %f\n", timeFinish);
        //printf("Total time: %f\n", timeFinish - timeStart);
        std::cout << "<- Process 0 recieved a confirmation message from process " << recvRank << std::endl;
    }
    else
    {
        MPI_Recv(arr, arrSize, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        //std::cout << "<- Process " << procRank << " recieved the message" << std::endl;
        arrSize--;
        int* newArr = new int[arrSize];
        for (int i = 0; i < arrSize; i++) 
        {
            newArr[i] = arr[i + 1];
        }
        int destProc = newArr[0];
        if (destProc == 0) //start of the info part 
        { //the process must report to the process 0
            MPI_Send(&procRank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            //std::cout << "-> Process " << procRank << " sent the confirmation message" << std::endl;
        }
        else 
        { //the process sends the data forward
            MPI_Send(newArr, arrSize, MPI_INT, destProc, 0, MPI_COMM_WORLD);
            //std::cout << "-> Process " << procRank << " sent the message to process " << destProc << std::endl;
        }
    }
    MPI_Finalize();
    return 0;
}