#include <iostream>
#include <stdio.h>
#include <mpi.h>

#define ARRAY_SIZE 3

#define ROOT 0

using namespace std;

int main(int argc, char** argv)
{
    MPI_Comm comm;
    double timeStart, timeFinish;
    int rank, color, procNum, *arrSend, *arrRecv;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Setting the color
    if (rank % 3 == 0) {
        color = 0;
    }
    else {
        color = MPI_UNDEFINED;
    }
    timeStart = MPI_Wtime();
    // Split - creating a new communicator
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    
    if (rank % 3 == 0) {
        arrSend = new int[ARRAY_SIZE];
        for (int i = 0; i < ARRAY_SIZE; i++) {
            arrSend[i] = i + ARRAY_SIZE * (rank / 3);
        }
        cout << "Process " << rank << " sends integers: ";
        for (int i = 0; i < ARRAY_SIZE; i++) {
            cout << arrSend[i] << " ";
        }
        int recvArrSize = (1 + (procNum - 1) / 3) * ARRAY_SIZE;
        arrRecv = new int[recvArrSize];
        // Sending to a root process
        MPI_Gather(arrSend, ARRAY_SIZE, MPI_INT, arrRecv, ARRAY_SIZE, MPI_INT, ROOT, comm);
        timeFinish = MPI_Wtime();
        //printf("%d.Start time: %f\n", rank, timeStart);
        //printf("%d.End time: %f\n", rank, timeFinish);
        if (rank == ROOT) {
            cout << "and process " << ROOT << " received integers: ";
            for (int i = 0; i < recvArrSize; i++) {
                cout << arrRecv[i] << " ";
            }
            cout << endl;
        }
    }
    MPI_Finalize();
}