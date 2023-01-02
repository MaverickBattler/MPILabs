#include <iostream>
#include <stdio.h>
#include <mpi.h>
#include <algorithm>

using namespace std;

int main(int argc, char** argv)
{
	double timeStart, timeFinish;
	int procNum, procRank;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	int arrSize = procNum * 3;
	int* arr = new int[arrSize];
	for (int i = 0; i < arrSize; i++) {
		arr[i] = i + arrSize * procRank;
	}
	int* recvArr = new int[arrSize];

	timeStart = MPI_Wtime();

	MPI_Alltoall(arr, 3, MPI_INT, recvArr, 3, MPI_INT, MPI_COMM_WORLD);

	timeFinish = MPI_Wtime();

	//sort(recvArr, recvArr + arrSize);
	for (int i = 0; i < arrSize; i++) {
		cout << "Process " << procRank << " recieved integer = " << recvArr[i] << ", sender rank = " << i / 3 << endl;
	}

	//printf("%d.Start time: %f\n", procRank, timeStart);
	//printf("%d.End time: %f\n", procRank, timeFinish);

	MPI_Finalize();

	return 0;
}