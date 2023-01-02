#include <iostream>
#include <stdio.h>
#include <mpi.h>

#define DIMENSIONS 3

using namespace std;

int main(int argc, char** argv)
{
    MPI_Comm matrix_Kdiv4x2x2,
        row_2;
    int procNum, originalRank, rankAfter;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &originalRank);

    int* dims = new int[DIMENSIONS];
    dims[0] = 2;
    dims[1] = 2;
    dims[2] = procNum / 4;
    int* periods = new int[DIMENSIONS];
    periods[0] = 0;
    periods[1] = 0; //1
    periods[2] = 0;
    
    MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, dims, periods, false, &matrix_Kdiv4x2x2);

    //finding coords of process in new comm (just to output this info)
    int* coords = new int[DIMENSIONS];
    MPI_Cart_coords(matrix_Kdiv4x2x2, originalRank, DIMENSIONS, coords);

    int* remainDims = new int[DIMENSIONS];
    remainDims[0] = 0;
    remainDims[1] = 1;
    remainDims[2] = 0;

    MPI_Cart_sub(matrix_Kdiv4x2x2, remainDims, &row_2);

    int* rootCoords = new int[1];
    rootCoords[0] = 0;
    int rootRank;

    //Finding a rank of the main process of every row_2 communicator
    MPI_Cart_rank(row_2, rootCoords, &rootRank);

    double number = originalRank * 1.8 + 0.3;
    double result;

    MPI_Reduce(&number, &result, 1, MPI_DOUBLE, MPI_PROD, rootRank, row_2);

    MPI_Comm_rank(row_2, &rankAfter);

    if (rankAfter == rootRank)
        cout << "Root process " << originalRank << " , coords: [" <<
        coords[0] << "," << coords[1] << "," << coords[2] << "], number = " << 
        number << ", result calculated: " << result << endl;
    else
        cout << "Process " << originalRank << ", coords: [" <<
        coords[0] << "," << coords[1] << "," << coords[2] << "], number = " << number << endl;

    MPI_Finalize();
}

//int main(int argc, char** argv)
//{
//    MPI_Comm matrix_Kdiv4x2x2, row_2;
//    double timeStart, timeFinish;
//    int procNum, originalRank;
//
//    MPI_Init(&argc, &argv);
//
//    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
//    MPI_Comm_rank(MPI_COMM_WORLD, &originalRank);
//
//    int* dims = new int[DIMENSIONS];
//    dims[0] = 2;
//    dims[1] = 2;
//    dims[2] = procNum / 4;
//
//    int* periods = new int[DIMENSIONS];
//    periods[0] = 0;
//    periods[1] = 0; //1
//    periods[2] = 0;
//
//    int* remainDims = new int[DIMENSIONS];
//    remainDims[0] = 0;
//    remainDims[1] = 1;
//    remainDims[2] = 0;
//
//    int* rootCoords = new int[1];
//    rootCoords[0] = 0;
//    int rootRank;
//
//    double number = originalRank * 1.8 + 0.3;
//    double result;
//
//    timeStart = MPI_Wtime();
//
//    MPI_Cart_create(MPI_COMM_WORLD, DIMENSIONS, dims, periods, false, &matrix_Kdiv4x2x2);
//
//    MPI_Cart_sub(matrix_Kdiv4x2x2, remainDims, &row_2);
//    
//    MPI_Cart_rank(row_2, rootCoords, &rootRank);
//
//    MPI_Reduce(&number, &result, 1, MPI_DOUBLE, MPI_PROD, rootRank, row_2);
//
//    timeFinish = MPI_Wtime();
//    printf("%d.Start time: %f\n", originalRank, timeStart);
//    printf("%d.End time: %f\n", originalRank, timeFinish);
//
//    MPI_Finalize();
//}
