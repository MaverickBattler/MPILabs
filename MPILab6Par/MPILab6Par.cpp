//#include <windows.h>
//#include <iostream>
//#include "mpi.h"
//int ProcCount = 0;
//int Rank = 0;
//int* A = 0, * B = 0, * C = 0, * A1 = 0, * B1 = 0, * C1 = 0, * T1 = 0;
//int ActualMatrSize = 0, MatrSize = 0, BlockSize = 0, GridSize = 0;
//int GridSize1 = 0, GridSize2 = 0;
//int GridCoords[2];
//double timeStart, timeFinish;
//MPI_Comm GridComm, ColComm, RowComm;
//MPI_Datatype MPI_BLOCK;
//
//using namespace std;
//
//void printMatrix(int** matrix, int matrixSize) {
//	for (int i = 0; i < matrixSize; i++) {
//		for (int j = 0; j < matrixSize; j++) {
//			cout << matrix[i][j] << " ";
//		}
//		cout << endl;
//	}
//}
//
//bool powerOfTwo(int n)
//{
//	if (n == 0) { return 0; }
//	while (n != 1)
//	{
//		n = n / 2;
//		if (n % 2 != 0 && n != 1) { return false; }
//	}
//	return true;
//}
//
//bool isPerfectSquare(int x)
//{
//	if (x >= 0) {
//		int sr = (int)sqrt(x);
//		return (sr * sr == x);
//	}
//	return false;
//}
//
//void CreateGridComm()
//{
//	int DimSize[2] = { GridSize1, GridSize2 };
//	int Periodic[2] = { 0,0 };
//	int SubDims[2];
//	MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 0, &GridComm);
//	MPI_Cart_coords(GridComm, Rank, 2, GridCoords);
//	SubDims[0] = 0;
//	SubDims[1] = 1;
//	MPI_Cart_sub(GridComm, SubDims, &RowComm);
//	SubDims[0] = 1;
//	SubDims[1] = 0;
//	MPI_Cart_sub(GridComm, SubDims, &ColComm);
//}
//void InitMatrices()
//{
//	if (Rank == 0)
//		ActualMatrSize = 10; // задание порядка матрицы
//	MPI_Bcast(&ActualMatrSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	
//	if (!powerOfTwo(ActualMatrSize)) {
//		int powerOfTwo = 1;
//		while (powerOfTwo < ActualMatrSize) {
//			powerOfTwo *= 2;
//		}
//		MatrSize = powerOfTwo;
//	}
//	else
//		MatrSize = ActualMatrSize;
//
//	BlockSize = MatrSize / GridSize;
//	A1 = new int[BlockSize * BlockSize];
//	B1 = new int[BlockSize * BlockSize];
//	C1 = new int[BlockSize * BlockSize];
//	T1 = new int[BlockSize * BlockSize];
//	for (int i = 0; i < BlockSize * BlockSize; i++)
//		C1[i] = 0;
//	if (Rank == 0)
//	{
//		A = new int[MatrSize * MatrSize];
//		B = new int[MatrSize * MatrSize];
//		C = new int[MatrSize * MatrSize];
//		for (int i = 0; i < MatrSize; i++)
//			for (int j = 0; j < MatrSize; j++)
//			{
//				if (i < ActualMatrSize && j < ActualMatrSize) {
//					A[i * MatrSize + j] = i + i + j;
//					B[i * MatrSize + j] = i + j + j;
//				}
//				else {
//					A[i * MatrSize + j] = 0;
//					B[i * MatrSize + j] = 0;
//				}
//				
//			}
//		cout << "Matrix A" << endl;
//		for (int i = 0; i < MatrSize; i++)
//		{
//			for (int j = 0; j < MatrSize; j++)
//				cout << A[i * MatrSize + j] << " ";
//			cout << endl;
//		}
//		cout << "Matrix B" << endl;
//		for (int i = 0; i < MatrSize; i++)
//		{
//			for (int j = 0; j < MatrSize; j++)
//				cout << B[i * MatrSize + j] << " ";
//			cout << endl;
//		}
//	}
//}
//void DataDistribution()
//{
//	MPI_Type_vector(BlockSize, BlockSize, MatrSize, MPI_INT, &MPI_BLOCK);
//	MPI_Type_commit(&MPI_BLOCK);
//	MPI_Request req1;
//	MPI_Request req2;
//
//	MPI_Irecv(A1, BlockSize * BlockSize, MPI_INT, 0, 0, GridComm, &req1);
//	MPI_Irecv(B1, BlockSize * BlockSize, MPI_INT, 0, 0, GridComm, &req2);
//
//	if (Rank == 0)
//		for (int r = 0; r < ProcCount; r++)
//		{
//			int c[2];
//			MPI_Cart_coords(GridComm, r, 2, c);
//			MPI_Send(A + c[0] * MatrSize * BlockSize + c[1] * BlockSize, 1,
//				MPI_BLOCK, r, 0, GridComm);
//			MPI_Send(B + c[0] * MatrSize * BlockSize + c[1] * BlockSize, 1,
//				MPI_BLOCK, r, 0, GridComm);
//		}
//	MPI_Status s;
//	MPI_Wait(&req1, &s);
//	MPI_Wait(&req2, &s);
//	cout << "Block A" << endl;
//	for (int i = 0; i < BlockSize; i++)
//	{
//		for (int j = 0; j < BlockSize; j++)
//			cout << Rank << ": " << A1[i * BlockSize + j] << " ";
//		cout << endl;
//	}
//	cout << "Block B" << endl;
//	for (int i = 0; i < BlockSize; i++)
//	{
//		for (int j = 0; j < BlockSize; j++)
//			cout << Rank << ": " << B1[i * BlockSize + j] << " ";
//		cout << endl;
//	}
//}
//void A1SendRow(int iter)
//{
//	int p = (GridCoords[0] + iter) % GridSize;
//	if (GridCoords[1] == p)
//		for (int i = 0; i < BlockSize * BlockSize; i++)
//			T1[i] = A1[i];
//	MPI_Bcast(T1, BlockSize * BlockSize, MPI_INT, p, RowComm);
//}
//void A1B1Mult()
//{
//	for (int i = 0; i < BlockSize; i++)
//		for (int j = 0; j < BlockSize; j++)
//		{
//			int t = 0;
//			for (int k = 0; k < BlockSize; k++)
//				t += T1[i * BlockSize + k] * B1[k * BlockSize + j];
//			C1[i * BlockSize + j] += t;
//		}
//}
//void B1SendCol()
//{
//	MPI_Status s;
//	int NextProc = GridCoords[0] + 1;
//	if (GridCoords[0] == GridSize - 1)
//		NextProc = 0;
//	int PrevProc = GridCoords[0] - 1;
//	if (GridCoords[0] == 0)
//		PrevProc = GridSize - 1;
//	MPI_Sendrecv_replace(B1, BlockSize * BlockSize, MPI_INT, PrevProc, 0,
//		NextProc, 0, ColComm, &s);
//}
//void ParallelCalc()
//{
//	for (int i = 0; i < GridSize; i++)
//	{
//		A1SendRow(i);
//		A1B1Mult();
//		B1SendCol();
//	}
//}
//void GatherResult()
//{
//	cout << "Block C" << endl;
//	for (int i = 0; i < BlockSize; i++)
//	{
//		for (int j = 0; j < BlockSize; j++)
//			cout << C1[i * BlockSize + j] << " ";
//		cout << endl;
//	}
//	MPI_Request req;
//	MPI_Isend(C1, BlockSize * BlockSize, MPI_INT, 0, 0, GridComm, &req);
//	if (Rank == 0)
//	{
//		MPI_Status s;
//		for (int r = 0; r < ProcCount; r++)
//		{
//			int c[2];
//			MPI_Cart_coords(GridComm, r, 2, c);
//			MPI_Recv(C + c[0] * MatrSize * BlockSize + c[1] * BlockSize, 1,
//				MPI_BLOCK, r, 0, GridComm, &s);
//			MPI_Wait(&req, &s);
//		}
//		cout << "Matrix C" << endl;
//		for (int i = 0; i < ActualMatrSize; i++)
//		{
//			for (int j = 0; j < ActualMatrSize; j++)
//				cout << C[i * MatrSize + j] << " ";
//			cout << endl;
//		}
//	}
//}
//
//void Termination() {
//
//}
//
//int main(int argc, char** argv) {
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &ProcCount);
//	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
//
//	if (powerOfTwo(ProcCount)) {
//
//		GridSize = (int)sqrt(ProcCount);
//		if (isPerfectSquare(ProcCount)) {
//			GridSize1 = GridSize2 = (int)sqrt(ProcCount);
//		}
//		else {
//			GridSize1 = (int)sqrt(ProcCount / 2) * 2;
//			GridSize2 = (int)sqrt(ProcCount / 2);
//		}
//
//		if (Rank == 0)
//			timeStart = MPI_Wtime();
//		CreateGridComm();
//		InitMatrices();
//		DataDistribution();
//		ParallelCalc();
//		GatherResult();
//		if (Rank == 0) {
//			timeFinish = MPI_Wtime();
//			cout << "Time spent: " << timeFinish - timeStart << endl;
//		}
//		// Termination(); 
//	}
//	else {
//		if (Rank == 0) {
//			cout << "Amount of processes is not a power of 2. Terminating the application" << endl;
//		}
//	}
//	MPI_Finalize();
//}

#include <windows.h>
#include <iostream>
#include "mpi.h"
int ProcCount = 0;
int Rank = 0;
int* A = 0, * B = 0, * C = 0, * A1 = 0, * B1 = 0, * C1 = 0, * T1 = 0;
int ActualMatrSize = 0, MatrSize = 0, BlockSize = 0, GridSize = 0;
int GridCoords[2];
double timeStart, timeFinish;
MPI_Comm GridComm, ColComm, RowComm;
MPI_Datatype MPI_BLOCK;

using namespace std;

void printMatrix(int** matrix, int matrixSize) {
	for (int i = 0; i < matrixSize; i++) {
		for (int j = 0; j < matrixSize; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

bool isDivisibleBy(int a, int b)
{
	if (a % b == 0)
		return true;
	else
		return false;
}

bool isPerfectSquare(int x)
{
	if (x >= 0) {
		int sr = (int)sqrt(x);
		return (sr * sr == x);
	}
	return false;
}

void CreateGridComm()
{
	int DimSize[2] = { GridSize, GridSize };
	int Periodic[2] = { 1,1 };
	int SubDims[2];
	MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 0, &GridComm);
	MPI_Cart_coords(GridComm, Rank, 2, GridCoords);
	SubDims[0] = 0;
	SubDims[1] = 1;
	MPI_Cart_sub(GridComm, SubDims, &RowComm);
	SubDims[0] = 1;
	SubDims[1] = 0;
	MPI_Cart_sub(GridComm, SubDims, &ColComm);
}
void InitMatrices()
{
	if (Rank == 0)
		ActualMatrSize = 200; // задание порядка матрицы
	MPI_Bcast(&ActualMatrSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MatrSize = ActualMatrSize;
	while (!isDivisibleBy(MatrSize, GridSize)) {
		MatrSize++;
	}
	BlockSize = MatrSize / GridSize;
	A1 = new int[BlockSize * BlockSize];
	B1 = new int[BlockSize * BlockSize];
	C1 = new int[BlockSize * BlockSize];
	T1 = new int[BlockSize * BlockSize];
	for (int i = 0; i < BlockSize * BlockSize; i++)
		C1[i] = 0;
	if (Rank == 0)
	{
		A = new int[MatrSize * MatrSize];
		B = new int[MatrSize * MatrSize];
		C = new int[MatrSize * MatrSize];
		for (int i = 0; i < MatrSize; i++)
			for (int j = 0; j < MatrSize; j++)
			{
				if (i < ActualMatrSize && j < ActualMatrSize) {
					A[i * MatrSize + j] = i + i + j;
					B[i * MatrSize + j] = i + j + j;
				}
				else {
					A[i * MatrSize + j] = 0;
					B[i * MatrSize + j] = 0;
				}

			}
		/*cout << "Matrix A" << endl;
		for (int i = 0; i < MatrSize; i++)
		{
			for (int j = 0; j < MatrSize; j++)
				cout << A[i * MatrSize + j] << " ";
			cout << endl;
		}
		cout << "Matrix B" << endl;
		for (int i = 0; i < MatrSize; i++)
		{
			for (int j = 0; j < MatrSize; j++)
				cout << B[i * MatrSize + j] << " ";
			cout << endl;
		}*/
	}
}
void DataDistribution()
{
	MPI_Type_vector(BlockSize, BlockSize, MatrSize, MPI_INT, &MPI_BLOCK);
	MPI_Type_commit(&MPI_BLOCK);
	MPI_Request req1;
	MPI_Request req2;

	MPI_Irecv(A1, BlockSize * BlockSize, MPI_INT, 0, 0, GridComm, &req1);
	MPI_Irecv(B1, BlockSize * BlockSize, MPI_INT, 0, 0, GridComm, &req2);

	if (Rank == 0)
		for (int r = 0; r < ProcCount; r++)
		{
			int c[2];
			MPI_Cart_coords(GridComm, r, 2, c);
			MPI_Send(A + c[0] * MatrSize * BlockSize + c[1] * BlockSize, 1,
				MPI_BLOCK, r, 0, GridComm);
			MPI_Send(B + c[0] * MatrSize * BlockSize + c[1] * BlockSize, 1,
				MPI_BLOCK, r, 0, GridComm);
		}
	MPI_Status s;
	MPI_Wait(&req1, &s);
	MPI_Wait(&req2, &s);
}
void A1SendRow(int iter)
{
	int p = (GridCoords[0] + iter) % GridSize;
	if (GridCoords[1] == p)
		for (int i = 0; i < BlockSize * BlockSize; i++)
			T1[i] = A1[i];
	MPI_Bcast(T1, BlockSize * BlockSize, MPI_INT, p, RowComm);
}
void A1B1Mult()
{
	for (int i = 0; i < BlockSize; i++)
		for (int j = 0; j < BlockSize; j++)
		{
			int t = 0;
			for (int k = 0; k < BlockSize; k++)
				t += T1[i * BlockSize + k] * B1[k * BlockSize + j];
			C1[i * BlockSize + j] += t;
		}
}
void B1SendCol()
{
	MPI_Status s;
	int NextProc = GridCoords[0] + 1;
	if (GridCoords[0] == GridSize - 1)
		NextProc = 0;
	int PrevProc = GridCoords[0] - 1;
	if (GridCoords[0] == 0)
		PrevProc = GridSize - 1;
	MPI_Sendrecv_replace(B1, BlockSize * BlockSize, MPI_INT, PrevProc, 0,
		NextProc, 0, ColComm, &s);
}
void ParallelCalc()
{
	for (int i = 0; i < GridSize; i++)
	{
		A1SendRow(i);
		A1B1Mult();
		B1SendCol();
	}
}
void GatherResult()
{
	MPI_Request req;
	MPI_Status s;

	MPI_Isend(C1, BlockSize * BlockSize, MPI_INT, 0, 0, GridComm, &req);
	if (Rank == 0)
	{
		for (int r = 0; r < ProcCount; r++)
		{
			int c[2];
			MPI_Cart_coords(GridComm, r, 2, c);
			MPI_Recv(C + c[0] * MatrSize * BlockSize + c[1] * BlockSize, 1,
				MPI_BLOCK, r, 0, GridComm, &s);
		}
		MPI_Wait(&req, &s);
		/*cout << "Matrix C" << endl;
		for (int i = 0; i < ActualMatrSize; i++)
		{
			for (int j = 0; j < ActualMatrSize; j++)
				cout << C[i * MatrSize + j] << " ";
			cout << endl;
		}*/
	} else
		MPI_Wait(&req, &s);
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

	if (isPerfectSquare(ProcCount)) {
		GridSize = (int)sqrt(ProcCount);
		InitMatrices();
		if (Rank == 0)
			timeStart = MPI_Wtime();
		CreateGridComm();
		DataDistribution();
		ParallelCalc();
		GatherResult();
		if (Rank == 0) {
			timeFinish = MPI_Wtime();
			cout << "Matrix size = " << ActualMatrSize << endl;
			cout << "Number of processes = " << ProcCount << endl;
			cout << "Time spent: " << timeFinish - timeStart << endl;
		}
	}
	else {
		if (Rank == 0) {
			cout << "Amount of processes is not a perfect square. Terminating the application." << endl;
		}
	}
	MPI_Finalize();
}