#include <iostream>
#include <chrono>

using namespace std;

int** MultiplyMatrices(int** A, int** B, int matrixSize) {

    int** C = new int* [matrixSize];
    for (int i = 0; i < matrixSize; i++) {
        C[i] = new int[matrixSize];
    }
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            C[i][j] = 0;
        }
    }

    for (int row = 0; row < matrixSize; row++) {
        for (int col = 0; col < matrixSize; col++) {
            // Multiply the row of A by the column of B to get the row, column of product.
            for (int inner = 0; inner < matrixSize; inner++) {
                C[row][col] += A[row][inner] * B[inner][col];
            }
        }
    }
    return C;
}

void printMatrix(int** matrix, int matrixSize) {
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char** argv) {
    int m;
    cout << "Matrix multiplication, sequential alrorithm." << endl;
    cout << "Please, enter the size of the matrix: ";
    cin >> m;
    int** A = new int* [m];
    for (int i = 0; i < m; i++) {
        A[i] = new int[m];
    }
    int** B = new int* [m];
    for (int i = 0; i < m; i++) {
        B[i] = new int[m];
    }
    if (m < 6) { //Manual input
        cout << "Please, enter the matrix A line by line:" << endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                cin >> A[i][j];
            }
        }
        cout << "Please, enter the matrix B line by line:" << endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                cin >> B[i][j];
            }
        }
    }
    else { //The contents of matrices are pre-defined
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = i + j;
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                B[i][j] = i + j;
            }
        }
    }
    
    int sumTime = 0;
    for (int i = 0; i < 3; i++) {
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        int** C = MultiplyMatrices(A, B, m);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        //printMatrix(A, m);
        //cout << "  *  " << endl;
        //printMatrix(B, m);
        //cout << "  =  " << endl;
        //printMatrix(C, m);
        int timeElapsed = chrono::duration_cast<chrono::nanoseconds> (end - begin).count();
        sumTime += timeElapsed;
    }
    int avgTime = sumTime / 3;

    cout << "Matrix size = " << m << endl;
    cout << "Average time = " << avgTime << "ns" << endl;

}