#include "pch.h"
#include "mpi.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>

struct Error
{
	double absolute, relative;
	size_t row, col;

	Error() : absolute(0), relative(0) {}
};

Error getError(double** accurateMatrix, double** approximatedMatrix, size_t rows, size_t columns);

int main(int argc, char *argv[])
{
	ApproximationObj approxObj(1, -0.75, 1);
	size_t xIterations, tIterations;

	std::cout << "Enter x iterations: ";
	std::cin >> xIterations;
	std::cout << "Enter t iterations: ";
	std::cin >> tIterations;

	double tauStep = 1.0 / (tIterations - 1);
	double hStep = 1.0 / (xIterations - 1);

	MPI_Init(&argc, &argv);

	int myid, numprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	double** accurateMatrix = new double*[tIterations];
	double** approximateMatrix = new double*[tIterations];

	for (size_t i = 0; i < tIterations; ++i) {
		accurateMatrix[i] = new double[xIterations];
		approximateMatrix[i] = new double[xIterations];
	}

	double currentX = 0, currentT = 0;

	std::ofstream AccurateFile("Accurate Matrix.txt");
	if (AccurateFile.is_open())
	{
		std::cout << std::endl << "Accurate Matrix is in file: Accurate Matrix.txt";
//#pragma omp parallel for private(currentX, currentT)

		for (size_t i = 0; i < tIterations; ++i) {
			currentT = tauStep * i;
			MPI_Bcast(accurateMatrix[i - 1], xIterations, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
			for (size_t j = 0; j < xIterations; ++j) {
				currentX = hStep * j;
				accurateMatrix[i][j] = approxObj.getAccurate(currentX, currentT);
				//std::cout << accurateMatrix[i][j] << " ";
				AccurateFile << "{" << currentX << ", " << currentT << ", " << accurateMatrix[i][j] << "}," << "\n";
				//AccurateFile << accurateMatrix[i][j] << " ";
			}
		}
		AccurateFile.close();
	}
	else std::cout << "Unable to open file";

	currentX = 0;
	currentT = 0;

	for (size_t i = 0; i < tIterations; ++i) {
		approximateMatrix[i][0] = approxObj.leftLimit(currentT);
		approximateMatrix[i][xIterations - 1] = approxObj.rightLimit(currentT);
		currentT += tauStep;
	}

	//currentX += hStep;
	for (size_t i = 1; i < xIterations - 1; ++i) {
		approximateMatrix[0][i] = approxObj.bottomLimit(currentX);
		currentX += hStep;
	}

	currentX = 0;
	currentT = 0;
	std::cout << std::endl << "Approximate Matrix is in file: Approximate Matrix.txt" << std::endl;
	for (size_t i = 1; i < tIterations; ++i) {
		currentX = 0;
//#pragma omp paralel for
		MPI_Bcast(approximateMatrix[i - 1], xIterations, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
		for (size_t j = 1; j < xIterations - 1; ++j) {
			approximateMatrix[i][j] = approxObj.getApproximation(approximateMatrix[i - 1][j - 1], approximateMatrix[i - 1][j + 1], approximateMatrix[i - 1][j], tauStep, hStep);
			//std::cout << approximateMatrix[i][j] << " ";
			currentX += hStep;
		}
		currentT += tauStep;
	}


	std::ofstream ApproximateFile("Approximate Matrix.txt");
	currentX = 0;
	currentT = 0;
	if (ApproximateFile.is_open())
	{
		std::cout << std::endl << "Approximate Matrix is in file: Approximate Matrix.txt" << std::endl;
		for (size_t i = 1; i < tIterations; ++i) {
			currentX = 0;

			for (size_t j = 0; j < xIterations; ++j) {

				//std::cout << approximateMatrix[i][j] << " ";
				ApproximateFile << "{" << currentX << ", " << currentT << ", " << approximateMatrix[i][j] << "}," << "\n";
				currentX += hStep;
			}
			currentT += tauStep;
		}
		ApproximateFile.close();
	}
	else std::cout << "Unable to open file";

	Error error = getError(accurateMatrix, approximateMatrix, tIterations, xIterations);
	std::cout << std::endl;
	std::cout << "Absolute error: " << " i[" << error.col << "] - j[" << error.row << "] - [" << error.absolute << "] " << std::endl;
	std::cout << "Relative error: " << " i[" << error.col << "] - j[" << error.row << "] - [" << error.relative << "] " << std::endl;

	MPI_Finalize();
	system("pause");
}

Error getError(double** accurateMatrix, double** approximatedMatrix, size_t rows, size_t columns)
{
	Error error;
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			if (fabs(approximatedMatrix[i][j] - accurateMatrix[i][j]) > error.absolute) {
				error.absolute = fabs(approximatedMatrix[i][j] - accurateMatrix[i][j]);
				error.row = i;
				error.col = j;
			}
		}
	}
	error.relative = error.absolute / accurateMatrix[error.row][error.col];
	return error;
}