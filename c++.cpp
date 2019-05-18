
#include "pch.h"
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
	ApproximationObj approxObj(1, -1, 1);
	size_t xIterations, tIterations;

	std::cout << "Enter x iterations: ";
	std::cin >> xIterations;
	std::cout << "Enter t iterations: ";
	std::cin >> tIterations;

	double tauStep = 1.0 / (tIterations);
	double hStep = 1.0 / (xIterations);

	double** accurateMatrix = new double*[tIterations];
	double** approximatedMatrix = new double*[tIterations];
	for (size_t i = 0; i < tIterations; ++i) {
		accurateMatrix[i] = new double[xIterations];
		approximatedMatrix[i] = new double[xIterations];
	}

	double currentX = 0, currentT = 0;

	std::cout << std::endl << "Accurate Matrix" << std::endl;
	#pragma omp parallel for private(currentX, currentT)
	for (size_t i = 0; i < tIterations; ++i) {
		currentT = tauStep * i;
		for (size_t j = 0; j < xIterations; ++j) {
			currentX = hStep * j;
			accurateMatrix[i][j] = approxObj.getAccurate(currentX, currentT);
			std::cout << accurateMatrix[i][j] << " ";
		}
	}
	   
	currentX = 0;
	currentT = 0;

	for (size_t i = 0; i < tIterations; ++i) {
		approximatedMatrix[i][0] = approxObj.leftLimit(currentT);
		approximatedMatrix[i][xIterations - 1] = approxObj.rightLimit(currentT);
		currentT += tauStep;
	}

	currentX += hStep;
	for (size_t i = 1; i < xIterations - 1; ++i) {
		approximatedMatrix[0][i] = approxObj.bottomLimit(currentX);
		currentX += hStep;
	}

	std::cout << std::endl << "Approximated Matrix" << std::endl;
	for (size_t i = 1; i < tIterations; ++i) {
	#pragma omp paralel for
		for (size_t j = 1; j < xIterations - 1; ++j) {
			approximatedMatrix[i][j] = approxObj.getApproximation(approximatedMatrix[i - 1][j - 1], approximatedMatrix[i - 1][j + 1], approximatedMatrix[i - 1][j], tauStep, hStep);
			std::cout << approximatedMatrix[i][j] << " ";
			
		}
	}

Error error = getError(accurateMatrix, approximatedMatrix, tIterations, xIterations);
	std::cout << std::endl;
	std::cout << "Absolute error: " << " i [" << error.col << "] - j [" << error.row << "] - [" << error.absolute << "] " << std::endl;
	std::cout << "Relative error: " << " i [" << error.col << "] - j [" << error.row << "] - [" << error.relative << "] " << std::endl;

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