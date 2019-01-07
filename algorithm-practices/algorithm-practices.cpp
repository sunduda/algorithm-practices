// algorithm-practices.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <iomanip>
#include "matrix2D.cpp"

#define ROWS 5
#define COLS 3

int main() {
	Matrix2D<double> m_a(ROWS, COLS);
	std::cout.setf(std::ios::right);
	for (int i = 1; i <= ROWS; i++) {
		for (int j = 1; j <= COLS; j++) {
			std::cout.width(4);
			m_a.matrix[i - 1][j - 1] = rand() % 100 + 1;
			std::cout << m_a.matrix[i - 1][j - 1];
			if (j == COLS) std::cout << std::endl;
			else std::cout << '\t';
		}
	}
	std::cout << std::endl;
	/*for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			m_a.matrix[i][j] = 25 - 5 * i - j;
			std::cout << m_a.matrix[i][j];
			if (j < N - 1) std::cout << ' ';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;*/
	m_a.LUDecomposition();

	m_a.DisplayLU(10);

	std::cout << std::endl << m_a.Determinant();

	/*Matrix2D<double> mtest((int)m_a.matrix.size());
	mtest = m_a.LUComposition();

	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			std::cout << mtest.matrix[i - 1][j - 1] << ' ';
		}
	}*/

	return 0;
}