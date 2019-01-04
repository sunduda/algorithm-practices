// algorithm-practices.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <iomanip>
#include "matrix2D.cpp"

#define N 5

int main() {
	Matrix2D<double> m_a(N);
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			m_a.matrix[i - 1][j - 1] = rand() % 100 + 1;
			std::cout << m_a.matrix[i - 1][j - 1] << ' ';
		}
	}
	std::cout << std::endl << std::endl;
	m_a.LUDecomposition();

	m_a.DisplayLU("combined");
	m_a.DisplayLU("separated");

	Matrix2D<double> mtest((int)m_a.matrix.size());
	mtest = m_a.LUComposition();

	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			std::cout << mtest.matrix[i - 1][j - 1] << ' ';
		}
	}

	return 0;
}