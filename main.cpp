// main.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include "include/xinda_math/tensor.h"
#include "include/xinda_math/tensor.cpp"

constexpr int ROWS = 5;
constexpr int COLS = 5;

int main() {
//	Matrix2D<double> m_a(ROWS, COLS);
//	for (int i = 0; i < ROWS; i++) {
//		for (int j = 0; j < COLS; j++) {
//			m_a.AssignElement(i, j, 25 - 5 * i - j);
//		}
//	}
//
//	std::cout.precision(2);
//	std::cout << m_a << std::endl;
//	m_a.LUDecomposition();
//	m_a.DisplayLU(2);
    std::size_t sizes[2] = {2, 2};
    xinda_math::tensor<2, int> mtest(sizes, 2);
	return 0;
}