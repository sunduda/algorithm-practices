#include "pch.h"
#ifndef _MATRIX2D_CPP
#define _MATRIX2D_CPP

#include <string>
#include <algorithm>
#include "matrix2D.h"

using namespace myla;


template <typename T>
void Matrix2D<T>::Pivot(int k) {
	T max_val = INT_MIN;
	int max_index = 0;
	int pk = 0;
	T x;
	for (int i = k; i < this->n_rows_; i++) {
		x = abs(this->matrix[this->matrix_p_[i]][k]);
		if (x > max_val) {
			max_val = x;
			max_index = i
		}
	}

	pk = this->matrix_p_[k];
	this->matrix_p_[k] = this->matrix_p_[max_index];
	this->matrix_p_[max_index] = pk;
}

template <typename T>
bool Matrix2D<T>::LUDecomposition() {
	int i = 0;
	int s = std::min(this->n_rows_, this->n_cols_)
	long double temp;

	for (int c = 0; c < this->n_cols_; c++) this->matrix_lu_[0][c] = this->matrix[0][c];

	while (i < s) {
		if (i < this->n_rows_ - 1) {
			for (int r = i + 1; r < this->n_rows_; r++) {
				if (i > 0) {
					temp = matrix[r][i] - matrix_lu_[r][0] * matrix_lu_[0][i];
					for (int j = 1; j < i; j++) {
						temp -= matrix_lu_[r][j] * matrix_lu_[j][i];
					}
				}
				else {
					temp = matrix[r][0];
				}
				matrix_lu_[r][i] = temp / matrix_lu_[i][i];
			}
		}

		for (int c = i + 1; c < this->n_cols_; c++) {
			temp = matrix[i + 1][c] - matrix_lu_[i+1][0] * matrix_lu_[0][c];
			for (int j = 1; j < i + 1; j++) {
				temp -= matrix_lu_[i+1][j] * matrix_lu_[j][i+1];
			}
			matrix_lu_[i+1][c] = temp;
		}

		++i;
	}
	return true;
}

template <typename T>
static Matrix2D<T> Matrix2D<T>::LUComposition(std::vector< std::vector<T> > matrix_lu_) {
	if (this->n_rows_ != this->n_cols_) return NULL;
	Matrix2D<T> ans(this->n_rows_);
	for (int i = 0; i < this->n_rows_; i++) {
		for (int j = 0; j < this->n_cols_; j++) {
			ans.matrix[i][j] = 0;
			for (int k = 0; k <= std::min(i, j); k++) {
				ans.matrix[i][j] += this->matrix_l_[(1 + i)*i / 2 + k] * this->matrix_u_[(1 + j)*j / 2 + k];
			}
		}
	}
	return ans;
}

template <typename T>
void Matrix2D<T>::DisplayLU(std::string mode) {
	if (mode.compare("combined") == 0) {
		for (int i = 0; i < n_rows_; i++) {
			for (int k = 0; k < i; k++) {
				std::cout << matrix_l_[(1 + i)*i / 2 + k] << "\t\t";
			}
			for (int k = i; k < n_cols_; k++) {
				std::cout << std::left << std::setw(25) << matrix_u_[(1 + k)*k / 2 + i];
				if ((k == n_cols_ - 1) && (i != n_rows_)) std::cout << std::endl << std::endl;
				else if (k == n_cols_ - 1) std::cout << std::endl;
			}
		}

	}
	else if (mode.compare("separated") == 0) {
		for (int i = 0; i < n_rows_; i++) {
			for (int k = i; k < n_cols_; k++) {
				std::cout << matrix_u_[(1 + k)*k / 2 + i];
				if (k != n_cols_ - 1) std::cout << "\t";
				else std::cout << std::endl;
			}
		}

		std::cout << std::endl;

		for (int i = 0; i < n_rows_; i++) {
			for (int k = 0; k <= i; k++) {
				std::cout << matrix_l_[(1 + i)*i / 2 + k];
				if (k != i) std::cout << "\t";
				else std::cout << std::endl;
			}
		}
	}
}

template <typename T>
Matrix2D<T> Matrix2D<T>::operator*(const Matrix2D &m) {
	if (this->n_cols_ != m.n_rows_) return NULL;
	Matrix2D<T> ans(this->n_rows_, m.n_cols_);
	for (int i = 0; i < this->n_rows_; i++) {
		for (int j = 0; j < m.n_cols_; j++) {
			ans.matrix[i][j] = 0;
			for (int k = 0; k < this->n_cols_; k++) {
				ans.matrix[i][j] += this->matrix[i][k] * m.matrix[k][j];
			}
		}
	}
	return ans;
}

#endif