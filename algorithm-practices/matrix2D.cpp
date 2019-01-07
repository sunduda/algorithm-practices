#include "pch.h"
#ifndef _MATRIX2D_CPP
#define _MATRIX2D_CPP

#include <string>
#include <algorithm>
#include "matrix2D.h"

#define PRECISION 10e-12

using namespace myla;

template <typename T>
void Matrix2D<T>::Pivot(int k) {
	std::vector< std::vector<double> >& mlu = this->matrix_lu_;
	std::vector<int>& pr = this->matrix_pr_;
	std::vector<int>& pc = this->matrix_pc_;
	T max_val = -128;
	int max_i = k;
	int max_j = k;
	int pk = k;
	for (int i = k; i < this->n_rows_; i++) {
		for (int j = k; j < this->n_cols_; j++) {
			if (abs(mlu[pr[i]][pc[j]]) > max_val) {
				max_val = abs(mlu[pr[i]][pc[j]]);
				max_i = i;
				max_j = j;
			}
		}
	}

	pk = pr[k];
	pr[k] = pr[max_i];
	pr[max_i] = pk;
	pk = pc[k];
	pc[k] = pc[max_j];
	pc[max_j] = pk;
}

template <typename T>
bool Matrix2D<T>::LUDecomposition() {
	std::vector< std::vector<double> >& mlu = this->matrix_lu_;
	std::vector<int>& pr = this->matrix_pr_;
	std::vector<int>& pc = this->matrix_pc_;
	int s = std::min(this->n_rows_, this->n_cols_);
	double divisor = 1.0;

	for (int i = 0; i < s; i++) {
		if (i <= 0) mlu = this->matrix;
		else {
			for (int r = i; r < this->n_rows_; r++) {
				for (int c = i; c < this->n_cols_; c++) {
					mlu[pr[r]][pc[c]] -= mlu[pr[r]][pc[i-1]] * mlu[pr[i-1]][pc[c]];
					if (isfinite(mlu[pr[r]][pc[c]]) && (abs(mlu[pr[r]][pc[c]]) < PRECISION)) mlu[pr[r]][pc[c]] = 0;
				}
			}
		}

		this->Pivot(i);
		divisor = 1 / mlu[pr[i]][pc[i]];
		for (int r = i + 1; r < this->n_rows_; r++) {
			mlu[pr[r]][pc[i]] *= divisor;
			if (isfinite(mlu[pr[r]][pc[i]]) && (abs(mlu[pr[r]][pc[i]]) < PRECISION)) mlu[pr[r]][pc[i]] = 0;
		}
	}
	return true;
}  

template <typename T>
static Matrix2D<T> Matrix2D<T>::LUComposition(std::vector< std::vector<double> > m_lu) {
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
void Matrix2D<T>::DisplayLU(unsigned pcs, int x) {
	std::vector< std::vector<double> >& mlu = this->matrix_lu_;
	std::vector<int>& pr = this->matrix_pr_;
	std::vector<int>& pc = this->matrix_pc_;
	std::cout.setf(std::ios::right);
	for (int i = 0; i < this->n_rows_; i++) {
		for (int j = 0; j < this->n_cols_; j++) {
			std::cout.width(pcs);
			if (x == 0) {
				std::cout << mlu[pr[i]][pc[j]];
			}
			else {
				std::cout << mlu[i][j];
			}
			if (j < this->n_cols_ - 1) std::cout << '\t';
		}
		std::cout << std::endl;
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

template <typename T>
double Matrix2D<T>::Determinant() {
	std::vector< std::vector<double> >& mlu = this->matrix_lu_;
	std::vector<int>& pr = this->matrix_pr_;
	std::vector<int>& pc = this->matrix_pc_;
	double det = 1;
	if (this->n_rows_ != this->n_cols_) return 0;
	for (int i = 0; i < this->n_rows_; i++) det *= mlu[pr[i]][pc[i]];
	return det;
}

#endif