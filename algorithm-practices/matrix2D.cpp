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
	std::vector< std::vector<T> >& m = this->matrix;
	std::vector<int>& pr = this->matrix_pr_;
	std::vector<int>& pc = this->matrix_pc_;
	T max_val = -128;
	int max_i = 0;
	int max_j = 0;
	int pk = 0;
	for (int i = k; i < this->n_rows_; i++) {
		for (int j = k; j < this->n_cols_; j++) {
			if (m[pr[i]][pc[j]] > max_val) {
				max_val = m[pr[i]][pc[j]];
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

	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			std::cout << m[pr[i]][pc[j]];
			if (j == 4) std::cout << std::endl;
			else std::cout << '\t';
		}
	}
	std::cout << std::endl;
}

template <typename T>
bool Matrix2D<T>::LUDecomposition() {
	std::vector< std::vector<T> >& m = this->matrix;
	std::vector< std::vector<double> >& mlu = this->matrix_lu_;
	std::vector<int>& pr = this->matrix_pr_;
	std::vector<int>& pc = this->matrix_pc_;
	int i = 0;
	int s = 0;
	double temp, divisor = 1.0;

	if (this->n_rows_ <= this->n_cols_) s = this->n_rows_ - 1;
	else s = this->n_cols_;
	while (i < s) {
		this->Pivot(i);
		if (i == 0) {
			for (int c = 0; c < this->n_cols_; c++) mlu[pr[0]][pc[c]] = m[pr[0]][pc[c]];
		}
		divisor = 1 / mlu[pr[i]][pc[i]];
		for (int r = i + 1; r < this->n_rows_; r++) {
			if (i > 0) {
				temp = m[pr[r]][pc[i]] - mlu[pr[r]][pc[0]] * mlu[pr[0]][pc[i]];
				for (int j = 1; j < i; j++) temp -= mlu[pr[r]][pc[j]] * mlu[pr[j]][pc[i]];
			}
			else {
				temp = m[pr[r]][pc[0]];
			}
			mlu[pr[r]][pc[i]] = temp * divisor;
			if (isfinite(mlu[pr[r]][pc[i]])) {
				mlu[pr[r]][pc[i]] = (abs(mlu[pr[r]][pc[i]]) >= PRECISION) ? mlu[pr[r]][pc[i]] : 0;
			}
		}
		if (!((this->n_rows_ > this->n_cols_) && (i == s - 1))) {
			for (int c = i + 1; c < this->n_cols_; c++) {
				temp = m[pr[i+1]][pc[c]] - mlu[pr[i + 1]][pc[0]] * mlu[pr[0]][pc[c]];
				for (int j = 1; j < i + 1; j++) temp -= mlu[pr[i + 1]][pc[j]] * mlu[pr[j]][pc[c]];
				if (isfinite(temp)) mlu[pr[i + 1]][pc[c]] = (abs(temp) >= PRECISION) ? temp : 0;
			}
		}
		++i;
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
void Matrix2D<T>::DisplayLU() {
	for (int i = 0; i < this->n_rows_; i++) {
		for (int j = 0; j < this->n_cols_; j++) {
			std::cout << this->matrix_lu_[i][j];
			if (j < this->n_cols_ - 1) std::cout << "\t\t";
		}
		std::cout << std::endl << std::endl;
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