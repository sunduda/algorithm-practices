#ifndef _MATRIX2D_H
#define _MATRIX2D_H

#include <string>
#include <vector>
namespace myla {
	template <typename T>
	class Matrix2D {
	public:
		Matrix2D(int r, int c = -1) {
			if (c < 0) c = r;
			this->matrix.resize(r, std::vector<T>(c));
			this->matrix_lu_.resize(r, std::vector<double>(c));
			this->n_rows_ = r;
			this->n_cols_ = c;
			for (int i = 0; i < r; i++) this->matrix_pr_.push_back(i);
			for (int j = 0; j < c; j++) this->matrix_pc_.push_back(j);
		}
		bool LUDecomposition();
		static Matrix2D<T> LUComposition(std::vector< std::vector<double> > m_lu = this->matrix_lu_);
		void DisplayLU();
		Matrix2D<T> operator*(const Matrix2D &m);

		std::vector< std::vector<T> > matrix;

		~Matrix2D() {};

	protected:

	private:
		void Pivot(int k);
		std::vector< std::vector<double> > matrix_lu_;
		std::vector<int> matrix_pr_, matrix_pc_;
		int n_rows_, n_cols_;
	};
}

#endif

