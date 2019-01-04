#ifndef _MATRIX2D_H
#define _MATRIX2D_H

#include <string>
#include <vector>
namespace myla {
	template <typename T>
	class Matrix2D {
	public:
		Matrix2D(int r, int c = r) {
			this->matrix.resize(r, std::vector<T>(c));
			this->matrix_lu_.resize(r, std::vector<long double>(c));
			this->n_rows_ = r;
			this->n_cols_ = c;
		}
		bool LUDecomposition();
		static Matrix2D<T> LUComposition(std::vector< std::vector<long double> > m_lu = this->matrix_lu_);
		void DisplayLU(std::string mode);
		Matrix2D<T> operator*(const Matrix2D &m);

		std::vector< std::vector<T> > matrix;

		~Matrix2D() {};

	private:
		void Pivot(int k);
		std::vector< std::vector<long double> > matrix_lu_;
		std::vector<int> matrix_p_;
		int n_rows_, n_cols_;
	};
}

#endif

