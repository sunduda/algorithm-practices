#ifndef _MATRIX2D_H
#define _MATRIX2D_H

#include <string>
#include <vector>

namespace xinda_linear_algebra {
	constexpr double PRECISION = 10e-12;

	template <typename T>
	class Matrix2D {
	public:
		Matrix2D(unsigned r, unsigned c = 0)
			:n_rows_(3)
			,n_cols_(3){
			Initialisation(r, c);
		}

		Matrix2D<T>& operator= (std::vector< std::vector<T> >& m);
		Matrix2D<T>& operator= (std::initializer_list< std::initializer_list<T> > il);
		bool LUDecomposition();
		static Matrix2D<T> LUComposition(std::vector< std::vector<double> > m_lu = this->matrix_lu_);
		void DisplayLU(unsigned pcs = 10);
		double Determinant();
		Matrix2D<T> operator*(const Matrix2D &m);

		~Matrix2D() {};

		friend std::ostream& operator<< (std::ostream& stream, const Matrix2D& matrix) {
			stream.setf(std::ios::right);
			for (long i = 0; i < matrix.n_rows_; i++) {
				for (long j = 0; j < matrix.n_cols_; j++) {
					stream.width(this->display_width_);
					stream << matrix[i][j];
					if (j < matrix.n_cols_ - 1) steam << '\t';
				}
				steam << std::endl;
			}
			return stream;
		}

	protected:
		bool Initialisation(const unsigned& r = 0, const unsigned& c = 0);

	private:
		void Pivot(unsigned k);
		std::vector< std::vector<T> > main_matrix_;
		std::vector< std::vector<double> > matrix_lu_;
		std::vector<unsigned> matrix_pr_, matrix_pc_;
		unsigned n_rows_, n_cols_;
		unsigned display_width_;
		bool is_decomposed_;
		bool is_invertible_;
	};
}

#endif

