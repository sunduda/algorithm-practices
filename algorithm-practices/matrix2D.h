#ifndef _MATRIX2D_H
#define _MATRIX2D_H

#include <string>
#include <vector>

namespace xinda_math {
	constexpr double PRECISION = 10e-12;

	template <typename T>
	class Matrix2D {
	public:
		// Constructor
		Matrix2D(int r = 0, int c = 0);
		Matrix2D<T>& operator= (std::vector< std::vector<T> >& m);
		Matrix2D<T>& operator= (std::initializer_list< std::initializer_list<T> >& il);
		std::vector<T>& operator[] (const int& index);
		void AssignElement(const int& r = 0, const int& c = 0, const T& val = 0);
		bool LUDecomposition();
		static Matrix2D<T>& LUComposition(const std::vector< std::vector<double> >& m_lu = matrix_lu_);
		void DisplayLU(int pcs = 10);
		double Determinant();
		Matrix2D<T> operator*(const Matrix2D &m);
		// Destructor
		virtual ~Matrix2D() {};

		// Overload operators
		friend std::ostream& operator<< (std::ostream& stream, Matrix2D<T>& matrix) {
			stream.setf(std::ios::right);
			stream.setf(std::ios::fixed);
			matrix.SetDisplayWidth(stream.precision());
			for (int i = 0; i < matrix.n_rows_; i++) {
				for (int j = 0; j < matrix.n_cols_; j++) {
					stream.width(matrix.display_width_);
					stream << matrix[i][j];
					if (j < matrix.n_cols_ - 1) stream << '\t';
				}
				stream << std::endl;
			}
			return stream;
		};

	protected:
		bool Initialisation(const int& r = 0, const int& c = 0);
		void FlagReset();
		void SetDisplayWidth(const int& pcs = 2);

	private:
		void Pivot(int k);
		std::vector< std::vector<T> > main_matrix_;
		std::vector< std::vector<double> > matrix_lu_;
		std::vector<int> matrix_pr_, matrix_pc_;
		int n_rows_, n_cols_;
		int display_width_;
		bool is_decomposed_;
		bool is_invertible_;
	};
};

#endif

