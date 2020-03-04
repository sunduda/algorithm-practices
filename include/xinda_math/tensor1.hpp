#ifndef XINDA_MATH_TENSOR_HPP
#define XINDA_MATH_TENSOR_HPP

#include <vector>

namespace xinda_math {
	constexpr double PRECISION = 10e-12;
	constexpr unsigned long GIGABYTE = 1073741824;

	template <size_t D, typename T>
	class tensor1 {
	public:
		// Constructor
		tensor1(size_t (&size)[D] = {0 });
		bool reset(size_t(&size)[D] = { 0 });
		tensor1<D, T>& operator= (std::vector< std::vector<T> >& m);
		tensor1<D, T>& operator= (std::initializer_list< std::initializer_list<T> >& il);
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
		template <size_t dims>
		multidimensional_vector<dims, T>::type& select_dimension(multidimensional_vector<dims, T>::type &mv);
		multidimensional_vector<D, T>::type mTensor;
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

	template <size_t D, typename T>
	struct multidimensional_vector
	{
		typedef std::vector< typename multidimensional_vector<D - 1, T>::type > type;
	};

	long long pow(int x, int p);
};

#endif

