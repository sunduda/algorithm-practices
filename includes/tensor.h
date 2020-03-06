//
// Created by Christopher on 3/4/20.
//

#ifndef TENSOR_TOOLKIT_TENSOR_H
#define TENSOR_TOOLKIT_TENSOR_H


#include <cstddef>
#include <vector>
#include <ostream>

namespace xinda_math {
    constexpr double PRECISION = 10e-12;
    constexpr unsigned long GIGABYTE = 1073741824;

    template<std::size_t D, typename T>
    class tensor {
    public:
        // Constructor
        template<typename cT>
        explicit tensor(cT (&sizes)[D], T &value = 0);

        template<typename cT>
        explicit tensor(cT (&sizes)[D], T value = 0);

        explicit tensor(const std::size_t sizes[D], T &value = 0);

        explicit tensor(const std::size_t sizes[D], T value = 0);

        // Destructor
        virtual ~tensor() = default;

    protected:
        std::vector<tensor<D - 1, T>> get_value() { return this->mTensor; };

    private:
        std::vector<tensor<D - 1, T>> mTensor;

        // Actual initialization method used in the constructor
        template<typename cT>
        void initializer(cT (&sizes)[D], T &value);

    };

    template<typename T>
    class tensor<1, T> {
    public:
        // Constructor
        template<typename cT>
        explicit tensor(cT (&sizes)[1], T &value = 0);

        template<typename cT>
        explicit tensor(cT (&sizes)[1], T value = 0);

        explicit tensor(const std::size_t sizes[1], T &value = 0);

        explicit tensor(const std::size_t sizes[1], T value = 0);

        // Destructor
        virtual ~tensor() = default;

    protected:
        std::vector<T> get_value() { return this->mTensor; };

    private:
        std::vector<T> mTensor;

        // Actual initialization method used in the constructor
        template<typename cT>
        void initializer(cT (&sizes)[1], T &value = 0);
    };

};


#endif //TENSOR_TOOLKIT_TENSOR_H

// TODO: Need to be defined methods
// tensor<D, T> &operator=(std::vector<std::vector<T> > &m);
// tensor<D, T> &operator=(std::initializer_list<std::initializer_list<T> > &il);
// std::vector<T> &operator[](const int &index);
// void set_element(const int &r = 0, const int &c = 0, const T &val = 0);
// tensor<D, T> operator*(const tensor &m);
// typename multidimensional_vector<D, T>::type &
// select_dimension(typename multidimensional_vector<mD, mT>::type &mv);

//// Overload operators <<
//friend std::ostream &operator<<(std::ostream &stream, tensor<D, T> &vtensor) {
//    stream.setf(std::ios::right);
//    stream.setf(std::ios::fixed);
//    vtensor.SetDisplayWidth(stream.precision());
//    for (int i = 0; i < vtensor.n_rows_; i++) {
//        for (int j = 0; j < vtensor.n_cols_; j++) {
//            stream.width(vtensor.display_width_);
//            stream << vtensor[i][j];
//            if (j < vtensor.n_cols_ - 1) stream << '\t';
//        }
//        stream << std::endl;
//    }
//    return stream;
//};
//
//// Overload operators >>
//friend std::ostream &operator>>(std::ostream &stream, tensor<D, T> &vtensor) {
//    stream.setf(std::ios::right);
//    stream.setf(std::ios::fixed);
//    vtensor.SetDisplayWidth(stream.precision());
//    for (int i = 0; i < vtensor.n_rows_; i++) {
//        for (int j = 0; j < vtensor.n_cols_; j++) {
//            stream.width(vtensor.display_width_);
//            stream << vtensor[i][j];
//            if (j < vtensor.n_cols_ - 1) stream << '\t';
//        }
//        stream << std::endl;
//    }
//    return stream;
//};


// TODO: Realize the following methods when define the child class 'matrix'
// bool LUDecomposition();
// static tensor<D, T> &LUComposition(const std::vector<std::vector<double> > &m_lu = matrix_lu_);
// void DisplayLU(int pcs = 10);
// double Determinant();
// void Pivot(int k);
// std::vector<std::vector<T> > main_matrix_;
// std::vector<std::vector<double> > matrix_lu_;
// std::vector<int> matrix_pr_, matrix_pc_;
// int n_rows_, n_cols_;
// int display_width_;
// bool is_decomposed_;
// bool is_invertible_;
// void FlagReset();
// void SetDisplayWidth(const int &pcs = 2);