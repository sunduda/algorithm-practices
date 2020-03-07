// main.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "tensor.h"

int main() {
    std::vector<std::size_t> dims({8, 9, 3, 14, 7});
    tensor_toolkit::tensor<int> mtest(dims);
    int val = 0;
    for (auto i4 = 0; i4 < dims[4]; i4++) {
        for (auto i3 = 0; i3 < dims[3]; i3++) {
            for (auto i2 = 0; i2 < dims[2]; i2++) {
                for (auto i1 = 0; i1 < dims[1]; i1++) {
                    for (auto i0 = 0; i0 < dims[0]; i0++) {
                        mtest[i0][i1][i2][i3][i4] = val;
                        val++;
                    }
                }
            }
        }
    }
    std::cout << mtest[0][1][1][5][3] << std::endl;
    mtest[0][1][1][5][3] = 687;
    std::cout << mtest[0][1][1][5][3] << std::endl;
    std::vector<double> dvec{0.0, 1.0, 2.0, 3.0};
    tensor_toolkit::tensor<int> ntest(dims, mtest.get_values());
    tensor_toolkit::tensor<double> xtest(dims, dvec);
    std::cout << ntest[0][1][1][5][3] << std::endl;
    return 0;
}