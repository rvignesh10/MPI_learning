#include <iostream>
#include "MatrixAlgebra.hpp"

int main(){
    Vector<int,5> a;
    Matrix<float,2,3> m;

    a.setValue(3,100);
    std::cout<< "The length of created vector is:" << a.getLength_() << std::endl;
    return 0;
}