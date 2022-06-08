#include <iostream>
#include "../Header/MatrixAlgebra.hpp"

int main(){
    const int N = 5;
    Vector<int,N> a;
    Matrix<float,2,3> m;

    a.setValue(6,100);
    a.displayVector();
    
    std::cout<< "The length of created vector is:" << a.getLength_() << std::endl;
    return 0;
}