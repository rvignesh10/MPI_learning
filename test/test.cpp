#include <iostream>
#include "../Header/MatrixAlgebra.hpp"

int main(){
    const int N = 5;
    double v[N] = {1.1,2.2,3.3,4.4,5.5};
    Vector<int,N> a1;
    Vector<double,N> a2(v); 
    Matrix<float,2,3> m;

    a1.setValue(6,100);
    a1.displayVector();
    a2.displayVector();
    
    return 0;
}