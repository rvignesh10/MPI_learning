#include <iostream>
#include "../Header/MatrixAlgebra.hpp"
#include "../Header/FiniteElement.hpp"

int main(){
    const int N = 20;
    Mesh<N> m(-10.0,10.0);
    m.gen_1D_mesh();
    
    const int order = 3;
    FiniteElementSpace<order,N> fes(m);
    fes.display_FE_Info(TWO);
    return 0;
}