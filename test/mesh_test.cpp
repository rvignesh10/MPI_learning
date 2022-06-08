#include <iostream>
#include "../Header/MatrixAlgebra.hpp"
#include "../Header/FiniteElement.hpp"

int main(){
    const int N = 20;
    Mesh<N> m(-1.0f,1.0f);
    m.gen_1D_mesh();
    
    const int order = 4;
    FiniteElementSpace<order,N> fes(m);
    fes.display_FE_Info(THREE);
    return 0;
}