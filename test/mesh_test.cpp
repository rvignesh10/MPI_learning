#include <iostream>
#include "../Header/MatrixAlgebra.hpp"
#include "../Header/FiniteElement.hpp"

int main(){
    const int N = 5;
    Mesh<N> m(0.0f,1.0f);
    m.gen_1D_mesh();
    m.disp_mesh();
    std::cout<< m.get_location(2)<<std::endl;

    const int order = 4;
    FiniteElementSpace<order,N> fes(m);
    fes.display_FE_Info(2);
    return 0;
}