#include <iostream>
#include <mpi.h>
#include <math.h>
//#include "../Header/IntegrationRule.hpp"
#include "../Header/FiniteElement.hpp"


int main(int arg, char *argv[]){

    // Initialize MPI environment
    MPI_Init(NULL,NULL);

    // Get the number of processes available
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    // Get rank of process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

    float g_xlim1 = 0.0f;
    float g_xlim2 = 10.0f;

    float DX = (g_xlim2-g_xlim1)/(float)world_size;

    float l_xlim1 = g_xlim1 + (float)world_rank*DX;
    float l_xlim2 = l_xlim1 + DX;

    const int Nelem = 100;
    Mesh<Nelem> m(l_xlim1,l_xlim2);
    FiniteElementSpace<FOUR,Nelem> fes(m);

    Vector<float,INT_PTS+TWO> points;
    Vector<float,INT_PTS+TWO> weight;
    Vector<float,TWO> l_bdr;
    float l_int = 0.0f;
    for (int i=0; i<Nelem; i++){
        fes.get_FE_Integration_Info(i,points,weight);
        fes.get_FE_Boundary_Info(i,l_bdr);
        Vector<float,INT_PTS+TWO> fval;
        for (int j=0; j<INT_PTS+TWO; j++){
            fval.setValue(j,powf(points.getValue(j),2)); // function is x^2
        }
        float s = fval.dotProduct(weight)*(l_bdr.getValue(1)-l_bdr.getValue(0))*0.5f;
        l_int += s;
    }
    
    float g_int;
    if (world_rank==0){
        g_int = l_int;
        for (int i=1; i<world_size; i++){
            MPI_Recv(&l_int,1,MPI_FLOAT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            g_int += l_int;
        }
        std::cout<< "The global numerical integration of the function x^2 from (0,1) is:" << g_int << std::endl;
        std::cout<< "The global exact integration of the function x^2 from (0,1) is:" <<(powf(g_xlim2,3.0f)/3.0f)<< std::endl;
    }
    else{
        //std::cout << "Local integration from world rank " << world_rank << "is:" << l_int << std::endl; 
        MPI_Send(&l_int, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

}