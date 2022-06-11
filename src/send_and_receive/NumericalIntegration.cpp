/* This program integrates the function f(x) = x^2 within a domain numerically. Then it compares it with the
exact integral. It uses simple MPI_Send() and MPI_Recv() functions and performs numerical integration. 
Author - Vignesh Ramakrishnan 
PhD student at Optimal Design Lab, Rensselaer Polytechnic Institute */
#include <iostream>
#include <mpi.h>
#include <math.h>
#include "../../Header/FiniteElement.hpp"


int main(int arg, char *argv[]){

    // Initialize MPI environment
    MPI_Init(NULL,NULL);

    // Get the number of processes available
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    // Get rank of process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

    double g_xlim1 = 3.3;
    double g_xlim2 = 8.45;

    double DX = (g_xlim2-g_xlim1)/(double)world_size;

    double l_xlim1 = g_xlim1 + (double)world_rank*DX;
    double l_xlim2 = l_xlim1 + DX;

    const int Nelem = 5;
    Mesh<Nelem> m(l_xlim1,l_xlim2);
    FiniteElementSpace<FOUR,Nelem> fes(m);
    
    Vector<double,INT_PTS+TWO> points;
    Vector<double,INT_PTS+TWO> weight;
    Vector<double,TWO> l_bdr;
    double l_int = 0.0;

    for (int i=0; i<Nelem; i++){
        fes.get_FE_Integration_Info(i,points,weight);
        fes.get_FE_Boundary_Info(i,l_bdr);
        Vector<double,INT_PTS+TWO> fval;
        for (int j=0; j<INT_PTS+TWO; j++){
            fval.setValue(j,pow(points.getValue(j),2.0)); // function is x^2
        }
        double s = fval.dotProduct(weight)*(l_bdr.getValue(1)-l_bdr.getValue(0))*0.5;
        l_int += s;
    }
    
    double g_int;
    if (world_rank==0){
        g_int = l_int;
        for (int i=1; i<world_size; i++){
            MPI_Recv(&l_int,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            g_int += l_int;
        }
        std::cout<< "The global numerical integration of the function x^2 within domain is:" << g_int << std::endl;
        std::cout<< "The global exact integration of the function x^2 within domain is:" <<((pow(g_xlim2,3.0)-pow(g_xlim1,3.))/3.0)<< std::endl;
    }
    else{
        //std::cout << "Local integration from world rank " << world_rank << "is:" << l_int << std::endl; 
        MPI_Send(&l_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

}