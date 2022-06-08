#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

#include <iostream>
#include "MatrixAlgebra.hpp"
#include "IntegrationRule.hpp"

template<const int Nelem>
class Mesh {
    private: 
        float xlim1;
        float xlim2;
        float dx; 
        Vector<float,Nelem+ONE> x;
    public:
        Mesh(float x1, float x2);
        void gen_1D_mesh();
        float get_location(int idx);
        void disp_mesh();
};

template<const int Nelem>
Mesh<Nelem>::Mesh(float x1, float x2): xlim1(0.0f), xlim2(1.0f) {
    xlim1 = x1;
    xlim2 = x2;
    
    dx    = (xlim2-xlim1)/(float)Nelem;
    gen_1D_mesh();
}

template<const int Nelem>
void Mesh<Nelem>::gen_1D_mesh(){
    x.setValue((int)0,xlim1);
    for (int i=1; i<Nelem+ONE; i++){
        x.setValue(i,xlim1+ (float)i*dx );
    }
}

template<const int Nelem>
float Mesh<Nelem>::get_location(int idx){
    return (float)x.getValue(idx);
}

template<const int Nelem>
void Mesh<Nelem>::disp_mesh(){
    x.displayVector();
}

/* -------------------------------------------------------------------------------- */
template<const int order>
class FiniteElement {

    private: 
        int ElemID;
        float DOF_loc;
        Vector<float,TWO> ElemBdr;
        Vector<float,order+ONE> InterpPts;
        Vector<int,order+ONE> InterpDOF;
        Vector<float,TWO+INT_PTS> IntPts;
        Vector<float,TWO+INT_PTS> IntWts;
    
    public:
        //FiniteElement();
        void set_Element_Info(int id, float current_loc, float next_loc);
        void print_Element_Info();
        void get_Integration_Info(Vector<float,INT_PTS+TWO> &p, Vector<float,INT_PTS+TWO> &w);
        void get_Boundary_Info(Vector<float,TWO> &b){b = ElemBdr;}
};

template<const int order>
void FiniteElement<order>::set_Element_Info(int id, float current_loc, float next_loc){
    ElemID = id+ONE;
    DOF_loc = (current_loc+next_loc)/2.0f ;
    ElemBdr.setValue(ZERO,current_loc);
    ElemBdr.setValue(ONE,next_loc);
    IntegrationRule<float,TWO+INT_PTS>(current_loc,next_loc,IntPts,IntWts);
}

template<const int order>
void FiniteElement<order>::print_Element_Info(){
    std::cout << "Element ID is: " << ElemID << std::endl;
    std::cout << "Element Center is: " << DOF_loc << std::endl << "Integration points are:";
    IntPts.displayVector();
    std::cout << std::endl << "Integration weights are:" << std::endl;
    IntWts.displayVector();
}

template<const int order>
void FiniteElement<order>::get_Integration_Info(Vector<float,INT_PTS+TWO> &p, Vector<float,INT_PTS+TWO> &w){
    p = IntPts;
    w = IntWts;
}

// /* -------------------------------------------------------------------------------- */
template<const int order, const int Nelem>
class FiniteElementSpace {
    private:
        FiniteElement<order> fe[Nelem];

    public:
        FiniteElementSpace(Mesh<Nelem> m);
        void display_FE_Info(int id){fe[id].print_Element_Info();}
        void get_FE_Integration_Info(int id, Vector<float,INT_PTS+TWO> &p, Vector<float,INT_PTS+TWO> &w);
        void get_FE_Boundary_Info(int id, Vector<float,TWO> &b);
        ~FiniteElementSpace(){std::cout<<"Object is deleted"<<std::endl;}
};

template<const int order, const int Nelem>
FiniteElementSpace<order,Nelem>::FiniteElementSpace(Mesh<Nelem> m){
    for (int i=0; i<Nelem; i++){
        fe[i].set_Element_Info(i, m.get_location(i), m.get_location(i+1));
    }
}

template<const int order, const int Nelem>
void FiniteElementSpace<order,Nelem>::get_FE_Integration_Info(int id, Vector<float,INT_PTS+TWO> &p, Vector<float,INT_PTS+TWO> &w){
    fe[id].get_Integration_Info(p,w);
}

template<const int order, const int Nelem>
void FiniteElementSpace<order,Nelem>::get_FE_Boundary_Info(int id, Vector<float,TWO> &b){
    fe[id].get_Boundary_Info(b);
}

#endif