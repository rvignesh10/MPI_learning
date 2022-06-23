#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

#include <iostream>
#include <string>
#include "MatrixAlgebra.hpp"
#include "IntegrationRule.hpp"
#include "LagrangePoly.hpp"

class AddDomainIntegrators;

template<const int Nelem>
class Mesh {
    private: 
        double xlim1;
        double xlim2;
        double dx; 
        Vector<double,Nelem+ONE> x;
    public:
        Mesh(double x1, double x2);
        void gen_1D_mesh();
        double get_location(int idx);
        void getX(Vector<double,Nelem+ONE> &x_out);
        void getDX(double &dx_out){dx_out = dx;}
        void disp_mesh();
};

template<const int Nelem>
Mesh<Nelem>::Mesh(double x1, double x2): xlim1(0.0), xlim2(1.0) {
    xlim1 = x1;
    xlim2 = x2;
    
    dx    = (xlim2-xlim1)/(double)Nelem;
    gen_1D_mesh();
}

template<const int Nelem>
void Mesh<Nelem>::gen_1D_mesh(){
    //x.setValue((int)0,xlim1);
    for (int i=0; i<Nelem+ONE; i++){
        x.setValue(i,xlim1+ (double)i*dx );
    }
}

template<const int Nelem>
double Mesh<Nelem>::get_location(int idx){
    return (double)x.getValue(idx);
}

template<const int Nelem>
void Mesh<Nelem>::getX(Vector<double,Nelem+ONE> &x_out){
    x_out = x;
}

template<const int Nelem>
void Mesh<Nelem>::disp_mesh(){
    x.displayVector();
}

/* --------------------------------------------------------------------------------------------------------------- */
template<const int order>
class FiniteElement {

    private: 
        int ElemID;
        Vector<double,TWO> ElemBdr;
        Vector<double,order+ONE> InterpPts;
        Vector<int,order+ONE> InterpDOF;
        Vector<double,TWO+INT_PTS> IntPts;
        Vector<double,TWO+INT_PTS> IntWts;
        Matrix<double,order+ONE,TWO+INT_PTS> ShapeFunc;
        Matrix<double,order+ONE,TWO+INT_PTS> DShapeFunc;
        double Jacobian;
    public:
        //FiniteElement();
        void set_Element_Info(int id, int &k, double current_loc, double next_loc);
        void set_Element_Info(int id, int &k, double current_loc, double next_loc, 
                              int DOF_prevWorld);
        void print_Element_Info();
        void get_Integration_Info(Vector<double,INT_PTS+TWO> &p, Vector<double,INT_PTS+TWO> &w);
        void get_Integration_Info(Vector<double,INT_PTS+TWO> &w, 
                                  Vector<int,order+ONE> &DOF,
                                  double &jac,
                                  Matrix<double,order+ONE,INT_PTS+TWO> *DShape = new Matrix<double,order+ONE,INT_PTS+TWO>,
                                  Matrix<double,order+ONE,INT_PTS+TWO> *Shape  = new Matrix<double,order+ONE,INT_PTS+TWO>);
        void get_Boundary_Info(Vector<double,TWO> &b){b = ElemBdr;}
        void H1_ContinuousShapeFn();
};

template<const int order>
void FiniteElement<order>::set_Element_Info(int id, int &k, double current_loc, double next_loc){
    ElemID = id+ONE;
    Jacobian= (next_loc - current_loc)/2.0 ;
    ElemBdr.setValue(ZERO,current_loc); ElemBdr.setValue(ONE,next_loc);
    Mesh<order> local_mesh(current_loc,next_loc);
    local_mesh.getX(InterpPts);
    IntegrationRule<double,TWO+INT_PTS>(current_loc,next_loc,IntPts,IntWts);
    H1_ContinuousShapeFn();
    for (int i=0; i<order+ONE; i++){
        InterpDOF.setValue(i,k);
        ++k;
    }
    k -= 1;
}

template<const int order>
void FiniteElement<order>::set_Element_Info(int id, int &k, double current_loc, double next_loc,
                                             int DOF_prevWorld){
    
    ElemID = id+ONE;
    Jacobian= (next_loc - current_loc)/2.0 ;
    ElemBdr.setValue(ZERO,current_loc); ElemBdr.setValue(ONE,next_loc);
    Mesh<order> local_mesh(current_loc,next_loc);
    local_mesh.getX(InterpPts);
    IntegrationRule<double,TWO+INT_PTS>(current_loc,next_loc,IntPts,IntWts);
    H1_ContinuousShapeFn();
    for (int i=0; i<order+ONE; i++){
        InterpDOF.setValue(i,(DOF_prevWorld-ONE)+k);
        ++k;
    }
    k -= 1;
}

template<const int order>
void FiniteElement<order>::print_Element_Info(){
    std::cout << "Element ID is: " << ElemID << std::endl;
    std::cout << "Integration points are:" << std::endl;
    IntPts.displayVector();
    std::cout << std::endl << "Integration weights are:" << std::endl;
    IntWts.displayVector();
    std::cout << std::endl << "Degrees of Freedom are located at: " << std::endl;
    InterpPts.displayVector();
    std::cout << std::endl << "Shape Function values are: " << std::endl;
    ShapeFunc.displayMatrix();
    std::cout << "DShape Function values are: " << std::endl;
    DShapeFunc.displayMatrix();
    std::cout << "Interpolation DOF are: " << std::endl;
    InterpDOF.displayVector();std::cout << std::endl;
}

template<const int order>
void FiniteElement<order>::get_Integration_Info(Vector<double,INT_PTS+TWO> &p, Vector<double,INT_PTS+TWO> &w){
    p = IntPts;
    w = IntWts;
}

template<const int order>
void FiniteElement<order>::get_Integration_Info(Vector<double,INT_PTS+TWO> &w,
                                                Vector<int,order+ONE> &DOF,
                                                double &jac, 
                                                Matrix<double,order+ONE,INT_PTS+TWO> *DShape,
                                                Matrix<double,order+ONE,INT_PTS+TWO> *Shape){
    w = IntWts;
    DOF = InterpDOF;
    *Shape =  ShapeFunc;
    *DShape = DShapeFunc;
    jac = Jacobian;
}

template<const int order>
void FiniteElement<order>::H1_ContinuousShapeFn(){
    for(int i=0; i<order+ONE; i++){
        Vector<double,TWO+INT_PTS> t1, t2;
        LagrangePolynomial(IntPts, InterpPts.getValue(i), InterpPts, t1, t2);
        ShapeFunc.setRow(i,t1);
        DShapeFunc.setRow(i,t2);
    }
}

// /* ------------------------------------------------------------------------------------------------------------ */
template<const int order, const int Nelem>
class FiniteElementSpace {
    private:
        FiniteElement<order> fe[Nelem];

    public:
        FiniteElementSpace(Mesh<Nelem> m);
        FiniteElementSpace(Mesh<Nelem> m, int world_rank);
        void display_FE_Info(int id){fe[id].print_Element_Info();}
        void get_FE_Integration_Info(int id, Vector<double,INT_PTS+TWO> &p, Vector<double,INT_PTS+TWO> &w);
        void get_FE_Boundary_Info(int id, Vector<double,TWO> &b);
        friend class AddDomainIntegrators;
};

template<const int order, const int Nelem>
FiniteElementSpace<order,Nelem>::FiniteElementSpace(Mesh<Nelem> m){
    int k = 1;
    for (int i=0; i<Nelem; i++){
        fe[i].set_Element_Info(i, k, m.get_location(i), m.get_location(i+ONE));
    }
}

template<const int order, const int Nelem>
FiniteElementSpace<order,Nelem>::FiniteElementSpace(Mesh<Nelem> m, int world_rank){
    // Total DOF in each batch
    int NelemTot = world_rank*Nelem;
    int DOF_prevWorld = (NelemTot+1) + NelemTot * (order-ONE); 
    int k = 1;
    for (int i=0; i<Nelem; i++){
        fe[i].set_Element_Info(i, k, m.get_location(i), m.get_location(i+1), DOF_prevWorld);
    }
}

template<const int order, const int Nelem>
void FiniteElementSpace<order,Nelem>::get_FE_Integration_Info(int id, Vector<double,INT_PTS+TWO> &p, Vector<double,INT_PTS+TWO> &w){
    fe[id].get_Integration_Info(p,w);
}

template<const int order, const int Nelem>
void FiniteElementSpace<order,Nelem>::get_FE_Boundary_Info(int id, Vector<double,TWO> &b){
    fe[id].get_Boundary_Info(b);
}

/* -------------------------------------------------------------------------------------------------------------------------------------------- */

class AddDomainIntegrators{
    private:
    public:
    template<const int order, const int Nelem>
    void DiffusionIntegrator(FiniteElementSpace<order,Nelem> fes, AppendList **Head_diff );
    template<const int order, const int Nelem>
    void MassIntegrator(FiniteElementSpace<order,Nelem> fes, AppendList **Head_mass);
};


template<const int order, const int Nelem>
void AddDomainIntegrators::DiffusionIntegrator(FiniteElementSpace<order,Nelem> fes,
                                               AppendList **Head_diff){
    Vector<double,INT_PTS+TWO> weights;
    Vector<int,order+ONE> DOF;
    Matrix<double,order+ONE,INT_PTS+TWO> DShape;
    double jac;
    Vector<double,INT_PTS+TWO> Nix, Njx,f;
    AppendList *Node = new AppendList;
    AppendList *temp = Node;
    *Head_diff = Node;
    for (int e=0; e<Nelem; e++){
        fes.fe[e].get_Integration_Info(weights,DOF,jac,&DShape);
        for (int i=0; i<order+ONE; i++){
            DShape.getRow(i,Nix);
            for (int j=0; j<order+ONE; j++){ 
                Node->i = DOF.getValue(i)-1;
                Node->j = DOF.getValue(j)-1;
                DShape.getRow(j,Njx);
                Nix.ElementMultiplication(Njx,f);
                Node->value = weights.dotProduct(f) * jac;
                Node = new AppendList;
                temp->next = Node;
                temp = Node;
            }
        }
    }
    Node->next = NULL;
}

template<const int order, const int Nelem>
void AddDomainIntegrators::MassIntegrator(FiniteElementSpace<order,Nelem> fes,
                                               AppendList **Head_mass){
    Vector<double,INT_PTS+TWO> weights;
    Vector<int,order+ONE> DOF;
    Matrix<double,order+ONE,INT_PTS+TWO> Shape;
    double jac;
    Vector<double,INT_PTS+TWO> Ni, Nj,f;
    AppendList *Node = new AppendList;
    AppendList *temp = Node;
    *Head_mass = Node;
    for (int e=0; e<Nelem; e++){
        fes.fe[e].get_Integration_Info(weights,DOF,jac,new Matrix<double,order+ONE,INT_PTS+TWO>,&Shape);
        for (int i=0; i<order+ONE; i++){
            Shape.getRow(i,Ni);
            for (int j=0; j<order+ONE; j++){ 
                Node->i = DOF.getValue(i)-1;
                Node->j = DOF.getValue(j)-1;
                Shape.getRow(j,Nj);
                Ni.ElementMultiplication(Nj,f);
                Node->value = weights.dotProduct(f) * jac;
                Node = new AppendList;
                temp->next = Node;
                temp = Node;
            }
        }
    }
    Node->next = NULL;
}


/* ---------------------------------------------------------------------------------------------------------------------------------------- */

#endif