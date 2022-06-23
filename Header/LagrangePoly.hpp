#ifndef LAGRANGE_POLY
#define LAGRANGE_POLY

#include <iostream>
#include <math.h>
#include "MatrixAlgebra.hpp"

template<const int m, const int n>
void LagrangePolynomial(Vector<double,n> IntPts, double CenterLoc, Vector<double,m> InterpPts,
                         Vector<double,n> &ShapeFunc, Vector<double,n> &DShapeFunc){
    
    double abs_tol = 1.0e-12;
    for(int i=0; i<n; i++){
        double x = IntPts.getValue(i);
        double p = 1.0;
        int count = 0;
        Vector<double,m-ONE> Numerator;
        Vector<double,m-ONE> Denominator;
        for(int j=0; j<m; j++){
            if( abs(CenterLoc-InterpPts.getValue(j)) <= abs_tol ){
                continue;
            }
            else{
                Numerator.setValue(count, x-InterpPts.getValue(j));
                Denominator.setValue(count,CenterLoc-InterpPts.getValue(j));
                p *= ( x-InterpPts.getValue(j) )/( CenterLoc-InterpPts.getValue(j) );
                count++;
            }
        }
        ShapeFunc.setValue(i,(p<1.0)?p:1.0);
        double sum = 0.0;
        for(int l=0; l<m-ONE; l++){
            double p2 = 1.0;
            for(int u=0; u<m-ONE; u++){
                if(l==u){
                    p2 *= (1.0/Denominator.getValue(u));
                }
                else{
                    p2 *= (Numerator.getValue(u)/Denominator.getValue(u));
                }
            }
            sum += p2;
        }
        DShapeFunc.setValue(i,sum);
    }

}

#endif