#ifndef INTEGRATION_RULE
#define INTEGRATION_RULE


#include <iostream>
#include <math.h>
#include "MatrixAlgebra.hpp"

// can only integrate polynomials upto degree 8 accurately.
// this is the limitation of this program

template<typename T, const int N>
void IntegrationRule(double xL, double xR, Vector<T,N> &v, Vector<T,N> &w){
    double w2 = (322.0-13.0*sqrt(70.0))/900.0;
    double w1 = (322.0+13.0*sqrt(70.0))/900.0;

    double p1 = (double)(1.0/3.0)*sqrt(5.0-2.0*sqrt((double)(10.0/7.0)));
    double p2 = (double)(1.0/3.0)*sqrt(5.0+2.0*sqrt((double)(10.0/7.0)));

    double xi[INT_PTS] = {-p2,-p1, 0.0, p1, p2};
    double wt[INT_PTS] = {w2,w1,(double)(128.0/225.0),w1,w2};
    for (int i=0; i<N; i++){
        if (N != INT_PTS+TWO){
            std::cerr<<"Incompatible vector size"<<std::endl;
        }
        else{
            if (i==0){v.setValue(i,xL);w.setValue(i,0);}
            else if(i==N-ONE){v.setValue(i,xR);w.setValue(i,0);}
            else{
                v.setValue(i,(xL+xR)/2.0 + (xR-xL)*xi[i-ONE]/2.0);
                w.setValue(i,wt[i-ONE]);
                }
        }
    }
}

// double pointsToLeft(){

// }

#endif