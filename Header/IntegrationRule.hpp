#ifndef INTEGRATION_RULE
#define INTEGRATION_RULE


#include <iostream>
#include <math.h>
#include "MatrixAlgebra.hpp"

// can only integrate polynomials upto degree 8 accurately.
// this is the limitation of this program

template<typename T, const int N>
void IntegrationRule(float xL, float xR, Vector<T,N> &v, Vector<T,N> &w){
    float w2 = (322.0f-13.0f*sqrtf(70.0f))/900.0f;
    float w1 = (322.0f+13.0f*sqrtf(70.0f))/900.0f;

    float p1 = (float)(1.0f/3.0f)*sqrtf(5.0f-2.0f*sqrtf((float)(10.0f/7.0f)));
    float p2 = (float)(1.0f/3.0f)*sqrtf(5.0f+2.0f*sqrtf((float)(10.0f/7.0f)));

    float xi[INT_PTS] = {-p2,-p1, 0.0f, p1, p2};
    float wt[INT_PTS] = {w2,w1,(float)(128.0f/125.0f),w1,w2};
    for (int i=0; i<N; i++){
        if (N != INT_PTS+TWO){
            std::cerr<<"Incompatible vector size"<<std::endl;
        }
        else{
            if (i==0){v.setValue(i,xL);w.setValue(i,0);}
            else if(i==N-ONE){v.setValue(i,xR);w.setValue(i,0);}
            else{
                v.setValue(i,(xL+xR)/2.0f + (xR-xL)*xi[i-ONE]/2.0f);
                w.setValue(i,wt[i-ONE]);
                }
        }
    }
}

// float pointsToLeft(){

// }

#endif