#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

//this is what links the c++ function into the fortran program
extern "C"
{
    void srcbes_(double* biz, double* gam0, double* gam1);
}

//this is the function
void srcbes_(double* biz, double* gam0, double* gam1)
{
    //initialize variables and derefernece biz for the loop
    double t1, t2;
    double biz_val = *biz;

    //check for the value of biz and opperates accordingly
    if (biz_val > 3.75)
    {
        t2 = 1.0 / sqrt(*biz);
        t1=3.75/(*biz);
        *gam0=(t2*((((((((.00392377*t1-.01647633)*t1+.02635537)*t1-.02057706)*t1+.00916281)*t1-.00157565)*t1+.00225319)*t1+.01328592)*t1+.39894228));
        *gam1=(t2*((((((((-.00420059*t1+.01787654)*t1-.02895312)*t1+.02282967)*t1-.01031555)*t1+.00163801)*t1-.00362018)*t1-.03988024)*t1+.39894228));
    }
    else
    {
        t1=pow((*biz/3.75),2.0);
        t2=exp(-*biz);
        *gam0=(t2*((((((.0045813*t1+.0360768)*t1+.2659732)*t1+1.2067492)*t1+3.0899424)*t1+3.5156229)*t1+1.0));
        *gam1=(t2*(*biz)*((((((.00032411*t1+.00301532)*t1+.02658733)*t1+.15084934)*t1+.51498869)*t1+.87890594)*t1+.5));
    }
}