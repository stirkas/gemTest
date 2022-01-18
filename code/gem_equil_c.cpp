#include <iostream>

#include "gem_util_c.h"
#include "gem_equil_c.h"

using namespace std;

double** dbdr_cptr   = nullptr; 
double** dbdth_cptr  = nullptr; 
double** grcgt_cptr  = nullptr; 
double** bfld_cptr   = nullptr;
double** radius_cptr = nullptr;
double** dydr_cptr   = nullptr;
double** qhat_cptr   = nullptr;
double** gr_cptr     = nullptr;
double** gxdgy_cptr  = nullptr; 
double** curvbz_cptr = nullptr;
double** bdcrvb_cptr = nullptr;
double** grdgt_cptr  = nullptr;
double** t0s_cptr    = nullptr;
double** capts_cptr  = nullptr; 
double** capns_cptr  = nullptr; 
double** xn0s_cptr   = nullptr;
double** vparsp_cptr = nullptr;


void new_gem_equil_c_(double s)
{
   //cout << "Allocating C arrays for gem equil." << endl;
   Allocate2dPointerArrays_gem_equil(s);
   Allocate3dPointerArrays_gem_equil();
   Allocate4dPointerArrays_gem_equil();
}

void Allocate2dPointerArrays_gem_equil(double s)
{
   //Allocate2dPointerArray<int>(array2dptr, array2dptr_c, ix, iy);
   //Allocate2dPointerArray<complex<float>>(complex2dptr, complex2dptr_c,ix,iy);
   Allocate2dPointerArray<double>(dbdr_ptr, dbdr_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(dbdth_ptr, dbdth_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(grcgt_ptr, grcgt_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(bfld_ptr, bfld_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(radius_ptr, radius_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(dydr_ptr, dydr_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(qhat_ptr,qhat_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(gr_ptr, gr_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(gxdgy_ptr, gxdgy_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(curvbz_ptr, curvbz_cptr, nr+1, ntheta+1);
   Allocate2dPointerArray<double>(grdgt_ptr, grdgt_cptr, nr+1, ntheta+1);

   Allocate2dPointerArray<double>(t0s_ptr, t0s_cptr, 5, nr+1);
   Allocate2dPointerArray<double>(capts_ptr, capts_cptr, 5, nr+1);
   Allocate2dPointerArray<double>(capns_ptr, capns_cptr, 5, nr+1);
   Allocate2dPointerArray<double>(xn0s_ptr, xn0s_cptr, 5, nr+1);
   Allocate2dPointerArray<double>(vparsp_ptr, vparsp_cptr, 5, nr+1);
}

void Allocate3dPointerArrays_gem_equil()
{
   //Allocate3dPointerArray<int>(array3dptr, array3dptr_c, ix, iy, iz);
}

void Allocate4dPointerArrays_gem_equil()
{
   //Allocate4dPointerArray<int>(array4dptr, array4dptr_c, ix, iy, iz, it);
}