#include <iostream>

#include "gem_com_c.h"
#include "gem_util_c.h"

using namespace std;

void new_gem_com_c_()
{
   //cout << "Allocating C arrays for gem com." << endl;
   Allocate2dPointerArrays_gem_com();
   Allocate3dPointerArrays_gem_com();
   Allocate4dPointerArrays_gem_com();
}

void Allocate2dPointerArrays_gem_com()
{
   //2d arrays for spec()
   Allocate2dPointerArray<double>(pfle_es_ptr, pfle_es_cptr,nsubd,nmx+1);//check y dim
   cout << "nsubd" << nsubd << "nmx" << nmx << endl;
   Allocate2dPointerArray<double>(avewi_ptr, avewi_cptr, 3, nmx+1);
   Allocate2dPointerArray<double>(efle_es_ptr, efle_es_cptr, nsubd, nmx+1);
   Allocate2dPointerArray<double>(yyre_ptr, yyre_cptr, jmx, 5);
   Allocate2dPointerArray<double>(yyim_ptr, yyim_cptr, jmx, 5);
   Allocate2dPointerArray<double>(yyamp_ptr, yyamp_cptr, jmx, 5);
   Allocate2dPointerArray<double>(efle_em_ptr, efle_em_cptr, nsubd, nmx+1);
   Allocate2dPointerArray<double>(pfle_em_ptr, pfle_em_cptr, nsubd, nmx+1);

   Allocate2dPointerArray<complex<double>>(phihis_ptr, phihis_cptr,7,jcnt-1+1);
   Allocate2dPointerArray<complex<double>>(aparhis_ptr, aparhis_cptr,7,jcnt-1+1);
}

void Allocate3dPointerArrays_gem_com()
{
   //3d arrays for spec()
   Allocate3dPointerArray<double>(pfl_es_ptr, pfl_es_cptr, 3, nsubd, nmx+1);
   Allocate3dPointerArray<double>(efl_es_ptr,efl_es_cptr, 3, nsubd, nmx+1);
   Allocate3dPointerArray<double>(pfl_em_ptr, pfl_em_cptr, 3, nsubd, nmx+1);
   Allocate3dPointerArray<double>(efl_em_ptr,efl_em_cptr, 3, nsubd, nmx+1);
}

void Allocate4dPointerArrays_gem_com()
{
   //Allocate4dPointerArray<int>(array4dptr, array4dptr_c, ix, iy, iz, it);
}
