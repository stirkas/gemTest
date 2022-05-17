#include <iostream>

#include "gem_com_c.h"
#include "gem_util_c.h"

using namespace std;

//Define global pointers for header access.
double** pfle_es_cptr = nullptr;
double** avewi_cptr   = nullptr;
double** efle_es_cptr = nullptr;
double** yyre_cptr    = nullptr;
double** yyim_cptr    = nullptr;
double** yyamp_cptr   = nullptr;
double** efle_em_cptr = nullptr;
double** pfle_em_cptr = nullptr;
double** x2_cptr      = nullptr;
double** z2_cptr      = nullptr;
double** u2_cptr      = nullptr;
double** mu_cptr      = nullptr;
double** y2_cptr      = nullptr;
double** w3_cptr      = nullptr;
double** w2_cptr      = nullptr;
double** pzi_cptr     = nullptr;
double** xii_cptr     = nullptr;
double** z0i_cptr     = nullptr;
double** uoi_cptr     = nullptr;
double** u3_cptr      = nullptr;
double** x3_cptr      = nullptr;
double** y3_cptr      = nullptr;
double** z3_cptr      = nullptr;
double** eki_cptr     = nullptr;
double** nos_cptr     = nullptr;
double** ke_cptr      = nullptr;


complex<double>** phihis_cptr  = nullptr;
complex<double>** aparhis_cptr = nullptr;

double*** pfl_es_cptr  = nullptr;
double*** efl_es_cptr  = nullptr;
double*** pfl_em_cptr  = nullptr;
double*** efl_em_cptr  = nullptr;

double*** ex_cptr   = nullptr;
double*** ey_cptr   = nullptr;
double*** ez_cptr   = nullptr;
double*** delbx_cptr= nullptr;
double*** delby_cptr= nullptr;
double*** dpdz_cptr = nullptr;
double*** dadz_cptr = nullptr;
double*** apar_cptr = nullptr;

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
   //cout << "nsubd" << nsubd << "nmx" << nmx << endl;
   Allocate2dPointerArray<double>(avewi_ptr, avewi_cptr, 3, nmx+1);
   Allocate2dPointerArray<double>(efle_es_ptr, efle_es_cptr, nsubd, nmx+1);
   Allocate2dPointerArray<double>(yyre_ptr, yyre_cptr, jmx, 5);
   Allocate2dPointerArray<double>(yyim_ptr, yyim_cptr, jmx, 5);
   Allocate2dPointerArray<double>(yyamp_ptr, yyamp_cptr, jmx, 5);
   Allocate2dPointerArray<double>(efle_em_ptr, efle_em_cptr, nsubd, nmx+1);
   Allocate2dPointerArray<double>(pfle_em_ptr, pfle_em_cptr, nsubd, nmx+1);

   Allocate2dPointerArray<complex<double>>(phihis_ptr, phihis_cptr,7,jcnt-1+1);
   Allocate2dPointerArray<complex<double>>(aparhis_ptr, aparhis_cptr,7,jcnt-1+1);

   Allocate2dPointerArray<double>(x2_ptr, x2_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(z2_ptr, z2_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(u2_ptr, u2_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(mu_ptr, mu_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(w3_ptr, w3_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(w2_ptr, w2_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(pzi_ptr, pzi_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(xii_ptr, xii_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(uoi_ptr, uoi_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(u3_ptr, u3_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(x3_ptr, x3_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(y3_ptr, y3_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(z3_ptr, z3_cptr,nsmx, mmx);
   Allocate2dPointerArray<double>(eki_ptr, eki_cptr,nsmx, mmx);

   Allocate2dPointerArray<double>(nos_ptr, nos_cptr,nsmx, nmx);
   Allocate2dPointerArray<double>(ke_ptr, ke_cptr,nsmx, nmx+1);

}

void Allocate3dPointerArrays_gem_com()
{
   //3d arrays for spec()
   Allocate3dPointerArray<double>(pfl_es_ptr, pfl_es_cptr, 3, nsubd, nmx+1);
   Allocate3dPointerArray<double>(efl_es_ptr,efl_es_cptr, 3, nsubd, nmx+1);
   Allocate3dPointerArray<double>(pfl_em_ptr, pfl_em_cptr, 3, nsubd, nmx+1);
   Allocate3dPointerArray<double>(efl_em_ptr,efl_em_cptr, 3, nsubd, nmx+1);

   Allocate3dPointerArray<double>(ex_ptr,ex_cptr, nxpp+1, jmx+1, 1+1);
   Allocate3dPointerArray<double>(ey_ptr,ey_cptr, nxpp+1, jmx+1, 1+1);
   Allocate3dPointerArray<double>(ez_ptr,ez_cptr, nxpp+1, jmx+1, +11);
   Allocate3dPointerArray<double>(delbx_ptr,delbx_cptr, nxpp+1, jmx+1, 1+1);
   Allocate3dPointerArray<double>(delby_ptr,delby_cptr, nxpp+1, jmx+1, 1+1);
   Allocate3dPointerArray<double>(dpdz_ptr,dpdz_cptr, nxpp+1, jmx+1, 1+1);
   Allocate3dPointerArray<double>(dadz_ptr,dadz_cptr, nxpp+1, jmx+1, 1+1);
   Allocate3dPointerArray<double>(apar_ptr,apar_cptr, nxpp+1, jmx+1, 1+1);
}

void Allocate4dPointerArrays_gem_com()
{
   //Allocate4dPointerArray<int>(array4dptr, array4dptr_c, ix, iy, iz, it);
}