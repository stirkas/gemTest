#include <iostream>

#include "gem_util_c.h"
#include "gem_equil_c.h"

using namespace std;

void new_gem_equil_c_()
{
   //cout << "Allocating C arrays for gem equil." << endl;
   Allocate2dPointerArrays_gem_equil();
   Allocate3dPointerArrays_gem_equil();
   Allocate4dPointerArrays_gem_equil();
}

void Allocate2dPointerArrays_gem_equil()
{
   //Allocate2dPointerArray<int>(array2dptr, array2dptr_c, ix, iy);
   //Allocate2dPointerArray<complex<float>>(complex2dptr, complex2dptr_c,ix,iy);
}

void Allocate3dPointerArrays_gem_equil()
{
   //Allocate3dPointerArray<int>(array3dptr, array3dptr_c, ix, iy, iz);
}

void Allocate4dPointerArrays_gem_equil()
{
   //Allocate4dPointerArray<int>(array4dptr, array4dptr_c, ix, iy, iz, it);
}