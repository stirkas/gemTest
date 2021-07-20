#pragma once

// this is going to hold all the externed global variavles that
// were held in gem_equil

extern double cn0e,rhoia;

extern double mimp;

extern int nr2;

//1D arrays
extern double *xn0e_ptr,*t0e_ptr;

extern "C"
{
    void new_gem_equil_c_();
}

void Allocate2dPointerArrays_gem_equil();
void Allocate3dPointerArrays_gem_equil();
void Allocate4dPointerArrays_gem_equil();