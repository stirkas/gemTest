#pragma once

#include <complex>
#include <cmath>

// this is going to hold all the externed global variavles that
// were held in gem_com

extern double tcurr,totvol,dt;

extern int myid,imx,nmx,jmx,jcnt,nsubd,master, gclr, tclr, kmx;

//1d pointers
extern double *vol_ptr,*rmsapa_ptr,*rmsphi_ptr,*avewe_ptr,*mdhis_ptr,*mdhisa_ptr;

extern double *mdhisb_ptr,*mdhisc_ptr,*mdhisd_ptr;

//2d pointers
extern double *pfle_es_ptr,*avewi_ptr,*efle_es_ptr,*yyre_ptr,*yyim_ptr;

extern double *yyamp_ptr, *efle_em_ptr, *pfle_em_ptr;

inline double **pfle_es_cptr = nullptr,**avewi_cptr = nullptr,**efle_es_cptr = nullptr,**yyre_cptr = nullptr ,**yyim_cptr = nullptr;

inline double **yyamp_cptr = nullptr, **efle_em_cptr = nullptr, **pfle_em_cptr = nullptr;

extern std::complex<double> *phihis_ptr,*aparhis_ptr;

inline std::complex<double> **phihis_cptr = nullptr,**aparhis_cptr = nullptr;

//3d pointers
extern double *pfl_es_ptr, *efl_es_ptr, *efl_em_ptr, *pfl_em_ptr;

inline double ***pfl_es_cptr = nullptr, ***efl_es_cptr = nullptr, ***pfl_em_cptr = nullptr, ***efl_em_cptr = nullptr;

//functions to allocate the pointers for arrays of dimension 2+
extern "C"
{
    void new_gem_com_c_();
}

void Allocate2dPointerArrays_gem_com();
void Allocate3dPointerArrays_gem_com();
void Allocate4dPointerArrays_gem_com();























