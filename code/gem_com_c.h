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

extern double **pfle_es_cptr,**avewi_cptr,**efle_es_cptr,**yyre_cptr,**yyim_cptr;

extern double **yyamp_cptr, **efle_em_cptr, **pfle_em_cptr;

extern std::complex<double> *phihis_ptr,*aparhis_ptr;

extern std::complex<double> **phihis_cptr,**aparhis_cptr;

//3d pointers
extern double *pfl_es_ptr, *efl_es_ptr, *efl_em_ptr, *pfl_em_ptr;

extern double ***pfl_es_cptr, ***efl_es_cptr, ***pfl_em_cptr, ***efl_em_cptr;

//functions to allocate the pointers for arrays of dimension 2+
extern "C"
{
    void new_gem_com_c_();
}

void Allocate2dPointerArrays_gem_com();
void Allocate3dPointerArrays_gem_com();
void Allocate4dPointerArrays_gem_com();





















