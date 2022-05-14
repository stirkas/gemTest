#pragma once

#include <complex>
#include <cmath>

// this is going to hold all the externed global variavles that
// were held in gem_com

extern double tcurr,totvol,dt, vwidth, lx, lr0, tor, br0, ly, vexbsw, vparsw, pzcrit[5], r0, lz, dx, dy, dz;

extern int myid,imx,nmx, mmx, nsmx,jmx,jcnt,nsubd,master, gclr, tclr, kmx, xnplt, ipred, icorr, jpred, jcorr, nopz, noen, nowe, idg, iput, ierr, nopi[5], iorb;

extern int iflr, ipara, nonlin[5], peritr, kcnt, nxpp;

extern double pi;

extern int iorb;


//1d pointers
extern double *vol_ptr,*rmsapa_ptr,*rmsphi_ptr,*avewe_ptr,*mdhis_ptr,*mdhisa_ptr;

extern double *mdhisb_ptr,*mdhisc_ptr,*mdhisd_ptr, *mm_ptr, *mims_ptr, *q_ptr;

extern int *lr_ptr;

extern double*q_ptr;

extern int *tmm_ptr;

//2d pointers
extern double *pfle_es_ptr,*avewi_ptr,*efle_es_ptr,*yyre_ptr,*yyim_ptr;

extern double *yyamp_ptr, *efle_em_ptr, *pfle_em_ptr;

extern double *x3_ptr, *y3_ptr, *z3_ptr, *eki_ptr;

extern double **pfle_es_cptr,**avewi_cptr,**efle_es_cptr,**yyre_cptr,**yyim_cptr;

extern double **yyamp_cptr, **efle_em_cptr, **pfle_em_cptr;

extern double **x3_cptr, **y3_cptr, **z3_cptr, **eki_cptr;

extern double *x2_ptr, *z2_ptr, *u2_ptr, *mu_ptr, *y2_ptr, *w3_ptr, *w2_ptr, *pzi_ptr, *xii_ptr, *z0i_ptr;

extern double **x2_cptr, **z2_cptr, **u2_cptr, **mu_cptr, **y2_cptr, **w3_cptr, **w2_cptr, **pzi_cptr, **xii_cptr, **z0i_cptr;

extern double *uoi_ptr, *u3_ptr;

extern double **uoi_cptr, **u3_cptr;

extern std::complex<double> *phihis_ptr,*aparhis_ptr;

extern std::complex<double> **phihis_cptr,**aparhis_cptr;

//3d pointers
extern double *pfl_es_ptr, *efl_es_ptr, *efl_em_ptr, *pfl_em_ptr;

extern double ***pfl_es_cptr, ***efl_es_cptr, ***pfl_em_cptr, ***efl_em_cptr, ***phi_ptr;

extern double *ex_ptr, *ey_ptr, *ez_ptr, *delbx_ptr, *delby_ptr, *dpdz_ptr, *dadz_ptr, *apar_ptr;

extern double ***ex_cptr, ***ey_cptr, ***ez_cptr, ***delbx_cptr, ***delby_cptr, ***dpdz_cptr, ***dadz_cptr, ***apar_cptr;

//functions to allocate the pointers for arrays of dimension 2+
extern "C"
{
    void new_gem_com_c_();
}

void Allocate2dPointerArrays_gem_com();
void Allocate3dPointerArrays_gem_com();
void Allocate4dPointerArrays_gem_com();





















