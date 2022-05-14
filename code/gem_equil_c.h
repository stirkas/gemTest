#pragma once

// this is going to hold all the externed global variavles that
// were held in gem_equil

extern double cn0e,rhoia, delz, rin, dr, dth, q0;

extern double mimp, s;

extern int nr2, ildu, itube, nr, iperidf, ntheta;

//1D arrays
extern double *xn0e_ptr,*t0e_ptr,*thfnz_ptr, *f_ptr, *jfn_ptr, *psip_ptr, *phincp_ptr, *psip2_ptr, *dipdr_ptr, *tgis_ptr, *sf_ptr;

extern double *psi_ptr, *cn0s_ptr;

//2d arrays
extern double *dbdr_ptr, *dbdth_ptr, *grcgt_ptr, *bfld_ptr, *radius_ptr, *dydr_ptr, *qhat_ptr, *gr_ptr, *gxdgy_ptr, *curvbz_ptr, *bdcrvb_ptr, *grdgt_ptr;
extern double **dbdr_cptr, **dbdth_cptr, **grcgt_cptr, **bfld_cptr, **radius_cptr, **dydr_cptr, **qhat_cptr, **gr_cptr, **gxdgy_cptr, **curvbz_cptr, **bdcrvb_cptr, **grdgt_cptr;

extern double *t0s_ptr, *capts_ptr, *capns_ptr, *xn0s_ptr, *vparsp_ptr;
extern double **t0s_cptr, **capts_cptr, **capns_cptr, **xn0s_cptr, **vparsp_cptr;

extern "C"
{
    void new_gem_equil_c_(double s);
}

void Allocate2dPointerArrays_gem_equil(double s);
void Allocate3dPointerArrays_gem_equil();
void Allocate4dPointerArrays_gem_equil();