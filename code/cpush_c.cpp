#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
#include <mpi.h>

#include "gem_com_c.h"
#include "gem_equil_c.h"
#include "gem_functions_c.h"

#include <unistd.h>

using namespace std;

void cpush_c_(int& n,int& ns)
{
    int n;
    double phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp;
    double wx0,wx1,wy0,wy1,wz0,wz1,w3old;
    int m,i,j,k,ns,l;
    int np_old,np_new;
    double vfac,vpar,vxdum,dum,xdot,ydot,zdot,pzdot,edot,pzd0,vp0;
    double rhog,xt,yt,zt,kap,xs,pidum,dum1,kaptp,kapnp,xnp;
    double b,th,r,enerb,cost,sint,qr,laps,sz,ter,bstar;
    double myke,mypfl_es[nsubd],mypfl_em[nsubd],myavewi;
    double myefl_es[nsubd],myefl_em[nsubd],mynos;
    double ketemp,pfltemp;
    double efltemp,nostemp;
    double sbuf[10],rbuf[10];
    double dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp;
    double grp,gxdgyp,rhox[4],rhoy[4],psp,pzp,vncp,vparspp,psip2p,bdcrvbp,curvbzp,dipdrp;

    sbuf[9] = 0.;
    rbuf[9] = 0.;
    myavewi = 0.;
    myke =0.;
    mypfl_es[nsubd] =0.;
    mypfl_em[nsubd] =0.;
    myefl_es[nsubd] =0.;
    myefl_em[nsubd] =0.;
    mynos=0.;
    ketemp=0.;
    pfltemp=0.;
    efltemp=0.;
    nostemp=0.;
    pidum = 1/(pi*2)^1.5*vwidth^3;

    for(m == 1; m < mm_ptr[ns]; m++)
    {
        r = x3_cptr[ns][m]-0.5*lx+lr0;        
        k = int(z3_cptr[ns][m]/delz);
        wz0 = ((k+1)*delz-z3_cptr[ns][m])/delz;
        wz1 = 1-wz0;
        th = wz0*thfnz_ptr[k]+wz1*thfnz_ptr[k+1];     
        i = int((r-rin)/dr);
        wx0 = (rin+(i+1)*dr-r)/dr;
        wx1 = 1.-wx0;
        k = int((th+pi)/dth);
        wz0 = (-pi+(k+1)*dth-th)/dth;
        wz1 = 1.-wz0;
        dbdrp   = wx0*wz0*dbdr_cptr[i][k]+wx0*wz1*dbdr_cptr[i][k+1]+wx1*wz0*dbdr_cptr[i+1][k]+wx1*wz1*dbdr_cptr[i+1][k+1];
        dbdtp   = wx0*wz0*dbdth_cptr[i][k]+wx0*wz1*dbdth_cptr[i][k+1]+wx1*wz0*dbdth_cptr[i+1][k]+wx1*wz1*dbdth_cptr[i+1][k+1];
        grcgtp  = wx0*wz0*grcgt_cptr[i][k]+wx0*wz1*grcgt_cptr[i][k+1]+wx1*wz0*grcgt_cptr[i+1][k]+wx1*wz1*grcgt_cptr[i+1][k+1]; 
        bfldp   = wx0*wz0*bfld_cptr[i][k]+wx0*wz1*bfld_cptr[i][k+1]+wx1*wz0*bfld_cptr[i+1][k]+wx1*wz1*bfld_cptr[i+1][k+1];
        radiusp = wx0*wz0*radius_cptr[i][k]+wx0*wz1*radius_cptr[i][k+1]+wx1*wz0*radius_cptr[i+1][k]+wx1*wz1*radius_cptr[i+1][k+1];
        dydrp   = wx0*wz0*dydr_cptr[i][k]+wx0*wz1*dydr_cptr[i][k+1]+wx1*wz0*dydr_cptr[i+1][k]+wx1*wz1*dydr_cptr[i+1][k+1];
        qhatp   = wx0*wz0*qhat_cptr[i][k]+wx0*wz1*qhat_cptr[i][k+1]+wx1*wz0*qhat_cptr[i+1][k]+wx1*wz1*qhat_cptr[i+1][k+1];
        grp     = wx0*wz0*gr_cptr[i][k]+wx0*wz1*gr_cptr[i][k+1]+wx1*wz0*gr_cptr[i+1][k]+wx1*wz1*gr_cptr[i+1][k+1];
        gxdgyp  = wx0*wz0*gxdgy_cptr[i][k]+wx0*wz1*gxdgy_cptr[i][k+1]+wx1*wz0*gxdgy_cptr[i+1][k]+wx1*wz1*gxdgy_cptr[i+1][k+1];       
        curvbzp = wx0*wz0*curvbz_cptr[i][k]+wx0*wz1*curvbz_cptr[i][k+1]+wx1*wz0*curvbz_cptr[i+1][k]+wx1*wz1*curvbz_cptr[i+1][k+1];
        bdcrvbp = wx0*wz0*bdcrvb_cptr[i][k]+wx0*wz1*bdcrvb_cptr[i][k+1]+wx1*wz0*bdcrvb_cptr[i+1][k]+wx1*wz1*bdcrvb_cptr[i+1][k+1];
        grdgtp  = wx0*wz0*grdgt_cptr[i][k]+wx0*wz1*grdgt_cptr[i][k+1]+wx1*wz0*grdgt_cptr[i+1][k]+wx1*wz1*grdgt_cptr[i+1][k+1];       
        fp = wx0*f_ptr[i]+wx1*f_ptr[i+1];        
        jfnp = wz0*jfn_ptr[k]+wz1*jfn_ptr[k+1];
        psipp = wx0*psip_ptr[i]+wx1*psip_ptr[i+1];        
        psp = wx0*psi_ptr[i]+wx1*psi_ptr[i+1];        
        ter = wx0*t0s_cptr[ns][i]+wx1*t0s_cptr[ns][i+1];        
        kaptp = wx0*capts_cptr[ns][i]+wx1*capts_cptr[ns][i+1];        
        kapnp = wx0*capns_cptr[ns][i]+wx1*capns_cptr[ns][i+1];        
        xnp = wx0*xn0s_cptr[ns][i]+wx1*xn0s_cptr[ns][i+1];        
        b=1.-tor+tor*bfldp;
        psip2p = wx0*psip2_ptr[i]+wx1*psip2_ptr[i+1];        
        dipdrp = wx0*dipdr_ptr[i]+wx1*dipdr_ptr[i+1];        
        pzp = mims_ptr[ns]*u3_cptr[ns][m]/b-q_ptr[ns]*psp/br0;
        vncp = wx0*phincp_ptr[i]+wx1*phincp_ptr[i+1];        
        vparspp = wx0*vparsp_cptr[ns][i]+wx1*vparsp_cptr[ns][i+1];               
        rhog=sqrt(2.*b*mu_cptr[ns][m]*mims_ptr[ns])/(q_ptr[ns]*b)*iflr;     
        rhox[1] = rhog*(1-tor)+rhog*grp*tor;
        rhoy[1] = rhog*gxdgyp/grp*tor;
        rhox[2] = -rhox[1];
        rhoy[2] = -rhoy[1];
        rhox[3] = 0;
        rhoy[3] = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor;
        rhox[4] = 0;
        rhoy[4] = -rhoy[3];
        //    calculate avg. e-field...
        //    do 1,2,4 point average, where lr is the no. of points..        
        phip=0.;
        exp1=0.;
        eyp=0.;
        ezp=0.;
        delbxp = 0.;
        delbyp = 0.;
        dpdzp = 0.;
        dadzp = 0.;
        aparp = 0;     
        //  4 pt. avg. written out explicitly for vectorization...
        for(l == 1; l < lr_ptr[1]; l++)
        {
           xs=x3_cptr[ns][m]+rhox_ptr[l]; //rwx(1,l)*rhog
           yt=y3_cptr[ns][m]+rhoy_ptr[l]; //(rwy(1,l)+sz*rwx(1,l))*rhog
           //   BOUNDARY
           xt=int(xs+800.*lx)% int(lx);
           yt=int(yt+800.*ly)%int(ly);
           xt = min(xt,lx-1.0e-8);
           yt = min(yt,ly-1.0e-8);        
           include "cpushngp.h";
        }
                
        exp1=exp1/4.;
        eyp=eyp/4.;
        ezp=ezp/4.;
        delbxp=delbxp/4.;
        delbyp=delbyp/4.;
        dpdzp = dpdzp/4.;
        dadzp = dadzp/4.;
        aparp = aparp/4     ;
        vfac = 0.5*(mims_ptr[ns]*u3_cptr[ns][m]^2 + 2.*mu_cptr[ns][m]*b);
        vp0 = 1./b^2*lr0/q0*qhatp*fp/radiusp*grcgtp;
        vp0 = vp0*vncp*vexbsw;
        vpar = u3_cptr[ns][m]-q_ptr[ns]/mims_ptr[ns]*aparp*nonlin_ptr[ns]*0.;
        bstar = b*(1+mims_ptr[ns]*vpar/(q_ptr[ns]*b)*bdcrvbp);
        enerb=(mu_cptr[ns][m]+mims_ptr[ns]*vpar*vpar/b)/q_ptr[ns]*b/bstar*tor;      
        kap = kapnp - (1.5-vfac/ter)*kaptp-vpar*mims_ptr[ns]/ter*vparspp*vparsw;
        dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp;
        vxdum = (eyp/b+vpar/b*delbxp)*dum1;
        xdot = vxdum*nonlin[ns] -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp;
        ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin[ns]+iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0+enerb/(bfldp^2)*psipp*lr0/q0/radiusp^2*(dbdrp*grp*2+dbdtp*grdgtp)-mims_ptr[ns]*vpar^2/(q_ptr[ns]*bstar*b)*(psip2p*grp^2/radiusp+curvbzp)*lr0/(radiusp*q0)-dipdrp/radiusp*mims_ptr[ns]*vpar^2/(q_ptr[ns]*bstar*b)*grcgtp*lr0/q0*qhatp;  
        zdot =  vpar*b/bstar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp+q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp-1./b^2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp-dipdrp/radiusp*mims_ptr[ns]*vpar^2/(q_ptr[ns]*bstar*b)*q0*br0*grcgtp/jfnp;     
        pzd0 = tor*(-mu_cptr[ns][m]/mims_ptr[ns]/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar+mu_cptr[ns][m]*vpar/(q_ptr[ns]*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp;
        pzdot = pzd0 + (q_ptr[ns]/mims_ptr[ns]*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp+q_ptr[ns]/mims_ptr[ns]*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara;        
        edot = q_ptr[ns]*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp)+q_ptr[ns]*pzdot*aparp*tor+q_ptr[ns]*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)-q_ptr[ns]*vpar*delbxp*vp0;      
        x3_cptr[ns][m] = x2_cptr[ns][m] + dt*xdot;
        y3_cptr[ns][m] = y2_cptr[ns][m] + dt*ydot;
        z3_cptr[ns][m] = z2_cptr[ns][m] + dt*zdot;
        u3_cptr[ns][m] = u2_cptr[ns][m] + dt*pzdot;      
        dum = 1-w3_cptr[ns][m]*nonlin_ptr[ns]*0.;

        if(ildu == 1)
        {
             dum = (tgis_ptr[ns]/ter)^1.5*exp(vfac*(1/tgis_ptr[ns]-1./ter));
        }

        //         vxdum = eyp+vpar/b*delbxp
        vxdum = (eyp/b+vpar/b*delbxp)*dum1;
        w3old = w3_cptr[ns][m];
        w3_cptr[ns][m] = w2_cptr[ns][m] + dt*(vxdum*kap+edot/ter)*dum*xnp;  

        if(abs(w3_cptr[ns][m])>1.0 && nonlin[ns]==1)
        {
           w3_cptr[ns][m] = 0.;
           w2_cptr[ns][m] = 0.;
        }     

        laps= int((z3_cptr[ns][m]/lz)-.5)*(1-peritr);
        r=x3_cptr[ns][m]-0.5*lx+lr0;
        i = int((r-rin)/dr);
        i = min(i,nr-1);
        i = max(i,0);
        wx0 = (rin+(i+1)*dr-r)/dr;
        wx1 = 1.-wx0;
        qr = wx0*sf_ptr[i]+wx1*sf_ptr[i+1];
        y3_cptr[ns][m]= int(y3_cptr[ns][m]-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly)%int(ly);

        if(x3_cptr[ns][m]>lx && iperidf==0)
        {
           x3_cptr[ns][m] = lx-1.e-8;
           z3_cptr[ns][m]=lz-z3_cptr[ns][m];
           x2_cptr[ns][m] = x3_cptr[ns][m];
           z2_cptr[ns][m] = z3_cptr[ns][m];
           w2_cptr[ns][m] = 0.;
           w3_cptr[ns][m] = 0.;
        }
        if(x3_cptr[ns][m]<0.&& iperidf==0)
        {
           x3_cptr[ns][m] = 1.e-8;
           z3_cptr[ns][m]=lz-z3_cptr[ns][m];
           x2_cptr[ns][m] = x3_cptr[ns][m];
           z2_cptr[ns][m] = z3_cptr[ns][m];
           w2_cptr[ns][m] = 0.;
           w3_cptr[ns][m] = 0.;
        }
           
        z3_cptr[ns][m]=int(z3_cptr[ns][m]+8.*lz)%int(lz);
        x3_cptr[ns][m]=int(x3_cptr[ns][m]+800.*lx)%int(lx);        
        x3_cptr[ns][m] = min(x3_cptr[ns][m],lx-1.0e-8);
        y3_cptr[ns][m] = min(y3_cptr[ns][m],ly-1.0e-8);
        z3_cptr[ns][m] = min(z3_cptr[ns][m],lz-1.0e-8);      
        //     particle diagnostics done here because info is available...
        k = int(x3_cptr[ns][m]/(lx/nsubd));
        k = min(k,nsubd-1);
        k = k+1;
        mypfl_es[k]=mypfl_es[k] + w3old*(eyp);
        mypfl_em[k]=mypfl_em[k] + w3old*(vpar*delbxp/b);
        myefl_es[k]=myefl_es[k] + vfac*w3old*(eyp);
        myefl_em[k]=myefl_em[k] + vfac*w3old*(vpar*delbxp/b);
        myke=myke + vfac*w3_cptr[ns][m];
        mynos=mynos + w3_cptr[ns][m];
        myavewi = myavewi+abs(w3_cptr[ns][m]);      
        //     xn+1 becomes xn...
        u2_cptr[ns][m]=u3_cptr[ns][m];
        x2_cptr[ns][m]=x3_cptr[ns][m];
        y2_cptr[ns][m]=y3_cptr[ns][m];
        z2_cptr[ns][m]=z3_cptr[ns][m];
        w2_cptr[ns][m]=w3_cptr[ns][m];        
        //     100     continue
     } 

   sbuf[1]=myke;
   sbuf[2]=myefl_es[nsubd/2];
   sbuf[3]=mypfl_es[nsubd/2];
   sbuf[4]=mynos;
   sbuf[5]=myavewi;
   ierr = MPI_Allreduce(sbuf,rbuf,10,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD);
   
   ketemp =rbuf[1];
   efltemp=rbuf[2];
   pfltemp=rbuf[3];
   nostemp=rbuf[4];
   avewi_cptr[ns][n] = rbuf[5]/( float(tmm_ptr[1]) );
   nos[1][n]=nostemp/( float(tmm_ptr[1]) );
   ke[1][n]=ketemp/( 2.*float(tmm_ptr[1])*mims_ptr[ns] );
   
   ierr = MPI_Barrier(MPI_COMM_WORLD);
   
   sbuf[nsubd] = myefl_es[nsubd];
   ierr = MPI_Allreduce(sbuf,rbuf,10,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD);
   for(k == 1; k < nsubd; k++)
   {
      efl_es_cptr[ns][k][n]=rbuf[k]/( float(tmm_ptr[1]) )*totvol/vol_ptr[k]*cn0s_ptr[ns];
   }
   
   sbuf[nsubd] = myefl_em[nsubd];
   ierr = MPI_Allreduce(sbuf,rbuf,10,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD);
   for(k ==1; k < nsubd; k++)
   {
      efl_em_cptr[ns][k][n]=rbuf[k]/( float(tmm_ptr[1]) )*totvol/vol_ptr[k]*cn0s_ptr[ns];
   }
   
   sbuf[nsubd] = mypfl_es[nsubd];
   ierr = MPI_Allreduce(sbuf,rbuf,10,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD);
   for(k == 1;k < nsubd; k++)
   {
      pfl_es_cptr[ns][k][n]=rbuf[k]/( float(tmm_ptr[1]) )*totvol/vol_ptr[k]*cn0s_ptr[ns];
   }
   
   sbuf[nsubd] = mypfl_em[nsubd];
   ierr = MPI_Allreduce(sbuf,rbuf,10,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD);
   for(k == 1;k < nsubd; k++)
   {
      pfl_em_cptr[ns][k][n]=rbuf[k]/( float(tmm_ptr[1]) )*totvol/vol_ptr[k]*cn0s_ptr[ns];
   }
   
   //      pfl(1,n)=pfltemp/( float(tmm(1)) )
   //      efl(1,n)=mims(ns)/tets(1)*efltemp/( float(tmm(1)) )
   
   np_old=mm_ptr[ns]; 
   ierr = MPI_Barrier(MPI_COMM_WORLD);

  /* call init_pmove(z3(ns,:),np_old,lz,ierr)
  !     
  call pmove(x2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(x3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(y2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(y3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w2(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w3(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(mu(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(xii(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z0i(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(pzi(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(eki(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(u0i(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call end_pmove(ierr)
  mm(ns)=np_new
  !     write(*,*)MyId,mm(ns)

  !      return*/
}    