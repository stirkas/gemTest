#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

#include "gem_com_c.h"
#include "gem_equil_c.h"

using namespace std;

extern "C"
{
    void spec_c_(int& n);
};

void spec_c_(int& n)
{
    int i,j,k,l,m;//,n;
    double pf, efe,efi,pfi,efc,pfc,efb,pfb,pflxgb,eflxgb,x;
    double pf_em, efe_em, efi_em, pfi_em,efc_em,pfc_em, efb_em, pfb_em;
    double tdum;
    double v[imx];

    eflxgb = xn0e_ptr[nr2]*cn0e*t0e_ptr[nr2]*sqrt(t0e_ptr[nr2]/mimp)*pow(rhoia,2);
    pflxgb = eflxgb/t0e_ptr[nr2];
    i = tcurr-dt;
    pf = 0.;
    efe = 0.;
    efi = 0.;
    pfi = 0.;
    efc = 0.;
    pfc = 0.;
    efb = 0.;
    pfb = 0.;
    pf_em = 0.;
    efe_em = 0.;
    efi_em = 0.;
    pfi_em = 0.;
    efc_em = 0.;
    pfc_em = 0.;
    efb_em = 0.;
    pfb_em = 0.;
    k = 2;
    x = float(nsubd)/float(nsubd-2*k);

    j=1+k;
    for(j;j < nsubd-k;j++)
    {
        pf = pf+pfle_es_cptr[n][j]*vol_ptr[j]/(totvol*x);
        efe = efe+efle_es_cptr[n][j]*vol_ptr[j]/(totvol*x);
        efi = efi+efl_es_cptr[n][j][1]*vol_ptr[j]/(totvol*x);
        pfi = pfi+pfl_es_cptr[n][j][1]*vol_ptr[j]/(totvol*x);
        efc = efc+efl_es_cptr[n][j][2]*vol_ptr[j]/(totvol*x);
        pfc = pfc+pfl_es_cptr[n][j][2]*vol_ptr[j]/(totvol*x);
        efb = efb+efl_es_cptr[n][j][3]*vol_ptr[j]/(totvol*x);
        pfb = pfb+pfl_es_cptr[n][j][3]*vol_ptr[j]/(totvol*x);
    }

    j = 1+k;
    for(j;j < nsubd-k;j++)
    {
        pf_em = pf_em+pfle_em_cptr[n][j]*vol_ptr[j]/(totvol*x);
        efe_em = efe_em+efle_em_cptr[n][j]*vol_ptr[j]/(totvol*x);
        efi_em = efi_em+efl_em_cptr[n][j][1]*vol_ptr[j]/(totvol*x);
        pfi_em = pfi_em+pfl_em_cptr[n][j][1]*vol_ptr[j]/(totvol*x);
        efc_em = efc_em+efl_em_cptr[n][j][2]*vol_ptr[j]/(totvol*x);
        pfc_em = pfc_em+pfl_em_cptr[n][j][2]*vol_ptr[j]/(totvol*x);
        efb_em = efb_em+efl_em_cptr[n][j][3]*vol_ptr[j]/(totvol*x);
        pfb_em = pfb_em+pfl_em_cptr[n][j][3]*vol_ptr[j]/(totvol*x);
    }

    tdum = tcurr-dt;
    if(myid == master)
    {
        //i/o file things that I need to figure out in fortran
        fstream plotFile;
        plotFile.open("plot",ios::app);

        fstream fluxFile;
        fluxFile.open("flux",ios::app);

        fstream yyreFile;
        yyreFile.open("yyre",ios::app);

        cout << i << scientific << setprecision(12) << rmsphi_ptr[n] << "   " << rmsapa_ptr[n] << "   " << pf << "   " << efe << "   " << pfi << "   " << efi;
        cout << " " << avewi_cptr[n][1] << "   " << avewe_ptr[n] << "   " << yyre_cptr[0][1];
        cout << " " << yyim_cptr[0][1] << "   " << yyamp_cptr[0][1] << endl;

        if(plotFile.is_open())
        {
            plotFile << i << scientific << setprecision(12) << "   " << rmsphi_ptr[n] << "   " << rmsapa_ptr[n] << "   " << pf << "   " << efe << "   " << pfi << "   " << efi;
            plotFile <<  "   " << avewi_cptr[n][1] << "   " << avewe_ptr[n] << "   " << yyre_cptr[0][1];
            plotFile << "   " << yyim_cptr[0][1] << "   " << yyamp_cptr[0][1]  << endl;  
        }
        
        if(fluxFile.is_open())
        {
            fluxFile << i << scientific << setprecision(12) << "   " << pf / pflxgb << "   " << pfi / pflxgb << "   " << pfc / pflxgb<< "   " ;
            fluxFile << efe / eflxgb << "   " << efi / eflxgb << "   "  << efc / eflxgb << "   " ;
            fluxFile << pf_em / pflxgb << "   "  << pfi_em / pflxgb << "   " << pfc_em / pflxgb << "   " ;
            fluxFile << efe_em / eflxgb << "   "  << efi_em / eflxgb << "   "  << efc_em / eflxgb << endl;
        }
        
        if(yyreFile.is_open())
        {
            yyreFile << i << scientific << setprecision(12) << "   " << yyre_cptr[0][1] << "   " << yyre_cptr[1][1] << "   " << yyre_cptr[3][1] << "   " << yyre_cptr[4][1] << endl;
        }
        

        plotFile.close();
        fluxFile.close();
        yyreFile.close();
    }

  /* This is a copy of code that is in GEM but is not utilized, thus it is saved here 
  for furture reference but it is not converted
  
    return
 if(gclr==kmx/2 .and. tclr==0)then
    open(22, file='yyre2', status='unknown',position='append')
    open(23, file='mdhis', status='unknown',position='append')
    open(24, file='mdhisa', status='unknown',position='append')
    open(25, file='stress', status='unknown',position='append')
    write(23,16)tdum,mdhis(0),mdhis(1),mdhis(2),mdhis(3),mdhis(4),&
         mdhisa(0),mdhisa(1),mdhisa(2),mdhisa(3),mdhisa(4)
    write(24,16)tdum,mdhisb(0),mdhisb(1),mdhisc(0),mdhisc(1),&
         mdhisd(0),mdhisd(1)
    write(25,17)tdum,(v(i),i = 0,imx-1)
    do  i = 0,6
       write(22,14)tdum,i,real(phihis(i,0)),(real(phihis(i,j)), &
            aimag(phihis(i,j)), j = 1,jcnt-2,2)
       write(22,14)tdum,i,real(aparhis(i,0)),(real(aparhis(i,j)), &
            aimag(aparhis(i,j)), j = 1,jcnt-2,2)
    end do
    close(22)
    close(23)
    close(24)
    close(25)
14   format(1x,f10.1,1x,i2,10(2x,e12.5))
16   format(1x,f10.1,1x,10(2x,e12.5))
17   format(1x,f10.1,256(2x,e12.5))
 end if*/

}