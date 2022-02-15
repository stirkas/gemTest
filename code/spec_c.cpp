#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>

#include "gem_com_c.h"
#include "gem_equil_c.h"
#include "gem_functions_c.h"

using namespace std;

void spec_c_(int& n)
{
    int i,k,l,m;//,n;
    double pf, efe,efi,pfi,efc,pfc,efb,pfb,pflxgb,eflxgb,x;
    double pf_em, efe_em, efi_em, pfi_em,efc_em,pfc_em, efb_em, pfb_em;
    double tdum;
    double v[imx];


    eflxgb = xn0e_ptr[nr2]*cn0e*t0e_ptr[nr2]*sqrt(t0e_ptr[nr2]/mimp)*pow(rhoia,2);
    pflxgb = eflxgb/t0e_ptr[nr2];
    i = tcurr-dt;
    pf = 0.0;
    efe = 0.0;
    efi = 0.0;
    pfi = 0.0;
    efc = 0.0;
    pfc = 0.0;
    efb = 0.0;
    pfb = 0.0;
    pf_em = 0.0;
    efe_em = 0.0;
    efi_em = 0.0;
    pfi_em = 0.0;
    efc_em = 0.0;
    pfc_em = 0.0;
    efb_em = 0.0;
    pfb_em = 0.0;
    k = 2;
    x = float(nsubd)/float(nsubd-2*k);

    for(int j = 1+k;j < nsubd-k;j++)
    {
        pf  += pfle_es_cptr[n][j]*vol_ptr[j-1]/totvol*x;
        efe += efle_es_cptr[n][j]*vol_ptr[j]/totvol*x;
        efi += efl_es_cptr[0][j][n]*vol_ptr[j]/totvol*x;
        pfi += pfl_es_cptr[n][j][0]*vol_ptr[j]/totvol*x;
        efc += efl_es_cptr[n][j][1]*vol_ptr[j]/totvol*x;
        pfc += pfl_es_cptr[n][j][1]*vol_ptr[j]/totvol*x;
        efb += efl_es_cptr[n][j][2]*vol_ptr[j]/totvol*x;
        pfb += pfl_es_cptr[n][j][2]*vol_ptr[j]/totvol*x;
    }

    for(int j=k;j < nsubd-k;j++)
    {
        pf_em += pfle_em_cptr[n][j]*vol_ptr[j]/totvol*x;
        efe_em += efle_em_cptr[n][j]*vol_ptr[j]/totvol*x;
        efi_em += efl_em_cptr[n][j][0]*vol_ptr[j]/totvol*x;
        pfi_em += pfl_em_cptr[n][j][0]*vol_ptr[j]/totvol*x;
        efc_em += efl_em_cptr[n][j][1]*vol_ptr[j]/totvol*x;
        pfc_em += pfl_em_cptr[n][j][1]*vol_ptr[j]/totvol*x;
        efb_em += efl_em_cptr[n][j][2]*vol_ptr[j]/totvol*x;
        pfb_em += pfl_em_cptr[n][j][2]*vol_ptr[j]/totvol*x;
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

        cout << scientific << setprecision(10) << setfill(' ') << uppercase;
        cout << setw(7) << i;
        cout << setw(19) << rmsphi_ptr[n] << setw(19) << rmsapa_ptr[n] << setw(19) << pf << setw(19) << efe << setw(19) << pfi << setw(19) << efi;
        cout << setw(19) << avewi_cptr[n][0] << setw(19) << avewe_ptr[n] << setw(19) << yyre_cptr[0][0] << setw(19) << yyim_cptr[0][0] << setw(19) << yyamp_cptr[0][0];
        cout << endl;

        if(plotFile.is_open())
        {
            plotFile << scientific << setprecision(3) << setfill(' ') << uppercase;
            plotFile << setw(7) << i;
            plotFile << setw(12) << rmsphi_ptr[n] << setw(12) << rmsapa_ptr[n] << setw(12) << pf << setw(12) << efe << setw(12) << pfi << setw(12) << efi;
            plotFile << setw(12) << avewi_cptr[n][0] << setw(12) << avewe_ptr[n] << setw(12) << yyre_cptr[0][0] << setw(12) << yyim_cptr[0][0] << setw(12) << yyamp_cptr[0][0];
            plotFile << endl;  
        }
        
        if(fluxFile.is_open())
        {
            fluxFile << scientific << setprecision(5) << setfill(' ') << uppercase;
            fluxFile << setw(7) << i;
            fluxFile << setw(14) << pf / pflxgb << setw(14) << pfi / pflxgb << setw(14) << pfc / pflxgb;
            fluxFile << setw(14) << efe / eflxgb << setw(14) << efi / eflxgb << setw(14) << efc / eflxgb;
            fluxFile << setw(14) << pf_em / pflxgb << setw(14) << pfi_em / pflxgb << setw(14) << pfc_em / pflxgb;
            fluxFile << setw(14) << efe_em / eflxgb << setw(14) << efi_em / eflxgb << setw(14) << efc_em / eflxgb;
            fluxFile << endl;
        }
        
        if(yyreFile.is_open())
        {
            yyreFile << scientific << setprecision(5) << setfill(' ') << uppercase;
            yyreFile << "what the fuck" << setw(7) << i;
            yyreFile << setw(14) << yyre_cptr[0][0] << setw(14) << yyre_cptr[1][0]<< setw(14) << yyre_cptr[2][0] << setw(14) << yyre_cptr[3][0] << setw(14) << yyre_cptr[4][0] << "why";
            yyreFile << endl;
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