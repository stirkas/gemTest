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

using namespace std;

void reporter_c_(int& n)
{
    int i,j,k,ip;

    if(n % xnplt == 0)
    {
        spec_c_(n);
    }

    if(myid == master)
    {
        fstream indicatorFile;
        indicatorFile.open("indicator", ios::app);

        indicatorFile << setfill(' ');
        indicatorFile << setw(7) << n;
        indicatorFile << setw(9) << ipred << setw(9) << icorr << setw(9) << jpred << setw(9) << jcorr;
        indicatorFile << setw(9) << nopz << setw(9) << noen << setw(9) << nowe << endl;

        indicatorFile.close();
    }

    ierr = MPI_Barrier(MPI_COMM_WORLD);

    if(iput == 1 && n+1 % 500 == 0)
    {
        restart_(2,n);
    }
    outd_(n);
    if(idg == 1)
    {
        cout << "pass outd" << endl;
    }
}

