#pragma once
#include <mpi.h>

extern "C"
{
    void spec_c_(int& n);

    void reporter_c_(int& n);

    void restart_(int iflag, int& n);

    void outd_(int& n);

    void ppush_c_(int& n, int& ns);

    //pputil functions

    void init_pmove_(double const &xp, int const &np, double const &lz, int &ierr);

    void pmove_(double * xp, int &np_old, int &np_new, int &ierr);

    void ppexit_();

    void end_pmove_(int &ierr);
};