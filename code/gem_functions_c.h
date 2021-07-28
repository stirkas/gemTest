#pragma once

extern "C"
{
    void spec_c_(int& n);

    void reporter_c_(int& n);

    void restart_(int iflag, int& n);

    void outd_(int& n);
};