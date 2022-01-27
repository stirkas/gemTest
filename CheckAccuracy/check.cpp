
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

void fileCheck(string fileName1, string fileName2)
{
    string line1;
    string line2;

    bool check = true;

    ifstream checkFile1;
    checkFile1.open(fileName1);

    ifstream checkFile2;
    checkFile2.open(fileName2);

    if(checkFile1.is_open() && checkFile2.is_open())
    {
        while(getline(checkFile1,line1) && getline(checkFile2,line2))
        {
            if(line1 != line2)
            {
                cout << "FAILED" << endl;
                cout << line1 << endl;
                cout << line2 << endl;
                check = false;
                break;
            }
        }

        if(check == true)
        {
            cout << "PASSED" << endl;
        }
    }
}

int main()
{
    cout << "run.out compare" << endl;
    fileCheck("runOriginal.out", "runPPUSHTEST.out");

    cout << "plot compare" << endl;
    fileCheck("plotOriginal", "plotPPUSH");

     cout << "flux compare" << endl;
     fileCheck("fluxOriginal", "fluxPPUSH");

     cout << "yyre compare" << endl;
     fileCheck("yyreOriginal","yyrePPUSH");
    return 0;
}