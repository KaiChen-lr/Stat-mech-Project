#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int main()
{
    fstream PEin("Potential Energy Test.txt", ios::in);
    vector<double> PETest;
    for (int i = 0; i < 50000; i++)
    {
        double PE;
        PEin >> PE;
        PETest.emplace_back(PE);
    }
    // Time correlation function
    fstream TCFout("Time correlation function.txt", ios::out);
    vector<double> CUU;
    for (int i = 0; i <= 20000; i++)
    {
        double CUUValue = 0.0;
        double count = 0.0;
        for (int j = 10000; j < 50000 - i; j += 1000)
        {
            CUUValue += PETest[j] * PETest[j + i];
            count++;
        }
        CUUValue /= count;
        CUU.emplace_back(CUUValue);
    }
    for (int i = 0; i < CUU.size(); i++)
        TCFout << CUU[i] << endl;
    return 0;
}
