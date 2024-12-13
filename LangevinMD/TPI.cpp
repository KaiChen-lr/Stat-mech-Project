#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

double calcDistance(vector<double> A, vector<double> B, double length)
{
    double r2 = 0.0;
    for (int i = 0; i < 3; i++)
    {
        double x = A[i] - B[i];
        if (x > length / 2)
            x -= length;
        else if (x < -length / 2)
            x += length;
        r2 += x * x;
    }
    double r = sqrt(r2);
    return r;
}

double calcLJEnergy(vector<double> A, vector<double> B, int iA, int iB, double length)
{
    double eps1 = 0.9977, eps2 = 0.3100; // kJ/mol
    double sig1 = 0.32, sig2 = 0.2782;   // nm
    double eps = 0.0, sig = 0.0;
    if (iA < 108 && iB < 108)
    {
        eps = eps1;
        sig = sig1;
    }
    else if (iA >= 108 && iB >= 108)
    {
        eps = eps2;
        sig = sig2;
    }
    else
    {
        eps = sqrt(eps1 * eps2);
        sig = 0.5 * (sig1 + sig2);
    }
    double sigDr = sig / calcDistance(A, B, length);
    double LJE = 4 * eps * (pow(sigDr, 12.0) - pow(sigDr, 6.0));
    return LJE;
}

double calcPE(vector<vector<double>> positions, double length)
{
    double PE = 0.0;
    for (int i = 0; i < 216; i++)
    {
        for (int j = i + 1; j < 216; j++)
            PE += calcLJEnergy(positions[i], positions[j], i, j, length);
    }
    return PE;
}

int main()
{
    vector<vector<vector<double>>> prodTrj;
    fstream trajIn1("productionTrj.txt", ios::in);
    for (int i = 0; i < 594; i++)
    {
        double N, x, y, z;
        trajIn1 >> N;
        vector<vector<double>> snapshot;
        for (int j = 0; j < 216; j++)
        {
            trajIn1 >> x >> y >> z;
            vector<double> rParticle = {x, y, z};
            snapshot.emplace_back(rParticle);
        }
        prodTrj.emplace_back(snapshot);
    }
    // Insert one Argon atom to the system
    vector<double> dEList;
    for (int i = 0; i < 594; i++)
    {
        double count = 6.0 * 6.0 * 6.0;
        double BFdE = 0.0;
        for (int x = 0; x < 6; x++)
        {
            for (int y = 0; y < 6; y++)
            {
                for (int z = 0; z < 6; z++)
                {
                    vector<double> rTest = {x + 0.5, y + 0.5, z + 0.5};
                    for (int j = 0; j < 216; j++)
                    {
                        double upperLimit = 25.89;
                        double ENow = calcLJEnergy(rTest, prodTrj[i][j], 0, j, 6.0);
                        if (ENow != ENow)
                            BFdE += 0.0;
                        else
                            BFdE += exp(-ENow * 1000.0 / (8.314 * 300.0));
                    }
                }
            }
        }
        BFdE /= count;
        dEList.emplace_back(BFdE);
    }
    fstream dEOut("Delta E of TPI.txt", ios::out);
    for (int i = 0; i < dEList.size(); i++)
        dEOut << dEList[i] << endl;
    return 0;
}