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

double Vsphere(double r)
{
    double V = (4.0 / 3.0) * M_PI * r * r * r;
    return V;
}

void calcRDF(vector<vector<vector<double>>> trj, vector<double> distance, vector<double> &gr, int N, double length)
{
    // Calculate the histogram
    for (int n = 1; n <= N; n++)
    {
        for (int i = 0; i < 108; i++)
        {
            for (int j = 108; j < 216; j++)
            {
                double r = calcDistance(trj[n][i], trj[n][j], length);
                int index = int(r / 0.01) + 1;
                gr[index] += 1.0;
            }
        }
    }
    // Normalization
    double rho = 108.0 / (6.0 * 6.0 * 6.0);
    double dr = 0.01;
    for (int i = 0; i < gr.size(); i++)
        gr[i] /= (Vsphere(distance[i]) - Vsphere(distance[i] - dr)) * rho * double(N * 108);
}

int main()
{
    vector<vector<vector<double>>> prodTrj;
    vector<vector<vector<double>>> testTrj;
    fstream trajIn1("productionTrj.txt", ios::in);
    fstream trajIn2("trj.txt", ios::in);
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
    for (int i = 0; i < 50000; i++)
    {
        double N, x, y, z;
        trajIn2 >> N;
        vector<vector<double>> snapshot;
        for (int j = 0; j < 216; j++)
        {
            trajIn2 >> x >> y >> z;
            vector<double> rParticle = {x, y, z};
            snapshot.emplace_back(rParticle);
        }
        testTrj.emplace_back(snapshot);
    }
    // Potential energy curve
    fstream PEout1("Potential Energy Prod.txt", ios::out);
    for (int i = 0; i < prodTrj.size(); i++)
    {
        double PENow = calcPE(prodTrj[i], 6.0);
        PEout1 << PENow << endl;
    }
    vector<double> PETest;
    fstream PEout2("Potential Energy Test.txt", ios::out);
    for (int i = 0; i < testTrj.size(); i++)
    {
        double PENow = calcPE(testTrj[i], 6.0);
        PETest.emplace_back(PENow);
        PEout2 << PENow << endl;
    }
    // Radial distribution function
    vector<double> r;
    vector<int> NList = {50, 100, 200, 250};
    double distance = 0.0;
    while (distance <= 2.0)
    {
        distance += 0.01;
        r.emplace_back(distance);
    }
    vector<vector<double>> grSet; // 4 g(r) with 50, 100, 200, 250 snapshots
    for (int n = 0; n < 4; n++)
    {
        vector<double> gr(200, 0.0);
        calcRDF(prodTrj, r, gr, NList[n], 6.0);
        grSet.emplace_back(gr);
    }
    fstream rdfOut("rdf.txt", ios::out);
    for (int i = 0; i < grSet[0].size(); i++)
        rdfOut << r[i] << " " << grSet[0][i] << " " << grSet[1][i] << " " << grSet[2][i] << " " << grSet[3][i] << " " << endl;

    return 0;
}