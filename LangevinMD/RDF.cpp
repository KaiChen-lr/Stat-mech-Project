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

double Vsphere(double r)
{
    double V = (4.0 / 3.0) * M_PI * r * r * r;
    return V;
}

void calcRDF(vector<vector<vector<double>>> trj, vector<double> distance, vector<double> &gr, int N, double length)
{
    // Calculate the histogram
    for (int n = 50; n < N + 50; n++)
    {
        for (int i = 0; i < 216; i++)
        {
            for (int j = i + 1; j < 216; j++)
            {
                double r = calcDistance(trj[n][i], trj[n][j], length);
                int index = (int)(r / 0.02) + 1;
                if (index < 100)
                    gr[index] += 1.0;
            }
        }
    }
    // Normalization
    double rho = 216.0 / (6.0 * 6.0 * 6.0);
    double dr = 0.02;
    for (int i = 0; i < gr.size(); i++)
        gr[i] /= (Vsphere(distance[i]) - Vsphere(distance[i] - dr)) * rho * (N * 108.0);
}

void calcRDF_Ar_Ar(vector<vector<vector<double>>> trj, vector<double> distance, vector<double> &gr, int N, double length)
{
    // Calculate the histogram
    for (int n = 50; n < N + 50; n++)
    {
        for (int i = 0; i < 108; i++)
        {
            for (int j = i + 1; j < 108; j++)
            {
                double r = calcDistance(trj[n][i], trj[n][j], length);
                int index = (int)(r / 0.02) + 1;
                if (index < 100)
                    gr[index] += 1.0;
            }
        }
    }
    // Normalization
    double rho = 108.0 / (6.0 * 6.0 * 6.0);
    double dr = 0.02;
    for (int i = 0; i < gr.size(); i++)
        gr[i] /= (Vsphere(distance[i]) - Vsphere(distance[i] - dr)) * rho * (N * 54.0);
}

void calcRDF_Ne_Ne(vector<vector<vector<double>>> trj, vector<double> distance, vector<double> &gr, int N, double length)
{
    // Calculate the histogram
    for (int n = 50; n < N + 50; n++)
    {
        for (int i = 108; i < 216; i++)
        {
            for (int j = i + 1; j < 216; j++)
            {
                double r = calcDistance(trj[n][i], trj[n][j], length);
                int index = (int)(r / 0.02) + 1;
                if (index < 100)
                    gr[index] += 1.0;
            }
        }
    }
    // Normalization
    double rho = 108.0 / (6.0 * 6.0 * 6.0);
    double dr = 0.02;
    for (int i = 0; i < gr.size(); i++)
        gr[i] /= (Vsphere(distance[i]) - Vsphere(distance[i] - dr)) * rho * (N * 54.0);
}

void calcRDF_Ar_Ne(vector<vector<vector<double>>> trj, vector<double> distance, vector<double> &gr, int N, double length)
{
    // Calculate the histogram
    for (int n = 50; n < N + 50; n++)
    {
        for (int i = 0; i < 108; i++)
        {
            for (int j = 108; j < 216; j++)
            {
                double r = calcDistance(trj[n][i], trj[n][j], length);
                int index = (int)(r / 0.02) + 1;
                if (index < 100)
                    gr[index] += 1.0;
            }
        }
    }
    // Normalization
    double rho = 108.0 / (6.0 * 6.0 * 6.0);
    double dr = 0.02;
    for (int i = 0; i < gr.size(); i++)
        gr[i] /= (Vsphere(distance[i]) - Vsphere(distance[i] - dr)) * rho * (N * 108.0);
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
    // Radial distribution function
    vector<double> r;
    vector<int> NList = {50, 100, 200, 500};
    double distance = 0.0;
    while (distance <= 2.0)
    {
        distance += 0.02;
        r.emplace_back(distance);
    }
    vector<vector<double>> grSet; // 4 g(r) with 50, 100, 200, 250 snapshots
    for (int n = 0; n < 4; n++)
    {
        vector<double> gr(100, 0.0);
        calcRDF(prodTrj, r, gr, NList[n], 6.0);
        grSet.emplace_back(gr);
    }
    fstream rdfOut("rdf different nSnapshots.txt", ios::out);
    for (int i = 0; i < grSet[0].size(); i++)
        rdfOut << r[i] << " " << grSet[0][i] << " " << grSet[1][i] << " " << grSet[2][i] << " " << grSet[3][i] << " " << endl;

    vector<double> grNe_Ne(100, 0.0);
    calcRDF_Ne_Ne(prodTrj, r, grNe_Ne, 500, 6.0);
    vector<double> grAr_Ar(100, 0.0);
    calcRDF_Ar_Ar(prodTrj, r, grAr_Ar, 500, 6.0);
    vector<double> grAr_Ne(100, 0.0);
    calcRDF_Ar_Ne(prodTrj, r, grAr_Ne, 500, 6.0);
    fstream rdfOut2("rdf different types.txt", ios::out);
    for (int i = 0; i < grNe_Ne.size(); i++)
        rdfOut2 << grNe_Ne[i] << " " << grAr_Ar[i] << " " << grAr_Ne[i] << endl;
    return 0;
}