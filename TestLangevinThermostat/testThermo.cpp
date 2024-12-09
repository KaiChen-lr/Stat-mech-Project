#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

double kB = 1.380649e-23;
double NA = 6.02214e23;
class particle
{
public:
    vector<double> p; // Momentum, kg*nm/fs
    vector<double> q; // Position, nm
    double m;         // Mass, kg
    double epsilon;   // Two parameters for LJ potential
    double sigma;
};

class vec
{
public:
    vector<double> vecComp;
    double vecMod;
};

// Function calculate the distance between A and B
void calcDistance(particle A, particle B, vec &rAB, double length)
{
    double r2 = 0.0;
    for (int i = 0; i < 3; i++)
    {
        double x = A.q[i] - B.q[i];
        if (x > length / 2)
            x -= length;
        else if (x < -length / 2)
            x += length;
        rAB.vecComp[i] = x;
        r2 += x * x;
    }
    rAB.vecMod = sqrt(r2);
}

// Function calculate the force on A because of B
void fLJ(vec rAB, vec &fA, vec &fB, particle A, particle B, double length)
{
    double eps = sqrt(A.epsilon * B.epsilon); // kg*nm^2/fs^2
    double sig = 0.5 * (A.sigma + B.sigma);   // nm
    double sr = sig / rAB.vecMod;
    double f = 0.0;
    if (rAB.vecMod < length / 2)
        f = (4.0 * eps / rAB.vecMod) * (12.0 * pow(sr, 12.0) - 6.0 * pow(sr, 6.0));
    for (int i = 0; i < 3; i++)
    {
        double fi = f * rAB.vecComp[i] / rAB.vecMod;
        fA.vecComp[i] += fi;
        fB.vecComp[i] -= fi;
    }
}

// Function to calculate the force for the whole system
void updateLJForce(vector<vec> &forces, vector<vector<vec>> rMatrix, vector<particle> sys, double length)
{
    for (int i = 0; i < forces.size(); i++)
    {
        forces[i].vecComp = {0.0, 0.0, 0.0};
        forces[i].vecMod = 0.0;
    }
    for (int i = 0; i < forces.size(); i++)
    {
        for (int j = i + 1; j < forces.size(); j++)
            fLJ(rMatrix[i][j], forces[i], forces[j], sys[i], sys[j], length);
    }
    for (int i = 0; i < forces.size(); i++)
        forces[i].vecMod = sqrt(forces[i].vecComp[0] * forces[i].vecComp[0] + forces[i].vecComp[1] * forces[i].vecComp[1] + forces[i].vecComp[2] * forces[i].vecComp[2]);
}

// Function to calculate the r matrix
void rMatrixUpdater(vector<particle> sys, vector<vector<vec>> &rMatrix, double length)
{
    for (int i = 0; i < sys.size(); i++)
    {
        for (int j = i + 1; j < sys.size(); j++)
        {
            calcDistance(sys[i], sys[j], rMatrix[i][j], length);
            for (int k = 0; k < 3; k++)
                rMatrix[j][i].vecComp[k] = -rMatrix[i][j].vecComp[k];
            rMatrix[j][i].vecMod = rMatrix[i][j].vecMod;
        }
    }
}

// Function to calculate the temperature using Virial theorem
double calcT(vector<particle> sys, double N)
{
    double EK = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            EK += sys[i].p[j] * sys[i].p[j] * 1.0e12;
        }
    }
    EK /= 2.0 * sys[0].m;
    double T = EK / (1.5 * N * kB);
    return T;
}

int main()
{
    vector<particle> system;
    // Construct the simulation cell
    double N = 216;
    double L = 6.0;
    // Add 108 Argon atoms and then 108 Neon atoms
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {
                particle temp;
                temp.q = {i * 1.0, j * 1.0, k * 1.0};
                temp.p = {0.0, 0.0, 0.0};
                if (i < 3)
                {
                    temp.m = 0.039948 / (6.02214e23);
                    temp.epsilon = 1.66e-33; // kg*nm^2/fs^2
                    temp.sigma = 0.34;       // nm
                }
                else if (i >= 3)
                {
                    temp.m = 0.0201797 / (6.02214e23);
                    temp.epsilon = 5.148605430000001e-34; // kg*nm^2/fs^2
                    temp.sigma = 0.2782;
                }
                system.emplace_back(temp);
            }
        }
    }

    // Create containers for the interatomic distances and particle forces
    vector<vector<vec>> rMatrix; // rMatrix[i][j] is the distance vector between particle i and j, rMatrix[i][i] will be a zero vector
    vector<vec> force;           // force[i] is the force vector for particle i
    for (int i = 0; i < N; i++)
    {
        vector<vec> temp;
        for (int j = 0; j < N; j++)
        {
            vec tempVec;
            tempVec.vecComp = {0.0, 0.0, 0.0};
            tempVec.vecMod = 0.0;
            temp.emplace_back(tempVec);
        }
        force.emplace_back(temp[0]);
        rMatrix.emplace_back(temp);
    }
    rMatrixUpdater(system, rMatrix, L);
    updateLJForce(force, rMatrix, system, L);

    // Start simulation
    fstream trjOut("trj.txt", ios::out);
    fstream TOut("Temperature.txt", ios::out);
    double dt = 1.0;      // fs
    double gamma = 0.001; // fs^-1
    double TB = 300.0;    // Bath temperature
    default_random_engine generator((unsigned int)(time(0)));
    normal_distribution<double> distribution(0.0, 1.0);

    for (int step = 1; step <= 10000; step++)
    {
        if (step % 100 == 0 || step == 1)
        {
            // Calculate the temperature and print to the file
            double T = calcT(system, N);
            TOut << step << " " << T << endl;
            // Print positions to the trajectory file
            trjOut << step << endl;
            for (int i = 0; i < N; i++)
                trjOut << system[i].q[0] << " " << system[i].q[1] << " " << system[i].q[2] << endl;
        }
        // Numerically integrate Langevin's EOM
        // Update the momentum to t+dt/2
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
                system[i].p[j] += 0.5 * dt * force[i].vecComp[j];
        }
        // Update the positions to t+dt/2
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                system[i].q[j] += 0.5 * dt * system[i].p[j] / system[i].m;
                if (system[i].q[j] < 0)
                    system[i].q[j] += L;
                else if (system[i].q[j] > L)
                    system[i].q[j] -= L;
            }
        }
        // Update momentum with friction and random force
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                // Gaussian random numbers
                double r = distribution(generator);
                double pLast = system[i].p[j];
                system[i].p[j] = pLast * exp(-gamma * dt) + 1.0e-6 * r * sqrt(kB * TB * system[i].m) * sqrt(1.0 - exp(-2 * gamma * dt));
            }
        }
        // Update positions to t+dt
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
                system[i].q[j] += 0.5 * dt * system[i].p[j] / system[i].m;
        }
        // Update forces and r matrix
        rMatrixUpdater(system, rMatrix, L);
        updateLJForce(force, rMatrix, system, L);
        // Update momentum to t+dt
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
                system[i].p[j] += 0.5 * dt * force[i].vecComp[j];
        }
    }
    return 0;
}