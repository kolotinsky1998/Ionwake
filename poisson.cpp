#include"poisson.h"
#include"converter.h"
#include "math.h"
#include "iostream"

using namespace std;

poisson::poisson(const converter &Converter) {


    hx = Converter.GetCoordinateStepX();
    hy = Converter.GetCoordinateStepY();
    hz = Converter.GetCoordinateStepZ();
    Nx = Converter.GetNx();
    Ny = Converter.GetNy();
    Nz = Converter.GetNz();
    n_i = Converter.GetDimensionlessIonConcentration();
    q = Converter.GetParticleCharge();
    Nx_0 = Converter.GetNx_0();
    Ny_0 = Converter.GetNy_0();
    Nz_0 = Converter.GetNz_0();

    potential = new double **[Nx];
    sfx = new double **[Nx];
    sfy = new double **[Nx];
    sfz = new double **[Nx];
    for (int i = 0; i < Nx; i++) {
        potential[i] = new double *[Ny];
        sfx[i] = new double *[Ny];
        sfy[i] = new double *[Ny];
        sfz[i] = new double *[Ny];
        for (int j = 0; j < Ny; j++) {
            potential[i][j] = new double[Nz];
            sfx[i][j] = new double[Nz];
            sfy[i][j] = new double[Nz];
            sfz[i][j] = new double[Nz];
        }
    }
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                potential[i][j][k] = 0;
                sfx[i][j][k] = 0;
                sfy[i][j][k] = 0;
                sfz[i][j][k] = 0;
            }
        }
    }

}


poisson::~poisson() {

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            delete[] potential[i][j];
            delete[] sfx[i][j];
            delete[] sfy[i][j];
            delete[] sfz[i][j];
        }
        delete[] potential[i];
        delete[] sfx[i];
        delete[] sfy[i];
        delete[] sfz[i];
    }
    delete[] potential;
    delete[] sfx;
    delete[] sfy;
    delete[] sfz;

}


void poisson::PoissonScheme(double ***density) {
    const double pi = 3.1415;
    double owx, owy, owz, tau, e;
    double owlx, owly, owlz;
    double Fi, Fj, Fk, F1;
    int f = 0;

    e = 0.0001;
    owx = pow(hx, 2);
    owy = pow(hy, 2);
    owz = pow(hz, 2);
    owlx = pow(hx * Nx, 2);
    owly = pow(hy * Ny, 2);
    owlz = pow(hz * Nz, 2);
    tau = 2. / (4. * (1. / owx + 1. / owy + 1. / owz) + pi * pi * (1. / owlx + 1. / owly + 1. / owlz));

    //cout << "*** Starting Poisson scheme ***" << endl;

    // scheme body
    do {
        f = 1;
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int k = 1; k < Nz - 1; k++) {
                    F1 = potential[i][j][k];
                    double d = 2.0 * F1;
                    Fi = (potential[i + 1][j][k] - d + potential[i - 1][j][k]) / owx;
                    Fj = (potential[i][j + 1][k] - d + potential[i][j - 1][k]) / owy;
                    Fk = (potential[i][j][k + 1] - d + potential[i][j][k - 1]) / owz;
                    potential[i][j][k] = potential[i][j][k]
                                         + (Fi + Fj + Fk + 4. * pi * (density[i][j][k] / n_i - 1. +
                                                                      q / (n_i * hx * hy * hz * Nz * Ny * Nz) -
                                                                      q * (i == Nx_0) * (j == Ny_0) * (k == Nz_0)) /
                                                           (n_i * hx * hy * hz)) * tau;
                    if (potential[i][j][k] != 0) {
                        if (abs((potential[i][j][k] - F1) / potential[i][j][k]) > e) f = 0;
                    }
                }
            }
        }

    } while (f == 0);

}


void poisson::GradientScheme() {

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                sfx[i][j][k] = -(potential[i + 1][j][k] - potential[i - 1][j][k]) / (2.0 * hx);
                sfy[i][j][k] = -(potential[i][j + 1][k] - potential[i][j - 1][k]) / (2.0 * hy);
                sfz[i][j][k] = -(potential[i][j][k + 1] - potential[i][j][k - 1]) / (2.0 * hz);
                sfx[i][j][k] = sfx[i][j][k] * n_i;
                sfy[i][j][k] = sfy[i][j][k] * n_i;
                sfz[i][j][k] = sfz[i][j][k] * n_i;
            }
        }
    }
}

double ***poisson::GetSelfConsistentForceFieldX() const {
    return sfx;
}

double ***poisson::GetSelfConsistentForceFieldY() const {
    return sfy;
}

double ***poisson::GetSelfConsistentForceFieldZ() const {
    return sfz;
}

double ***poisson::GetPotential() const {
    return potential;
}
