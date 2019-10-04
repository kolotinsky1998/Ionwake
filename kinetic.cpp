#include <cmath>
#include <iostream>

#include "constants.h"
#include "kinetic.h"
#include "poisson.h"
#include "omp.h"

using namespace std;

kinetic::kinetic(const converter &Converter_, const poisson &Poisson) : Converter(&Converter_) {

    hx = Converter->GetCoordinateStepX();
    hy = Converter->GetCoordinateStepY();
    hz = Converter->GetCoordinateStepZ();

    nx = Converter->GetNx();
    ny = Converter->GetNy();
    nz = Converter->GetNz();

    nx_0 = Converter->GetNx_0();
    ny_0 = Converter->GetNy_0();
    nz_0 = Converter->GetNz_0();

    Lx = Converter->GetLx();
    Ly = Converter->GetLy();
    Lz = Converter->GetLz();

    hvx = 0.2;
    hvy = 0.2;
    hvz = 0.2;

    vxcut = max(5., 8 * Converter->GetDimensionlessFlowVelocity());
    vyzcut = 5;
    nvx = int((vxcut + vyzcut) / hvx);
    nvy = int(2 * vyzcut / hvy);
    nvz = int(2 * vyzcut / hvz);

    nxi = int(MAX_XI / DELTA_XI);

    al = Converter->GetLAxelerationCofficient();

    deltaT = 0.5 * min(0.2 / (al / hvx), 1.0 / (vxcut / hx + vyzcut / hy + vyzcut / hz));
    CurT = 0;
    tau = Converter->GetDimensionlessTau();

    f = new double *****[nx];
    for (int i = 0; i < nx; ++i) {
        f[i] = new double ****[ny];
        for (int j = 0; j < ny; ++j) {
            f[i][j] = new double ***[nz];
            for (int k = 0; k < nz; ++k) {
                f[i][j][k] = new double **[nvx];
                for (int a = 0; a < nvx; ++a) {
                    f[i][j][k][a] = new double *[nvy];
                    for (int b = 0; b < nvy; ++b) {
                        f[i][j][k][a][b] = new double[nvz];
                    }
                }
            }
        }
    }

    f_time = new double *****[nx];
    for (int i = 0; i < nx; ++i) {
        f_time[i] = new double ****[ny];
        for (int j = 0; j < ny; ++j) {
            f_time[i][j] = new double ***[nz];
            for (int k = 0; k < nz; ++k) {
                f_time[i][j][k] = new double **[nvx];
                for (int a = 0; a < nvx; ++a) {
                    f_time[i][j][k][a] = new double *[nvy];
                    for (int b = 0; b < nvy; ++b) {
                        f_time[i][j][k][a][b] = new double[nvz];
                    }
                }
            }
        }
    }

    vx = new double[nvx];
    vy = new double[nvy];
    vz = new double[nvz];

    for (int a = 0; a < nvx; ++a) {
        vx[a] = -vyzcut + hvx * a;
    }

    for (int b = 0; b < nvy; ++b) {
        vy[b] = -vyzcut + hvy * b;
    }

    for (int c = 0; c < nvz; ++c) {
        vz[c] = -vyzcut + hvz * c;
    }

    x = new double[nx];
    y = new double[ny];
    z = new double[nz];

    for (int i = 0; i < nx; ++i) {
        x[i] = hx * i;
    }

    for (int j = 0; j < ny; ++j) {
        y[j] = hy * j;
    }

    for (int k = 0; k < nz; ++k) {
        z[k] = hz * k;
    }

    f_n = new double **[nvx];
    for (int a = 0; a < nvx; ++a) {
        f_n[a] = new double *[nvy];
        for (int b = 0; b < nvy; ++b) {
            f_n[a][b] = new double[nvz];
        }
    }

    for (int a = 0; a < nvx; ++a) {
        for (int b = 0; b < nvy; ++b) {
            for (int c = 0; c < nvz; ++c) {
                f_n[a][b][c] = exp(-0.5 * (vx[a] * vx[a] + vy[b] * vy[b] + vz[c] * vz[c])) / pow(2 * M_PI, 1.5);
            }
        }
    }

    InitialConditions();

    profile_x = new double[nvx];
    profile_y = new double[nvy];
    profile_z = new double[nvz];
    analytical_profile_x = new double[nvx];

    for (int a = 0; a < nvx; ++a) {
        profile_x[a] = 0;
    }

    for (int b = 0; b < nvy; ++b) {
        profile_y[b] = 0;
    }

    for (int c = 0; c < nvz; ++c) {
        profile_z[c] = 0;
    }

    vfl_x = new double **[nx];
    vfl_y = new double **[nx];
    vfl_z = new double **[nx];
    acx = new double **[nx];
    acy = new double **[nx];
    acz = new double **[nx];
    density = new double **[nx];
    for (int i = 0; i < nx; ++i) {
        vfl_x[i] = new double *[ny];
        vfl_y[i] = new double *[ny];
        vfl_z[i] = new double *[ny];
        acx[i] = new double *[ny];
        acy[i] = new double *[ny];
        acz[i] = new double *[ny];
        density[i] = new double *[ny];
        for (int j = 0; j < ny; ++j) {
            vfl_x[i][j] = new double[nz];
            vfl_y[i][j] = new double[nz];
            vfl_z[i][j] = new double[nz];
            acx[i][j] = new double[nz];
            acy[i][j] = new double[nz];
            acz[i][j] = new double[nz];
            density[i][j] = new double[nz];
        }
    }


    asx = Poisson.GetSelfConsistentForceFieldX();
    asy = Poisson.GetSelfConsistentForceFieldY();
    asz = Poisson.GetSelfConsistentForceFieldZ();

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                density[i][j][k] = Converter->GetDimensionlessIonConcentration();

                acx[i][j][k] = -Converter->GetCAxelerationCofficient() * (x[i] - x[nx_0] + 0.5 * hx)
                                                                         / pow(pow((x[i] - x[nx_0] + 0.5 * hx), 2.) + pow((y[j] - y[ny_0] + 0.5 * hy), 2.) +
                                                                                       pow((z[k] - z[nz_0] + 0.5 * hz), 2.), 1.5);

                acy[i][j][k] = -Converter->GetCAxelerationCofficient() * (y[j] - y[ny_0] + 0.5 * hy)
                                                                         / pow(pow((x[i] - x[nx_0] + 0.5 * hx), 2.) + pow((y[j] - y[ny_0] + 0.5 * hy), 2.) +
                                                                                       pow((z[k] - z[nz_0] + 0.5 * hz), 2.), 1.5);

                acz[i][j][k] = -Converter->GetCAxelerationCofficient() * (z[k] - z[nz_0] + 0.5 * hz)
                               / pow(pow((x[i] - x[nx_0] + 0.5 * hx), 2.) + pow((y[j] - y[ny_0] + 0.5 * hy), 2.) +
                                     pow((z[k] - z[nz_0] + 0.5 * hz), 2.), 1.5);

                vfl_x[i][j][k] = 0;
                vfl_y[i][j][k] = 0;
                vfl_z[i][j][k] = 0;
            }
        }
    }

    v_fl_an = al * tau;
    SetAnalyticalProfileX();
}

kinetic::~kinetic() {

    delete[] vx;
    delete[] vy;
    delete[] vz;

    delete[] x;
    delete[] y;
    delete[] z;

    delete[] profile_x;
    delete[] profile_y;
    delete[] profile_z;
    delete[] analytical_profile_x;


    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            delete[] density[i][j];
            delete[] vfl_x[i][j];
            delete[] vfl_y[i][j];
            delete[] vfl_z[i][j];
            delete[] acx[i][j];
            delete[] acy[i][j];
            delete[] acz[i][j];
        }
        delete[] density[i];
        delete[] vfl_x[i];
        delete[] vfl_y[i];
        delete[] vfl_z[i];
        delete[] acx[i];
        delete[] acy[i];
        delete[] acz[i];
    }

    delete[] density;
    delete[] vfl_x;
    delete[] vfl_y;
    delete[] vfl_z;
    delete[] acx;
    delete[] acy;
    delete[] acz;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int a = 0; a < nvx; ++a) {
                    for (int b = 0; b < nvy; ++b) {
                        delete[] f[i][j][k][a][b];
                    }
                    delete[] f[i][j][k][a];
                }
                delete[] f[i][j][k];
            }
            delete[] f[i][j];
        }
        delete[] f[i];
    }

    delete[] f;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int a = 0; a < nvx; ++a) {
                    for (int b = 0; b < nvy; ++b) {
                        delete[] f_time[i][j][k][a][b];
                    }
                    delete[] f_time[i][j][k][a];
                }
                delete[] f_time[i][j][k];
            }
            delete[] f_time[i][j];
        }
        delete[] f_time[i];
    }

    delete[] f_time;


    for (int a = 0; a < nvx; ++a) {
        for (int b = 0; b < nvy; ++b) {
            delete[] f_n[a][b];
        }
        delete[] f_n[a];
    }

    delete[] f_n;


}


void kinetic::CoordinatePart() {
    const double deltaThx = deltaT / hx;
    const double deltaThy = deltaT / hy;
    const double deltaThz = deltaT / hz;

    //cout << "***  Starting coordinate part ***" << endl;
    SaveState();
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int a = 1; a < nvx - 1; ++a) {
                    double **f_ijka = f[i][j][k][a];
                    double **f_time_ijka = f_time[i][j][k][a];
                    if (vx[a] > 0) {
                        if (abs(vx[a] * deltaThx) > 1) cout << "x" << vx[a] * deltaThx << endl;
                        for (int b = 1; b < nvy - 1; ++b) {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijka[b][c] = f_time_ijka[b][c] - (deltaThx * vx[a])
                                                                   * (f_time_ijka[b][c] -
                                                                      f_time[(nx + i - 1) %
                                                                             nx][j][k][a][b][c]);
                            }
                        }
                    } else {
                        for (int b = 1; b < nvy - 1; ++b) {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijka[b][c] = f_time_ijka[b][c] - (deltaThx * vx[a]) *
                                                                   (f_time[(i + 1) % nx][j][k][a][b][c] -
                                                                    f_time_ijka[b][c]);
                            }
                        }
                    }
                }
            }
        }
    }

    SaveState();
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int a = 1; a < nvx - 1; ++a) {
                    for (int b = 1; b < nvy - 1; ++b) {
                        double *f_ijkab = f[i][j][k][a][b];
                        double *f_time_ijkab = f_time[i][j][k][a][b];
                        if (vy[b] > 0) {
                            if (abs(vy[b] * deltaThy) > 1) cout << "y" << vy[b] * deltaThy << endl;
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijkab[c] = f_time_ijkab[c] - (deltaThy * vy[b]) *
                                                               (f_time_ijkab[c] -
                                                                f_time[i][(ny + j - 1) % ny][k][a][b][c]);
                            }
                        } else {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijkab[c] = f_time_ijkab[c] - (deltaThy * vy[b]) *
                                                               (f_time[i][(j + 1) % ny][k][a][b][c] - f_time_ijkab[c]);
                            }
                        }
                    }
                }
            }
        }
    }
    SaveState();
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int a = 1; a < nvx - 1; ++a) {
                    for (int b = 1; b < nvy - 1; ++b) {
                        for (int c = 1; c < nvz - 1; ++c) {
                            if (vz[c] > 0) {
                                if (abs(vz[c] * deltaThz) > 1) cout << "z" << vz[c] * deltaThz << endl;
                                f[i][j][k][a][b][c] =
                                        f_time[i][j][k][a][b][c] - (deltaThz * vz[c]) * (f_time[i][j][k][a][b][c] -
                                                                                         f_time[i][j][(nz + k - 1) %
                                                                                                      nz][a][b][c]);
                            } else {
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] - (deltaThz * vz[c]) *
                                                                                 (f_time[i][j][(k + 1) % nz][a][b][c] -
                                                                                  f_time[i][j][k][a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }

    //cout << "***  Finalizing coordinate part***" << endl;

}


void kinetic::VelocityPart() {

    const double sAxelerationCofficient = Converter->GetSAxelerationCofficient();
    const double deltaThvx = deltaT / hvx;
    const double deltaThvy = deltaT / hvy;
    const double deltaThvz = deltaT / hvz;

    SaveState();
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                const double fx_ = sAxelerationCofficient * asx[i][j][k] + al + acx[i][j][k];
                double ***f_ijk = f[i][j][k];
                double ***f_time_ijk = f_time[i][j][k];
                if (fx_ > 0) {
                    if (abs(fx_ * deltaThvx) > 1) cout << "vx " << fx_ * deltaThvx << endl;
                    for (int a = 1; a < nvx - 1; ++a) {
                        for (int b = 1; b < nvy - 1; ++b) {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (deltaThvx * fx_) *
                                                                       (f_time_ijk[a][b][c] -
                                                                        f_time_ijk[a - 1][b][c]);
                            }
                        }
                    }
                } else {
                    for (int a = 1; a < nvx - 1; ++a) {
                        for (int b = 1; b < nvy - 1; ++b) {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (deltaThvx * fx_) *
                                                                       (f_time_ijk[a + 1][b][c] -
                                                                        f_time_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }

    SaveState();
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                const double fy_ = sAxelerationCofficient * asy[i][j][k] + acy[i][j][k];
                double ***f_ijk = f[i][j][k];
                double ***f_time_ijk = f_time[i][j][k];
                if (fy_ > 0) {
                    if (abs(fy_ * deltaThvy) > 1) cout << "vy " << fy_ * deltaThvy << endl;
                    for (int a = 1; a < nvx - 1; ++a) {
                        for (int b = 1; b < nvy - 1; ++b) {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (deltaThvy * fy_) *
                                                                       (f_time_ijk[a][b][c] -
                                                                        f_time_ijk[a][b - 1][c]);
                            }
                        }
                    }
                } else {
                    for (int a = 1; a < nvx - 1; ++a) {
                        for (int b = 1; b < nvy - 1; ++b) {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_ijk[a][b][c] - (deltaThvy * fy_) *
                                                                  (f_ijk[a][b + 1][c] -
                                                                   f_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }

    SaveState();
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                const double fz_ = sAxelerationCofficient * asz[i][j][k] + acz[i][j][k];
                double ***f_ijk = f[i][j][k];
                double ***f_time_ijk = f_time[i][j][k];
                if (fz_ > 0) {
                    if (abs(fz_ * deltaThvz) > 1) cout << "vz " << fz_ * deltaThvz << endl;
                    for (int a = 1; a < nvx - 1; ++a) {
                        for (int b = 1; b < nvy - 1; ++b) {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (deltaThvz * fz_) *
                                                                       (f_time_ijk[a][b][c] -
                                                                        f_time_ijk[a][b][c - 1]);
                            }
                        }
                    }
                } else {
                    for (int a = 1; a < nvx - 1; ++a) {
                        for (int b = 1; b < nvy - 1; ++b) {
                            for (int c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (deltaThvz * fz_) *
                                                                       (f_time_ijk[a][b][c + 1] -
                                                                        f_time_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }

}

double kinetic::Relaxation(int i, int j, int k, int a, int b, int c) {
    return f[i][j][k][a][b][c] + deltaT * (density[i][j][k] * f_n[a][b][c] - f[i][j][k][a][b][c]) / tau;
}


void kinetic::IntegrateAll() {
    double deltaByTau = deltaT / tau;

    CoordinatePart();

    VelocityPart();
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                double ***f_ijk = f[i][j][k];
                double density_ijk = density[i][j][k];
                for (int a = 1; a < nvx - 1; ++a) {
                    for (int b = 1; b < nvy - 1; ++b) {
                        for (int c = 1; c < nvy - 1; ++c) {
                            f_ijk[a][b][c] =
                                    f_ijk[a][b][c] + deltaByTau * (density_ijk * f_n[a][b][c] - f_ijk[a][b][c]);
                        }
                    }
                }
            }
        }
    }

    CurT += deltaT;

    ComputeDensity();

    ComputeFlowVelocity();

}

void kinetic::ComputeDensity() {
    const double m = hvx * hvy * hvz;
    //cout << "*** Compute density ***" << endl;
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                double ***f_ijk = f[i][j][k];
                double sum = 0.0;
                for (int a = 0; a < nvx; ++a) {
                    for (int b = 0; b < nvy; ++b) {
                        for (int c = 0; c < nvz; ++c) {
                            sum += f_ijk[a][b][c];
                        }
                    }
                }
                density[i][j][k] = sum * m;
            }
        }
    }


}


void kinetic::ComputeFlowVelocity() {
    const double m = hvx * hvy * hvz;
    //cout << "*** Compute flow velocity ***" << endl;
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                double ***f_ijk = f[i][j][k];
                double sum_x = 0.0;
                double sum_y = 0.0;
                double sum_z = 0.0;
                for (int a = 0; a < nvx; ++a) {
                    for (int b = 0; b < nvy; ++b) {
                        for (int c = 0; c < nvz; ++c) {
                            sum_x += vx[a] * f_ijk[a][b][c];
                            sum_y += vy[b] * f_ijk[a][b][c];
                            sum_z += vz[c] * f_ijk[a][b][c];
                        }
                    }
                }
                vfl_x[i][j][k] = sum_x * m;
                vfl_y[i][j][k] = sum_y * m;
                vfl_z[i][j][k] = sum_z * m;
            }
        }
    }

}


void kinetic::InitialConditions() {
    const double dimensionlessIonConcentration = Converter->GetDimensionlessIonConcentration();
    for (int i = 0; i < 1; ++i) {
        for (int j = 0; j < 1; ++j) {
            for (int k = 0; k < 1; ++k) {
                for (int a = 0; a < nvx; ++a) {
                    for (int b = 0; b < nvy; ++b) {
                        for (int c = 0; c < nvz; ++c) {
                            double temp = 0.0;
                            double xi = 0;
                            while (xi < MAX_XI) {
                                double vx_a = vx[a] - xi * v_fl_an;
                                double vy_b = vy[b];
                                double vz_c = vz[c];
                                temp += exp(-xi) * exp(-0.5 * (vx_a * vx_a + vy_b * vy_b + vz_c * vz_c));
                                xi += DELTA_XI;
                            }
                            f[i][j][k][a][b][c] = dimensionlessIonConcentration * temp * DELTA_XI
                                                  / 15.749609945722419;//pow(2 * M_PI, 1.5);
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int a = 0; a < nvx; ++a) {
                    for (int b = 0; b < nvy; ++b) {
                        copy(f[0][0][0][a][b], f[0][0][0][a][b] + nvz, f[i][j][k][a][b]);
                    }
                }
            }
        }
    }
/*
    for (int a=0; a<nvx; ++a){
        for (int b=0; b<nvy; ++b){
            for (int c=0; c<nvz; ++c){
	        f[nx/2][ny/2][nz/2][a][b][c] += 0.1 * Converter->GetDimensionlessIonConcentration() * f_n[a][b][c];
	    }
        }
    }
*/
}


double ***kinetic::GetDensity() const {
    return density;
}

int kinetic::GetNvx() const {
    return nvx;
}

int kinetic::GetNvy() const {
    return nvy;
}

int kinetic::GetNvz() const {
    return nvz;
}

double kinetic::Gethvx() const {
    return hvx;
}

double kinetic::Gethvy() const {
    return hvy;
}

double kinetic::Gethvz() const {
    return hvz;
}

double kinetic::GetCutVelocityAlong() const {
    return vxcut;
}

double kinetic::GetCutVelocityNormal() const {
    return vyzcut;
}

double kinetic::GetCurrentTime() const {
    return CurT;
}


double ***kinetic::GetFlowVelocityX() const {
    return vfl_x;
}


double ***kinetic::GetFlowVelocityY() const {
    return vfl_y;
}


double ***kinetic::GetFlowVelocityZ() const {
    return vfl_z;
}


double *kinetic::GetProfileX() const {
    return profile_x;
}


double *kinetic::GetProfileY() const {
    return profile_y;
}


double *kinetic::GetProfileZ() const {
    return profile_z;
}

double *kinetic::GetAnalyticalProfileX() const {
    return analytical_profile_x;
}


void kinetic::SetProfileX(int i, int j, int k) {
    for (int a = 0; a < nvx; ++a) {
        profile_x[a] = f[i][j][k][a][nvy / 2][nvz / 2];
    }

}


void kinetic::SetProfileY(int i, int j, int k) {
    for (int b = 0; b < nvy; ++b) {
        profile_y[b] = f[i][j][k][nvx / 2][b][nvz / 2];
    }
}


void kinetic::SetProfileZ(int i, int j, int k) {
    for (int c = 0; c < nvz; ++c) {
        profile_z[c] = f[i][j][k][nvx / 2][nvy / 2][c];
    }
}

void kinetic::SetAnalyticalProfileX() {
    double dimensionlessIonConcentration = Converter->GetDimensionlessIonConcentration();
    for (int a = 0; a < nvx; ++a) {
        analytical_profile_x[a] = 0;
        double xi = 0;
        while (xi < MAX_XI) {
            analytical_profile_x[a] = analytical_profile_x[a]
                                      + exp(-xi) * exp(-0.5 * ((vx[a] - xi * v_fl_an) * (vx[a] - xi * v_fl_an) +
                                                               vy[nvy / 2] * vy[nvy / 2] + vz[nvz / 2] * vz[nvz / 2]));
            xi += DELTA_XI;
        }
        analytical_profile_x[a] = dimensionlessIonConcentration * analytical_profile_x[a] * DELTA_XI
                                  / 15.749609945722419; //pow(2 * M_PI, 1.5);
    }
}


double kinetic::GetDeltaT() const {
    return deltaT;
}

void kinetic::SaveState() {
#pragma omp parallel for collapse(3)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                for (int a = 0; a < nvx; ++a) {
                    for (int b = 0; b < nvy; ++b) {
                        std::copy(f[i][j][k][a][b], f[i][j][k][a][b] + nvz, f_time[i][j][k][a][b]);
                    }
                }
            }
        }
    }
}


double *kinetic::GetVelocitySetX() const {
    return vx;
}

double *kinetic::GetVelocitySetY() const {
    return vx;
}

double *kinetic::GetVelocitySetZ() const {
    return vx;
}


void kinetic::DefineDensity() {

    double n_i = Converter->GetDimensionlessIonConcentration();

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                density[i][j][k] = ((M_PI / Lx) * (M_PI / Lx) + (M_PI / Lx) * (M_PI / Lx) + (M_PI / Lz) * (M_PI / Lz))
                                   * sin(M_PI * x[i] / Lx) * sin(M_PI * y[j] / Ly) * sin(M_PI * z[k] / Lz) /
                                   (4. * M_PI) + n_i;
            }
        }
    }

}
