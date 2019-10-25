//
// Created by Igor on 25.10.2019.
//
#include "IonWake.h"

#include <vector>
#include <cmath>
#include <algorithm>

using namespace ionwake;

IonWake::IonWake(size_t nx, size_t ny, size_t nz, size_t nx_0, size_t ny_0, size_t nz_0, double lx, double ly,
                 double lz, double particle_charge, double t, double n_0, double tau, double el)
        : nx(nx), ny(ny), nz(nz), nx_0(nx_0), ny_0(ny_0), nz_0(nz_0), lx(lx), ly(ly), lz(lz),
          particleCharge(particle_charge), temperature(t), physicalIonConcentration(n_0), physicalTau(tau),
          strength(el), currentTime(0.0), hvx(0.2), hvy(0.2), hvz(0.2) {

    thermalVelocity = std::sqrt(t * BOLTZMANN_K / ION_MASS);
    debayRadius = sqrt(BOLTZMANN_K * t / (4.0 * M_PI * ELECTRON_CHARGE * ELECTRON_CHARGE * n_0));
    plasmasFrequency = thermalVelocity / debayRadius;
    convertPotential = ELECTRON_CHARGE / debayRadius;

    accelerationCoefficientS =
            ELECTRON_CHARGE * ELECTRON_CHARGE / (ION_MASS * thermalVelocity * thermalVelocity * debayRadius);
    accelerationCoefficientL = el * ELECTRON_CHARGE * debayRadius / (ION_MASS * thermalVelocity * thermalVelocity);
    accelerationCoefficientC =
            particle_charge * ELECTRON_CHARGE * ELECTRON_CHARGE /
            (ION_MASS * thermalVelocity * thermalVelocity * debayRadius);

    coordinateStepX = lx / nx;
    coordinateStepY = ly / ny;
    coordinateStepZ = lz / nz;

    dimensionlessFlowVelocity = ELECTRON_CHARGE * el * tau / (ION_MASS * thermalVelocity);
    dimensionlessIonConcentration = n_0 * (debayRadius * debayRadius * debayRadius);
    dimensionlessTau = tau * thermalVelocity / debayRadius;
    analyticalFlowVelocity = accelerationCoefficientL * dimensionlessTau;

    vxcut = std::max(5.0, 8.0 * dimensionlessFlowVelocity);
    vyzcut = 5.0;
    nvx = ((vxcut + vyzcut) / hvx);
    nvy = (2.0 * vyzcut / hvy);
    nvz = (2.0 * vyzcut / hvz);

    deltaT = 0.5 * std::min(0.2 / (accelerationCoefficientL / hvx), 1.0 / (vxcut / coordinateStepX +
                                                                           vyzcut / coordinateStepY +
                                                                           vyzcut / coordinateStepZ));
    vx = new double[nvx];
    vy = new double[nvy];
    vz = new double[nvz];
    for (size_t i = 0; i < nvx; ++i) { vx[i] = i * hvx - vyzcut; } // vyzcut???
    for (size_t i = 0; i < nvy; ++i) { vy[i] = i * hvy - vyzcut; }
    for (size_t i = 0; i < nvz; ++i) { vz[i] = i * hvz - vyzcut; }

    x = new double[nx];
    y = new double[ny];
    z = new double[nz];
    for (size_t i = 0; i < nx; ++i) { x[i] = i * coordinateStepX; }
    for (size_t i = 0; i < ny; ++i) { y[i] = i * coordinateStepY; }
    for (size_t i = 0; i < nz; ++i) { z[i] = i * coordinateStepZ; }

    f_n = new double **[nvx];
    for (size_t i = 0; i < nvx; ++i) {
        f_n[i] = new double *[nvy];
        for (size_t j = 0; j < nvy; ++j) {
            f_n[i][j] = new double[nvz];
            for (size_t k = 0; k < nvz; ++k) {
                const double x = vx[i] * vx[i];
                const double y = vy[j] * vy[j];
                const double z = vz[k] * vz[k];
                f_n[i][j][k] = std::exp(-0.5 * (x + y + z)) / std::pow(2.0 * M_PI, 1.5);
            }
        }
    }

    f = new double *****[nx];
    for (size_t i = 0; i < nx; ++i) {
        f[i] = new double ****[ny];
        for (size_t j = 0; j < ny; ++j) {
            f[i][j] = new double ***[nz];
            for (size_t k = 0; k < nz; ++k) {
                f[i][j][k] = new double **[nvx];
                for (size_t a = 0; a < nvx; ++a) {
                    f[i][j][k][a] = new double *[nvy];
                    for (size_t b = 0; b < nvy; ++b) {
                        f[i][j][k][a][b] = new double[nvz];
                    }
                }
            }
        }
    }

    f_time = new double *****[nx];
    for (size_t i = 0; i < nx; ++i) {
        f_time[i] = new double ****[ny];
        for (size_t j = 0; j < ny; ++j) {
            f_time[i][j] = new double ***[nz];
            for (size_t k = 0; k < nz; ++k) {
                f_time[i][j][k] = new double **[nvx];
                for (size_t a = 0; a < nvx; ++a) {
                    f_time[i][j][k][a] = new double *[nvy];
                    for (size_t b = 0; b < nvy; ++b) {
                        f_time[i][j][k][a][b] = new double[nvz];
                    }
                }
            }
        }
    }

    //Initial conditions
    for (size_t i = 0; i < 1; ++i) {
        for (size_t j = 0; j < 1; ++j) {
            for (size_t k = 0; k < 1; ++k) {
                for (size_t a = 0; a < nvx; ++a) {
                    for (size_t b = 0; b < nvy; ++b) {
                        for (size_t c = 0; c < nvz; ++c) {
                            double sum = 0.0;
                            double xi = 0;
                            while (xi < MAX_XI) {
                                double vx_a = vx[a] - xi * analyticalFlowVelocity;
                                double vy_b = vy[b];
                                double vz_c = vz[c];

                                sum += exp(-xi) * exp(-0.5 * (vx_a * vx_a + vy_b * vy_b + vz_c * vz_c));
                                xi += DELTA_XI;
                            }
                            f[i][j][k][a][b][c] = dimensionlessIonConcentration * sum * DELTA_XI / pow(2.0 * M_PI, 1.5);
                        }
                    }
                }
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                for (size_t a = 0; a < nvx; ++a) {
                    for (size_t b = 0; b < nvy; ++b) {
                        for (size_t c = 0; c < nvz; ++c) {
                            f[i][j][k][a][b][c] = f[0][0][0][a][b][c];
                        }
                    }
                }
            }
        }
    }

    acx = new double **[nx];
    acy = new double **[nx];
    acz = new double **[nx];
    flowVelocityX = new double **[nx];
    flowVelocityY = new double **[nx];
    flowVelocityZ = new double **[nx];
    selfConsistentForceFieldX = new double **[nx];
    selfConsistentForceFieldY = new double **[nx];
    selfConsistentForceFieldZ = new double **[nx];
    potential = new double **[nx];
    density = new double **[nx];
    for (size_t i = 0; i < nx; ++i) {
        acx[i] = new double *[ny];
        acy[i] = new double *[ny];
        acz[i] = new double *[ny];
        flowVelocityX[i] = new double *[ny];
        flowVelocityY[i] = new double *[ny];
        flowVelocityZ[i] = new double *[ny];
        selfConsistentForceFieldX[i] = new double *[ny];
        selfConsistentForceFieldY[i] = new double *[ny];
        selfConsistentForceFieldZ[i] = new double *[ny];
        potential[i] = new double *[ny];
        density[i] = new double *[ny];
        for (size_t j = 0; j < ny; ++j) {
            acx[i][j] = new double[nz]{};
            acy[i][j] = new double[nz]{};
            acz[i][j] = new double[nz]{};
            flowVelocityX[i][j] = new double[nz]{};
            flowVelocityY[i][j] = new double[nz]{};
            flowVelocityZ[i][j] = new double[nz]{};
            selfConsistentForceFieldX[i][j] = new double[ny];
            selfConsistentForceFieldY[i][j] = new double[ny];
            selfConsistentForceFieldZ[i][j] = new double[ny];
            potential[i][j] = new double[ny]{};

            density[i][j] = new double[nz];
            std::fill(density[i][j], density[i][j] + nz, dimensionlessIonConcentration);
        }
    }
}

void IonWake::nextStep() {
    poissonScheme();
    gradientScheme();

    saveStep();
    coordinatePart();

    saveStep();
    velocityPart();

#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                for (size_t a = 1; a < nvx - 1; ++a) {
                    for (size_t b = 1; b < nvy - 1; ++b) {
                        for (size_t c = 1; c < nvz - 1; ++c) {
                            f_ijk[a][b][c] +=
                                    deltaT * (density[i][j][k] * f_n[a][b][c] - f_ijk[a][b][c]) / dimensionlessTau;
                        }
                    }
                }
            }
        }
    }

    currentTime += deltaT;

    computeDensity();
    computeFlowVelocity();
}

void IonWake::poissonScheme() {
    const double hx = coordinateStepX;
    const double hy = coordinateStepY;
    const double hz = coordinateStepZ;
    const double owx = hx * hx;
    const double owy = hy * hy;
    const double owz = hz * hz;

    const double owlx = owx * nx * nx;
    const double owly = owy * ny * ny;
    const double owlz = owz * nz * nz;

    const double tau = 2. / (4. * (1. / owx + 1. / owy + 1. / owz) + M_PI * M_PI * (1. / owlx + 1. / owly + 1. / owlz));

    bool proceed = true;
    while (proceed) {
        proceed = false;
        for (size_t i = 1; i < nx - 1; ++i) {
            for (size_t j = 1; j < ny - 1; ++j) {
                for (size_t k = 1; k < nz - 1; ++k) {
                    const double f1 = potential[i][j][k];

                    const double double_f1 = 2.0 * f1;
                    const double fi = (potential[i + 1][j][k] - double_f1 + potential[i - 1][j][k]) / owx;
                    const double fj = (potential[i][j + 1][k] - double_f1 + potential[i][j - 1][k]) / owy;
                    const double fk = (potential[i][j][k + 1] - double_f1 + potential[i][j][k - 1]) / owz;

                    const double first = density[i][j][k] / dimensionlessIonConcentration - 1.0;
                    const double second =
                            particleCharge / (dimensionlessIonConcentration * hx * hy * hz * nz * ny * nz);
                    const double third = (i == nx_0 && j == ny_0 && k == nz_0 ? particleCharge : 0.0);
                    const double strange =
                            4. * M_PI * (first + second - third) / (dimensionlessIonConcentration * hx * hy * hz);

                    potential[i][j][k] += (fi + fj + fk + strange) * tau;

                    if (potential[i][j][k] == 0.0 && std::abs((potential[i][j][k] - f1) / potential[i][j][k]) > E) {
                        proceed = true;
                    }
                }
            }
        }
    }
}

void IonWake::gradientScheme() {
    const double x = 2.0 * coordinateStepX;
    const double y = 2.0 * coordinateStepY;
    const double z = 2.0 * coordinateStepZ;

#pragma omp parallel for collapse(3)
    for (size_t i = 1; i < nx - 1; ++i) {
        for (size_t j = 1; j < ny - 1; ++j) {
            for (size_t k = 1; k < nz - 1; ++k) {
                selfConsistentForceFieldX[i][j][k] =
                        -(potential[i + 1][j][k] - potential[i - 1][j][k]) / x * dimensionlessIonConcentration;
                selfConsistentForceFieldY[i][j][k] =
                        -(potential[i][j + 1][k] - potential[i][j - 1][k]) / y * dimensionlessIonConcentration;
                selfConsistentForceFieldZ[i][j][k] =
                        -(potential[i][j][k + 1] - potential[i][j][k - 1]) / z * dimensionlessIonConcentration;
            }
        }
    }
}

void IonWake::coordinatePart() {
    const double delta_x_step = deltaT / coordinateStepX;
    const double delta_y_step = deltaT / coordinateStepY;
    const double delta_z_step = deltaT / coordinateStepZ;
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                double ***const f_time_ijk = f_time[i][j][k];
                for (size_t a = 1; a < nvx - 1; ++a) {
                    const double vx_a = vx[a];
                    if (vx_a > 0.0) {
                        for (size_t b = 1; b < nvy - 1; ++b) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (delta_x_step * vx_a) * (f_time_ijk[a][b][c] -
                                                                                                f_time[(nx + i - 1) %
                                                                                                       nx][j][k][a][b][c]);
                            }
                        }
                    } else {
                        for (size_t b = 1; b < nvy - 1; ++b) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (delta_x_step * vx_a) *
                                                                       (f_time[(i + 1) % nx][j][k][a][b][c] -
                                                                        f_time_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                double ***const f_time_ijk = f_time[i][j][k];
                for (size_t a = 1; a < nvx - 1; ++a) {
                    for (size_t b = 1; b < nvy - 1; ++b) {
                        const double vy_b = vy[b];
                        if (vy_b > 0.0) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (delta_y_step * vy_b) * (f_time_ijk[a][b][c] -
                                                                                                f_time[i][(ny + j - 1) %
                                                                                                          ny][k][a][b][c]);
                            }
                        } else {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (delta_y_step * vy_b) *
                                                                       (f_time[i][(j + 1) % ny][k][a][b][c] -
                                                                        f_time_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                double ***const f_time_ijk = f_time[i][j][k];
                for (size_t a = 1; a < nvx - 1; ++a) {
                    for (size_t b = 1; b < nvy - 1; ++b) {
                        for (size_t c = 1; c < nvz - 1; ++c) {
                            const double vz_c = vz[c];
                            if (vz_c > 0.0) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (delta_z_step * vz_c) * (f_time_ijk[a][b][c] -
                                                                                                f_time[i][j][
                                                                                                        (nz + k - 1) %
                                                                                                        nz][a][b][c]);
                            } else {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] - (delta_z_step * vz_c) *
                                                                       (f_time[i][j][(k + 1) % nz][a][b][c] -
                                                                        f_time_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void IonWake::velocityPart() {
    const double delta_t_hvx = deltaT / hvx;
    const double delta_t_hvy = deltaT / hvy;
    const double delta_t_hvz = deltaT / hvz;

#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                double ***const f_time_ijk = f_time[i][j][k];
                const double fx = accelerationCoefficientS * selfConsistentForceFieldX[i][j][k] +
                                  accelerationCoefficientL + acx[i][j][k];
                if (fx > 0.0) {
                    for (size_t a = 1; a < nvx - 1; ++a) {
                        for (size_t b = 1; b < nvy - 1; ++b) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] -
                                                 (delta_t_hvx * fx) * (f_time_ijk[a][b][c] - f_time_ijk[a - 1][b][c]);
                            }
                        }
                    }
                } else {
                    for (size_t a = 1; a < nvx - 1; ++a) {
                        for (size_t b = 1; b < nvy - 1; ++b) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] -
                                                 (delta_t_hvx * fx) * (f_time_ijk[a + 1][b][c] - f_time_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                double ***const f_time_ijk = f_time[i][j][k];
                const double fy = accelerationCoefficientS * selfConsistentForceFieldY[i][j][k] + acy[i][j][k];
                if (fy > 0.0) {
                    for (size_t a = 1; a < nvx - 1; ++a) {
                        for (size_t b = 1; b < nvy - 1; ++b) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] -
                                                 (delta_t_hvy * fy) * (f_time_ijk[a][b][c] - f_time_ijk[a][b - 1][c]);
                            }
                        }
                    }
                } else {
                    for (size_t a = 1; a < nvx - 1; ++a) {
                        for (size_t b = 1; b < nvy - 1; ++b) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] -
                                                 (delta_t_hvy * fy) * (f_time_ijk[a][b + 1][c] - f_time_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                double ***const f_time_ijk = f_time[i][j][k];
                const double fz = accelerationCoefficientS * selfConsistentForceFieldZ[i][j][k] + acz[i][j][k];
                if (fz > 0.0) {
                    for (size_t a = 1; a < nvx - 1; ++a) {
                        for (size_t b = 1; b < nvy - 1; ++b) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] -
                                                 (delta_t_hvz * fz) * (f_time_ijk[a][b][c] - f_time_ijk[a][b][c - 1]);
                            }
                        }
                    }
                } else {
                    for (size_t a = 1; a < nvx - 1; ++a) {
                        for (size_t b = 1; b < nvy - 1; ++b) {
                            for (size_t c = 1; c < nvz - 1; ++c) {
                                f_ijk[a][b][c] = f_time_ijk[a][b][c] -
                                                 (delta_t_hvz * fz) * (f_time_ijk[a][b][c + 1] - f_time_ijk[a][b][c]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void IonWake::computeDensity() {
    const double m = hvx * hvy * hvz;
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                double sum = 0.0;
                for (size_t a = 0; a < nvx; ++a) {
                    for (size_t b = 0; b < nvy; ++b) {
                        for (size_t c = 0; c < nvz; ++c) {
                            sum += f_ijk[a][b][c];
                        }
                    }
                }
                density[i][j][k] = sum * m;
            }
        }
    }
}

void IonWake::computeFlowVelocity() {
    const double m = hvx * hvy * hvz;
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                double ***const f_ijk = f[i][j][k];
                double sum_x = 0.0;
                double sum_y = 0.0;
                double sum_z = 0.0;
                for (size_t a = 0; a < nvx; ++a) {
                    for (size_t b = 0; b < nvy; ++b) {
                        for (size_t c = 0; c < nvz; ++c) {
                            sum_x += vx[a] * f_ijk[a][b][c];
                            sum_y += vy[b] * f_ijk[a][b][c];
                            sum_z += vz[c] * f_ijk[a][b][c];
                        }
                    }
                }
                flowVelocityX[i][j][k] = sum_x * m;
                flowVelocityY[i][j][k] = sum_y * m;
                flowVelocityZ[i][j][k] = sum_z * m;
            }
        }
    }
}


void IonWake::saveStep() {
    auto temp = f;
    f = f_time;
    f_time = temp;
}


void IonWake::writeInitialPotential(std::ofstream &of) const {
    of << nx << "\t" << ny << "\t" << nz << "\n";
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                const double x1 = (double) (i - nx_0);
                const double y1 = (double) (j - ny_0);
                const double z1 = (double) (k - nz_0);
                const double x2 = lx * debayRadius / double(nx);
                const double y2 = ly * debayRadius / double(ny);
                const double z2 = lz * debayRadius / double(nz);

                double r = sqrt(x1 * x1 * x2 * x2 + y1 * y1 * y2 * y2 + z1 * z1 * z2 * z2);
                of << -particleCharge * ELECTRON_CHARGE / r - strength * (k * lz * debayRadius / double(nz)) << "\n";
            }
        }
    }
}

void IonWake::writeDensity(std::ofstream &of) const {
    of << nx << "\t" << ny << "\t" << nz << "\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                of << density[i][j][k] * pow(1.0 / debayRadius, 3) << "\n";
            }
        }
    }
}

void IonWake::writePotential(std::ofstream &of) const {
    of << nx << "\t" << ny << "\t" << nz << "\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                of << dimensionlessIonConcentration * convertPotential * potential[i][j][k] << "\n";
            }
        }
    }
}

void IonWake::writeVelocity(std::ofstream &of) const {
    of << nx << "\t" << ny << "\t" << nz << "\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                of << flowVelocityX[i][j][k] << "\t" << flowVelocityY[i][j][k] << "\t" << flowVelocityZ[i][j][k]
                   << "\n";
            }
        }
    }
}

void IonWake::plotDistributionFunctionX(size_t i, size_t j, size_t k, std::ofstream &of) const {
    for (size_t a = 0; a < nvx; ++a) {
        of << vx[a] << "\t" << f[i][j][k][a][nvy / 2][nvz / 2] << "\n";
    }
}

void IonWake::plotDistributionFunctionY(size_t i, size_t j, size_t k, std::ofstream &of) const {
    for (size_t b = 0; b < nvx; ++b) {
        of << vy[b] << "\t" << f[i][j][k][nvx / 2][b][nvz / 2] << "\n";
    }
}

void IonWake::plotDistributionFunctionZ(size_t i, size_t j, size_t k, std::ofstream &of) const {
    for (size_t c = 0; c < nvx; ++c) {
        of << vz[c] << "\t" << f[i][j][k][nvx / 2][nvy / 2][c] << "\n";
    }
}


double IonWake::getTemperature() const {
    return temperature;
}

double IonWake::getThermalVelocity() const {
    return thermalVelocity;
}

double IonWake::getDebayRadius() const {
    return debayRadius;
}

double IonWake::getDeltaT() const {
    return deltaT;
}

double IonWake::getParticleCharge() const {
    return particleCharge;
}

double IonWake::getConvertPotential() const {
    return convertPotential;
}

double IonWake::getPlasmasFrequency() const {
    return plasmasFrequency;
}

double IonWake::getPhysicalIonConcentration() const {
    return physicalIonConcentration;
}

double IonWake::getPhysicalTau() const {
    return physicalTau;
}

double IonWake::getStrength() const {
    return strength;
}

double IonWake::getCoordinateStepX() const {
    return coordinateStepX;
}

double IonWake::getCoordinateStepY() const {
    return coordinateStepY;
}

double IonWake::getCoordinateStepZ() const {
    return coordinateStepZ;
}

size_t IonWake::getNx() const {
    return nx;
}

size_t IonWake::getNy() const {
    return ny;
}

size_t IonWake::getNz() const {
    return nz;
}

size_t IonWake::getNvx() const {
    return nvx;
}

size_t IonWake::getNvy() const {
    return nvy;
}

size_t IonWake::getNvz() const {
    return nvz;
}

double IonWake::getHvx() const {
    return hvx;
}

double IonWake::getHvy() const {
    return hvy;
}

double IonWake::getHvz() const {
    return hvz;
}

double IonWake::getCutVelocityAlong() const {
    return vxcut;
}

double IonWake::getCutVelocityNormal() const {
    return vyzcut;
}

double IonWake::getDimensionlessIonConcentration() const {
    return dimensionlessIonConcentration;
}

double IonWake::getDimensionlessTau() const {
    return dimensionlessTau;
}

double IonWake::getAnalyticalFlowVelocity() const {
    return analyticalFlowVelocity;
}

double IonWake::getDimensionlessFlowVelocity() const {
    return dimensionlessFlowVelocity;
}

double IonWake::getCurrentTime() const {
    return currentTime;
}

IonWake::~IonWake() {
    delete[] vx;
    delete[] vy;
    delete[] vz;

    delete[] x;
    delete[] y;
    delete[] z;

    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            delete[] density[i][j];
            delete[] potential[i][j];
            delete[] flowVelocityX[i][j];
            delete[] flowVelocityY[i][j];
            delete[] flowVelocityZ[i][j];
            delete[] selfConsistentForceFieldX[i][j];
            delete[] selfConsistentForceFieldY[i][j];
            delete[] selfConsistentForceFieldZ[i][j];
            delete[] acx[i][j];
            delete[] acy[i][j];
            delete[] acz[i][j];
        }
        delete[] density[i];
        delete[] potential[i];
        delete[] flowVelocityX[i];
        delete[] flowVelocityY[i];
        delete[] flowVelocityZ[i];
        delete[] selfConsistentForceFieldX[i];
        delete[] selfConsistentForceFieldY[i];
        delete[] selfConsistentForceFieldZ[i];
        delete[] acx[i];
        delete[] acy[i];
        delete[] acz[i];
    }

    delete[] density;
    delete[] potential;
    delete[] flowVelocityX;
    delete[] flowVelocityY;
    delete[] flowVelocityZ;
    delete[] selfConsistentForceFieldX;
    delete[] selfConsistentForceFieldY;
    delete[] selfConsistentForceFieldZ;
    delete[] acx;
    delete[] acy;
    delete[] acz;

    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                for (size_t a = 0; a < nvx; ++a) {
                    for (size_t b = 0; b < nvy; ++b) {
                        delete[] f[i][j][k][a][b];
                        delete[] f_time[i][j][k][a][b];
                    }
                    delete[] f[i][j][k][a];
                    delete[] f_time[i][j][k][a];
                }
                delete[] f[i][j][k];
                delete[] f_time[i][j][k];
            }
            delete[] f[i][j];
            delete[] f_time[i][j];
        }
        delete[] f[i];
        delete[] f_time[i];
    }
    delete[] f;
    delete[] f_time;


    for (size_t a = 0; a < nvx; ++a) {
        for (size_t b = 0; b < nvy; ++b) {
            delete[] f_n[a][b];
        }
        delete[] f_n[a];
    }
    delete[] f_n;
}



//Debay radious


