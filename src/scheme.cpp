#define _USE_MATH_DEFINES

#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "converter.hpp"
#include "scheme.hpp"
#include "omp.h"

using namespace std;

Scheme::Scheme(
        const Converter &converter_, size_t nx_, size_t ny_, size_t nz_,
        double lx_, double ly_, double lz_,
        int nt_, double r0x_, double r0y_, double r0z_,
        double hvx_ = 0.2, double hvy_ = 0.2, double hvz_ = 0.2,
        double vminyz_ = -5, double vmaxyz_ = 5,
        double hxi_ = 0.002, double ximax_ = 13
) :
        converter(converter_),
        nx(nx_), ny(ny_), nz(nz_),
        lx(lx_), ly(ly_), lz(lz_),
        nt(nt_), r0x(r0x_), r0y(r0y_), r0z(r0z_),
        hvx(hvx_), hvy(hvy_), hvz(hvz_),
        vminyz(vminyz_),
        vmaxyz(vmaxyz_),
        hxi(hxi_), ximax(ximax_) {
    Eext = converter.externalField();
    wc = converter.collisionFrequency();
    vfl = converter.flowVelocity();
    rde = converter.electronDebyeRadious();
    q = converter.dustParticleCharge();

    hx = lx / double(nx);
    hy = ly / double(ny);
    hz = lz / double(nz);

    vminx = vminyz;
    vmaxx = vmaxxCompute();

    nvx = size_t((vmaxx - vminx) / hvx);
    nvy = size_t((vmaxyz - vminyz) / hvy);
    nvz = size_t((vmaxyz - vminyz) / hvz);

    dt = 0.5 * min(hvx / (Eext + (q / (hx * hx))), hx / vmaxx);
    t = 0;

    rx = new double[nx];
    for (size_t i = 0; i < nx; i++) {
        rx[i] = i * hx;
    }
    ry = new double[ny];
    for (size_t j = 0; j < ny; j++) {
        ry[j] = j * hy;
    }
    rz = new double[nz];
    for (size_t k = 0; k < nz; k++) {
        rz[k] = k * hz;
    }
    vx = new double[nvx];
    for (size_t a = 0; a < nvx; a++) {
        vx[a] = vminx + a * hvx;
    }
    vy = new double[nvy];
    for (size_t b = 0; b < nvy; b++) {
        vy[b] = vminyz + b * hvy;
    }
    vz = new double[nvz];
    for (size_t c = 0; c < nvz; c++) {
        vz[c] = vminyz + c * hvz;
    }

    fullCharge = 0;


    f = TArray<double, 6>({nx, ny, nz, nvx, nvy, nvz});
    ftime = TArray<double, 6>({nx, ny, nz, nvx, nvy, nvz});
    maxvell = TArray<double, 3>({nvx, nvy, nvz});

    n = TArray<double, 3>({nx, ny, nz});
    n.fill(1.0);
    fi = TArray<double, 3>({nx, ny, nz});
    ax = TArray<double, 3>({nx, ny, nz});
    ay = TArray<double, 3>({nx, ny, nz});
    az = TArray<double, 3>({nx, ny, nz});
#pragma omp parallel for collapse(3)
    for (size_t a = 0; a < nvx; a++) {
        for (size_t b = 0; b < nvy; b++) {
            for (size_t c = 0; c < nvz; c++) {
                f[{0, 0, 0, a, b, c}] = initialDisrtibutionFunction(vx[a], vy[b], vz[c]);
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                for (size_t a = 0; a < nvx; a++) {
                    for (size_t b = 0; b < nvy; b++) {
                        for (size_t c = 0; c < nvz; c++) {
                            f[{i, j, k, a, b, c}] = f[{0, 0, 0, a, b, c}];
                            ftime[{i, j, k, a, b, c}] = f[{0, 0, 0, a, b, c}];
                        }
                    }
                }
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (size_t a = 0; a < nvx; a++) {
        for (size_t b = 0; b < nvy; b++) {
            for (size_t c = 0; c < nvz; c++) {
                maxvell[{a, b, c}] = maxwellianDistribution(vx[a], vy[b], vz[c]);
            }
        }
    }

}

Scheme::~Scheme() {
    delete[] rx;
    delete[] ry;
    delete[] rz;
    delete[] vx;
    delete[] vy;
    delete[] vz;
}

double Scheme::vmaxxCompute() {
    double distribution;
    const double maxvell = exp(-vmaxyz * vmaxyz * 0.5);
    double vmax = 0;
    do {
        vmax = vmax + hvx;
        distribution = 0;
        double xi = 0;
        while (xi < ximax) {
            distribution += exp(-xi) * exp(-0.5 * (vmax - xi * vfl) * (vmax - xi * vfl)) * hxi;
            xi = xi + hxi;
        }
    } while (distribution > maxvell);
    return vmax;
}

double Scheme::initialDisrtibutionFunction(double vx, double vy, double vz) {
    double distribution = 0;
    double xi = 0;
    while (xi < ximax) {
        distribution += exp(-xi) * exp(-0.5 * (vx - xi * vfl) * (vx - xi * vfl)) * exp(-0.5 * vy * vy) *
                        exp(-0.5 * vz * vz) * hxi;
        xi = xi + hxi;
    }
    distribution = distribution / pow(2 * M_PI, 1.5);
    return distribution;
}

double Scheme::maxwellianDistribution(double vx, double vy, double vz) {
    return exp(-0.5 * vx * vx) * exp(-0.5 * vy * vy) * exp(-0.5 * vz * vz) / pow(2 * M_PI, 1.5);
}

void Scheme::density() {
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                double density = 0;
                for (size_t a = 0; a < nvx; a++) {
                    for (size_t b = 0; b < nvy; b++) {
                        for (size_t c = 0; c < nvz; c++) {
                            density += f[{i, j, k, a, b, c}];
                        }
                    }
                }
                n[{i, j, k}] = density * hvx * hvy * hvx;
            }
        }
    }
}

double Scheme::potentialDebye(double rx, double ry, double rz) {
    const double r = distance(rx, ry, rz, r0x, r0y, r0z);
    const double t = sqrt(hx * hx + hy * hy + hz * hz);
    const double distance = r < t ? t : r;
    return -q * exp(-distance / rde) / distance;
}

inline double Scheme::distance(double r1x, double r1y, double r1z, double r2x, double r2y, double r2z) {
    return sqrt((r1x - r2x) * (r1x - r2x) + (r1y - r2y) * (r1y - r2y) + (r1z - r2z) * (r1z - r2z));
}

void Scheme::force() {
#pragma omp parallel for collapse(3)
    for (size_t i = 1; i < nx - 1; i++) {
        for (size_t j = 1; j < ny - 1; j++) {
            for (size_t k = 1; k < nz - 1; k++) {
                ax[{i, j, k}] = -(fi[{i + 1, j, k}] - fi[{i - 1, j, k}]) / (2.0 * hx) + Eext;
                ay[{i, j, k}] = -(fi[{i, j + 1, k}] - fi[{i, j - 1, k}]) / (2.0 * hy);
                az[{i, j, k}] = -(fi[{i, j, k + 1}] - fi[{i, j, k - 1}]) / (2.0 * hz);
            }
        }
    }
}

void Scheme::kinetic() {
    ///coordinate propagation part
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                const size_t i_prev = i == 0 ? nx - 1 : i - 1;
                const size_t i_next = i == nx - 1 ? 0 : i + 1;

                for (size_t a = 1; a < nvx - 1; a++) {
                    if (vx[a] > 0) {
                        for (size_t b = 1; b < nvy - 1; b++) {
                            for (size_t c = 1; c < nvz - 1; c++) {
                                const double diff = ftime[{i, j, k, a, b, c}] - ftime[{i_prev, j, k, a, b, c}];
                                f[{i, j, k, a, b, c}] = ftime[{i, j, k, a, b, c}] - vx[a] * diff * dt / hx;
                            }
                        }
                    } else {
                        for (size_t b = 1; b < nvy - 1; b++) {
                            for (size_t c = 1; c < nvz - 1; c++) {
                                const double diff = ftime[{i_next, j, k, a, b, c}] - ftime[{i, j, k, a, b, c}];
                                f[{i, j, k, a, b, c}] = ftime[{i, j, k, a, b, c}] - vx[a] * diff * dt / hx;
                            }
                        }
                    }
                }
            }
        }
    }

#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                for (size_t a = 1; a < nvx - 1; a++) {
                    const size_t j_prev = j == 0 ? ny - 1 : j - 1;
                    const size_t j_next = j == ny - 1 ? 0 : j + 1;

                    for (size_t b = 1; b < nvy - 1; b++) {
                        if (vy[b] > 0) {
                            for (size_t c = 1; c < nvz - 1; c++) {
                                const double diff = f[{i, j, k, a, b, c}] - f[{i, j_prev, k, a, b, c}];
                                ftime[{i, j, k, a, b, c}] = f[{i, j, k, a, b, c}] - vy[b] * diff * dt / hy;
                            }
                        } else {
                            for (size_t c = 1; c < nvz - 1; c++) {
                                const double diff = f[{i, j_next, k, a, b, c}] - f[{i, j, k, a, b, c}];
                                ftime[{i, j, k, a, b, c}] = f[{i, j, k, a, b, c}] - vy[b] * diff * dt / hy;
                            }
                        }
                    }
                }
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                const size_t k_prev = k == 0 ? nz - 1 : k - 1;
                const size_t k_next = k == nz - 1 ? 0 : k + 1;

                for (size_t a = 1; a < nvx - 1; a++) {
                    for (size_t b = 1; b < nvy - 1; b++) {
                        for (size_t c = 1; c < nvz - 1; c++) {
                            if (vz[c] > 0) {
                                const double diff = ftime[{i, j, k, a, b, c}] - ftime[{i, j, k_prev, a, b, c}];
                                f[{i, j, k, a, b, c}] = ftime[{i, j, k, a, b, c}] - vz[c] * diff * dt / hz;
                            } else {
                                const double diff = ftime[{i, j, k_next, a, b, c}] - ftime[{i, j, k, a, b, c}];
                                f[{i, j, k, a, b, c}] = ftime[{i, j, k, a, b, c}] - vz[c] * diff * dt / hz;
                            }
                        }
                    }
                }
            }
        }
    }
    ///velocity propagation part
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                const double ax_ijk = ax[{i, j, k}];
                const size_t a_shift_1 = ax_ijk > 0 ? 0 : 1;
                const size_t a_shift_2 = ax_ijk > 0 ? 1 : 0;
                
                for (size_t a = 1; a < nvx - 1; a++) {
                    for (size_t b = 1; b < nvy - 1; b++) {
                        for (size_t c = 1; c < nvz - 1; c++) {
                            const double diff = f[{i, j, k, a + a_shift_1, b, c}] - f[{i, j, k, a - a_shift_2, b, c}];
                            ftime[{i, j, k, a, b, c}] = f[{i, j, k, a, b, c}] - ax_ijk * diff * dt / hvx;
                        }
                    }
                }
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                const double ay_ijk = ay[{i, j, k}];
                const size_t b_shift_1 = ay_ijk > 0 ? 0 : 1;
                const size_t b_shift_2 = ay_ijk > 0 ? 1 : 0;

                for (size_t a = 1; a < nvx - 1; a++) {
                    for (size_t b = 1; b < nvy - 1; b++) {
                        for (size_t c = 1; c < nvz - 1; c++) {
                            const double diff = ftime[{i, j, k, a, b + b_shift_1, c}] - ftime[{i, j, k, a, b - b_shift_2, c}];
                            f[{i, j, k, a, b, c}] = ftime[{i, j, k, a, b, c}] - ay_ijk * diff * dt / hvy;
                        }
                    }
                }
            }
        }
    }
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                double az_ijk = az[{i, j, k}];
                const size_t c_shift_1 = az_ijk > 0 ? 0 : 1;
                const size_t c_shift_2 = az_ijk > 0 ? 1 : 0;

                for (size_t a = 1; a < nvx - 1; a++) {
                    for (size_t b = 1; b < nvy - 1; b++) {
                        for (size_t c = 1; c < nvz - 1; c++) {
                            const double diff = f[{i, j, k, a, b, c + c_shift_1}] - f[{i, j, k, a, b, c - c_shift_2}];
                            ftime[{i, j, k, a, b, c}] = f[{i, j, k, a, b, c}] - az_ijk * diff * dt / hvz;
                        }
                    }
                }
            }
        }
    }
    ///Collision part
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                for (size_t a = 1; a < nvx - 1; a++) {
                    for (size_t b = 1; b < nvy - 1; b++) {
                        for (size_t c = 1; c < nvz - 1; c++) {
                            const double distribution = maxwellianDistribution(vx[a], vy[b], vz[c]);
                            const double value = ftime[{i, j, k, a, b, c}];
                            f[{i, j, k, a, b, c}] = value + wc * (distribution * n[{i, j, k}] - value) * dt;
                        }
                    }
                }
            }
        }
    }
}

void Scheme::poisson() {
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                double potential = potentialDebye(rx[i], ry[j], rz[k]);
                for (size_t i_ = 0; i_ < nx; i_++) {
                    for (size_t j_ = 0; j_ < ny; j_++) {
                        for (size_t k_ = 0; k_ < nz; k_++) {
                            double r = distance(rx[i], ry[j], rz[k], rx[i_], ry[j_], rz[k_]);
                            if (r != 0) {
                                potential += exp(-r / rde) * (n[{i_, j_, k_}] - 1.) * hx * hy * hz / (r * 4 * M_PI);
                            }
                        }
                    }
                }
                fi[{i, j, k}] = potential;
            }
        }
    }
}

void Scheme::schemeStep() {
    poisson();
    force();
    kinetic();
    density();
    calculateFullCharge();
    t = t + dt;
}

void Scheme::calculateFullCharge() {
    fullCharge = 0;
#pragma omp parallel for collapse(3)
    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nz; k++) {
                fullCharge += (n[{i, j, k}] - 1.0) * hx * hy * hz;
            }
        }
    }
}

void Scheme::printCurrentTime() {
    cout << "Current dimmensionless time: " << t << endl;
}

void Scheme::printFullCharge() {
    cout << "Current full charge in the computational box: " << fullCharge << endl;
}

void Scheme::writeDensityFile(const string &data = "./") {
    stringstream filename;
    filename << data << "/density_t" << t << ".dat";
    ofstream file;
    file.open(filename.str().c_str());
    file << nx << "\t" << ny << "\t" << nz << "\n";
    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                file << n[{i, j, k}] - 1.0 << "\n";
            }
        }
    }
    file.close();
}

void Scheme::writePotentialFile(const string &data = "./") {
    stringstream filename;
    filename << data << "/potential_t" << t << ".dat";
    ofstream file;
    file.open(filename.str().c_str());
    file << nx << "\t" << ny << "\t" << nz << "\n";
    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                file << fi[{i, j, k}] << "\n";
            }
        }
    }
    file.close();
}

void Scheme::InitialLogOut() {
    cout << "#########################################################" << endl;
    cout << "##******* Output provided by numerical scheme *********##" << endl;
    cout << "#########################################################" << endl;
    cout << "Particle charge: " << q << endl;
    cout << "External electricity field: " << Eext << endl;
    cout << "Ion flow velocity: " << vfl << endl;
    cout << "Ion-neutral collision frequency: " << wc << endl;
    cout << "Debye electron radious: " << rde << endl;
    cout << "Size of the computaional box:" << endl;
    cout << "nx:  " << nx << endl;
    cout << "ny:  " << ny << endl;
    cout << "nz:  " << nz << endl;
    cout << "nvx: " << nvx << endl;
    cout << "nvy: " << nvy << endl;
    cout << "nvz: " << nvz << endl;
    cout << "Sheme integration steps:" << endl;
    cout << "hx:  " << hx << endl;
    cout << "hy:  " << hy << endl;
    cout << "hz:  " << hz << endl;
    cout << "hvx: " << hvx << endl;
    cout << "hvy: " << hvy << endl;
    cout << "hvz: " << hvz << endl;
    cout << "dt:  " << dt << endl;
    cout << "Cut velocities in velocity space:" << endl;
    cout << "vmaxyz: " << vmaxyz << endl;
    cout << "vmaxx:  " << vmaxx << endl;
    cout << "vminyz: " << vminyz << endl;
    cout << "vminz: " << vminx << endl;
    cout << endl;
}
