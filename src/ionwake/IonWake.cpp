#include "IonWake.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iostream>

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
    for (size_t i = 0; i < nvx; ++i) { vx[i] = i * hvx - vyzcut; }
    for (size_t i = 0; i < nvy; ++i) { vy[i] = i * hvy - vyzcut; }
    for (size_t i = 0; i < nvz; ++i) { vz[i] = i * hvz - vyzcut; }

    coordinate_total_size = nx * ny * nz;
    velocity_total_size = nvx * nvy * nvz;
    total_size = coordinate_total_size * velocity_total_size;

    f_n = new double[velocity_total_size];
    for (size_t total_i = 0; total_i < velocity_total_size; ++total_i) {
        const size_t c = total_i % nvz;
        const size_t b = total_i / nvz % nvy;
        const size_t a = total_i / nvz / nvy;

        const double x = vx[a] * vx[a];
        const double y = vy[b] * vy[b];
        const double z = vz[c] * vz[c];
        f_n[total_i] = std::exp(-0.5 * (x + y + z)) / std::pow(2.0 * M_PI, 1.5);
    }

    f = new double[total_size];
    f_time = new double[total_size];

    for (size_t total_i = 0; total_i < total_size; ++total_i) {
        const size_t c = total_i % nvz;
        const size_t b = total_i / nvz % nvy;
        const size_t a = total_i / nvz / nvy % nvx;
        const size_t k = total_i / nvz / nvy / nvx % nz;
        const size_t j = total_i / nvz / nvy / nvx / nz % ny;
        const size_t i = total_i / nvz / nvy / nvx / nz / ny;

        if (i == 0 && j == 0 && k == 0) {
            double sum = 0.0;
            double xi = 0;
            while (xi < MAX_XI) {
                double vx_a = vx[a] - xi * analyticalFlowVelocity;
                double vy_b = vy[b];
                double vz_c = vz[c];

                sum += exp(-xi) * exp(-0.5 * (vx_a * vx_a + vy_b * vy_b + vz_c * vz_c));
                xi += DELTA_XI;
            }
            f[total_i] = dimensionlessIonConcentration * sum * DELTA_XI / pow(2.0 * M_PI, 1.5);
        }
    }

    for (size_t total_i = 0; total_i < total_size; ++total_i) {
        const size_t c = total_i % nvz;
        const size_t b = total_i / nvz % nvy;
        const size_t a = total_i / nvz / nvy % nvx;

        const size_t velocity_i = a * nvy * nvz + b * nvz + c;
        f[total_i] = f[velocity_i];
    }

    flowVelocityX = new double[coordinate_total_size]{0};
    flowVelocityY = new double[coordinate_total_size]{0};
    flowVelocityZ = new double[coordinate_total_size]{0};
    selfConsistentForceFieldX = new double[coordinate_total_size]{0};
    selfConsistentForceFieldY = new double[coordinate_total_size]{0};
    selfConsistentForceFieldZ = new double[coordinate_total_size]{0};
    potential = new double[coordinate_total_size]{0};
    density = new double[coordinate_total_size];
    std::fill(density, density + coordinate_total_size, dimensionlessIonConcentration);

    acx = new double[coordinate_total_size]{0};
    acy = new double[coordinate_total_size]{0};
    acz = new double[coordinate_total_size]{0};
    for (size_t total_i = 0; total_i < coordinate_total_size; ++total_i) {
        const size_t k = total_i % nz;
        const size_t j = total_i / nz % ny;
        const size_t i = total_i / nz / ny;

        double x_step = i * coordinateStepX - nx_0 * coordinateStepX;
        double y_step = j * coordinateStepY - ny_0 * coordinateStepY;
        double z_step = k * coordinateStepZ - nz_0 * coordinateStepZ;
        double down = pow(pow((x_step + 0.5 * coordinateStepX), 2.) +
                          pow((y_step + 0.5 * coordinateStepY), 2.) +
                          pow((z_step + 0.5 * coordinateStepZ), 2.), 1.5);

        acx[total_i] = -accelerationCoefficientC * (x_step + 0.5 * coordinateStepX) / down;
        acy[total_i] = -accelerationCoefficientC * (y_step + 0.5 * coordinateStepY) / down;
        acz[total_i] = -accelerationCoefficientC * (z_step + 0.5 * coordinateStepZ) / down;
    }
}

void IonWake::nextStep() {

    auto start_2 = std::chrono::high_resolution_clock::now();
    poissonScheme();
    std::cout << "Time taken by poissonScheme: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::high_resolution_clock::now() - start_2).count() / 1000000.0 << " seconds"
              << std::endl;

    auto start_1 = std::chrono::high_resolution_clock::now();
    gradientScheme();
    std::cout << "Time taken by gradientScheme: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::high_resolution_clock::now() - start_1).count() / 1000000.0 << " seconds"
              << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    coordinatePart();
    std::cout << "Time taken by coordinatePart: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::high_resolution_clock::now() - start).count() / 1000000.0 << " seconds" << std::endl;
    auto start2 = std::chrono::high_resolution_clock::now();
    velocityPart();
    std::cout << "Time taken by numerical velocityPart: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::high_resolution_clock::now() - start2).count() / 1000000.0 << " seconds"
              << std::endl;

    auto start3 = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
    for (size_t total_i = 0; total_i < total_size; ++total_i) {
        const size_t c = total_i % nvz;
        const size_t b = total_i / nvz % nvy;
        const size_t a = total_i / nvz / nvy % nvx;

        if (a != 0 && b != 0 && c != 0 && a != nvx - 1 && b != nvy - 1 && c != nvz - 1) {
            const size_t k = total_i / nvz / nvy / nvx % nz;
            const size_t j = total_i / nvz / nvy / nvx / nz % ny;
            const size_t i = total_i / nvz / nvy / nvx / nz / ny;

            const size_t coordinate_i = i * ny * nz + j * nz + k;
            const size_t velocity_i = a * nvy * nvz + b * nvz + c;
            f[total_i] += deltaT * (density[coordinate_i] * f_n[velocity_i] - f[total_i]) / dimensionlessTau;
        }
    }
    std::cout << "Time taken by ceter: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::high_resolution_clock::now() - start3).count() / 1000000.0 << " seconds"
              << std::endl;

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
        for (size_t total_i = 0; total_i < coordinate_total_size; ++total_i) {
            const size_t k = total_i % nz;
            const size_t j = total_i / nz % ny;
            const size_t i = total_i / nz / ny;

            if (i != 0 && j != 0 && k != 0 && i != nx - 1 && j != ny - 1 && k != nz - 1) {
                const double f1 = potential[total_i];

                const double double_f1 = 2.0 * f1;
                const double fi = (potential[total_i + nz * ny] - double_f1 + potential[total_i - nz * ny]) / owx;
                const double fj = (potential[total_i + nz] - double_f1 + potential[total_i - nz]) / owy;
                const double fk = (potential[total_i + 1] - double_f1 + potential[total_i - 1]) / owz;

                const double first = density[total_i] / dimensionlessIonConcentration - 1.0;
                const double second = particleCharge / (dimensionlessIonConcentration * hx * hy * hz * nz * ny * nz);
                const double third = (i == nx_0 && j == ny_0 && k == nz_0 ? particleCharge : 0.0);
                const double strange =
                        4. * M_PI * (first + second - third) / (dimensionlessIonConcentration * hx * hy * hz);

                potential[total_i] += (fi + fj + fk + strange) * tau;

                if (potential[total_i] == 0.0 && std::abs((potential[total_i] - f1) / potential[total_i]) > E) {
                    proceed = true;
                }
            }
        }
    }
}

void IonWake::gradientScheme() {
    const double x = 2.0 * coordinateStepX;
    const double y = 2.0 * coordinateStepY;
    const double z = 2.0 * coordinateStepZ;
    double x_concentration = dimensionlessIonConcentration / x;
    double y_concentration = dimensionlessIonConcentration / y;
    double z_concentration = dimensionlessIonConcentration / z;

    const size_t nx_nz = nz * ny;

#pragma omp parallel for
    for (size_t total_i = 0; total_i < coordinate_total_size; ++total_i) {
        const size_t k = total_i % nz;
        const size_t j = total_i / nz % ny;
        const size_t i = total_i / nz / ny;

        if (i != 0 && j != 0 && k != 0 && i != nx - 1 && j != ny - 1 && k != nz - 1) {
            selfConsistentForceFieldX[total_i] =
                    -(potential[total_i + nx_nz] - potential[total_i - nx_nz]) * x_concentration;
            selfConsistentForceFieldY[total_i] =
                    -(potential[total_i + nz] - potential[total_i - nz]) * y_concentration;
            selfConsistentForceFieldZ[total_i] =
                    -(potential[total_i + 1] - potential[total_i - 1]) * z_concentration;
        }
    }

}

void IonWake::coordinatePart() {
    const double delta_x_step = deltaT / coordinateStepX;
    const double delta_y_step = deltaT / coordinateStepY;
    const double delta_z_step = deltaT / coordinateStepZ;

    const size_t z_shift = nvz * nvy * nvx;
    const size_t y_shift = z_shift * nz;
    const size_t x_shift = y_shift * ny;

    saveStep();
#pragma omp parallel for
    for (size_t velocity_total_i = 0; velocity_total_i < velocity_total_size; ++velocity_total_i) {
        const size_t c = velocity_total_i % nvz;
        const size_t b = velocity_total_i / nvz % nvy;
        const size_t a = velocity_total_i / nvz / nvy;

        if (a != 0 && b != 0 && c != 0 && a != nvx - 1 && b != nvy - 1 && c != nvz - 1) {
            for (size_t coordinate_total_i = 0; coordinate_total_i < coordinate_total_size; ++coordinate_total_i) {
                const size_t total_shift = coordinate_total_i * velocity_total_size + velocity_total_i;

                const double vx_a = vx[a];
                const double shift = vx_a > 0.0
                                     ? f_time[total_shift] -
                                       f_time[(total_shift + (total_shift - x_shift)) % total_shift]
                                     : f_time[(total_shift + x_shift) % total_shift] - f_time[total_shift];

                f[total_shift] = f_time[total_shift] - (delta_x_step * vx_a) * shift;
            }
        }
    }

    saveStep();
#pragma omp parallel for
    for (size_t velocity_total_i = 0; velocity_total_i < velocity_total_size; ++velocity_total_i) {
        const size_t c = velocity_total_i % nvz;
        const size_t b = velocity_total_i / nvz % nvy;
        const size_t a = velocity_total_i / nvz / nvy;

        if (a != 0 && b != 0 && c != 0 && a != nvx - 1 && b != nvy - 1 && c != nvz - 1) {
            for (size_t coordinate_total_i = 0; coordinate_total_i < coordinate_total_size; ++coordinate_total_i) {
                const size_t total_shift = coordinate_total_i * velocity_total_size + velocity_total_i;

                const double vy_b = vy[b];
                const double shift = vy_b > 0.0
                                     ? f_time[total_shift] - f_time[(total_shift + (total_size - y_shift)) % total_size]
                                     : f_time[(total_shift + y_shift) % total_size] - f_time[total_shift];

                f[total_shift] = f_time[total_shift] - (delta_y_step * vy_b) * shift;
            }
        }
    }

    saveStep();
#pragma omp parallel for
    for (size_t velocity_total_i = 0; velocity_total_i < velocity_total_size; ++velocity_total_i) {
        const size_t c = velocity_total_i % nvz;
        const size_t b = velocity_total_i / nvz % nvy;
        const size_t a = velocity_total_i / nvz / nvy;

        if (a != 0 && b != 0 && c != 0 && a != nvx - 1 && b != nvy - 1 && c != nvz - 1) {
            for (size_t coordinate_total_i = 0; coordinate_total_i < coordinate_total_size; ++coordinate_total_i) {
                const size_t total_shift = coordinate_total_i * velocity_total_size + velocity_total_i;

                const double vz_c = vz[c];
                const double shift = vz_c > 0.0
                                     ? f_time[total_shift] - f_time[(total_shift + (total_size - z_shift)) % total_size]
                                     : f_time[(total_shift + z_shift) % total_size] - f_time[total_shift];

                f[total_shift] = f_time[total_shift] - (delta_z_step * vz_c) * shift;
            }
        }
    }
}

void IonWake::velocityPart() {
    const double delta_t_hvx = deltaT / hvx;
    const double delta_t_hvy = deltaT / hvy;
    const double delta_t_hvz = deltaT / hvz;

    const size_t ny_nz = ny * nz;
    const size_t nvx_nvy = nvz * nvy;

    saveStep();
#pragma omp parallel for
    for (size_t velocity_total_i = 0; velocity_total_i < velocity_total_size; ++velocity_total_i) {
        const size_t c = velocity_total_i % nvz;
        const size_t b = velocity_total_i / nvz % nvy;
        const size_t a = velocity_total_i / nvz / nvy;

        if (a != 0 && b != 0 && c != 0 && a != nvx - 1 && b != nvy - 1 && c != nvz - 1) {
            for (size_t coordinate_total_i = 0; coordinate_total_i < coordinate_total_size; ++coordinate_total_i) {
                size_t total = coordinate_total_i * velocity_total_size + velocity_total_i;

                const double fx = accelerationCoefficientS * selfConsistentForceFieldX[coordinate_total_i] +
                                  accelerationCoefficientL + acx[coordinate_total_i];
                const double shift = fx > 0.0
                                     ? f_time[total] - f_time[total - nvx_nvy]
                                     : f_time[total + nvx_nvy] - f_time[total];

                f[total] = f_time[total] - (delta_t_hvx * fx) * shift;
            }
        }
    }

    saveStep();
#pragma omp parallel for
    for (size_t velocity_total_i = 0; velocity_total_i < velocity_total_size; ++velocity_total_i) {
        const size_t c = velocity_total_i % nvz;
        const size_t b = velocity_total_i / nvz % nvy;
        const size_t a = velocity_total_i / nvz / nvy;

        if (a != 0 && b != 0 && c != 0 && a != nvx - 1 && b != nvy - 1 && c != nvz - 1) {
            for (size_t coordinate_total_i = 0; coordinate_total_i < coordinate_total_size; ++coordinate_total_i) {
                size_t total = coordinate_total_i * velocity_total_size + velocity_total_i;

                const double fy = accelerationCoefficientS * selfConsistentForceFieldY[coordinate_total_i] +
                                  acy[coordinate_total_i];
                const double shift = fy > 0.0
                                     ? f_time[total] - f_time[total - nvz]
                                     : f_time[total + nvz] - f_time[total];

                f[total] = f_time[total] - (delta_t_hvy * fy) * shift;
            }
        }
    }

    saveStep();
#pragma omp parallel for
    for (size_t velocity_total_i = 0; velocity_total_i < velocity_total_size; ++velocity_total_i) {
        const size_t c = velocity_total_i % nvz;
        const size_t b = velocity_total_i / nvz % nvy;
        const size_t a = velocity_total_i / nvz / nvy;

        if (a != 0 && b != 0 && c != 0 && a != nvx - 1 && b != nvy - 1 && c != nvz - 1) {
            for (size_t coordinate_total_i = 0; coordinate_total_i < coordinate_total_size; ++coordinate_total_i) {
                size_t total = coordinate_total_i * velocity_total_size + velocity_total_i;

                const double fz = accelerationCoefficientS * selfConsistentForceFieldZ[coordinate_total_i] +
                                  acz[coordinate_total_i];
                const double shift = fz > 0.0
                                     ? f_time[total] - f_time[total - 1]
                                     : f_time[total + 1] - f_time[total];

                f[total] = f_time[total] - (delta_t_hvz * fz) * shift;
            }
        }
    }
}

void IonWake::computeDensity() {
    const double m = hvx * hvy * hvz;
#pragma omp parallel for
    for (size_t total_i = 0; total_i < coordinate_total_size; ++total_i) {
        const size_t velocity_shift = total_i * velocity_total_size;

        double sum = 0.0;
        for (size_t total_j = 0; total_j < velocity_total_size; ++total_j) {
            sum += f[velocity_shift + total_j];
        }
        density[total_i] = sum * m;
    }
}

void IonWake::computeFlowVelocity() {
    const double m = hvx * hvy * hvz;
#pragma omp parallel for
    for (size_t total_i = 0; total_i < coordinate_total_size; ++total_i) {
        size_t velocity_shift = total_i * velocity_total_size;
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;
        for (size_t total_j = 0; total_j < velocity_total_size; ++total_j) {
            const size_t c = total_i % nvz;
            const size_t b = total_i / nvz % nvy;
            const size_t a = total_i / nvz / nvy;

            sum_x += vx[a] * f[velocity_shift + total_j];
            sum_y += vy[b] * f[velocity_shift + total_j];
            sum_z += vz[c] * f[velocity_shift + total_j];
        }

        flowVelocityX[total_i] = sum_x * m;
        flowVelocityY[total_i] = sum_y * m;
        flowVelocityZ[total_i] = sum_z * m;
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
                const size_t coordinate_i = i * ny * nz + j * nz + k;
                of << density[coordinate_i] * pow(1.0 / debayRadius, 3) << "\n";
            }
        }
    }
}

void IonWake::writePotential(std::ofstream &of) const {
    of << nx << "\t" << ny << "\t" << nz << "\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                const size_t coordinate_i = i * ny * nz + j * nz + k;
                of << dimensionlessIonConcentration * convertPotential * potential[coordinate_i] << "\n";
            }
        }
    }
}

void IonWake::writeVelocity(std::ofstream &of) const {
    of << nx << "\t" << ny << "\t" << nz << "\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                const size_t coordinate_i = i * ny * nz + j * nz + k;
                of << flowVelocityX[coordinate_i] << "\t" << flowVelocityY[coordinate_i] << "\t"
                   << flowVelocityZ[coordinate_i] << "\n";
            }
        }
    }
}

void IonWake::plotDistributionFunctionX(size_t i, size_t j, size_t k, std::ofstream &of) const {
    const size_t coordinate_i = i * ny * nz + j * nz + k;
    const size_t shift = coordinate_i * velocity_total_size;

    for (size_t a = 0; a < nvx; ++a) {
        const size_t velocity_i = a * nvy * nvz + (nvy / 2) * nvz + (nvz / 2);
        of << vx[a] << "\t" << f[shift + velocity_i] << "\n";
    }
}

void IonWake::plotDistributionFunctionY(size_t i, size_t j, size_t k, std::ofstream &of) const {
    const size_t coordinate_i = i * ny * nz + j * nz + k;
    const size_t shift = coordinate_i * velocity_total_size;

    for (size_t b = 0; b < nvy; ++b) {
        const size_t velocity_i = (nvx / 2) * nvy * nvz + b * nvz + (nvz / 2);
        of << vy[b] << "\t" << f[shift + velocity_i] << "\n";
    }
}

void IonWake::plotDistributionFunctionZ(size_t i, size_t j, size_t k, std::ofstream &of) const {
    const size_t coordinate_i = i * ny * nz + j * nz + k;
    const size_t shift = coordinate_i * velocity_total_size;

    for (size_t c = 0; c < nvz; ++c) {
        const size_t velocity_i = (nvx / 2) * nvy * nvz + (nvy / 2) * nvz + c;
        of << vz[c] << "\t" << f[shift + velocity_i] << "\n";
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
    delete[] f;
    delete[] f_time;
    delete[] f_n;
}
