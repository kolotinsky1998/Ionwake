#ifndef OPTIMIZATION_TSCHEME_H
#define OPTIMIZATION_TSCHEME_H


#include <cstddef>
#include <cmath>
#include <algorithm>
#include <ostream>

#include "omp.h"

///Fundamental constants
constexpr double E = 4.8032e-10; ///electron charge
constexpr double K_B = 1.3806e-16; ///Boltzmann constant

class TScheme final {

private:

    const size_t x;
    const size_t y;
    const size_t z;

    const size_t vx;
    const size_t vy;
    const size_t vz;

    double *const f;
    double *const f_time;

    const double hx;
    const double hy;
    const double hz;

    const double hvx;
    const double hvy;
    const double hvz;

    double *const n;
    double *const fi;
    double *const precomputed_potential_debye;

    double *const ax;
    double *const ay;
    double *const az;

    const double q;
    const double rde;
    const double Eext;
    const double dt;
    const double vminx;
    const double vminyz;
    const double wc;

    TScheme(
            const size_t x, const size_t y, const size_t z,
            const size_t vx, const size_t vy, const size_t vz,
            const double hx, const double hy, const double hz,
            const double hvx, const double hvy, const double hvz,
            const double rde, const double q, const double eext, const double dt, const double vminx,
            const double vminyz, const double wc, double *const f, double *const f_time, double *const n,
            double *const precomputed_potential_debye
    ) :
            x(x),
            y(y),
            z(z),
            vx(vx),
            vy(vy),
            vz(vz),
            hvx(hvx),
            hvy(hvy),
            hvz(hvz),
            f(f),
            f_time(f_time),
            precomputed_potential_debye(precomputed_potential_debye),
            ax(new double[x * y * z]{0}),
            ay(new double[x * y * z]{0}),
            az(new double[x * y * z]{0}),
            n(n),
            fi(new double[x * y * z]{0}), hy(hy), hx(hx), hz(hz), rde(rde), q(q), Eext(eext),
            dt(dt), vminx(vminx), vminyz(vminyz), wc(wc) {}

    static inline double distance(double x1, double y1, double z1, double x2, double y2, double z2) noexcept {
        const double x = x1 - x2;
        const double y = y1 - y2;
        const double z = z1 - z2;
        return sqrt(x * x + y * y + z * z);
    }

    static inline void compute_poisson(
            const size_t x, const size_t y, const size_t z,
            const double hx, const double hy, const double hz,
            const double rde,
            const double *const precomputed_potential_debye,
            const double *const n,
            double *const fi
    ) noexcept {
        std::copy_n(precomputed_potential_debye, x * y * z, fi);
#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    const double x1 = i * hx;
                    const double y1 = j * hy;
                    const double z1 = k * hz;
                    double potential = 0;

                    for (size_t i2 = 0; i2 < x; ++i2) {
                        for (size_t j2 = 0; j2 < y; ++j2) {
                            for (size_t k2 = 0; k2 < z; ++k2) {
                                if (i == i2 && j == j2 && k == k2) continue;

                                const double r = distance(x1, y1, z1, i2 * hx, j2 * hy, k2 * hz);
                                potential += exp(-r / rde) * (n[i2 * y * z + j2 * z + k2] - 1.) * hx * hy * hz /
                                             (r * 4 * M_PI);
                            }
                        }
                    }

                    fi[i * y * z + j * z + k] += potential;
                }
            }
        }
    }

    static inline void compute_force(
            const size_t x, const size_t y, const size_t z,
            const double hx, const double hy, const double hz,
            const double Eext,
            const double *const fi,
            double *const ax,
            double *const ay,
            double *const az
    ) noexcept {
#pragma omp parallel for collapse(3)
        for (size_t i = 1; i < x - 1; ++i) {
            for (size_t j = 1; j < y - 1; ++j) {
                for (size_t k = 1; k < z - 1; ++k) {
                    const size_t ind = i * y * z + j * z + k;
                    ax[ind] = -(fi[(i + 1) * y * z + j * z + k] - fi[(i - 1) * y * z + j * z + k]) / (2.0 * hx) + Eext;
                    ay[ind] = -(fi[i * y * z + (j + 1) * z + k] - fi[i * y * z + (j - 1) * z + k]) / (2.0 * hy);
                    az[ind] = -(fi[i * y * z + j * z + (k + 1)] - fi[i * y * z + j * z + (k - 1)]) / (2.0 * hz);
                }
            }
        }
    }

    static inline double maxwellianDistribution(double vx, double vy, double vz) noexcept {
        return exp(-0.5 * vx * vx) * exp(-0.5 * vy * vy) * exp(-0.5 * vz * vz) / pow(2 * M_PI, 1.5);
    }

    void compute_kinetic() noexcept {
        const size_t d5 = vz;
        const size_t d4 = vy * d5;
        const size_t d3 = vx * d4;
        const size_t d2 = z * d3;
        const size_t d1 = y * d2;

#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    const size_t i_prev = i == 0 ? x - 1 : i - 1;
                    const size_t i_next = i == x - 1 ? 0 : i + 1;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        const double vx_a = vminx + a * hvx;
                        if (vx_a > 0) {
                            for (size_t b = 1; b < vy - 1; ++b) {
                                for (size_t c = 1; c < vz - 1; ++c) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_prev = i_prev * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const double diff = f_time[index] - f_time[index_prev];
                                    f[index] = f_time[index] - vx_a * diff * dt / hx;
                                }
                            }
                        } else {
                            for (size_t b = 1; b < vy - 1; ++b) {
                                for (size_t c = 1; c < vz - 1; ++c) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_next = i_next * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const double diff = f_time[index_next] - f_time[index];
                                    f[index] = f_time[index] - vx_a * diff * dt / hx;
                                }
                            }
                        }
                    }
                }
            }
        }

#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    const size_t j_prev = j == 0 ? y - 1 : j - 1;
                    const size_t j_next = j == y - 1 ? 0 : j + 1;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            const double vy_b = vminyz + b * hvy;
                            if (vy_b > 0) {
                                for (size_t c = 1; c < vz - 1; ++c) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_prev = i * d1 + j_prev * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const double diff = f[index] - f[index_prev];
                                    f_time[index] = f[index] - vy_b * diff * dt / hy;
                                }
                            } else {
                                for (size_t c = 1; c < vz - 1; ++c) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_next = i * d1 + j_next * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const double diff = f[index_next] - f[index];
                                    f_time[index] = f[index] - vy_b * diff * dt / hy;
                                }
                            }
                        }
                    }
                }
            }
        }

#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    const size_t k_prev = k == 0 ? z - 1 : k - 1;
                    const size_t k_next = k == z - 1 ? 0 : k + 1;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const double vz_c = vminyz + c * hvz;
                                if (vz_c > 0) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_prev = i * d1 + j * d2 + k_prev * d3 + a * d4 + b * d5 + c;
                                    const double diff = f_time[index] - f_time[index_prev];
                                    f[index] = f_time[index] - vz_c * diff * dt / hz;
                                } else {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_next = i * d1 + j * d2 + k_next * d3 + a * d4 + b * d5 + c;
                                    const double diff = f_time[index_next] - f_time[index];
                                    f[index] = f_time[index] - vz_c * diff * dt / hz;
                                }
                            }
                        }
                    }
                }
            }
        }


#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    const double ax_ijk = ax[i * y * z + j * z + k];
                    const size_t a_shift_1 = ax_ijk > 0 ? 0 : 1;
                    const size_t a_shift_2 = ax_ijk > 0 ? 1 : 0;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                const size_t index1 = index + a_shift_1 * d4;
                                const size_t index2 = index - a_shift_2 * d4;
                                const double diff = f[index1] - f[index2];
                                f_time[index] = f[index] - ax_ijk * diff * (dt / hvx);
                            }
                        }
                    }
                }
            }
        }
#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    const double ay_ijk = ay[i * y * z + j * z + k];
                    const size_t b_shift_1 = ay_ijk > 0 ? 0 : 1;
                    const size_t b_shift_2 = ay_ijk > 0 ? 1 : 0;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                const size_t index1 = index + b_shift_1 * d5;
                                const size_t index2 = index - b_shift_2 * d5;
                                const double diff = f_time[index1] - f_time[index2];
                                f[index] = f_time[index] - ay_ijk * diff * (dt / hvy);
                            }
                        }
                    }
                }
            }
        }
#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    double az_ijk = az[i * y * z + j * z + k];
                    const size_t c_shift_1 = az_ijk > 0 ? 0 : 1;
                    const size_t c_shift_2 = az_ijk > 0 ? 1 : 0;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                const size_t index1 = index + c_shift_1;
                                const size_t index2 = index - c_shift_2;
                                const double diff = f[index1] - f[index2];
                                f_time[index] = f[index] - az_ijk * diff * (dt / hvz);
                            }
                        }
                    }
                }
            }
        }

#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    double n_ijk = n[i * y * z + j * z + k];

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const double vx_a = vminx + a * hvx;
                                const double vy_b = vminyz + b * hvy;
                                const double vz_c = vminyz + c * hvz;
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                const double distribution = maxwellianDistribution(vx_a, vy_b, vz_c);
                                const double value = f_time[index];
                                f[index] = value + wc * dt * (distribution * n_ijk - value);
                            }
                        }
                    }
                }
            }
        }
    }

    void compute_density() noexcept {
        const size_t d5 = vz;
        const size_t d4 = vy * d5;
        const size_t d3 = vx * d4;
        const size_t d2 = z * d3;
        const size_t d1 = y * d2;

#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    double density = 0.0;
                    for (size_t a = 0; a < vx; ++a) {
                        for (size_t b = 0; b < vy; ++b) {
                            for (size_t c = 0; c < vz; ++c) {
                                density += f[i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c];
                            }
                        }
                    }
                    n[i * y * z + j * z + k] = density * (hvx * hvy * hvz);
                }
            }
        }
    }

public:

    class TBuilder;

    void next_step() noexcept {
        compute_poisson(x, y, z, hx, hy, hz, rde, precomputed_potential_debye, n, fi);
        compute_force(x, y, z, hx, hy, hz, Eext, fi, ax, ay, az);
        compute_kinetic();
        compute_density();
    }

    void write_density(std::ostream &of) const noexcept {
        of << x << "\t" << y << "\t" << z << "\n";
        for (size_t k = 0; k < z; ++k) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t i = 0; i < x; ++i) {
                    of << n[i * y * z + j * z + k] - 1.0 << "\n";
                }
            }
        }
    }

    void write_potential(std::ostream &of) const noexcept {
        of << x << "\t" << y << "\t" << z << "\n";
        for (size_t k = 0; k < z; ++k) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t i = 0; i < x; ++i) {
                    of << fi[i * y * z + j * z + k] << "\n";
                }
            }
        }
    }

    double get_full_charge() const noexcept {
        double full_charge = 0;
#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; i++) {
            for (size_t j = 0; j < y; j++) {
                for (size_t k = 0; k < z; k++) {
                    full_charge += (n[i * y * z + j * z + k] - 1.0) * hx * hy * hz;
                }
            }
        }
        return full_charge;
    }

    double get_dt() const noexcept {
        return dt;
    }

    ~TScheme() {
        delete[] f;
        delete[] f_time;
        delete[] ax;
        delete[] ay;
        delete[] az;
        delete[] n;
        delete[] fi;
        delete[] precomputed_potential_debye;
    }

};

class TScheme::TBuilder final {

    ///Physical parameters of the system in physical units
    double Te; ///electron temperature
    double Ti; ///ion temperature
    double ni; ///ion concentration
    double q; ///charge of the dust particle
    double Eext; ///external electricity field
    double wc; ///ion-neutral collisions frequency
    double mi; ///ion mass

    /// Calculation of initial distribution function
    double hxi; ///integration parameter
    double ximax; ///high limit of integration

    double vminyz, vmaxyz;
    double x0, y0, z0;
    double hvx, hvy, hvz;
    double lx, ly, lz;
    size_t nx, ny, nz;

    inline static double initialDisrtibutionFunction(
            double vx,
            double vy,
            double vz,
            double hxi,
            double ximax,
            double vfl
    ) {
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

    inline static double vmaxxCompute(double vmaxyz, double hvx, double hxi, double ximax, double vfl) noexcept {
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

    inline double
    potential_debye(
            double hx, double hy, double hz, double rx, double ry, double rz, double q, double rde
    ) const noexcept {
        const double r = distance(rx, ry, rz, x0, y0, z0);
        const double t = sqrt(hx * hx + hy * hy + hz * hz);
        const double distance = r < t ? t : r;
        return -q * exp(-distance / rde) / distance;
    }

public:

    TBuilder &set_integration_parameters(double ximax, double step) noexcept {
        this->ximax = ximax;
        this->hxi = step;
        return *this;
    }

    TBuilder &set_electron_temperature(double te) noexcept {
        this->Te = te;
        return *this;
    }

    TBuilder &set_ion_temperature(double ti) noexcept {
        this->Ti = ti;
        return *this;
    }

    TBuilder &set_ion_concentration(double ni) noexcept {
        this->ni = ni;
        return *this;
    }

    TBuilder &set_dust_particle_charge(double q) noexcept {
        this->q = q;
        return *this;
    }

    TBuilder &set_external_electricity_field(double Eext) noexcept {
        this->Eext = Eext;
        return *this;
    }

    TBuilder &set_ion_neutral_collisions_frequency(double wc) noexcept {
        this->wc = wc;
        return *this;
    }

    TBuilder &set_ion_mass(double mi) noexcept {
        this->mi = mi;
        return *this;
    }

    TBuilder &set_particle_position(double x0, double y0, double z0) noexcept {
        this->x0 = x0;
        this->y0 = y0;
        this->z0 = z0;
        return *this;
    }

    TBuilder &set_velocity_step(double hvx, double hvy, double hvz) noexcept {
        this->hvx = hvx;
        this->hvy = hvy;
        this->hvz = hvz;
        return *this;
    }

    TBuilder &set_box_size(double lx, double ly, double lz) noexcept {
        this->lx = lx;
        this->ly = ly;
        this->lz = lz;
        return *this;
    }

    TBuilder &set_box_split_step_count(size_t nx, size_t ny, size_t nz) noexcept {
        this->nx = nx;
        this->ny = ny;
        this->nz = nz;
        return *this;
    }

    TBuilder &set_velocity_bound_for_yz(double min_yz, double max_yz) noexcept {
        this->vminyz = min_yz;
        this->vmaxyz = max_yz;
        return *this;
    }

    TScheme build() const {
        const double rdi = sqrt(K_B * Ti / (4.0 * M_PI * E * E * ni));
        const double rde = sqrt(K_B * Te / (4.0 * M_PI * E * E * ni));
        const double wpi = sqrt(4.0 * M_PI * E * E * ni / mi);
        const double vT = sqrt(K_B * Ti / mi);
        const double vfl = E * Eext / (mi * wc);

        const double Eext_d = Eext * E * rdi / (K_B * Ti);
        const double wc_d = wc / wpi;
        const double vfl_d = vfl / vT;
        const double rde_d = rde / rdi;
        const double q_d = q / (4.0 * M_PI * ni * rdi * rdi * rdi);

        const double vminx = vminyz;
        const double vmaxx = vmaxxCompute(vmaxyz, hvx, hxi, ximax, vfl_d);
        const auto nvx = size_t(std::round((vmaxx - vminx) / hvx));
        const auto nvy = size_t(std::round((vmaxyz - vminyz) / hvy));
        const auto nvz = size_t(std::round((vmaxyz - vminyz) / hvz));

        const double hx = lx / nx;
        const double hy = ly / ny;
        const double hz = lz / nz;

        const double dt = 0.5 * std::fmin(hvx / (Eext_d + (q_d / (hx * hx))), hx / vmaxx);

        auto *const f = new double[nx * ny * nz * nvx * nvy * nvz]{0};
        auto *const f_time = new double[nx * ny * nz * nvx * nvy * nvz]{0};
        auto *const n = new double[nx * ny * nz];
        std::fill_n(n, nx * ny * nz, 1.0);

        const size_t d5 = nvz;
        const size_t d4 = nvy * d5;
        const size_t d3 = nvx * d4;
        const size_t d2 = nz * d3;
        const size_t d1 = ny * d2;

#pragma omp parallel for collapse(3)
        for (size_t a = 0; a < nvx; a++) {
            for (size_t b = 0; b < nvy; b++) {
                for (size_t c = 0; c < nvz; c++) {
                    const double vx_a = vminx + a * hvx;
                    const double vy_b = vminyz + b * hvy;
                    const double vz_c = vminyz + c * hvz;
                    f[a * d4 + b * d5 + c] = initialDisrtibutionFunction(vx_a, vy_b, vz_c, hxi, ximax, vfl_d);
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
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                f[index] = f[a * d4 + b * d5 + c];
                                f_time[index] = f[a * d4 + b * d5 + c];
                            }
                        }
                    }
                }
            }
        }

        auto *const precomputed_potential_debye = new double[nx * ny * nz];
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t k = 0; k < nz; ++k) {
                    const double x1 = i * hx;
                    const double y1 = j * hy;
                    const double z1 = k * hz;
                    const size_t index = i * ny * nz + j * nz + k;
                    precomputed_potential_debye[index] = potential_debye(hx, hy, hz, x1, y1, z1, q_d, rde_d);
                }
            }
        }

        return TScheme(
                nx, ny, nz,
                nvx, nvy, nvz,
                hx, hy, hz,
                hvx, hvy, hvz,
                rde_d, q_d, Eext_d, dt, vminx, vminyz, wc_d, f, f_time, n, precomputed_potential_debye
        );
    }
};

#endif //OPTIMIZATION_TSCHEME_H
