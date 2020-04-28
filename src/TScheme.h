#ifndef OPTIMIZATION_TSCHEME_H
#define OPTIMIZATION_TSCHEME_H


#include <cstddef>
#include <cmath>
#include <algorithm>
#include <ostream>

#include "omp.h"

template<class T>
class TScheme {

private:

    const size_t x;
    const size_t y;
    const size_t z;

    const size_t vx;
    const size_t vy;
    const size_t vz;

    T *const f;
    T *const f_time;

    const T hx;
    const T hy;
    const T hz;

    const T hvx;
    const T hvy;
    const T hvz;

    T *const n;
    T *const fi;

    const T x0;
    const T y0;
    const T z0;

    T *const ax;
    T *const ay;
    T *const az;

    const T q;
    const T rde;
    const T Eext;
    const T dt;
    const T vminx;
    const T vminyz;
    const T wc;

    TScheme(const size_t x, const size_t y, const size_t z,
            const size_t vx, const size_t vy, const size_t vz,
            const T hx, const T hy, const T hz,
            const T hvx, const T hvy, const T hvz,
            const T x0, const T y0, const T z0,
            const T rde, const T q, const T eext, const T dt, const T vminx,
            const T vminyz, const T wc, T *const f, T *const f_time, T *const n) :
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
            ax(new T[x * y * z]{0}),
            ay(new T[x * y * z]{0}),
            az(new T[x * y * z]{0}),
            n(n),
            fi(new T[x * y * z]{0}), hy(hy), hx(hx), hz(hz), x0(x0), y0(y0), z0(z0), rde(rde), q(q), Eext(eext),
            dt(dt), vminx(vminx), vminyz(vminyz), wc(wc) {}

    inline T potential_debye(T rx, T ry, T rz) const noexcept {
        const T r = distance(rx, ry, rz, x0, y0, z0);
        const T t = sqrt(hx * hx + hy * hy + hz * hz);
        const T distance = r < t ? t : r;
        return -q * exp(-distance / rde) / distance;
    }

    static inline T distance(T x1, T y1, T z1, T x2, T y2, T z2) noexcept {
        const T x = x1 - x2;
        const T y = y1 - y2;
        const T z = z1 - z2;
        return sqrt(x * x + y * y + z * z);
    }

    void compute_poisson() noexcept {
#pragma omp parallel for collapse(3)
        for (size_t i = 0; i < x; ++i) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t k = 0; k < z; ++k) {
                    const T x1 = i * hx;
                    const T y1 = j * hy;
                    const T z1 = k * hz;
                    T potential = potential_debye(x1, y1, z1);

                    for (size_t i2 = 0; i2 < x; ++i2) {
                        for (size_t j2 = 0; j2 < y; ++j2) {
                            for (size_t k2 = 0; k2 < z; ++k2) {
                                if (i != i2 || j != j2 || k != k2) {
                                    const T r = distance(x1, y1, z1, i2 * hx, j2 * hy, k2 * hz);
                                    potential += exp(-r / rde) * (n[i2 * y * z + j2 * z + k2] - 1.) * hx * hy * hz /
                                                 (r * 4 * M_PI);
                                }
                            }
                        }
                    }

                    fi[i * y * z + j * z + k] = potential;
                }
            }
        }
    }

    void compute_force() noexcept {
#pragma omp parallel for collapse(3)
        for (size_t i = 1; i < x - 1; ++i) {
            for (size_t j = 1; j < y - 1; ++j) {
                for (size_t k = 1; k < z - 1; ++k) {
                    ax[i * y * z + j * z + k] =
                            -(fi[(i + 1) * y * z + j * z + k] - fi[(i - 1) * y * z + j * z + k]) / (2.0 * hx) + Eext;
                    ay[i * y * z + j * z + k] =
                            -(fi[i * y * z + (j + 1) * z + k] - fi[i * y * z + (j - 1) * z + k]) / (2.0 * hy);
                    az[i * y * z + j * z + k] =
                            -(fi[i * y * z + j * z + (k + 1)] - fi[i * y * z + j * z + (k - 1)]) / (2.0 * hz);
                }
            }
        }
    }

    static inline T maxwellianDistribution(T vx, T vy, T vz) {
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
                        const T vx_a = vminx + a * hvx;
                        if (vx_a > 0) {
                            for (size_t b = 1; b < vy - 1; ++b) {
                                for (size_t c = 1; c < vz - 1; ++c) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_prev = i_prev * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const T diff = f_time[index] - f_time[index_prev];
                                    f[index] = f_time[index] - vx_a * diff * dt / hx;
                                }
                            }
                        } else {
                            for (size_t b = 1; b < vy - 1; ++b) {
                                for (size_t c = 1; c < vz - 1; ++c) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_next = i_next * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const T diff = f_time[index_next] - f_time[index];
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
                            const T vy_b = vminyz + b * hvy;
                            if (vy_b > 0) {
                                for (size_t c = 1; c < vz - 1; ++c) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_prev = i * d1 + j_prev * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const T diff = f[index] - f[index_prev];
                                    f_time[index] = f[index] - vy_b * diff * dt / hy;
                                }
                            } else {
                                for (size_t c = 1; c < vz - 1; ++c) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_next = i * d1 + j_next * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const T diff = f[index_next] - f[index];
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
                                const T vz_c = vminyz + c * hvz;
                                if (vz_c > 0) {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_prev = i * d1 + j * d2 + k_prev * d3 + a * d4 + b * d5 + c;
                                    const T diff = f_time[index] - f_time[index_prev];
                                    f[index] = f_time[index] - vz_c * diff * dt / hz;
                                } else {
                                    const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                    const size_t index_next = i * d1 + j * d2 + k_next * d3 + a * d4 + b * d5 + c;
                                    const T diff = f_time[index_next] - f_time[index];
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
                    const T ax_ijk = ax[i * y * z + j * z + k];
                    const size_t a_shift_1 = ax_ijk > 0 ? 0 : 1;
                    const size_t a_shift_2 = ax_ijk > 0 ? 1 : 0;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                const size_t index1 = index + a_shift_1 * d4;
                                const size_t index2 = index - a_shift_2 * d4;
                                const T diff = f[index1] - f[index2];
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
                    const T ay_ijk = ay[i * y * z + j * z + k];
                    const size_t b_shift_1 = ay_ijk > 0 ? 0 : 1;
                    const size_t b_shift_2 = ay_ijk > 0 ? 1 : 0;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                const size_t index1 = index + b_shift_1 * d5;
                                const size_t index2 = index - b_shift_2 * d5;
                                const T diff = f_time[index1] - f_time[index2];
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
                    T az_ijk = az[i * y * z + j * z + k];
                    const size_t c_shift_1 = az_ijk > 0 ? 0 : 1;
                    const size_t c_shift_2 = az_ijk > 0 ? 1 : 0;

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                const size_t index1 = index + c_shift_1;
                                const size_t index2 = index - c_shift_2;
                                const T diff = f[index1] - f[index2];
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
                    T n_ijk = n[i * y * z + j * z + k];

                    for (size_t a = 1; a < vx - 1; ++a) {
                        for (size_t b = 1; b < vy - 1; ++b) {
                            for (size_t c = 1; c < vz - 1; ++c) {
                                const T vx_a = vminx + a * hvx;
                                const T vy_b = vminyz + b * hvy;
                                const T vz_c = vminyz + c * hvz;
                                const size_t index = i * d1 + j * d2 + k * d3 + a * d4 + b * d5 + c;
                                const T distribution = maxwellianDistribution(vx_a, vy_b, vz_c);
                                const T value = f_time[index];
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
                    T density = 0.0;
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
        compute_poisson();
        compute_force();
        compute_kinetic();
        compute_density();
    }

    void write_density(std::ostream &of) const {
        of << x << "\t" << y << "\t" << z << "\n";
        for (size_t k = 0; k < z; ++k) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t i = 0; i < x; ++i) {
                    of << n[i * y * z + j * z + k] - 1.0 << "\n";
                }
            }
        }
    }

    void write_potential(std::ostream &of) const {
        of << x << "\t" << y << "\t" << z << "\n";
        for (size_t k = 0; k < z; ++k) {
            for (size_t j = 0; j < y; ++j) {
                for (size_t i = 0; i < x; ++i) {
                    of << fi[i * y * z + j * z + k] << "\n";
                }
            }
        }
    }

    T get_full_charge() const {
        T full_charge = 0;
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

    T get_dt() const {
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
    }

};


template<class T>
class TScheme<T>::TBuilder {

    ///Physical parameters of the system in physical units
    T Te; ///electron temperature
    T Ti; ///ion temperature
    T ni; ///ion concentration
    T q; ///charge of the dust particle
    T Eext; ///external electricity field
    T wc; ///ion-neutral collisions frequency
    T mi; ///ion mass

    /// Calculation of initial distribution function
    T hxi; ///integration parameter
    T ximax; ///high limit of integration

    T vminyz, vmaxyz;
    T x0, y0, z0;
    T hvx, hvy, hvz;
    T lx, ly, lz;
    size_t nx, ny, nz;

    inline static T initialDisrtibutionFunction(
            T vx,
            T vy,
            T vz,
            T hxi,
            T ximax,
            T vfl
    ) {
        T distribution = 0;
        T xi = 0;
        while (xi < ximax) {
            distribution += exp(-xi) * exp(-0.5 * (vx - xi * vfl) * (vx - xi * vfl)) * exp(-0.5 * vy * vy) *
                            exp(-0.5 * vz * vz) * hxi;
            xi = xi + hxi;
        }
        distribution = distribution / pow(2 * M_PI, 1.5);
        return distribution;
    }

    inline static T vmaxxCompute(T vmaxyz, T hvx, T hxi, T ximax, T vfl) {
        T distribution;
        const T maxvell = exp(-vmaxyz * vmaxyz * 0.5);
        T vmax = 0;
        do {
            vmax = vmax + hvx;
            distribution = 0;
            T xi = 0;
            while (xi < ximax) {
                distribution += exp(-xi) * exp(-0.5 * (vmax - xi * vfl) * (vmax - xi * vfl)) * hxi;
                xi = xi + hxi;
            }
        } while (distribution > maxvell);
        return vmax;
    }

public:

    TBuilder &set_integration_parameters(T ximax, T step) {
        this->ximax = ximax;
        this->hxi = step;
        return *this;
    }

    TBuilder &set_electron_temperature(T te) {
        this->Te = te;
        return *this;
    }

    TBuilder &set_ion_temperature(T ti) {
        this->Ti = ti;
        return *this;
    }

    TBuilder &set_ion_concentration(T ni) {
        this->ni = ni;
        return *this;
    }

    TBuilder &set_dust_particle_charge(T q) {
        this->q = q;
        return *this;
    }

    TBuilder &set_external_electricity_field(T Eext) {
        this->Eext = Eext;
        return *this;
    }

    TBuilder &set_ion_neutral_collisions_frequency(T wc) {
        this->wc = wc;
        return *this;
    }

    TBuilder &set_ion_mass(T mi) {
        this->mi = mi;
        return *this;
    }

    TBuilder &set_particle_position(T x0, T y0, T z0) noexcept {
        this->x0 = x0;
        this->y0 = y0;
        this->z0 = z0;
        return *this;
    }

    TBuilder &set_velocity_step(T hvx, T hvy, T hvz) noexcept {
        this->hvx = hvx;
        this->hvy = hvy;
        this->hvz = hvz;
        return *this;
    }

    TBuilder &set_box_size(T lx, T ly, T lz) noexcept {
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

    TBuilder &set_velocity_bound_for_yz(T min_yz, T max_yz) noexcept {
        this->vminyz = min_yz;
        this->vmaxyz = max_yz;
        return *this;
    }

    TScheme build() const {
        ///Fundamental constants
        constexpr T E = 4.8032e-10; ///electron charge
        constexpr T K_B = 1.3806e-16; ///Boltzmann constant

        const T rdi = sqrt(K_B * Ti / (4.0 * M_PI * E * E * ni));
        const T rde = sqrt(K_B * Te / (4.0 * M_PI * E * E * ni));
        const T wpi = sqrt(4.0 * M_PI * E * E * ni / mi);
        const T vT = sqrt(K_B * Ti / mi);
        const T vfl = E * Eext / (mi * wc);

        const T Eext_d = Eext * E * rdi / (K_B * Ti);
        const T wc_d = wc / wpi;
        const T vfl_d = vfl / vT;
        const T rde_d = rde / rdi;
        const T q_d = q / (4.0 * M_PI * ni * rdi * rdi * rdi);

        const T vminx = vminyz;
        const T vmaxx = vmaxxCompute(vmaxyz, hvx, hxi, ximax, vfl_d);
        const auto nvx = size_t(std::round((vmaxx - vminx) / hvx));
        const auto nvy = size_t(std::round((vmaxyz - vminyz) / hvy));
        const auto nvz = size_t(std::round((vmaxyz - vminyz) / hvz));

        const T hx = lx / nx;
        const T hy = ly / ny;
        const T hz = lz / nz;

        const T dt = 0.5 * std::fmin(hvx / (Eext_d + (q_d / (hx * hx))), hx / vmaxx);

        auto *const f = new T[nx * ny * nz * nvx * nvy * nvz]{0};
        auto *const f_time = new T[nx * ny * nz * nvx * nvy * nvz]{0};
        auto *const n = new T[nx * ny * nz];
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
                    const T vx_a = vminx + a * hvx;
                    const T vy_b = vminyz + b * hvy;
                    const T vz_c = vminyz + c * hvz;
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

        return TScheme(nx, ny, nz,
                       nvx, nvy, nvz,
                       hx, hy, hz,
                       hvx, hvy, hvz,
                       x0, y0, z0,
                       rde_d, q_d, Eext_d, dt, vminx, vminyz, wc_d, f, f_time, n);
    }
};

#endif //OPTIMIZATION_TSCHEME_H
