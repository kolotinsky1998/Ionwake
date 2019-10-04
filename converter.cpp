#include"converter.h"

#include <cmath>

using namespace std;

converter::converter(double T_, double tau_, double n_0_, double q_,
                     double El_, double Lx_, double Ly_, double Lz_, int Nx_, int Ny_, int Nz_, int Nx_0_, int Ny_0_,
                     int Nz_0_) :
        T(T_),
        tau(tau_),
        n_0(n_0_),
        q(q_),
        El(El_),
        Lx(Lx_),
        Ly(Ly_),
        Lz(Lz_),
        Nx(Nx_),
        Ny(Ny_),
        Nz(Nz_),
        Nx_0(Nx_0_),
        Ny_0(Ny_0_),
        Nz_0(Nz_0_) {
    pi = 3.1415926535897932;
    e = 4.8032 * pow(10., -10.);
    m = 6.6464 * pow(10., -23.);
    k = 1.3806 * pow(10., -16.);
    vt = sqrt(T * k / m);
    rd = sqrt(k * T / (4. * pi * e * e * n_0));
    wp = vt / rd;

    //this parameters should be visible in the outside of the converter
    as = e * e / (m * vt * vt * rd);
    al = El * e * rd / (m * vt * vt);
    ac = q * e * e / (m * vt * vt * rd);
    potential = e / rd;
    tau_d = tau * vt / rd;
    v_fl_an_d = e * El * tau / (m * vt);
    n_0_d = n_0 * (rd * rd * rd);
    hx = Lx / double(Nx);
    hy = Ly / double(Ny);
    hz = Lz / double(Nz);


}

converter::~converter() = default;


double converter::GetSAxelerationCofficient() const {
    return as;
}


double converter::GetLAxelerationCofficient() const {
    return al;
}


double converter::GetCAxelerationCofficient() const {
    return ac;
}


double converter::ConvertPotential() const {
    return potential;
}


double converter::GetDimensionlessTau() const {
    return tau_d;
}


double converter::GetDimensionlessFlowVelocity() const {
    return v_fl_an_d;
}


double converter::GetDimensionlessIonConcentration() const {
    return n_0_d;
}


double converter::GetCoordinateStepX() const {
    return hx;
}


double converter::GetCoordinateStepY() const {
    return hy;
}


double converter::GetCoordinateStepZ() const {
    return hz;
}


int converter::GetNx() const {
    return Nx;
}


int converter::GetNy() const {
    return Ny;
}


int converter::GetNz() const {
    return Nz;
}

int converter::GetNx_0() const {
    return Nx_0;
}

int converter::GetNy_0() const {
    return Ny_0;
}

int converter::GetNz_0() const {
    return Nz_0;
}

double converter::GetLx() const {
    return Lx;
}

double converter::GetLy() const {
    return Ly;
}

double converter::GetLz() const {
    return Lz;
}


double converter::GetTemperature() const {
    return T;
}


double converter::GetThermalVelocity() const {
    return vt;
}


double converter::GetPlasmasFrequency() const {
    return wp;
}


double converter::GetPhysicalIonConcentration() const {
    return n_0;
}


double converter::GetParticleCharge() const {
    return q;
}


double converter::GetStrength() const {
    return El;
}


double converter::GetPhysicalTau() const {
    return tau;
}


double converter::GetRd() const {
    return rd;
}

double converter::GetElementaryCharge() const {
    return e;
}


double converter::GetEl() const {
    return El;
}

double converter::GetPi() const {
    return pi;
}
