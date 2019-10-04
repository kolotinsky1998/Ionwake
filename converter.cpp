#include"converter.h"

#include <cmath>
#include "constants.h"

using namespace std;

converter::converter(double T_, double tau_, double n_0_, double q_,
                     double El_, double Lx_, double Ly_, double Lz_, int Nx_0_, int Ny_0_, int Nz_0_) :
        T(T_),
        tau(tau_),
        n_0(n_0_),
        q(q_),
        El(El_),
        Lx(Lx_),
        Ly(Ly_),
        Lz(Lz_),
        Nx_0(Nx_0_),
        Ny_0(Ny_0_),
        Nz_0(Nz_0_) {

    vt = sqrt(T * K / M);
    rd = sqrt(K * T / (4.0 * M_PI * E * E * n_0));
    wp = vt / rd;
    //this parameters should be visible in the outside of the converter
    as = E * E / (M * vt * vt * rd);
    al = El * E * rd / (M * vt * vt);
    ac = q * E * E / (M * vt * vt * rd);
    potential = E / rd;
    tau_d = tau * vt / rd;
    v_fl_an_d = E * El * tau / (M * vt);
    n_0_d = n_0 * (rd * rd * rd);
    hx = Lx / (double) Nx;
    hy = Ly / (double) Ny;
    hz = Lz / (double) Nz;
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

double converter::GetEl() const {
    return El;
}
