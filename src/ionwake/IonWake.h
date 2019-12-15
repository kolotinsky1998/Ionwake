//
// Created by Igor on 25.10.2019.
//

#include <cstddef>
#include <fstream>

#ifndef OPTIMIZATION_IONWAKE_H
#define OPTIMIZATION_IONWAKE_H

namespace ionwake {
    ///electron charge
    const double ELECTRON_CHARGE = 4.8032e-10;

    ///ion mass
    const double ION_MASS = 6.6464e-23;

    ///Boltzmann constant
    const double BOLTZMANN_K = 1.3806e-16;

    ///maximum value of integration patameter in velocity space
    const int MAX_XI = 15;

    ///discrete step for integration parameter
    const double DELTA_XI = 0.005;

    const double E = 0.0001;

    class IonWake {
    private:
        double *f;
        double *f_time;
        double *f_n;
        
        double *potential;
        double *density;

        const size_t nx;
        const size_t ny;
        const size_t nz;

        size_t coordinate_total_size;
        size_t velocity_total_size;
        size_t total_size;

        const double lx;
        const double ly;
        const double lz;

        size_t nvx;
        size_t nvy;
        size_t nvz;

        double hvx;
        double hvy;
        double hvz;

        double *vx;
        double *vy;
        double *vz;

        const size_t nx_0;
        const size_t ny_0;
        const size_t nz_0;

        double accelerationCoefficientS;
        double accelerationCoefficientL;
        double accelerationCoefficientC;

        double coordinateStepX;
        double coordinateStepY;
        double coordinateStepZ;

        double *acx;
        double *acy;
        double *acz;

        double *flowVelocityX;
        double *flowVelocityY;
        double *flowVelocityZ;

        double *selfConsistentForceFieldX;
        double *selfConsistentForceFieldY;
        double *selfConsistentForceFieldZ;

        double convertPotential;
        double plasmasFrequency;

        double deltaT;
        double currentTime;
        double debayRadius;
        double thermalVelocity;

        double dimensionlessFlowVelocity;
        double dimensionlessIonConcentration;
        double dimensionlessTau;
        double analyticalFlowVelocity;

        double vxcut;
        double vyzcut;

        const double particleCharge;
        const double temperature;
        const double physicalIonConcentration;
        const double physicalTau;
        const double strength;

        void poissonScheme();

        void gradientScheme();

        void coordinatePart();

        void velocityPart();

        void computeDensity();

        void computeFlowVelocity();

        void saveStep();

    public:

        IonWake(size_t nx, size_t ny, size_t nz, size_t nx_0, size_t ny_0, size_t nz_0, double lx, double ly, double lz,
                double particle_charge, double t, double n_0, double tau, double el);

        ~IonWake();

        void nextStep();

        void writeInitialPotential(std::ofstream &of) const;

        void writeDensity(std::ofstream &of) const;

        void writePotential(std::ofstream &of) const;

        void writeVelocity(std::ofstream &of) const;

        void plotDistributionFunctionX(size_t i, size_t j, size_t k, std::ofstream &of) const;

        void plotDistributionFunctionY(size_t i, size_t j, size_t k, std::ofstream &of) const;

        void plotDistributionFunctionZ(size_t i, size_t j, size_t k, std::ofstream &of) const;

        double getTemperature() const;

        double getThermalVelocity() const;

        double getDebayRadius() const;

        double getDeltaT() const;

        double getParticleCharge() const;

        double getConvertPotential() const;

        double getPlasmasFrequency() const;

        double getPhysicalIonConcentration() const;

        double getPhysicalTau() const;

        double getStrength() const;

        double getCoordinateStepX() const;

        double getCoordinateStepY() const;

        double getCoordinateStepZ() const;

        size_t getNx() const;

        size_t getNy() const;

        size_t getNz() const;

        size_t getNvx() const;

        size_t getNvy() const;

        size_t getNvz() const;

        double getHvx() const;

        double getHvy() const;

        double getHvz() const;

        double getCutVelocityAlong() const;

        double getCutVelocityNormal() const;

        double getDimensionlessIonConcentration() const;

        double getDimensionlessTau() const;

        double getAnalyticalFlowVelocity() const;

        double getDimensionlessFlowVelocity() const;

        double getCurrentTime() const;
    };
}


#endif //OPTIMIZATION_IONWAKE_H
