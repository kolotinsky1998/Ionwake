#ifndef scheme_h
#define scheme_h

#include <string>

#include "converter.hpp"
#include "TArray.h"

using namespace std;

class Scheme {
public:
    Scheme(
            const Converter &converter_, size_t nx_, size_t ny_, size_t nz_, double lx_, double ly_, double lz_,
            int nt_,
            double r0x_, double r0y_, double r0z_, double hvx_, double hvy_, double hvz_, double vminyz_,
            double vmaxyz_, double hxi_, double ximax_
    );

    ~Scheme();

    void schemeStep();

    void writeDensityFile(const string &data);

    void writePotentialFile(const string &data);

    void printCurrentTime() const;

    void printFullCharge() const;

    void InitialLogOut() const; ///write simulation parameters into the log file
private:
    Converter converter;
    /// Coordinate space
    size_t nx, ny, nz; ///amount of cells in coordinate space
    size_t nvx, nvy, nvz; ///amount of cells in velocity space

    double lx, ly, lz; ///size of the computational box
    double hx, hy, hz; ///integration step in coordinate space
    /// Velocity space
    double hvx, hvy, hvz; ///integration step in velocity space
    double vminx, vminyz; ///minimum velocity in velocity space
    double vmaxx, vmaxyz; ///maximum velocity in velocity space
    /// Calculation of initial distribution function
    double hxi; ///integration parameter
    double ximax; ///high limit of integration
    ///Time evolution
    double dt; ///time step
    int nt; ///number of time steps

    ///Computational lattice in coordinate space
    double *rx, *ry, *rz;
    ///Computational lattice in velocity space
    double *vx, *vy, *vz;
    ///current time
    double t;

    ///Unknown for which the sceme was written
    TArray<double, 6> f; ///distribution function of ions
    TArray<double, 6> ftime; ///helping distribution function
    TArray<double, 3> n; ///ions' concentration
    TArray<double, 3> fi; ///potential
    TArray<double, 3> ax, ay, az; ///accelaration of ions

    ///Helping variables
    TArray<double, 3> maxvell; ///save Maxwellian lattice function

    ///Physical parameters of the system
    double r0x, r0y, r0z; ///dust particle position in computational box
    double Eext; ///external electricity field
    double wc; ///collision frequency
    double vfl; ///velocity of the ion flow
    double rde; ///electrons' Debye radious
    double q; ///charge of the dust particle

    ///Control parameters
    double fullCharge; ///full charge of the system

    ///3 functions containing numerical scheme
    void kinetic();

    void poisson();

    void force();

    ///Helping functions
    static double maxwellianDistribution(double vx, double vy, double vz); ///compute Maxwell distribution
    void density(); ///compute density
    ///compute distance between two points in coordinate space
    inline static double distance(double r1x, double r1y, double r1z, double r2x, double r2y, double r2z);

    double vmaxxCompute() const; ///compute the maximum velocity in velocity space in x direction
    ///set the initial distribution function
    double initialDisrtibutionFunction(double vx, double vy, double vz) const;

    ///calculate Debye potential due to the dust particle
    double potentialDebye(double rx, double ry, double rz) const;

    void calculateFullCharge(); ///calculate full charge in computational box
};

#endif
