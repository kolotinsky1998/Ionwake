#pragma once

#include"converter.h"

///In this class numerical scheme to solve Poisson equation is implemented
class poisson {
public:

    ///constructor

    poisson(const converter &Converter);

    ///destructor
    ~poisson();

    ///using this you can acces to self-consistent force field x-component
    double ***GetSelfConsistentForceFieldX() const;

    ///using this you can acces to self-consistent force field y-component
    double ***GetSelfConsistentForceFieldY() const;

    ///using this you can acces to self-consistent force field z-component
    double ***GetSelfConsistentForceFieldZ() const;

    ///using this you can acces to self-consistent potential
    double ***GetPotential() const;

    ///numerical scheme for solving the Poisson equation
    void PoissonScheme(double ***density);

    ///numerical scheme for taking gradient
    void GradientScheme();

private:

    ///self consistent force field x-component
    double ***sfx;
    ///self consistent force field y-component
    double ***sfy;
    ///self consistent force field z-component
    double ***sfz;
    ///self-consistent potential
    double ***potential;
    ///coordinate discretization step x-direction
    double hx;
    ///coordinate discretization step y-direction
    double hy;
    ///coordinate discretization step z-direction
    double hz;
    ///grid size X
    int Nx;
    ///grid size Y
    int Ny;
    ///grid size Z
    int Nz;
    //equilibrium ion concentration
    double n_i;
    ///particle charge
    double q;
    //position at the grid of the dust particle X
    double Nx_0;
    //position at the grid of the dust particle Y
    double Ny_0;
    //position at the grid of the dust particle Z
    double Nz_0;


};