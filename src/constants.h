#ifndef OPTIMIZATION_CONSTANTS_H
#define OPTIMIZATION_CONSTANTS_H

///electron charge
const double E = 4.8032e-10;

///ion mass
const double M = 6.6464e-23;

///Boltzmann constant
const double K = 1.3806e-16;

/// maximum value of integration patameter in velocity space
const int MAX_XI = 15;

///discrete step for integration parameter
const double DELTA_XI = 0.005;

const double hvx = 0.2;
const double hvy = 0.2;
const double hvz = 0.2;

////grid size X
//const int Nx;
////grid size Y
//const int Ny;
////grid size Z
//const int Nz;

//number of sells in cordinate space x
const int Nx = 10;
//number of sells in cordinate space y
const int Ny = 10;
//number of sells in cordinate space z
const int Nz = 10;

#endif //OPTIMIZATION_CONSTANTS_H
