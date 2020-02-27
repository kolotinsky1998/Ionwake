#include <cmath>
#include <iostream>
#include <chrono>

#include "scheme.hpp"
#include "converter.hpp"
#include "omp.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]) {

    //#######################################//
    //## Physical parameters of the system ##//
    //#######################################//

    ///////////////////////////////////////////
    //ion temperature  (kelvin)
    const double Ti = 300;
    //electron temperature (kelvin)
    const double Te = 9000;
    // time between ion-neutral collisions in seconds
    const double wc = 361745;
    //equilibrium concentaration of ions [1/cm^3]
    const double ni = 3. * pow(10., 7.);
    //particle charge in electron units
    const double q = 10000;
    //electricity field in cgs units
    const double Eext = 0.257e-2;
    //ion mass
    const double mi = 6.6464 * pow(10., -23.);
    //In debay radious units
    //x-dimension of the computational box in Debye radious
    const double lx = 15;
    //x-dimension of the computational box in Debye radious
    const double ly = 6;
    //x-dimension of the computational box in Debye radious
    const double lz = 6;
    //coordinate-x of the dust particle
    const double r0x = 3;
    //coordinate-y of the dust particle
    const double r0y = 3;
    //coordinate-z of the dust particle
    const double r0z = 3;
    ///////////////////////////////////////////

    //########################################//
    //## Numerical parameters of the system ##//
    //########################################//

    ///////////////////////////////////////////
    //number of cells in cordinate space x
    const int nx = 10;
    //number of cells in cordinate space y
    const int ny = 10;
    //number of ells in cordinate space z
    const int nz = 10;
    //integration step in velocity space x
    const double hvx = 0.2;
    //integration step in velocity space y
    const double hvy = 0.2;
    //integration step in velocity space z
    const double hvz = 0.2;
    //maximum velocity yz
    const double vmaxyz = 1;
    //minimum velocity yz
    const double vminyz = -1;
    //integration step for calculation initial conditions
    const double hxi = 0.002;
    //integration limit for calculation initial conditions
    const double ximax = 5;
    ///////////////////////////////////////////


    //#######################################//
    //##     Define output frequency       ##//
    //#######################################//

    ///////////////////////////////////////////
    int ITmax = 100;
    int T_output = 20;
    ///////////////////////////////////////////


    //#######################################//
    //##         Inizialization            ##//
    //#######################################//

//    omp_set_num_threads(44);

    auto start = high_resolution_clock::now();
    Converter converter(Te, Ti, ni, q, Eext, wc, mi);
    Scheme scheme(
            converter, nx, ny, nz, lx, ly, lz,
            ITmax, r0x, r0y, r0z, hvx, hvy, hvz, vminyz, vmaxyz, hxi, ximax
    );
    scheme.InitialLogOut();
    converter.InitialLogOut();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by inizialisation: "
         << duration.count() / 1000000.0 << " seconds" << endl;

    //####################################//
    //##         Computation            ##//
    //####################################//
    cout << "##############################################" << endl;
    cout << "##******* Comutations start ****************##" << endl;
    cout << "##############################################" << endl;
    auto total_start = high_resolution_clock::now();
    for (int t = 0; t < ITmax; t++) {
        start = high_resolution_clock::now();
        cout << "# Step №" << t << endl;
        scheme.schemeStep();
        scheme.printCurrentTime();
        scheme.printFullCharge();
        if (t % T_output == 0) {
            scheme.writeDensityFile("./");
            scheme.writePotentialFile("./");
        }

        stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by one step: "
             << duration.count() / 1000000.0 << " seconds" << endl;
    }
    auto total_stop = high_resolution_clock::now();
    auto total_duration = duration_cast<microseconds>(total_stop - total_start);
    cout << "Time taken by scheme: "
         << total_duration.count() / 1000000.0 << " seconds" << endl;

    return 0;
}
