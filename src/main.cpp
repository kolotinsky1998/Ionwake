#include <cmath>
#include <iostream>
#include <chrono>
#include <sstream>
#include <fstream>

#include "scheme.hpp"
#include "TScheme.h"
#include "converter.hpp"
#include "omp.h"

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
    const size_t nx = 10;
    //number of cells in cordinate space y
    const size_t ny = 10;
    //number of ells in cordinate space z
    const size_t nz = 10;
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
    size_t ITmax = 3;
    size_t T_output = 20;
    ///////////////////////////////////////////


    //#######################################//
    //##         Inizialization            ##//
    //#######################################//

//    omp_set_num_threads(44);


    const auto start = std::chrono::high_resolution_clock::now();
    const time_t start_time_t = std::chrono::system_clock::to_time_t(start);
    std::cout << "Start scheme crating at " << std::ctime(&start_time_t);

    TScheme<double> scheme = TScheme<double>::TBuilder()
            .set_electron_temperature(Te)
            .set_ion_temperature(Ti)
            .set_ion_concentration(ni)
            .set_dust_particle_charge(q)
            .set_external_electricity_field(Eext)
            .set_ion_neutral_collisions_frequency(wc)
            .set_ion_mass(mi)
            .set_box_split_step_count(nx, ny, nz)
            .set_box_size(lx, ly, lz)
            .set_particle_position(r0x, r0y, r0z)
            .set_velocity_step(hvx, hvy, hvz)
            .set_velocity_bound_for_yz(vminyz, vmaxyz)
            .set_integration_parameters(ximax, hxi)
            .build();

    const auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken by inizialisation: "
              << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0
              << " seconds" << endl;

    //####################################//
    //##         Computation            ##//
    //####################################//
    std::cout << "##############################################\n"
                 "##******* Comutations start ****************##\n"
                 "##############################################\n";
    const auto total_start = std::chrono::high_resolution_clock::now();
    const time_t total_start_time_t = std::chrono::system_clock::to_time_t(start);
    std::cout << "Start time: " << std::ctime(&total_start_time_t);
    for (size_t t = 0; t < ITmax; ++t) {
        std::cout << "# Step " << t << endl;
        const auto step_start = std::chrono::high_resolution_clock::now();
        scheme.next_step();

        std::cout << "Current dimmensionless time: " << t * scheme.get_dt() << endl;
        std::cout << "Current full charge in the computational box: " << scheme.get_full_charge() << endl;
        if (t % T_output == 0) {
            ofstream file;

            stringstream filename;
            filename << "./density_t" << t << ".dat";
            file.open(filename.str().c_str());
            scheme.write_density(file);
            file.close();

            stringstream filename2;
            filename2 << "./potential_t" << t << ".dat";
            file.open(filename2.str().c_str());
            scheme.write_potential(file);
            file.close();
        }

        const auto step_stop = std::chrono::high_resolution_clock::now();
        std::cout << "Time taken by one step: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(step_stop - step_start).count() / 1000000.0
                  << " seconds" << endl;
    }
    const auto total_stop = std::chrono::high_resolution_clock::now();
    const time_t total_stop_time_t = std::chrono::system_clock::to_time_t(start);
    std::cout << "Stop time: " << std::ctime(&total_stop_time_t);
    std::cout << "Time taken by scheme: "
              << std::chrono::duration_cast<std::chrono::microseconds>(total_stop - total_start).count() / 1000000.0
              << " seconds" << endl;

    return 0;
}
