#include "kinetic.h"
#include "converter.h"
#include "poisson.h"
#include <cmath>
#include <iostream>
#include <chrono>
#include <algorithm>
#include "output.h"
#include "omp.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]) {

    //#######################################//
    //## Physical parameters of the system ##//
    //#######################################//

    ///////////////////////////////////////////
    //ion temperature  (kelvin)                      
    const double T_i = 300;
    // time between ion-neutral collisions in seconds                  
    const double tau = 0.4 * pow(10., -3.);
    //equilibrium concentaration of ions [1/cm^3]  
    const double n_0 = 2. * pow(10., 7.);
    //particle charge in electron units                      
    const double q = 2000;
    //electricity field in cgs units                     
    const double El = 0 * pow(10., -5.);
    //In debay radious units                 
    //x-dimension of the computational box in Debye radious  
    const double Lx = 6;
    //x-dimension of the computational box in Debye radious  
    const double Ly = 6;
    //x-dimension of the computational box in Debye radious   
    const double Lz = 6;
    //number of sells in cordinate space x   
    const int Nx = 10;
    //number of sells in cordinate space y   
    const int Ny = 10;
    //number of sells in cordinate space z   
    const int Nz = 10;
    //position of dust particle at the grid
    const int Nx_0 = 5;
    //position of dust particle at the grid
    const int Ny_0 = 5;
    //position of dust particle at the grid
    const int Nz_0 = 5;
    ///////////////////////////////////////////


    //#######################################//
    //##     Define output frequency       ##//
    //#######################################//

    ///////////////////////////////////////////
    int ITmax = 10;
    int T_term = 1;
    int T_output = 1;
    ///////////////////////////////////////////


    //#######################################//
    //##         Inizialization            ##//
    //#######################################//
//    omp_set_num_threads(44);
    ///////////////////////////////////////////
    auto start = high_resolution_clock::now();

    converter Converter(T_i, tau, n_0, q, El, Lx, Ly, Lz, Nx, Ny, Nz, Nx_0, Ny_0, Nz_0);
    poisson Poisson(Converter);
    kinetic Kinetic(Converter, Poisson);
    output Output(Converter, Poisson, Kinetic);
    Output.StartOutput();
    Output.CreateOutputDirectory();
    Output.WriteInitialPotential();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by inizialisation: "
         << duration.count() / 1000000.0 << " seconds" << endl;
    ///////////////////////////////////////////

    //#######################################//
    //##             Time loop             ##//
    //#######################################//

    start = high_resolution_clock::now();
    for (int IT = 0; IT < ITmax; IT++) {
        Output.TerminalOutput(T_term);
        Poisson.PoissonScheme(Kinetic.GetDensity());
        Poisson.GradientScheme();
        Kinetic.IntegrateAll();
        Output.VTKoutput(T_output);
        Output.PlotDistributionFunction(1, 1, 1, T_output);
        Output.UpdateTime();
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by numerical scheme: "
         << duration.count() / 1000000.0 << " seconds" << endl;

    return 0;

}
