#include <cmath>
#include <iostream>
#include <chrono>
#include <algorithm>

#include "kinetic.h"
#include "converter.h"
#include "poisson.h"
#include "output.h"
#include "omp.h"
#include "ionwake/IonWake.h"

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
    //cout << omp_get_num_threads() << endl;

    int numThreads;
    cout << "Enter numThreads: ";
    cin >> numThreads;
    omp_set_num_threads(numThreads);
    ///////////////////////////////////////////
    auto start = high_resolution_clock::now();
    system("mkdir gnuplot");
    system("mkdir density");
    system("mkdir potential");
    system("mkdir velocity");

    ionwake::IonWake ionWake(10, 10, 10, Nx_0, Ny_0, Nz_0, Lx, Ly, Lz, q, T_i, n_0, tau, El);
    cout << "######### Physical parameters of the system #############\n";
    cout << "#########################################################\n";
    cout << "Temperature       : " << ionWake.getTemperature() << " K\n";
    cout << "Thermal velocity  : " << ionWake.getThermalVelocity() << " cm/s\n";
    cout << "Debay radious     : " << ionWake.getDebayRadius() << " cm\n";
    cout << "Plasmas frequency : " << ionWake.getPlasmasFrequency() << " 1/s\n";
    cout << "Relaxation time   : " << ionWake.getPhysicalTau() << " s\n";
    cout << "Ion concenteration: " << ionWake.getPhysicalIonConcentration() << " 1/cm**3\n";
    cout << "Particle charge   : " << ionWake.getParticleCharge() << " (electron charges)\n";
    cout << "External force    : " << ionWake.getStrength() << " (cgs units)\n";
    cout << "Physical delta t  : " << ionWake.getDeltaT() * ionWake.getDebayRadius() / ionWake.getThermalVelocity()
         << "s\n";
    cout << "\n";


    cout << "########## Scheme parameters ############################\n";
    cout << "#########################################################\n";
    cout << "                        X       Y       Z     \n\n";
    cout << "coordinate step  " << "\t" << ionWake.getCoordinateStepX() << "\t" << ionWake.getCoordinateStepY()
         << "\t" << ionWake.getCoordinateStepZ() << "\n";
    cout << "velocity step    " << "\t" << ionWake.getHvx() << "\t" << ionWake.getHvy() << "\t" << ionWake.getHvz()
         << "\n";
    cout << "coordinate N     " << "\t" << ionWake.getNx() << "\t" << ionWake.getNy() << "\t"
         << ionWake.getNz() << "\n";
    cout << "velocity N       " << "\t" << ionWake.getNvx() << "\t" << ionWake.getNvy() << "\t" << ionWake.getNvz()
         << "\n";
    cout << "cut velocity     " << "\t" << ionWake.getCutVelocityAlong() << "\t" << ionWake.getCutVelocityNormal()
         << "\t" << ionWake.getCutVelocityNormal() << "\n";
    cout << "length of the box" << "\t" << ionWake.getCoordinateStepX() * ionWake.getNx() << "\t"
         << ionWake.getCoordinateStepY() * ionWake.getNy() << "\t"
         << ionWake.getCoordinateStepZ() * ionWake.getNz() << "\n";
    cout << "#########################################################\n\n";

    cout << "### Some important parameters for numerical scheme ######\n";
    cout << "#########################################################\n";
    cout << "Flow velocity    :" << ionWake.getDimensionlessFlowVelocity() << "\n";
    cout << "Relaxation time  :" << ionWake.getDimensionlessTau() << "\n";
    cout << "Ion concentration:" << ionWake.getDimensionlessIonConcentration() << "\n";
    cout << "Delta T          :" << ionWake.getDeltaT() << "\n";
    cout << endl;
//    converter Converter(T_i, tau, n_0, q, El, Lx, Ly, Lz, Nx_0, Ny_0, Nz_0);
//    poisson Poisson(Converter);
//    kinetic Kinetic(Converter, Poisson);
//    output Output(Converter, Poisson, Kinetic);
//    Output.StartOutput();
//    Output.CreateOutputDirectory();
//    Output.WriteInitialPotential();


    std::ofstream initialPotential("./InitialPotential.dat", std::ofstream::out);
    ionWake.writeInitialPotential(initialPotential);
    initialPotential.close();

    auto duration = duration_cast<microseconds>(high_resolution_clock::now() - start);
    cout << "Time taken by inizialisation: "
         << duration.count() / 1000000.0 << " seconds" << endl;
    ///////////////////////////////////////////

    //#######################################//
    //##             Time loop             ##//
    //#######################################//

    start = high_resolution_clock::now();
    for (int i = 0; i < ITmax; i++) {
        cout << "####### Simulation time = " << ionWake.getCurrentTime() << " #######" << endl << endl;

        ionWake.nextStep();

        std::ofstream density("./density/density_x" + std::to_string(i) + ".dat", std::ofstream::out);
        ionWake.writeDensity(density);
        density.close();

        std::ofstream potential("./potential/potential_x" + std::to_string(i) + ".dat", std::ofstream::out);
        ionWake.writePotential(potential);
        potential.close();

        std::ofstream velocity("./velocity/velocity_x" + std::to_string(i) + ".dat", std::ofstream::out);
        ionWake.writeVelocity(velocity);
        velocity.close();

        std::ofstream profile_x("./gnuplot/profile_x_t" + std::to_string(i) + ".dat", std::ofstream::out);
        ionWake.plotDistributionFunctionX(1, 1, 1, profile_x);
        profile_x.close();

        std::ofstream profile_y("./gnuplot/profile_y_t" + std::to_string(i) + ".dat", std::ofstream::out);
        ionWake.plotDistributionFunctionY(1, 1, 1, profile_y);
        profile_x.close();

        std::ofstream profile_z("./gnuplot/profile_z_t" + std::to_string(i) + ".dat", std::ofstream::out);
        ionWake.plotDistributionFunctionZ(1, 1, 1, profile_z);
        profile_z.close();

//        Output.TerminalOutput(T_term);
//        Poisson.PoissonScheme(Kinetic.GetDensity());
//        Poisson.GradientScheme();
//        Kinetic.IntegrateAll();
//        Output.VTKoutput(T_output);
//        Output.PlotDistributionFunction(1, 1, 1, T_output);
//        Output.UpdateTime();
    }
    duration = duration_cast<microseconds>(high_resolution_clock::now() - start);
    cout << "Time taken by numerical scheme: "
         << duration.count() / 1000000.0 << " seconds" << endl;

    return 0;

}
