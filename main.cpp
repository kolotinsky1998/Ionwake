#include "kinetic.h"
#include "converter.h"
#include "poisson.h"
#include "math.h"
#include "output.h"
#include "iostream"
#include <chrono>
#include <algorithm>  
using namespace std;
using namespace std::chrono;
int main(){

    //#######################################//
    //## Physical parameters of the system ##//
    //#######################################//

    ///////////////////////////////////////////
    //ion temperature  (kelvin)                      
    const double T_i = 300;                                             
    // relaxation time                       
    const double tau = 0.4 * pow(10.,-3.);                       
    //equilibrium concentaration of ions     
    const double n_0 = 2. * pow(10.,4.);                        
    //particle charge in electron units                      
    const double q = 0*100;                               
    //electricity force                      
    const double El= 2*pow(10.,-5.);                        
    //In debay radious units                 
    //x-dimension of the computational box   
    const double Lx=6;                       
    //x-dimension of the computational box   
    const double Ly=6;                        
    //x-dimension of the computational box   
    const double Lz=6;                        
    //number of sells in cordinate space x   
    const int Nx=40;                           
    //number of sells in cordinate space y   
    const int Ny=40;                           
    //number of sells in cordinate space z   
    const int Nz=40;                           
    ///////////////////////////////////////////


    //#######################################//
    //##     Define output frequency       ##//
    //#######################################//

    ///////////////////////////////////////////
    int ITmax = 1000;
    int T_term = 1;
    int T_output = 10;
    ///////////////////////////////////////////


    //#######################################//
    //##         Inizialization            ##//
    //#######################################//

    ///////////////////////////////////////////
    auto start = high_resolution_clock::now();

    converter Converter(T_i, tau, n_0, q, El, Lx, Ly, Lz, Nx, Ny, Nz);
    poisson Poisson(Converter);
    kinetic Kinetic(Converter, Poisson);
    output Output(Converter, Poisson, Kinetic);
    Output.StartOutput();
    Output.CleanData();
    Output.WriteColoumbForceField();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by inizialisation: "
         << duration.count() << " microseconds" << endl;
    ///////////////////////////////////////////

    //#######################################//
    //##             Time loop             ##//
    //#######################################//
    
    start = high_resolution_clock::now(); 
    for(int IT=0; IT<ITmax; IT ++){
        Output.TerminalOutput(T_term);
        //Kinetic.DefineDensity();
        Poisson.PoissonScheme(Kinetic.GetDensity());
        Poisson.GradientScheme();
        Kinetic.IntegrateAll();
        Output.VTKoutput(T_output);
        Output.PlotDistributionFunction(1,1,1,T_output);
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by numerical scheme: "
         << duration.count() << " microseconds" << endl;
    
    return 0;

}
