#include "kinetic.h"
#include "poisson.h"
#include "converter.h"
#include "output.h"
#include "iostream"
#include <iomanip>
#include <sstream> 
#include <fstream> 
#include "math.h"


using namespace std;

output::output(const converter & Converter_, const poisson & Poisson_,  kinetic & Kinetic_):
Converter(&Converter_),
Poisson(&Poisson_),
Kinetic(&Kinetic_)
{
    Time = 0;
}


output::~output(){

}


void output::StartOutput(){

    cout << "######### Physical parameters of the system #############" << endl;
    cout << "#########################################################" << endl;
    cout << "Temperature       : " << Converter->GetTemperature() << " K" << endl;
    cout << "Thermal velocity  : " << Converter->GetThermalVelocity() << " cm/s" << endl;
    cout << "Debay radious     : " << Converter->GetRd() << " cm" << endl;
    cout << "Plasmas frequency : " << Converter->GetPlasmasFrequency() << " 1/s" << endl;
    cout << "Relaxation time   : " << Converter->GetPhysicalTau() << " s" << endl;
    cout << "Ion concenteration: " << Converter->GetPhysicalIonConcentration() << " 1/cm**3" << endl;
    cout << "Particle charge   : " << Converter->GetParticleCharge() << " (electron charges)" << endl;
    cout << "External force    : " << Converter->GetStrength() << " (cgs units)" << endl;
    cout << "Physical delta t  : " << Kinetic->GetDeltaT() * Converter->GetRd() / Converter->GetThermalVelocity() << "s" << endl;
    cout << endl;


    cout << "########## Scheme parameters ###########################" << endl;
    cout << "#########################################################" << endl;
    cout << "                        X       Y       Z     " << endl;
    cout << endl;
    cout << "coordinate step  " << "\t"<< Converter->GetCoordinateStepX() << "\t" << Converter->GetCoordinateStepY() << "\t" << Converter->GetCoordinateStepZ() << endl;
    cout << "velocity step    " << "\t"<< Kinetic->Gethvx()<< "\t" << Kinetic->Gethvy()<< "\t" << Kinetic->Gethvz() << endl;
    cout << "coordinate N     " << "\t"<< Converter->GetNx()<< "\t" << Converter->GetNy()<< "\t" << Converter->GetNz() << endl;
    cout << "velocity N       " << "\t"<< Kinetic->GetNvx() << "\t"<< Kinetic->GetNvy() << "\t"<< Kinetic->GetNvz() << endl;
    cout << "cut velocity     " << "\t"<< Kinetic->GetCutVelocityAlong() << "\t"<< Kinetic->GetCutVelocityNormal() << "\t"<< Kinetic->GetCutVelocityNormal() << endl;
    cout << "length of the box" << "\t"<< Converter->GetCoordinateStepX() * Converter->GetNx() << "\t" << Converter->GetCoordinateStepY() * Converter->GetNy() << "\t" << Converter->GetCoordinateStepZ() * Converter->GetNz() << endl;
    cout << "#########################################################" << endl;
    cout << endl;

    cout << "### Some important parameters for numerical scheme ######" << endl;
    cout << "#########################################################" << endl;
    cout << "Flow velocity    :" << Converter->GetDimensionlessFlowVelocity() << endl;
    cout << "Relaxation time  :" << Converter->GetDimensionlessTau() << endl;
    cout << "Ion concentration:" << Converter->GetDimensionlessIonConcentration() << endl;
    cout << "Delta T          :" << Kinetic->GetDeltaT() << endl;
    cout << endl;

}

void output::TerminalOutput(int term_time){
    if(Time % term_time == 0){
        cout << endl;
        cout << "####### Simulation time = " << Kinetic->GetCurrentTime() << " " << " #######" << endl; 
    }  
}


void output::VTKoutput(int output_time){

    int nx, ny, nz;


    nx = Converter->GetNx();
    ny = Converter->GetNy();
    nz = Converter->GetNz();

    double ***density = Kinetic->GetDensity();
    double ***potential = Poisson->GetPotential();
    double ***vfl_x = Kinetic->GetFlowVelocityX();
    double ***vfl_y = Kinetic->GetFlowVelocityY();
    double ***vfl_z = Kinetic->GetFlowVelocityZ();

    /* This part for testing Poisson solver*/
    //double pi = 3.14159265359;
    /*****************************************/
    stringstream density_filename;
    stringstream potential_filename;
    stringstream velocity_filename;
    density_filename << "data/density_t" << Time << ".dat";
    potential_filename << "data/potential_t" << Time << ".dat";
    velocity_filename << "data/velocity_t" << Time << ".dat";

    ofstream output_density;
    ofstream output_potential;
    ofstream output_velocity;

    //cout << "*** Starting VTK output ***" << endl;

    if(Time % output_time == 0){
        output_density.open(density_filename.str().c_str());
   
        output_density << Converter->GetNx()<< "\n";
        for(int k = 0; k < nz; k++){
            for(int j = 0; j < ny; j++) {
                for(int i = 0; i < nx; i++) {
                    output_density <<  density[i][j][k] << "\n";
                }
            }
        }
        output_density.close();

        output_velocity.open(velocity_filename.str().c_str());

        output_velocity << Converter->GetNx()<< "\n";
        for(int k = 0; k < nz; k++){
            for(int j = 0; j < ny; j++) {
                for(int i = 0; i < nx; i++) {
                    output_velocity << vfl_x[i][j][k] << "\t" << vfl_y[i][j][k] << "\t" << vfl_z[i][j][k] <<"\n";
                }
            }
        }
        output_velocity.close();

        output_potential.open(potential_filename.str().c_str());

        output_potential << Converter->GetNx()<< "\n";
        for(int k = 0; k < nz; k++){
            for(int j = 0; j < ny; j++) { 
                for(int i = 0; i < nx; i++) {  
                    output_potential <<  potential[i][j][k] << "\n";
                    /* This part for testing Poisson solver*/
                    //output_potential <<  potential[i][j][k] - sin(pi*i/double(nx-1))*sin(pi*j/double(ny-1))*sin(pi*k/double(nz-1))<< "\n";
                    /*****************************************/
                }
            }
        }
        output_potential.close();
    }

    Time ++;

    //cout << "*** Finalizing VTK output ***" << endl;

}


void output::CleanData(){
    
    system("rm data/density_t*");
    system("rm data/potential_t*");
    system("rm data/velocity_t*");
    system("rm data/gnuplot/profile_x*");
    system("rm data/gnuplot/profile_y*");
    system("rm data/gnuplot/profile_z*");
    system("rm data/gnuplot/analytial_profile.dat");


}


void output::PlotDistributionFunction(int i, int j, int k, int output_time){

    Kinetic->SetProfileX(i,j,k);
    Kinetic->SetProfileY(i,j,k);
    Kinetic->SetProfileZ(i,j,k);

    int nvx, nvy, nvz;
    double *profile_x;
    double *profile_y;
    double *profile_z;
    double *velocity_x;
    double *velocity_y;
    double *velocity_z;

    profile_x = Kinetic->GetProfileX();
    profile_y = Kinetic->GetProfileY();
    profile_z = Kinetic->GetProfileZ();


    velocity_x = Kinetic->GetVelocitySetX();
    velocity_y = Kinetic->GetVelocitySetY();
    velocity_z = Kinetic->GetVelocitySetZ();

    double time;

    time = Kinetic->GetCurrentTime();

    nvx = Kinetic->GetNvx();
    nvy = Kinetic->GetNvy();
    nvz = Kinetic->GetNvz();

    ofstream output_profile_x;
    ofstream output_profile_y;
    ofstream output_profile_z;

    stringstream profile_x_filename;
    stringstream profile_y_filename;
    stringstream profile_z_filename;

    stringstream gnuplot_command_x;
    stringstream gnuplot_command_y;
    stringstream gnuplot_command_z;

    profile_x_filename << "data/gnuplot/profile_x_t" << Time << ".dat";
    profile_y_filename << "data/gnuplot/profile_y_t" << Time << ".dat";
    profile_z_filename << "data/gnuplot/profile_z_t" << Time << ".dat";
    
        ofstream analytical_profile("data/gnuplot/analytial_profile.dat");
        for (int a=0; a<nvx; a++){
            analytical_profile << velocity_x[a] << "\t" <<  Kinetic->GetAnalyticalProfileX()[a] << endl;
        }
        analytical_profile.close();
    

    if (Time % output_time == 0){
        output_profile_x.open(profile_x_filename.str().c_str());
        for (int a=0; a<nvx; a++){
            output_profile_x << velocity_x[a] << "\t" <<  profile_x[a] << endl;
        }
        output_profile_x.close();

        output_profile_y.open(profile_y_filename.str().c_str());
        for (int b=0; b<nvy; b++){
            output_profile_y << velocity_y[b] << "\t" << profile_y[b] << endl;
        }
        output_profile_y.close();

        output_profile_z.open(profile_z_filename.str().c_str());
        for (int c=0; c<nvz; c++){
            output_profile_z << velocity_z[c] << "\t" << profile_z[c] << endl;
        }
        output_profile_z.close();

    }
    
    

}


void output::WriteColoumbForceField(){


    int nx, ny, nz;


    nx = Converter->GetNx();
    ny = Converter->GetNy();
    nz = Converter->GetNz();

    ofstream output_field;
    
    output_field.open("data/ForceField.dat");

    double*** force_field_x = Kinetic->GetColoumbForceFieldX();
    double*** force_field_y = Kinetic->GetColoumbForceFieldY();
    double*** force_field_z = Kinetic->GetColoumbForceFieldZ();

    output_field << Converter->GetNx()<< "\n";
    for(int k = 0; k < nz; k++){
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                output_field <<  force_field_x[i][j][k] << "\t" <<  force_field_y[i][j][k] << "\t" << force_field_z[i][j][k] <<"\n";
            }
        }
    }
    output_field.close();
    
}
