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

    double Density;
    double Potential;

    double ***density = Kinetic->GetDensity();
    double ***potential = Poisson->GetPotential();
    double ***vfl_x = Kinetic->GetFlowVelocityX();
    double ***vfl_y = Kinetic->GetFlowVelocityY();
    double ***vfl_z = Kinetic->GetFlowVelocityZ();

    stringstream density_filename;
    stringstream potential_filename;
    stringstream velocity_filename;

    density_filename << data << "/density_t" << Time << ".dat";
    potential_filename << data << "/potential_t" << Time << ".dat";
    velocity_filename << data << "/velocity_t" << Time << ".dat";

    ofstream output_density;
    ofstream output_potential;
    ofstream output_velocity;


    if(Time % output_time == 0){
        output_density.open(density_filename.str().c_str());
   
        output_density << nx << "\t" << ny << "\t" << nz <<  "\n";
        for(int k = 0; k < nz; k++){
            for(int j = 0; j < ny; j++) {
                for(int i = 0; i < nx; i++) {
                    Density = density[i][j][k] * pow(Converter->GetRd()/Converter->GetThermalVelocity(),3.);
                    output_density <<  Density << "\n";
                }
            }
        }
        output_density.close();

        output_velocity.open(velocity_filename.str().c_str());

        output_velocity << nx << "\t" << ny << "\t" << nz <<  "\n";
        for(int k = 0; k < nz; k++){
            for(int j = 0; j < ny; j++) {
                for(int i = 0; i < nx; i++) {
                    output_velocity << vfl_x[i][j][k] << "\t" << vfl_y[i][j][k] << "\t" << vfl_z[i][j][k] <<"\n";
                }
            }
        }
        output_velocity.close();

        output_potential.open(potential_filename.str().c_str());

        output_potential << nx << "\t" << ny << "\t" << nz <<  "\n";
        for(int k = 0; k < nz; k++){
            for(int j = 0; j < ny; j++) { 
                for(int i = 0; i < nx; i++) {  
                    Potential = Converter->GetDimensionlessIonConcentration() * Converter->ConvertPotential() * potential[i][j][k];
                    output_potential <<  Potential << "\n";
                }
            }
        }
        output_potential.close();
    }

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

    nvx = Kinetic->GetNvx();
    nvy = Kinetic->GetNvy();
    nvz = Kinetic->GetNvz();

    ofstream output_profile_x;
    ofstream output_profile_y;
    ofstream output_profile_z;

    stringstream profile_x_filename;
    stringstream profile_y_filename;
    stringstream profile_z_filename;

    profile_x_filename << data << "/gnuplot/profile_x_t" << Time << ".dat";
    profile_y_filename << data << "/gnuplot/profile_y_t" << Time << ".dat";
    profile_z_filename << data << "/gnuplot/profile_z_t" << Time << ".dat";

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


void output::WriteInitialPotential(){


    int nx, ny, nz;
    double lx, ly, lz;

    double potential;
    double r;
    double rd;

    nx = Converter->GetNx();
    ny = Converter->GetNy();
    nz = Converter->GetNz();

    lx = Converter->GetLx();
    ly = Converter->GetLy();
    lz = Converter->GetLz();

    rd = Converter->GetRd();

    stringstream output_potential_filename;
    ofstream output_potential;
    
    output_potential_filename << data << "/InitialPotential.dat";
    output_potential.open(output_potential_filename.str().c_str());


    output_potential << Converter->GetNx()<< "\t" << Converter->GetNy() << "\t" << Converter->GetNz() << "\n";
    for(int k = 0; k < nz; k++){
        for(int j = 0; j < ny; j++) {
            for(int i = 0; i < nx; i++) {
                r = (i - Converter->GetNx_0())*(i - Converter->GetNx_0()) * (lx * rd /double(nx)) * (lx * rd /double(nx));  
                r = r + (j - Converter->GetNy_0())*(j - Converter->GetNy_0()) * (ly * rd /double(ny)) * (ly * rd /double(ny));
                r = r + (k - Converter->GetNz_0())*(k - Converter->GetNz_0()) * (lz * rd /double(nz)) * (lz * rd /double(nz));
                r = sqrt(r);
                potential = - Converter->GetParticleCharge() * Converter->GetElementaryCharge() / r;
                potential = potential - Converter->GetEl() * (k * lz * rd /double(nz)); 
                output_potential << potential << "\n";
            }
        }
    }
    output_potential.close();
    
}


void output::CreateOutputDirectory(string data_="data"){
    stringstream command, command1;
    data = data_;
    command << "mkdir " << data;
    command1 << "mkdir " << data << "/gnuplot";
    system(command.str().c_str());
    system(command1.str().c_str());
}

void output::UpdateTime(){

	Time ++;

}
