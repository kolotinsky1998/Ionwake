#include "kinetic.h"
#include "poisson.h"
#include "math.h"
#include <iostream>
#include "omp.h"
using namespace std;

kinetic::kinetic(const converter &Converter_, const poisson &Poisson):
Converter(&Converter_)
{
    

    pi = Converter->GetPi();

    hx = Converter->GetCoordinateStepX();
    hy = Converter->GetCoordinateStepY();
    hz = Converter->GetCoordinateStepZ();

    nx = Converter->GetNx();
    ny = Converter->GetNy();
    nz = Converter->GetNz();
    
    nx_0 = Converter->GetNx_0();
    ny_0 = Converter->GetNy_0();
    nz_0 = Converter->GetNz_0();

    Lx = Converter->GetLx();
    Ly = Converter->GetLy();
    Lz = Converter->GetLz();

    hvx = 0.2;
    hvy = 0.2;
    hvz = 0.2;

    vxcut = max(5.,8 * Converter->GetDimensionlessFlowVelocity());
    vyzcut = 5;
    nvx = int ((vxcut + vyzcut) / hvx);
    nvy = int (2 * vyzcut / hvy);
    nvz = int (2 * vyzcut / hvz);

    maxXi = 15;
    deltaXi = 0.005;
    nxi = int (maxXi / deltaXi);

    al = Converter->GetLAxelerationCofficient();

    deltaT = 0.5 * min(0.2 / (al / hvx ), 1.0 / (vxcut / hx + vyzcut / hy + vyzcut / hz ));
    CurT = 0;
    tau = Converter->GetDimensionlessTau();

    f = new double *****[nx];
    for (int i=0; i<nx; i++){
        f[i] = new double**** [ny];
        for (int j=0; j<ny; j++){
            f[i][j] = new double*** [nz];
            for (int k=0; k<nz; k++){
                f[i][j][k] = new double** [nvx];
                for (int a=0; a<nvx; a++){
                    f[i][j][k][a] = new double* [nvy];
                    for (int b=0; b<nvy; b++){
                        f[i][j][k][a][b] = new double [nvz];
                    }
                } 
            }
        }
    }

    f_time = new double *****[nx];
    for (int i=0; i<nx; i++){
        f_time[i] = new double**** [ny];
        for (int j=0; j<ny; j++){
            f_time[i][j] = new double*** [nz];
            for (int k=0; k<nz; k++){
                f_time[i][j][k] = new double** [nvx];
                for (int a=0; a<nvx; a++){
                    f_time[i][j][k][a] = new double* [nvy];
                    for (int b=0; b<nvy; b++){
                        f_time[i][j][k][a][b] = new double [nvz];
                    }
                } 
            }
        }
    }

    vx = new double [nvx];
    vy = new double [nvy];
    vz = new double [nvz];

    for (int a=0; a<nvx; a++){
        vx[a] = -vyzcut + hvx * a; 
    }

    for (int b=0; b<nvy; b++){
        vy[b] = -vyzcut + hvy * b; 
    }

    for (int c=0; c<nvz; c++){
        vz[c] = -vyzcut + hvz * c; 
    }

    x = new double [nx];
    y = new double [ny];
    z = new double [nz];

    for (int i=0; i<nx; i++){
        x[i] = hx * i;
    }

    for (int j=0; j<ny; j++){
        y[j] = hy * j;
    }

    for (int k=0; k<nz; k++){
        z[k] = hz * k;
    }

    f_n = new double **[nvx];
    for (int a=0; a<nvx; a++){
        f_n[a] = new double* [nvy];
        for (int b=0; b<nvy; b++){
            f_n[a][b] = new double [nvz];
        }
    }

    for (int a=0; a<nvx; a++){
        for (int b=0; b<nvy; b++){
            for (int c=0; c<nvz; c++){
                f_n[a][b][c] = exp (-0.5 * (vx[a] * vx[a] + vy[b] * vy[b] + vz[c] * vz[c])) / pow(2*pi,1.5);
            }
        }
    }

    InitialConditions();

    profile_x = new double [nvx];
    profile_y = new double [nvy];
    profile_z = new double [nvz];
    analytical_profile_x = new double [nvx];

    for(int a=0; a<nvx; a++){
        profile_x[a] = 0;
    }

    for(int b=0; b<nvy; b++){
        profile_y[b] = 0;
    }

    for(int c=0; c<nvz; c++){
        profile_z[c] = 0;
    }

    vfl_x = new double **[nx];
    vfl_y = new double **[nx];
    vfl_z = new double **[nx];
    acx = new double **[nx];
    acy = new double **[nx];
    acz = new double **[nx];
    density = new double **[nx];
    for (int i=0; i<nx; i++){
        vfl_x[i] = new double *[ny];
        vfl_y[i] = new double *[ny];
        vfl_z[i] = new double *[ny];
        acx[i] = new double * [ny];
        acy[i] = new double * [ny];
        acz[i] = new double * [ny];
        density[i] = new double* [ny];
        for (int j=0; j<ny; j++){
            vfl_x[i][j] = new double [nz];
            vfl_y[i][j] = new double [nz];
            vfl_z[i][j] = new double [nz];
            acx[i][j] = new double [nz];
            acy[i][j] = new double [nz];
            acz[i][j] = new double [nz];
            density[i][j] = new double [nz];
        }
    } 


    asx = Poisson.GetSelfConsistentForceFieldX();
    asy = Poisson.GetSelfConsistentForceFieldY();
    asz = Poisson.GetSelfConsistentForceFieldZ();

    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                density[i][j][k] = Converter->GetDimensionlessIonConcentration();

                acx[i][j][k] =  - Converter->GetCAxelerationCofficient() * (x[i] - x[nx_0] + 0.5 * hx)
                / pow(pow((x[i] - x[nx_0] + 0.5 * hx),2.) + pow(( y[j] - y[ny_0] + 0.5 * hy),2.) + pow(( z[k] - z[nz_0] + 0.5 * hz),2.),1.5);

                acy[i][j][k] = - Converter->GetCAxelerationCofficient() * (y[j] - y[ny_0] + 0.5 * hy)
                / pow(pow((x[i] - x[nx_0] + 0.5 * hx),2.) + pow(( y[j] - y[ny_0] + 0.5 * hy),2.) + pow(( z[k] - z[nz_0] + 0.5 * hz),2.),1.5);

                acz[i][j][k] = - Converter->GetCAxelerationCofficient() * (z[k] - z[nz_0] + 0.5 * hz)
                / pow(pow((x[i] - x[nx_0] + 0.5 * hx),2.) + pow(( y[j] - y[ny_0] + 0.5 * hy),2.) + pow(( z[k] - z[nz_0] + 0.5 * hz),2.),1.5);

                vfl_x[i][j][k] = 0;
                vfl_y[i][j][k] = 0;
                vfl_z[i][j][k] = 0;
            }
        }
    }

    ComputeAnalyticalFlowVelocity();
    SetAnalyticalProfileX();
    
}

kinetic::~kinetic(){

    delete [] vx;
    delete [] vy;
    delete [] vz;

    delete [] x;
    delete [] y;
    delete [] z;

    delete [] profile_x;
    delete [] profile_y;
    delete [] profile_z;
    delete [] analytical_profile_x;
   

    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            delete [] density[i][j];
            delete [] vfl_x[i][j];
            delete [] vfl_y[i][j];
            delete [] vfl_z[i][j];
            delete [] acx[i][j];
            delete [] acy[i][j];
            delete [] acz[i][j];
        }
        delete [] density[i];
        delete [] vfl_x[i];
        delete [] vfl_y[i];
        delete [] vfl_z[i];
        delete [] acx[i];
        delete [] acy[i];
        delete [] acz[i];
    }

    delete [] density;
    delete [] vfl_x;
    delete [] vfl_y;
    delete [] vfl_z;
    delete [] acx;
    delete [] acy;
    delete [] acz;

    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                for (int a=0; a<nvx; a++){
                    for (int b=0; b<nvy; b++){
                        delete[] f[i][j][k][a][b];
                    }
                    delete[] f[i][j][k][a];
                }
                delete[] f[i][j][k]; 
            }
            delete[] f[i][j];
        }
        delete[] f[i];
    }

    delete[] f;

    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                for (int a=0; a<nvx; a++){
                    for (int b=0; b<nvy; b++){
                        delete[] f_time[i][j][k][a][b];
                    }
                    delete[] f_time[i][j][k][a];
                }
                delete[] f_time[i][j][k]; 
            }
            delete[] f_time[i][j];
        }
        delete[] f_time[i];
    }

    delete[] f_time;
    

    for (int a=0; a<nvx; a++){
        for (int b=0; b<nvy; b++){
            delete[] f_n[a][b];
        }
        delete[] f_n[a];
    }

    delete[] f_n;


}


void kinetic::CoordinatePart(){

  //cout << "***  Starting coordinate part ***" << endl; 
  SaveState();
#pragma omp parallel
{
#pragma omp for
  for (int i=0; i<nx; i++){    
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                for (int a=1; a<nvx-1; a++){
                    if (vx[a] > 0){
                        if(abs(vx[a]*deltaT/hx)>1) cout << "x"<< vx[a]*deltaT/hx << endl;
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] - (deltaT * vx[a] / hx)
                                 * (f_time[i][j][k][a][b][c] - f_time[(nx+i-1)%nx][j][k][a][b][c]);
                            } 
                        }
                    } 
                    else{
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] - (deltaT * vx[a] / hx)
                                 * (f_time[(i+1)%nx][j][k][a][b][c] - f_time[i][j][k][a][b][c]);
                            }   
                        }
                    }
                }
            }
        }
    }
}
    SaveState();
    for (int i=0; i<nx; i++){    
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                for (int a=1; a<nvx-1; a++){
                    for (int b=1; b<nvy-1; b++){
                        if (vy[b] > 0){
                            if(abs(vy[b]*deltaT/hy)>1) cout << "y"<< vy[b]*deltaT/hy << endl;
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] - (deltaT * vy[b] / hy)
                                 * (f_time[i][j][k][a][b][c] - f_time[i][(ny+j-1)%ny][k][a][b][c]);
                            }   
                        }
                        else{
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] - (deltaT * vy[b] / hy)
                                 * (f_time[i][(j+1)%ny][k][a][b][c] - f_time[i][j][k][a][b][c]);
                            }   
                        }
                    }
                }
            }
        }
    } 
    SaveState();
    for (int i=0; i<nx; i++){    
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                for (int a=1; a<nvx-1; a++){
                    for (int b=1; b<nvy-1; b++){
                        for (int c=1; c<nvz-1; c++){
                            if (vz[c] > 0){
                                if(abs(vz[c]*deltaT/hz)>1) cout << "z"<< vz[c]*deltaT/hz << endl;
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] - (deltaT * vz[c] / hz)
                                 * (f_time[i][j][k][a][b][c] - f_time[i][j][(nz+k-1)%nz][a][b][c]);
                            }   
                            else{
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] - (deltaT * vz[c] / hz)
                                 * (f_time[i][j][(k+1)%nz][a][b][c] - f_time[i][j][k][a][b][c]);
                            }   
                        }
                    }
                }
            }
        }
    }

    //cout << "***  Finalizing coordinate part***" << endl;

}



void kinetic::VelocityPart(){

   double fx_, fy_, fz_;

   SaveState();

   for (int i=0; i<nx; i++){    
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                fx_ = Converter->GetSAxelerationCofficient() * asx[i][j][k] + al + acx[i][j][k];
                if(fx_ > 0){
                    if (abs(fx_*deltaT/hvx)>1) cout << "vx "<< fx_*deltaT/hvx << endl;
                    for (int a=1; a<nvx-1; a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] 
                                 - (deltaT * fx_ / hvx) * (f_time[i][j][k][a][b][c] - f_time[i][j][k][a-1][b][c]);
                            }
                        }
                    }
                }
                else
                {
                    for (int a=1; a<nvx-1; a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] 
                                 - (deltaT * fx_ / hvx) * (f_time[i][j][k][a+1][b][c] - f_time[i][j][k][a][b][c]);
                            }
                        }
                    }
                }  
            }
        }
    }  

    SaveState();

    for (int i=0; i<nx; i++){    
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                fy_ = Converter->GetSAxelerationCofficient() * asy[i][j][k] + acy[i][j][k];
                if(fy_ > 0){
                    if (abs(fy_*deltaT/hvy)>1) cout << "vy "<< fy_*deltaT/hvy << endl;
                    for (int a=1; a<nvx-1; a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] 
                                 - (deltaT * fy_ / hvy) * (f_time[i][j][k][a][b][c] - f_time[i][j][k][a][b-1][c]);
                            }
                        }
                    }
                }
                else
                {
                    for (int a=1; a<nvx-1; a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f[i][j][k][a][b][c] 
                                 - (deltaT * fy_ / hvy) * (f[i][j][k][a][b+1][c] - f[i][j][k][a][b][c]);
                            }
                        }
                    }
                }  
            }
        }
    }

    SaveState();

    for (int i=0; i<nx; i++){    
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                fz_ = Converter->GetSAxelerationCofficient() * asz[i][j][k] + acz[i][j][k];
                if(fz_ > 0){
                    if (abs(fz_*deltaT/hvz)>1) cout << "vz "<< fz_*deltaT/hvz << endl;
                    for (int a=1; a<nvx-1; a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] 
                                 - (deltaT * fz_ / hvz) * (f_time[i][j][k][a][b][c] - f_time[i][j][k][a][b][c-1]);
                            }
                        }
                    }
                }
                else
                {
                    for (int a=1; a<nvx-1; a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f[i][j][k][a][b][c] = f_time[i][j][k][a][b][c] 
                                 - (deltaT * fz_ / hvz) * (f_time[i][j][k][a][b][c+1] - f_time[i][j][k][a][b][c]);
                            }
                        }
                    }
                }  
            }
        }
    }

}

double kinetic::Relaxation(int i, int j, int k, int a, int b, int c){
    return f[i][j][k][a][b][c] + deltaT * (density[i][j][k] * f_n[a][b][c] - f[i][j][k][a][b][c]) / tau;
}


void kinetic::IntegrateAll(){

    CoordinatePart();

    VelocityPart();
    
    for (int i=0; i<nx; i++){    
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                for (int a=1; a<nvx-1; a++){
                    for (int b=1; b<nvy-1; b++){
                        for (int c=1; c<nvy-1; c++){
                            f[i][j][k][a][b][c]=Relaxation(i,j,k,a,b,c); 
                        }
                    }
                }
            }
        }
    } 

    CurT+=deltaT;
    
    ComputeDensity();

    ComputeFlowVelocity();

}

void kinetic::ComputeDensity(){

    //cout << "*** Compute density ***" << endl;

    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                density[i][j][k] = 0;
                for (int a=0; a<nvx; a++){
                    for (int b=0; b<nvy; b++){
                        for (int c=0; c<nvz; c++){
                            density[i][j][k] = density[i][j][k] + f[i][j][k][a][b][c];
                        }
                    }
                } 
                density[i][j][k] = density[i][j][k] * hvx * hvy * hvz;
            }
        }
    }

    


}


void kinetic::ComputeFlowVelocity(){

    //cout << "*** Compute flow velocity ***" << endl;

    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                vfl_x[i][j][k] = 0;
                vfl_y[i][j][k] = 0;
                vfl_z[i][j][k] = 0;
                for (int a=0; a<nvx; a++){
                    for (int b=0; b<nvy; b++){
                        for (int c=0; c<nvz; c++){
                            vfl_x[i][j][k] = vfl_x[i][j][k] + vx[a] * f[i][j][k][a][b][c];
                            vfl_y[i][j][k] = vfl_y[i][j][k] + vy[b] * f[i][j][k][a][b][c];
                            vfl_z[i][j][k] = vfl_z[i][j][k] + vz[c] * f[i][j][k][a][b][c];
                        }
                    }
                } 
                vfl_x[i][j][k] = vfl_x[i][j][k] * hvx * hvy * hvz;
                vfl_y[i][j][k] = vfl_y[i][j][k] * hvx * hvy * hvz;
                vfl_z[i][j][k] = vfl_z[i][j][k] * hvx * hvy * hvz;
            }
        }
    }

}


void kinetic::InitialConditions(){
ComputeAnalyticalFlowVelocity();

    for (int i=0; i<1; i++){
        for (int j=0; j<1; j++){
            for (int k=0; k<1; k++){
                for (int a=0; a<nvx; a++){
                    for (int b=0; b<nvy; b++){
                        for (int c=0; c<nvz; c++){
			    float Xi;
        	            Xi = 0;
         		    f[i][j][k][a][b][c] = 0;
        		    while(Xi < maxXi){
         		        f[i][j][k][a][b][c] = f[i][j][k][a][b][c]
                                +  exp(-Xi) * exp(- 0.5 * (vx[a] - Xi * v_fl_an) * (vx[a] - Xi * v_fl_an)
                                - 0.5 * vy[b] * vy[b] - 0.5 * vz[c] * vz[c]);
                                Xi = Xi + deltaXi;
                            }
        		    f[i][j][k][a][b][c] = Converter->GetDimensionlessIonConcentration() * f[i][j][k][a][b][c] * deltaXi
             		     /pow(2*pi,1.5);
                        }
                    }
                }
            } 
        }
    }

    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                for (int a=0; a<nvx; a++){
                    for (int b=0; b<nvy; b++){
                        for (int c=0; c<nvz; c++){
         		    f[i][j][k][a][b][c] = f[0][0][0][a][b][c];
                        }
                    }
                }
            } 
        }
    }
/*
    for (int a=0; a<nvx; a++){
        for (int b=0; b<nvy; b++){
            for (int c=0; c<nvz; c++){
	        f[nx/2][ny/2][nz/2][a][b][c] += 0.1 * Converter->GetDimensionlessIonConcentration() * f_n[a][b][c];
	    }
        }
    }
*/
}

void kinetic::ComputeAnalyticalFlowVelocity(){
    v_fl_an = al * tau ;
}


double*** kinetic::GetDensity() const{
    return density;
}

int kinetic::GetNvx() const{
    return nvx;
}

int kinetic::GetNvy() const{
    return nvy;
}

int kinetic::GetNvz() const{
    return nvz;
}

double kinetic::Gethvx() const{
    return hvx;
}

double kinetic::Gethvy() const{
    return hvy;
}

double kinetic::Gethvz() const{
    return hvz;
}

double kinetic::GetCutVelocityAlong() const{
    return vxcut;
}

double kinetic::GetCutVelocityNormal() const{
    return vyzcut;
}

double kinetic::GetCurrentTime() const{
    return CurT;
}


double*** kinetic::GetFlowVelocityX() const{
    return vfl_x;
}
            

double*** kinetic::GetFlowVelocityY() const{
    return vfl_y;
}


double*** kinetic::GetFlowVelocityZ() const{
    return vfl_z;
}


double* kinetic::GetProfileX() const{
    return profile_x;
}


double* kinetic::GetProfileY() const{
    return profile_y;
}


double* kinetic::GetProfileZ() const{
    return profile_z;
}

double* kinetic::GetAnalyticalProfileX() const{
    return analytical_profile_x;
}


void kinetic::SetProfileX(int i, int j, int k){
    for(int a=0; a<nvx; a++){
        profile_x[a] = f[i][j][k][a][nvy/2][nvz/2];
    }
    
}


void kinetic::SetProfileY(int i, int j, int k){
    for(int b=0; b<nvy; b++){
        profile_y[b] = f[i][j][k][nvx/2][b][nvz/2];
    }
}


void kinetic::SetProfileZ(int i, int j, int k){
    for(int c=0; c<nvz; c++){
        profile_z[c] = f[i][j][k][nvx/2][nvy/2][c];
    }
}

void kinetic::SetAnalyticalProfileX(){
    float Xi;
    for(int a=0; a<nvx; a++){
        analytical_profile_x[a] = 0;
        Xi = 0;
        while(Xi < maxXi){
            analytical_profile_x[a] = analytical_profile_x[a]
             +  exp(-Xi) * exp(- 0.5 * (vx[a] - Xi * v_fl_an) * (vx[a] - Xi * v_fl_an)
             - 0.5 * vy[nvy/2] * vy[nvy/2] - 0.5 * vz[nvz/2] * vz[nvz/2]);
            Xi = Xi + deltaXi;
        }
        analytical_profile_x[a] = Converter->GetDimensionlessIonConcentration() * analytical_profile_x[a] * deltaXi
        / pow(2*pi,1.5);
    }
}


double kinetic::GetDeltaT() const{
    return deltaT;
}

void kinetic::SaveState(){

    for (int i=0; i<nx; i++){    
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                for (int a=0; a<nvx; a++){
                    for (int b=0; b<nvy; b++){
                        for (int c=0; c<nvz; c++){
                            f_time[i][j][k][a][b][c]=f[i][j][k][a][b][c]; 
                        }
                    }
                }
            }
        }
    } 

}


double* kinetic::GetVelocitySetX() const{
    return vx;
}
           
double* kinetic::GetVelocitySetY() const{
    return vx;
}
            
double* kinetic::GetVelocitySetZ() const{
    return vx;
}


void kinetic::DefineDensity() {

    double n_i = Converter->GetDimensionlessIonConcentration();

    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                density[i][j][k] = ((pi/Lx)*(pi/Lx)+(pi/Lx)*(pi/Lx)+(pi/Lz)*(pi/Lz)) 
                 * sin(pi*x[i]/Lx)*sin(pi*y[j]/Ly)*sin(pi*z[k]/Lz) / (4. * pi) + n_i ;
            }
        }
    }

}
