#define _USE_MATH_DEFINES
#include "converter.hpp"
#include "scheme.hpp" 
#include <cmath>
#include <algorithm>
#include "nik-mdarray-hpp/nik-mdarray6.hpp"
#include <iostream>
#include <sstream> 
#include <fstream> 
#include "string"
#include "omp.h"
using namespace std;
Scheme::Scheme(Converter converter_, int nx_, int ny_, int nz_, 
                double lx_, double ly_, double lz_,
                int nt_, double r0x_, double r0y_, double r0z_,
                double hvx_=0.2, double hvy_=0.2, double hvz_=0.2,
                double vminyz_=-5, double vmaxyz_=5,
                 double hxi_=0.002, double ximax_=13):
    converter(converter_),
    nx(nx_), ny(ny_), nz(nz_),
    lx(lx_), ly(ly_), lz(lz_),
    nt(nt_), r0x(r0x_), r0y(r0y_), r0z(r0z_),
    hvx(hvx_), hvy(hvy_), hvz(hvz_),
    vminyz(vminyz_),
    vmaxyz(vmaxyz_),
    hxi(hxi_),ximax(ximax_)
    {
        Eext = converter.externalField(); 
        wc = converter.collisionFrequency(); 
        vfl = converter.flowVelocity(); 
        rde = converter.electronDebyeRadious(); 
        q = converter.dustParticleCharge(); 

        hx = lx / double(nx);
        hy = ly / double(ny);
        hz = lz / double(nz);

        vminx = vminyz;
        vmaxx = vmaxxCompute();
        
        nvx = int((vmaxx-vminx) / hvx);
        nvy = int((vmaxyz-vminyz) / hvy);
        nvz = int((vmaxyz-vminyz) / hvz);

        dt = 0.5 * min(hvx/(Eext+(q/(hx*hx))), hx/vmaxx);
        t = 0;

        rx = new double[nx];
        for (int i=0; i<nx; i++){
            rx[i] = i*hx;
        }
        ry = new double[ny];
        for (int j=0; j<ny; j++){
            ry[j] = j*hy;
        }
        rz = new double[nz];
        for (int k=0; k<nz; k++){
            rz[k] = k*hz;
        }
        vx = new double [nvx];
        for (int a=0; a<nvx; a++){
            vx[a] = vminx + a*hvx;
        }
        vy = new double [nvy];
        for (int b=0; b<nvy; b++){
            vy[b] = vminyz + b*hvy;
        }
        vz = new double [nvz];
        for (int c=0; c<nvz; c++){
            vz[c] = vminyz + c*hvz;
        }
        
        fullCharge = 0;

        f.setdimensions(nx,ny,nz,nvx,nvy,nvz);
        ftime.setdimensions(nx,ny,nz,nvx,nvy,nvz);
        n.setdimensions(nx,ny,nz);
        fi.setdimensions(nx,ny,nz);
        ax.setdimensions(nx,ny,nz);
        ay.setdimensions(nx,ny,nz);
        az.setdimensions(nx,ny,nz);
        maxvell.setdimensions(nvx,nvy,nvz);
        
        n.fill(1.0);

        #pragma omp parallel for collapse(3)
        for (int a=0; a<nvx; a++){
            for (int b=0; b<nvy; b++){
                for (int c=0; c<nvz; c++){
                    f(0,0,0,a,b,c) = initialDisrtibutionFunction(vx[a], vy[b], vz[c]);
                }
            }
        }

        #pragma omp parallel for collapse(3)
        for (int i=0; i<nx; i++){
            for (int j=0; j<ny; j++){
                for (int k=0; k<nz; k++){
                    for (int a=0; a<nvx; a++){
                        for (int b=0; b<nvy; b++){
                            for (int c=0; c<nvz; c++){
                                f(i,j,k,a,b,c) = f(0,0,0,a,b,c);
                                ftime(i,j,k,a,b,c) = f(0,0,0,a,b,c);
                            }
                        }
                    }
                }
            }
        }

        #pragma omp parallel for collapse(3)
        for (int a=0; a<nvx; a++){
            for (int b=0; b<nvy; b++){
                for (int c=0; c<nvz; c++){
                    maxvell(a,b,c) = maxwellianDistribution(vx[a], vy[b], vz[c]);
                }
            }
        }
        
    }

Scheme::~Scheme(){
    delete [] rx;
    delete [] ry;
    delete [] rz;
    delete [] vx;
    delete [] vy;
    delete [] vz;
}

double Scheme::vmaxxCompute() {
    double vmax;
    double distribution;
    double maxvell;
    double xi;
    maxvell = exp(-vmaxyz*vmaxyz*0.5);
    vmax = 0;
    while(true){
	    vmax = vmax + hvx;
        distribution = 0;
        xi = 0;
        while(xi < ximax){
            distribution = distribution
             +  exp(-xi) * exp(- 0.5 * (vmax - xi * vfl) * (vmax - xi * vfl))*hxi;
            xi = xi + hxi;
        }
	    if(distribution<maxvell){
	        break;
	    }
    }
    return vmax;
}

double Scheme::initialDisrtibutionFunction(double vx, double vy, double vz){
    double distribution;
    double xi;
    distribution = 0;
    xi = 0;
    while(xi < ximax){
        distribution = distribution
        + exp(-xi) * exp(- 0.5 * (vx - xi * vfl) * (vx - xi * vfl)) * exp(-0.5*vy*vy) *exp(-0.5*vz*vz)*hxi;
        xi = xi + hxi;
    }
    distribution = distribution / pow(2*M_PI,1.5);
    return distribution;
}

double Scheme::maxwellianDistribution(double vx, double vy, double vz){
    double maxvell;
    maxvell = exp(-0.5*vx*vx) * exp(-0.5*vy*vy) * exp(-0.5*vz*vz) / pow(2*M_PI,1.5);
    return maxvell;
}

void Scheme::density(){
    double density;
    #pragma omp parallel for collapse(3)
    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                density = 0;
                for (int a=0; a<nvx; a++){
                    for (int b=0; b<nvy; b++){
                        for (int c=0; c<nvz; c++){
                            density = density + f(i,j,k,a,b,c);
                        }
                    }
                }
                n(i,j,k) = density*hvx*hvy*hvx;
            }
        }
    }
}

double Scheme::potentialDebye(double rx, double ry, double rz){
    double r;
    r = distance(rx, ry, rz, r0x, r0y, r0z);
    if (r<sqrt(hx*hx+hy*hy+hz*hz)){
        return -q * exp(-sqrt(hx*hx+hy*hy+hz*hz)/rde)/sqrt(hx*hx+hy*hy+hz*hz);
    } else {
        return -q * exp(-r/rde)/r;
    }
}

double Scheme::distance(double r1x, double r1y, double r1z, double r2x, double r2y, double r2z){
    return sqrt((r1x-r2x)*(r1x-r2x) + (r1y-r2y)*(r1y-r2y) + (r1z-r2z)*(r1z-r2z));
}

void Scheme::force(){
    #pragma omp parallel for collapse(3)
    for ( int i=1; i<nx-1; i++){
        for ( int j=1; j<ny-1; j++){
            for (int k=1; k<nz-1; k++){
                ax(i,j,k) = -(fi(i+1,j,k)-fi(i-1,j,k))/(2.0*hx) + Eext;
                ay(i,j,k) = -(fi(i,j+1,k)-fi(i,j-1,k))/(2.0*hy);
                az(i,j,k) = -(fi(i,j,k+1)-fi(i,j,k-1))/(2.0*hz);
            }
        }
    }
}

void Scheme::kinetic(){
    ///coordinate propagation part
    #pragma omp parallel for collapse(3)
    for ( int i=1; i<nx-1; i++){
        for ( int j=1; j<ny-1; j++){
            for (int k=1; k<nz-1; k++){
                for (int a=1; a<nvx-1; a++){
                    if(vx[a]>0){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f(i,j,k,a,b,c) = ftime(i,j,k,a,b,c) - vx[a]*(ftime(i,j,k,a,b,c) - ftime(i-1,j,k,a,b,c))*dt/hx;
                            }
                        }
                    } else{
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f(i,j,k,a,b,c) = ftime(i,j,k,a,b,c) - vx[a]*(ftime(i+1,j,k,a,b,c) - ftime(i,j,k,a,b,c))*dt/hx;
                            }
                        }
                    }
                }
            }
        }
    }
    #pragma omp parallel for collapse(3)
    for ( int i=1; i<nx-1; i++){
        for ( int j=1; j<ny-1; j++){
            for (int k=1; k<nz-1; k++){
                for (int a=1; a<nvx-1; a++){
                    for (int b=1; b<nvy-1; b++){
                        if(vy[b]>0){
                            for (int c=1; c<nvz-1; c++){
                                ftime(i,j,k,a,b,c) = f(i,j,k,a,b,c) - vy[b]*(f(i,j,k,a,b,c) - f(i,j-1,k,a,b,c))*dt/hy;
                            }
                        } else{
                            for (int c=1; c<nvz-1; c++){
                                ftime(i,j,k,a,b,c) = f(i,j,k,a,b,c) - vy[b]*(f(i,j+1,k,a,b,c) - f(i,j,k,a,b,c))*dt/hy;
                            }
                        }
                    }
                }
            }
        }
    }
    #pragma omp parallel for collapse(3)
    for ( int i=1; i<nx-1; i++){
        for ( int j=1; j<ny-1; j++){
            for (int k=1; k<nz-1; k++){
                for (int a=1; a<nvx-1; a++){
                    for (int b=1; b<nvy-1; b++){
                        for (int c=1; c<nvz-1; c++){
                            if(vz[c]>0){
                                f(i,j,k,a,b,c) = ftime(i,j,k,a,b,c) - vz[c]*(ftime(i,j,k,a,b,c) - ftime(i,j,k-1,a,b,c))*dt/hz;
                            } else {
                                f(i,j,k,a,b,c) = ftime(i,j,k,a,b,c) - vz[c]*(ftime(i,j,k+1,a,b,c) - ftime(i,j,k,a,b,c))*dt/hz;
                            }
                        }
                    }
                }
            }
        }
    }
    ///velocity propagation part
    #pragma omp parallel for collapse(3)
    for ( int i=1; i<nx-1; i++){
        for ( int j=1; j<ny-1; j++){
            for (int k=1; k<nz-1; k++){
                if (ax(i,j,k)>0){
                    for (int a=1; a<nvx-1;a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                ftime(i,j,k,a,b,c) = f(i,j,k,a,b,c) - ax(i,j,k)*(f(i,j,k,a,b,c) - f(i,j,k,a-1,b,c))*dt/hvx;
                            }
                        }
                    } 
                } else {
                    for (int a=1; a<nvx-1;a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                ftime(i,j,k,a,b,c) = f(i,j,k,a,b,c) - ax(i,j,k)*(f(i,j,k,a+1,b,c) - f(i,j,k,a,b,c))*dt/hvx;
                            }
                        }
                    } 
                }
            }
        }
    }
    #pragma omp parallel for collapse(3)
    for ( int i=1; i<nx-1; i++){
        for ( int j=1; j<ny-1; j++){
            for (int k=1; k<nz-1; k++){
                if (ay(i,j,k)>0){
                    for (int a=1; a<nvx-1;a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f(i,j,k,a,b,c) = ftime(i,j,k,a,b,c) - ay(i,j,k)*(ftime(i,j,k,a,b,c) - ftime(i,j,k,a,b-1,c))*dt/hvy;
                            }
                        }
                    } 
                } else {
                    for (int a=1; a<nvx-1;a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                f(i,j,k,a,b,c) = ftime(i,j,k,a,b,c) - ay(i,j,k)*(ftime(i,j,k,a,b+1,c) - ftime(i,j,k,a,b,c))*dt/hvy;
                            }
                        }
                    } 
                }
            }
        }
    }
    #pragma omp parallel for collapse(3)
    for ( int i=1; i<nx-1; i++){
        for ( int j=1; j<ny-1; j++){
            for (int k=1; k<nz-1; k++){
                if (az(i,j,k)>0){
                    for (int a=1; a<nvx-1;a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                ftime(i,j,k,a,b,c) = f(i,j,k,a,b,c) - az(i,j,k)*(f(i,j,k,a,b,c) - f(i,j,k,a,b,c-1))*dt/hvz;
                            }
                        }
                    } 
                } else {
                    for (int a=1; a<nvx-1;a++){
                        for (int b=1; b<nvy-1; b++){
                            for (int c=1; c<nvz-1; c++){
                                ftime(i,j,k,a,b,c) = f(i,j,k,a,b,c) - az(i,j,k)*(f(i,j,k,a,b,c+1) - f(i,j,k,a,b,c))*dt/hvz;
                            }
                        }
                    } 
                }
            }
        }
    }
    ///Collision part
    #pragma omp parallel for collapse(3)
    for ( int i=1; i<nx-1; i++){
        for ( int j=1; j<ny-1; j++){
            for (int k=1; k<nz-1; k++){
                for (int a=1; a<nvx-1;a++){
                    for (int b=1; b<nvy-1; b++){
                        for (int c=1; c<nvz-1; c++){
                            f(i,j,k,a,b,c) = ftime(i,j,k,a,b,c) + wc*(maxwellianDistribution(vx[a],vy[b],vz[c])*n(i,j,k) - ftime(i,j,k,a,b,c))*dt;
                        }
                    }
                }
            }
        }
    }
}

void Scheme::poisson(){
    double potential;
    double r;
    #pragma omp parallel for collapse(3)
    for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                potential = potentialDebye(rx[i], ry[j], rz[k]);
                for (int i_=0; i_<nx; i_++){
                    for (int j_=0; j_<ny; j_++){
                        for (int k_=0; k_<nz; k_++){
                            r = distance(rx[i],ry[j],rz[k],rx[i_],ry[j_],rz[k_]);
                            if(r!=0){
                                potential = potential + exp(-r/rde)*(n(i_,j_,k_)-1.)*hx*hy*hz/(r*4*M_PI);
                            }
                        }
                    }
                }
                fi(i,j,k) = potential;
            }
        }
    }
}

void Scheme::schemeStep(){
    poisson();
    force();
    kinetic();
    density();
    calculateFullCharge();
    t = t + dt;
}

void Scheme::calculateFullCharge(){
    fullCharge = 0;
    #pragma omp parallel for collapse(3)
    for ( int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            for (int k=0; k<nz; k++){
                fullCharge = fullCharge + (n(i,j,k) - 1.0)*hx*hy*hz;
            }
        }
    }
}

void Scheme::printCurrentTime(){
    cout << "Current dimmensionless time: " << t << endl;
}

void Scheme::printFullCharge(){
    cout << "Current full charge in the computational box: " << fullCharge << endl;
}

void Scheme::writeDensityFile(string data="./"){
    stringstream filename;
    filename << data << "/density_t" << t << ".dat";
    ofstream file;
    file.open(filename.str().c_str());
    file << nx << "\t" << ny << "\t" << nz << "\n";
    for (int k=0; k<nz; k++){
        for (int j=0; j<ny; j++){
            for (int i=0; i<nx; i++){
                file << n(i,j,k) - 1.0 << "\n";
            }
        }
    }
    file.close();
}

void Scheme::writePotentialFile(string data="./"){
    stringstream filename;
    filename << data << "/potential_t" << t << ".dat";
    ofstream file;
    file.open(filename.str().c_str());
    file << nx << "\t" << ny << "\t" << nz << "\n";
    for (int k=0; k<nz; k++){
        for (int j=0; j<ny; j++){
            for (int i=0; i<nx; i++){
                file << fi(i,j,k) << "\n";
            }
        }
    }
    file.close();
}

void Scheme::InitialLogOut(){
    cout << "#########################################################" << endl;
    cout << "##******* Output provided by numerical scheme *********##" << endl;
    cout << "#########################################################" << endl;
    cout << "Particle charge: " << q <<endl;
    cout << "External electricity field: " << Eext << endl;
    cout << "Ion flow velocity: " << vfl << endl;
    cout << "Ion-neutral collision frequency: " << wc << endl;
    cout << "Debye electron radious: " << rde << endl;
    cout << "Size of the computaional box:" << endl;
    cout << "nx:  " << nx << endl;
    cout << "ny:  " << ny << endl;
    cout << "nz:  " << nz << endl;
    cout << "nvx: " << nvx << endl;
    cout << "nvy: " << nvy << endl;
    cout << "nvz: " << nvz << endl;
    cout << "Sheme integration steps:" << endl;
    cout << "hx:  " << hx << endl;
    cout << "hy:  " << hy << endl;
    cout << "hz:  " << hz << endl;
    cout << "hvx: " << hvx << endl;
    cout << "hvy: " << hvy << endl;
    cout << "hvz: " << hvz << endl;
    cout << "dt:  " << dt << endl;
    cout << "Cut velocities in velocity space:" << endl;
    cout << "vmaxyz: " << vmaxyz << endl;
    cout << "vmaxx:  " << vmaxx << endl;
    cout << "vminyz: " << vminyz << endl;
    cout << "vminz: " << vminx << endl;
    cout << endl;
}
