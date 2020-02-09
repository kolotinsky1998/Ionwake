#define _USE_MATH_DEFINES

#include "converter.hpp"
#include <cmath>
#include <iostream>

using namespace std;

Converter::Converter(double Te_, double Ti_, double ni_, double q_, double Eext_, double wc_, double mi_) :
        Te(Te_), Ti(Ti_), ni(ni_), q(q_), Eext(Eext_), wc(wc_), mi(mi_) {
    e = 4.8032 * pow(10., -10.);
    kB = 1.3806 * pow(10., -16.);

    rdi = sqrt(kB * Ti / (4.0 * M_PI * ni * e * e));
    rde = sqrt(kB * Te / (4.0 * M_PI * ni * e * e));
    wpi = sqrt(4.0 * M_PI * ni * e * e / mi);
    vT = sqrt(kB * Ti / mi);
    vfl = e * Eext / (mi * wc);

    Eext_d = Eext * e * rdi / (kB * Ti);
    wc_d = wc / wpi;
    vfl_d = vfl / vT;
    rde_d = rde / rdi;
    q_d = q / (4.0 * M_PI * ni * rdi * rdi * rdi);

}

Converter::~Converter() {

}

double Converter::externalField() {
    return Eext_d;
}

double Converter::collisionFrequency() {
    return wc_d;
}

double Converter::flowVelocity() {
    return vfl_d;
}

double Converter::electronDebyeRadious() {
    return rde_d;
}

double Converter::dustParticleCharge() {
    return q_d;
}

void Converter::InitialLogOut() {
    cout << "#########################################################" << endl;
    cout << "##******* Output provided by converter ****************##" << endl;
    cout << "#########################################################" << endl;
    cout << "Temperature of electrons: " << Te << " [K]" << endl;
    cout << "Temperature of ions: " << Ti << " [K]" << endl;
    cout << "Ion concentration: " << ni << " [1/cm^3]" << endl;
    cout << "Charge of the dust particle: " << q << " [electron charges]" << endl;
    cout << "External electricity field: " << Eext << " [cgs units]" << endl;
    cout << "Ion-neutral collision frequency: " << wc << " [1/s]" << endl;
    cout << "Mass of ions: " << mi << " [g]" << endl;
    cout << "Debye ion radious: " << rdi << " [cm]" << endl;
    cout << "Debye electron radious: " << rde << " [cm]" << endl;
    cout << "Plasma ion frequency: " << wpi << " [1/s]" << endl;
    cout << "Thermal velocity of neutrals: " << vT << " [K]" << endl;
    cout << "Velocity of the ion flow: " << vfl << " [cm/s]" << endl;
    cout << endl;
}
