#ifndef converter_h
#define converter_h

class Converter {
public:
    Converter(double Te_, double Ti_, double ni_, double q_, double Eext, double wc_, double mi_);

    ~Converter();

    double externalField();

    double collisionFrequency();

    double flowVelocity();

    double electronDebyeRadious();

    double dustParticleCharge();

    void InitialLogOut(); ///write simulation parameters into the log file

private:
    ///Physical parameters of the system in dimmensional units
    double Eext_d; ///external electricity field
    double wc_d; ///collision frequency
    double vfl_d; ///velocity of the ion flow
    double rde_d; ///electrons' Debye radious
    double q_d; ///charge of the dust particle

    ///Physical parameters of the system in physical units
    double Te; ///electron temperature
    double Ti; ///ion temperature
    double ni; ///ion concentration
    double q; ///charge of the dust particle
    double Eext; ///external electricity field
    double wc; ///ion-neutral collisions frequency
    double mi; ///ion mass
    double rdi; ///ion Debye radious
    double rde; ///electron Debye radious
    double wpi; ///ion plasma frequency
    double vT; ///termal velocity of neutrals
    double vfl; ///velocity of the ion flow

    ///Fundamental constants
    double kB; ///Boltzmann constant
    double e; ///electron charge
};

#endif