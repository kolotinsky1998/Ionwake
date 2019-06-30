#pragma once

///This class is used to save and convert physical values  
class converter{
    public:
            ///constuctor
            converter();
            ///constuctor
            converter(double T_, double tau_, double n_0_, double q_, 
                      double El_, double Lx_, double Ly_, double Lz_, int Nx_, int Ny_, int Nz_);
            ///destructor
            ~converter();
            ///get cofficient near the self-consistent axeleration term
            double GetSAxelerationCofficient()const ;
            ///get cofficient near the longitudinal axeleration term
            double GetLAxelerationCofficient()const ;
            ///get cofficient near the Coloumb axeleration term
            double GetCAxelerationCofficient()const ;
            ///convertion cofficent for dimensionless potential
            double ConvertPotential()const ;
            ///convert physical mean free time in dimensionless units
            double GetDimensionlessTau()const ;
            ///get dimensionless flow velocity
            double GetDimensionlessFlowVelocity()const ;
            ///get dimensionless ion concentration
            double GetDimensionlessIonConcentration()const ;
            ///get coordinate step X
            double GetCoordinateStepX()const ;
            ///get coordinate step Y
            double GetCoordinateStepY()const ;
            ///get coordinate step Z
            double GetCoordinateStepZ()const ;
            ///get grid size X
            int GetNx()const ;
            ///get grid size Y
            int GetNy()const ;
            ///get grid size Z
            int GetNz()const ;
            ///get temperature
            double GetTemperature()const ;
            ///get thermal velocity
            double GetThermalVelocity()const ;
            ///get plasmas frequency
            double GetPlasmasFrequency() const;
            ///get Physical equilibrium ion concentration
            double GetPhysicalIonConcentration() const ;
            ///get particle charge
            double GetParticleCharge()const ;
            ///get external electricity strength
            double GetStrength()const ;
            ///get physical relaxation time
            double GetPhysicalTau() const;
            ///get Debay radious 
            double GetRd() const ;


            



    private:
    /******* some fundamental quantities ******/
            //electron charge 
            double e;
            //ion mass
            double m;
            //Boltzmann constant
            double k;
            //pi number
            double pi;

    /******* system quantities in physical units *******/
            //temperature of neutrals
            double T;
            //debay radious
            double rd;
            ///plasmas frequency
            double wp;
            //mean free time
            double tau;
            //ion equillibrium concentration
            double n_0;
            //dust particle charge
            double q;
            //External electricity force
            double El;
    /******* derivative quantities in physical units *******/
            //thermal velocity of neutrals
            double vt;
            //cofficient near the self-consistent axeleration term
            double as;
            //cofficient near the longitudinal axeleration term
            double al;
            ///cofficient near the Coloumb axeleration term
            double ac;
            //cofficent to сonvert potential from dimensionless units to physical
            double potential;
    /******* dimensionless parameters *******/
            //size of computational box in debay units X
            double Lx;
            //size of computational box in debay units Y
            double Ly;
            //size of computational box in debay units Z
            double Lz;
            //grid size X
            int Nx;
            //grid size Y
            int Ny;
            //grid size Z
            int Nz;
            //dimensionless relaxation time
            double tau_d;
            //dimensionless ion concentration
            double n_0_d;
            //flow velocity
            double v_fl_an_d;
            //coordinate discretization step X
            double hx;
            //coordinate discretization step Y
            double hy;
            //coordinate discretization step Z
            double hz;


};