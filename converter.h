#pragma once

///This class is used to save and convert physical values  
class converter{
    public:
            ///constuctor
            converter();
            ///constuctor
            converter(double T_, double tau_, double n_0_, double q_, 
                      double El_, double Lx_, double Ly_, double Lz_, int Nx_, int Ny_, int Nz_, int Nx_0, int Ny_0, int Nz_0);
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
            ///get position at the grid of the duct particle X
            int GetNx_0()const ;
            ///get position at the grid of the duct particle Y
            int GetNy_0()const ;
            ///get position at the grid of the duct particle Z
            int GetNz_0()const ;
            ///get length of computational box X
            double GetLx()const ;
            ///get length of computational box Y
            double GetLy()const ;
            ///get length of computational box Z
            double GetLz()const ;
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
            ///get elementary charge
            double GetElementaryCharge() const;
            ///get longitudinal electricity field
            double GetEl() const;
            ///get value of pi constant
            double GetPi() const;


            



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
            //temperature of neutrals in Kelvin
            double T;
            //debay radious [cm]
            double rd;
            ///plasmas frequency  [Hz]
            double wp;
            //mean free time [s]
            double tau;
            //ion equillibrium concentration [1/cm^3] 
            double n_0;
            //External electricity field [gauss units]
            double El;
            //thermal velocity of neutrals [cm/s]
            double vt;
  /******* The coefficient between physical and dimmension units *******/
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
            //position at the grid of the dust particle X
            int Nx_0;
            //position at the grid of the dust particle Y
            int Ny_0;
            //position at the grid of the dust particle Z
            int Nz_0;
            //grid size X
            int Nx;
            //grid size Y
            int Ny;
            //grid size Z
            int Nz;
            // dust particle charge
            double q;
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