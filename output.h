#pragma once
#include "kinetic.h"
#include "poisson.h"
#include "converter.h"


///This class provides output options
class output{
    public:
            ///constructor
            output(const converter &, const poisson &, kinetic &);
            ///destructor
            ~output();
            ///output initial simulation parameter
            void StartOutput();
            ///write down current simulation parameters in the output
            void TerminalOutput(int term_time);
            ///write down in the VTK file denstity and potential
            void VTKoutput(int output_time);
            /// Clean directory with data files
            void CleanData();
            /// Plot distribution function inside the  (i, j, k)-cell
            void PlotDistributionFunction(int i, int j, int k, int output_time);
            /// Write Coloumb force field
            void WriteColoumbForceField();

    private:
            // access to converter members
            const converter *Converter;
            // access to poisson members
            const poisson *Poisson;
            // access to kinetic members
            kinetic *Kinetic;
            //current time
            int Time;
            
};
