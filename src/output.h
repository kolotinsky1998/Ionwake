#pragma once

#include "kinetic.h"
#include "poisson.h"
#include "converter.h"
#include "string"

using namespace std;


///This class provides output options
class output {
public:
    ///constructor
    output(const converter &, const poisson &, kinetic &);

    ///output initial simulation parameter
    void StartOutput();

    ///write down current simulation parameters in the output
    void TerminalOutput(int term_time);
    ///This function writes down data.
    /*!
    The output files contain ion concentration, self-consistent potential and momenta flow.
    All these files will be avaible in created data directory.
     */
    void VTKoutput(int output_time);

    /// Plot distribution function inside the  (i, j, k)-cell
    void PlotDistributionFunction(int i, int j, int k, int output_time);
    /// \brief This function creates data file with initial potential around the dust particle.
    /*!
    The initial potential is an external homogenious field potential and coloumb potential due to the dust particle.
    The potential is written in sgc units.
    \f[
        \phi = - E_{ext}z - \frac{Q}{r}
    \f]
    */
    void WriteInitialPotential();
    ///\brief This founction creates a directory in which all output files will be saved.
    /*!
    The name of the directory define the user.
     */
    void CreateOutputDirectory();

    ///\brief This function updates class member "Time" which is reasponseble for output
    void UpdateTime();
    ///\brief This function is responsible for output of distribution function
    /*!
     The output is nessesary to continue calcalution after the reaching the permitted time for the job on queue.
    The distribution from the output will be used as initial conditions for further calculation.
     */
    void WriteDistribution();

private:
    // access to converter members
    const converter *Converter;
    // access to poisson members
    const poisson *Poisson;
    // access to kinetic members
    kinetic *Kinetic;
    //current time
    int Time;
    /// directory for the output data
    string data;
};
