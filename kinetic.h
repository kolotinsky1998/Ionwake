#pragma once

#include "converter.h"
#include "poisson.h"


///\brief This class is used to solve kinetic equation for ions numericaly
/*!
\f[
\frac{\partial f}{\partial t} + \overline{v}\bigtriangledown f+\frac{\overline{F}}{m}\frac{\partial f}{\partial \overline{v}}=St[f]
\f]
\f[
St[f] = \nu (f - n_i \Phi_M)
\f]
\f[
\Phi_M = \frac{1}{(2 \pi v_T^2)^{\frac{3}{2}}} \exp(-\frac{v^2}{2v_T^2}) 
\f]
\f[
n_i = \int f dv
\f]
*/
class kinetic {
public:
    /// \brief kinetic constructor
    /*!
    \param &Converter the link on the object which store the information about system parameters
    \param &Poisson the link on the object which solve Poisson equation
     */
    kinetic(const converter &Converter, const poisson &Poisson);

    ///destructor
    ~kinetic();

    ///you can get density using this function
    double ***GetDensity() const;

    /// integrate the hole system
    void IntegrateAll();

    ///return number of sells in velocity space X
    int GetNvx() const;

    ///return number of sells in velocity space Y
    int GetNvy() const;

    ///return number of sells in velocity space Z
    int GetNvz() const;

    ///return velocity integration step X
    double Gethvx() const;

    ///return velocity integration step Y
    double Gethvy() const;

    ///return velocity integration step Z
    double Gethvz() const;

    ///return cut velocity along the flow
    double GetCutVelocityAlong() const;

    ///return cut velocity normal the flow
    double GetCutVelocityNormal() const;

    ///return current time
    double GetCurrentTime() const;

    ///return velocity field x-component
    double ***GetFlowVelocityX() const;

    ///return velocity field y-component
    double ***GetFlowVelocityY() const;

    ///return velocity field z-component
    double ***GetFlowVelocityZ() const;

    ///return profile of the distribution function x-component
    double *GetProfileX() const;

    ///return profile of the distribution function y-component
    double *GetProfileY() const;

    ///return profile of the distribution function z-component
    double *GetProfileZ() const;

    ///return analytical profile of the distribution function z-component
    double *GetAnalyticalProfileX() const;

    ///set profile of the distribution function x-component
    void SetProfileX(int i, int j, int k);

    ///set profile of the distribution function y-component
    void SetProfileY(int i, int j, int k);

    ///set profile of the distribution function z-component
    void SetProfileZ(int i, int j, int k);

    ///return set delta t
    double GetDeltaT() const;

    //get velocity set x-component
    double *GetVelocitySetX() const;

    //get velocity set y-component
    double *GetVelocitySetY() const;

    //get velocity set z-component
    double *GetVelocitySetZ() const;

    ///define density as ((pi/Lx)^2+(pi/Lx)^2+(pi/Lz)^2 * sin(pi*x/Lx)*sin(pi*y/Ly)*sin(pi*z/Lz) / 4pi + n_i/4pi
    void DefineDensity();


private:
    ///distribution function
    double ******f;
    ///temporal distribution function
    double ******f_time;
    ///distribution function of neutrals
    double ***f_n;
    ///velocity field (x component)
    double *vx;
    ///velocity field (y component)
    double *vy;
    ///velocity field (z component)
    double *vz;
    ///self-consistent force field (x component axelaration term)
    double ***asx;
    ///self-consistent force field (y component axelaration term)
    double ***asy;
    ///self-consistent force field (z component axelaration term)
    double ***asz;
    ///Coloumb force field (x component axelaration term)
    double ***acx;
    ///Coloumb force field (y component axelaration term)
    double ***acy;
    ///Coloumb force field (z component axelaration term)
    double ***acz;
    ///longitudinal force field (axelaration term directed along the x axis)
    double al;
    ///coordinate space X
    double *x;
    ///coordinate space Y
    double *y;
    ///coordinate space Z
    double *z;
    ///ion density
    double ***density;
    ///flow velocity numerically calculated x-component
    double ***vfl_x;
    ///flow velocity numerically calculated y-component
    double ***vfl_y;
    ///flow velocity numerically calculated z-component
    double ***vfl_z;
    ///profile of the distribution function
    double *profile_x;
    ///profile of the distribution function
    double *profile_y;
    ///profile of the distribution function
    double *profile_z;
    ///analytical profile of the distribution function
    double *analytical_profile_x;
    ///number of coordinate grids x-direction
    int nx;
    ///number of coordinate grids x-direction
    int ny;
    ///number of coordinate grids x-direction
    int nz;
    ///number of position of the particle X
    int nx_0;
    ///number of position of the particle Y
    int ny_0;
    ///number of position of the particle Z
    int nz_0;
    ///number of velocity grids x-direction
    int nvx;
    ///number of velocity grids y-direction
    int nvy;
    ///number of velocity grids z-direction
    int nvz;
    ///discretisation step in coordinate space X
    double hx;
    ///discretisation step in coordinate space X
    double hy;
    ///discretisation step in coordinate space X
    double hz;
    ///X-length of the computation box
    double Lx;
    ///Y-length of the computation box
    double Ly;
    ///z-length of the computation box
    double Lz;
    ///discretisation step in velocity space X
    double hvx;
    ///discretisation step in velocity space Y
    double hvy;
    ///discretisation step in velocity space Z
    double hvz;
    ///maximum cut velocity in x direction in dimensionless units
    double vxcut;
    ///maximum cut velocity in yz directions in dimensionless units
    double vyzcut;
    ///discretisation step in time space
    double deltaT;
    ///current time
    double CurT;
    /// relaxation time
    double tau;
    /// maximum value of integration patameter in velocity space
    double maxXi;
    ///discrete step for integration parameter
    double deltaXi;
    ///grid size for integration parameter
    int nxi;
    ///analytically calculated velocity of the ion flow
    double v_fl_an;
    ///converter
    const converter *Converter;

    ///implementation of the coordinate part of the integration step
    void CoordinatePart();

    ///implementation of the velocity part of the integration step
    void VelocityPart();

    ///relaxation due to collisions
    double Relaxation(int i, int j, int k, int a, int b, int c);

    ///set initial distribution function
    void InitialConditions();

    ///compute density
    void ComputeDensity();

    ///compute analytical velocity of the ion flow
    void ComputeAnalyticalFlowVelocity();

    ///compute velocity of the ion flow
    void ComputeFlowVelocity();

    ///save previous disctribution function set
    void SaveState();

    ///set analytical profile of the distribution function z-component
    void SetAnalyticalProfileX();


};
