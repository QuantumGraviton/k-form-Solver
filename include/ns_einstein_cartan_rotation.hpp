#pragma once

#include <utility>  // for std::swap

#include "vector.hpp"
#include "eos.hpp"
#include "integrator.hpp"
#include "nsmodel.hpp"
#include "ns_einstein_cartan.hpp"

namespace FBS {

// this is a class modeling a neutron star in Eistein Cartan gravity
// constructor: EOS (ptr), EOS2 (ptr)
class NSEinsteinCartanRotation : public NSEinsteinCartan {
protected:
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

public:

	double rho_0;	// initial condition, central density of the NS fluid
    double M_T, R_NS, R_99; // total mass M_T; radius of neutron star where pressure is zero R_NS; radius where 99% of restmass is included R_99 
	double M_rest, C; // total restmass of the NS fluid; compactness C:=M/R


	NSEinsteinCartanRotation(std::shared_ptr<EquationOfState> EOS, double rho_0_in, double beta_eff_in)
        : NSEinsteinCartan(EOS,rho_0_in,beta_eff_in,0.0), rho_0(rho_0_in), M_T(0.), R_NS(0.), R_99(0.), M_rest(0.), C(0.) {}

    //vector dy_dr(const double r, const vector& vars);  // holds the system of ODEs
	/* The differential equations describing the neutron star in Einstein Cartan gravity. The quantities are a, alpha, P */
    vector dy_dr(const double r, const vector& vars) const;
    vector get_initial_conditions(const double r_init=R_INIT) const; // holds the FBS init conditions

	void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();

    friend std::ostream& operator<<(std::ostream&, const NSEinsteinCartanRotation&);
    static std::vector<std::string> labels();

    static const integrator::Event Pressure_zero;
    static const integrator::Event Pressure_diverging;
};

std::ostream& operator<<(std::ostream&, const NSEinsteinCartanRotation&);

}
