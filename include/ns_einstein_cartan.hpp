#pragma once

#include <utility>  // for std::swap

#include "vector.hpp"
#include "eos.hpp"
#include "integrator.hpp"
#include "nsmodel.hpp"

namespace FBS {

// this is a class modeling a neutron star in Eistein Cartan gravity
// constructor: EOS (ptr), EOS2 (ptr)
class NSEinsteinCartan : public NSmodel {
protected:
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);

public:

	double beta, gamma;	// parameter in the "polytropic" ansatz for the spin density s^2 = beta*P^gamma. Corresponds to model b) in our paper
    double rho_0;	// initial condition, central density of the NS fluid
    double M_T, R_NS, R_99; // total mass M_T; radius of neutron star where pressure is zero R_NS; radius where 99% of restmass is included R_99 
	double M_rest, C; // total restmass of the NS fluid; compactness C:=M/R


	NSEinsteinCartan(std::shared_ptr<EquationOfState> EOS, double rho_0_in, double beta_in, double gamma_in)
        : NSmodel(EOS), beta(beta_in), gamma(gamma_in), rho_0(rho_0_in), M_T(0.), R_NS(0.), R_99(0.), M_rest(0.), C(0.) {}

    //vector dy_dr(const double r, const vector& vars);  // holds the system of ODEs
	/* The differential equations describing the neutron star in Einstein Cartan gravity. The quantities are a, alpha, P */
    vector dy_dr(const double r, const vector& vars) const;

    vector get_initial_conditions(const double r_init=R_INIT) const; // holds the FBS init conditions
    void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();

    // optimizes the central density to find a star with a specific mass
    void shooting_constant_Mass(double wanted_mass, std::string quantity_label, double accuracy=1e-6, int max_steps=200);

    friend std::ostream& operator<<(std::ostream&, const NSEinsteinCartan&);
    static std::vector<std::string> labels();
    double get_quantity(std::string quantity_label);

    static const integrator::Event Pressure_zero;
    static const integrator::Event Pressure_diverging;
};

std::ostream& operator<<(std::ostream&, const NSEinsteinCartan&);

}
