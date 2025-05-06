#pragma once

#include <utility>  // for std::swap

#include "vector.hpp"
//#include "eos.hpp"
#include "integrator.hpp"
#include "nsmodel.hpp"

namespace FBS {

// this is a class modeling a neutron star in Eistein Cartan gravity
// constructor: EOS (ptr), EOS2 (ptr)
class KFormStar : public NSmodel {
protected:
    void calculate_star_parameters(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events);
    
    // placeholder EOS so that the cod edoed not break:
    //std::shared_ptr<EquationOfState> EOS; // TODO: FIX THIS

public:
    
    double M_T; // total ADM mass M_T

    double theta = 1.5; // polar angle in radians
    double alpha_Tmunu = 1.0; // alpha-parameter in energy momentum tensor
    double beta_Tmunu = 1.0; // beta-parameter in energy momentum tensor
    // initial field values:
    double phi1_0 = 0.01;
    double phi2_0 = 0.01;
    double phi3_0 = 0.01;


	KFormStar(double theta_in, double alpha_Tmunu_in, double beta_Tmunu_in, double phi1_0_in, double phi2_0_in, double phi3_0_in)
        : NSmodel(), theta(theta_in), alpha_Tmunu(alpha_Tmunu_in), beta_Tmunu(beta_Tmunu_in), phi1_0(phi1_0_in), phi2_0(phi2_0_in), phi3_0(phi3_0_in), M_T(0.) {}

    //vector dy_dr(const double r, const vector& vars);  // holds the system of ODEs
	/* The differential equations describing the neutron star in Einstein Cartan gravity. The quantities are a, alpha, P */
    vector dy_dr(const double r, const vector& vars) const;

    vector get_initial_conditions(const double r_init=R_INIT) const; // holds the FBS init conditions
    void evaluate_model(std::vector<integrator::step>& results, std::string filename="");
    void evaluate_model();

    // optimizes the central density to find a star with a specific mass
    void shooting_constant_Mass(double wanted_mass, std::string quantity_label, double accuracy=1e-6, int max_steps=200);

    friend std::ostream& operator<<(std::ostream&, const KFormStar&);
    static std::vector<std::string> labels();
    double get_quantity(std::string quantity_label);

    static const integrator::Event Field_diverging;
};

std::ostream& operator<<(std::ostream&, const KFormStar&);

}
