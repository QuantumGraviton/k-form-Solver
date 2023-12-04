#include "ns_einstein_cartan_rotation.hpp"

using namespace FBS;

/****************************
 * NSEinsteinCartanRotation *
 ***************************/

// Event to stop the integration when the pressure reaches zero
const integrator::Event NSEinsteinCartanRotation::Pressure_zero = integrator::Event([](const double r, const double dr, const vector &y, const vector &dy, const void *params)
                                                                           { return ((y[2] <= 0.1*P_ns_min) ); }, true);
// Event to stop the integration when the pressure diverges (derivative gets positive)
const integrator::Event NSEinsteinCartanRotation::Pressure_diverging = integrator::Event([](const double r, const double dr, const vector &y, const vector &dy, const void *params)
                                                                           { return ((dy[2] > 0.0) ); }, true);


vector NSEinsteinCartanRotation::dy_dr(const double r, const vector &vars) const {

	// rename variables for convenience
    const double a = vars[0], alpha = vars[1]; double P = vars[2];

    // load the EOS for the NS fluid:
    EquationOfState &myEOS = *(this->EOS);

    // call the EOS and compute the wanted values:
    double etot=0., dP_de =0., de_dP =0.;	// total energy density of the fluid
    if (P <= 0. || P < myEOS.min_P()) {
        P = 0.;
        etot = 0.;
		dP_de = 0.;
    }
    else {
        etot = myEOS.get_e_from_P(P);
		dP_de = etot > myEOS.min_e() ? myEOS.dP_de(etot) : 0.;
        de_dP = dP_de > 0. ? 1./dP_de : 0.;
    }

	// define helper variables:
	double Gamma = 1./ std::sqrt(1. - r*r* this->Omega_rot*this->Omega_rot);
	double s2 = std::pow(r,4.) * this->Omega_rot*this->Omega_rot * std::pow(Gamma,4.) * std::pow( etot + P,2.); // rotation ansatz the spin fluid: s^2 = ...
	double division_term = 1. - 16.*M_PI*std::pow(r,4.) * this->Omega_rot*this->Omega_rot * std::pow(Gamma,4.) * (etot + P) *(1. + de_dP); // derivative of s^2 with respect to P

	// compute the ODE:
	double da_dr = 0.5* a *      ( (1.-a*a) / r + 8.*M_PI*r*a*a*( etot - 8.*M_PI*s2 ) );
    double dalpha_dr = 0.5* alpha * ( (a*a-1.) / r + 8.*M_PI*r*a*a*( P - 8.*M_PI*s2 ) );
	double dP_dr = ( ( 32.*M_PI*s2 *(1./r + r*this->Omega_rot*this->Omega_rot* Gamma*Gamma) ) - (etot + P - 16.*M_PI*s2) * dalpha_dr/alpha ) / division_term;

    return vector({da_dr, dalpha_dr, dP_dr});
}


void NSEinsteinCartanRotation::evaluate_model(std::vector<integrator::step> &results, std::string filename) {

    // define variables used in the integrator and events during integration:
    integrator::IntegrationOptions intOpts;
    intOpts.save_intermediate = true;
    // stop integration if pressure is zero:
    std::vector<integrator::Event> events = {Pressure_zero, Pressure_diverging};
    results.clear();

    this->integrate(results, events, this->get_initial_conditions(), intOpts); // integrate the star

    // option to save all the radial profiles into a txt file:
    if (!filename.empty())
    {
		// add two columns for the energy density and restmass density to the results array:
		for(unsigned i=0; i< results.size(); i++) {
			results[i].second = vector({results[i].second[0], results[i].second[1], results[i].second[2], this->EOS->get_e_from_P(results[i].second[2]), this->EOS->get_rho_from_P(results[i].second[2])});
		}
        plotting::save_integration_data(results, {0, 1, 2, 3, 4}, {"a", "alpha", "P", "e", "rho"}, filename);
    }

    NSEinsteinCartan::calculate_star_parameters(results, events);
}

std::ostream& FBS::operator<<(std::ostream &os, const NSEinsteinCartanRotation &fbs) {

    return os << fbs.M_T << " "					// total gravitational mass in [M_sun]
              << fbs.rho_0 << " "				// central density of NS fluid
			  << fbs.R_NS * 1.476625061 << " "  // radius where P(r)=0 in [km]
              << fbs.R_99 * 1.476625061 << " "	// radius (99% matter included) of NS fluid in [km]
              << fbs.M_rest << " "				// total particle number of NS fluid (total restmass) in [M_sun]
              << fbs.C << " "					// Compactness = M_T / R_NS
			  << fbs.Omega_rot					// anguar rotation velocity
			  ;
}

std::vector<std::string> NSEinsteinCartanRotation::labels() {
	// labels for the Einstein-Cartan NS case:
    return std::vector<std::string>({"M_T", "rho_0", "R_NS", "R_99", "M_rest", "C", "Omega_rot"});
}