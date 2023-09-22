#include "ns_einstein_cartan.hpp"

using namespace FBS;

/***********************
 * NSEinsteinCartan *
 ***********************/

// Event to stop the integration when the pressure reaches zero
const integrator::Event NSEinsteinCartan::Pressure_zero = integrator::Event([](const double r, const double dr, const vector &y, const vector &dy, const void *params)
                                                                           { return ((y[2] <= 0.1*P_ns_min) ); },
                                                                           true);

// initial conditions for two arbitrary fluids
vector NSEinsteinCartan::get_initial_conditions(const double r_init) const {

    return vector({1.0, 1.0, rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(rho_0, 0.) : 0.});
}

vector NSEinsteinCartan::dy_dr(const double r, const vector &vars) const {

	// rename variables for convenience
    const double a = vars[0], alpha = vars[1]; double P = vars[2];

    // load the EOS for the NS fluid:
    EquationOfState &myEOS = *(this->EOS);

	// define helper variables:
	double s2 = this->beta * std::pow(P, this->gamma); // polytropic ansatz for the spin fluid: s^2 = beta * P^gamma
	double s2_prime = this->beta*this->gamma*std::pow(P, this->gamma - 1.); // derivative of s^2 with respect to P

    // call the EOS and compute the wanted values:
    double etot=0.;	// total energy density of the fluid
    if (P <= 0. || P < myEOS.min_P() || (P-8.*M_PI*s2) <= 0.) {
        P = 0.;
        etot = 0.;
		s2 = 0.; s2_prime=0.;
    }
    else {
        etot = myEOS.get_e_from_P(P);
		// helper variables
		s2 = this->beta * std::pow(P, this->gamma);
		s2_prime = this->beta*this->gamma*std::pow(P, this->gamma - 1.);
    }

	// compute the ODE:
	double da_dr = 0.5* a *      ( (1.-a*a) / r + 8.*M_PI*r*a*a*( etot - 8.*M_PI*s2 ) );
    double dalpha_dr = 0.5* alpha * ( (a*a-1.) / r + 8.*M_PI*r*a*a*( P - 8.*M_PI*s2 ) );
	double dP_dr = -(etot + P - 16.*M_PI*s2)/(1. - 8.*M_PI* s2_prime) * dalpha_dr/alpha;

    return vector({da_dr, dalpha_dr, dP_dr});
}

void NSEinsteinCartan::calculate_star_parameters(const std::vector<integrator::step> &results, const std::vector<integrator::Event> &events) {

    // compute all values of the star, including mass and radius:
    int last_index = results.size() - 1; // last index in the integration

	// we stop the integration when the pressure reaches zero. This corresponds to the radius of the NS.
	// we extract the total gravitational mass M_T also at this point (outside of the NS is just Schwarzschild):
	//R_NS = results[last_index].first;
	for (unsigned i = 1; i < results.size(); i++) {
        if (results[i].second[2] < P_ns_min || results[i].second[2] < this->EOS->min_P() || results[i].second[2] < 1.01*results[last_index].second[2]) {
            R_NS = results[i].first; // radius uf NS
            break;
        }
    }
	M_T = results[last_index].first / 2. * (1. - 1./pow(results[last_index].second[0], 2));	// M = R/2 * (1 - 1/ a(R)^2 )
	double C = M_T / R_NS; // compactness of the NS

	// compute the total restmass of the NS using the conserved Noether current (conservation of fluid flow: J^mu = rho u^mu)
    std::vector<double> r(results.size()), N_F_integrand(results.size());
    vector v;
    double rho, eps;
	// compute restmass at each r:
    for(unsigned int i = 0; i < results.size(); i++) {
        r[i] = results[i].first;
        v = results[i].second;
        // P = v[2]
        if (v[2] < P_ns_min || v[2] < this->EOS->min_P() || v[2] < 1.01*results[last_index].second[2]) {rho = 0.;}
        else {this->EOS->callEOS(rho, eps, v[2]);}
        N_F_integrand[i] = 4.*M_PI * v[0] * rho * r[i] * r[i];
	}
    // integrate:
    std::vector<double> N_F_integrated;
    integrator::cumtrapz(r, N_F_integrand, N_F_integrated);
    // compute the restmass density:
    double N_F =  N_F_integrated[last_index];

    // compute radius where 99% of the restmass is contained:
	// find index where this happens
	int i_F = 0;
    for(int i = 1; i < last_index; i++) {
        if(N_F_integrated[i] < 0.99*N_F) {i_F++;}       
    }
    // obtain radius from corresponding index:
    double R_99 = results[i_F].first;

    // ---------------------------------------------------------------------------------

    // update all the global star values:
    this->M_T = M_T;	// total gravitational mass
    this->R_99 = R_99;	// radius where 99% of restmass is contained (also calles 'tidal radius')
    this->R_NS = R_NS;	// radius of NS
    this->C = C;		// compactness
    this->M_rest = N_F;	// total restmass
}

void NSEinsteinCartan::evaluate_model() {

    std::vector<integrator::step> results;
    this->evaluate_model(results);
}

void NSEinsteinCartan::evaluate_model(std::vector<integrator::step> &results, std::string filename) {

    // define variables used in the integrator and events during integration:
    integrator::IntegrationOptions intOpts;
    intOpts.save_intermediate = true;
    // stop integration if pressure is zero:
    std::vector<integrator::Event> events = {NSEinsteinCartan::Pressure_zero};
    results.clear();

    this->integrate(results, events, this->get_initial_conditions(), intOpts); // integrate the star

    // option to save all the radial profiles into a txt file:
    if (!filename.empty())
    {
        plotting::save_integration_data(results, {0, 1, 2}, {"a", "alpha", "P"}, filename);
    }

    this->calculate_star_parameters(results, events);
}

std::ostream& FBS::operator<<(std::ostream &os, const NSEinsteinCartan &fbs) {

    return os << fbs.M_T << " "					// total gravitational mass in [M_sun]
              << fbs.rho_0 << " "				// central density of NS fluid
			  << fbs.R_NS * 1.476625061 << " "  // radius where P(r)=0 in [km]
              << fbs.R_99 * 1.476625061 << " "	// radius (99% matter included) of NS fluid in [km]
              << fbs.M_rest << " "				// total particle number of NS fluid (total restmass) in [M_sun]
              << fbs.C << " "					// Compactness = M_T / R_NS
			  << fbs.beta << " "				// beta-parameter for spin fluid model
			  << fbs.gamma 						// gamma-parameter for spin fluid model
			  ;
}

std::vector<std::string> NSEinsteinCartan::labels() {
	// labels for the Einstein-Cartan NS case:
    return std::vector<std::string>({"M_T", "rho_0", "R_NS", "R_99", "M_rest", "C", "beta", "gamma"});
}