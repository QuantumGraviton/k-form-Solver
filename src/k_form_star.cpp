#include "k_form_star.hpp"

using namespace FBS;

/***********************
 * KFormStar *
 ***********************/

// Event to stop the integration when the field diverges
const integrator::Event KFormStar::Field_diverging = integrator::Event([](const double r, const double dr, const vector &y, const vector &dy, const void *params)
                                                                           { return ((dy[2] > 10.0) || (dy[3] > 10.0) || (dy[4] > 10.0)); }, false);

// initial conditions for two arbitrary fluids
vector KFormStar::get_initial_conditions(const double r_init) const {
    //gtt, grr, phi1, phi2, phi3, Phi1, Phi2, Phi3
    return vector({-1.0, 1.0, this->phi1_0 *pow(r_init,3), this->phi2_0 *pow(r_init,3), this->phi3_0 *pow(r_init,3), 
                    this->phi1_0* 3.*pow(r_init,2), this->phi2_0* 3.*pow(r_init,2), this->phi3_0* 3.*pow(r_init,2)});
}

vector KFormStar::dy_dr(const double r, const vector &vars) const {

	// rename variables for convenience
    const double gtt = vars[0], grr = vars[1], phi1 = vars[2], phi2 = vars[3], phi3 = vars[4], Phi1 = vars[5], Phi2 = vars[6], Phi3 = vars[7];

    // helper variables:
    double th = this->theta;
    double P2 = 0.5* (3.* pow(std::cos(th),2)  - 1. );
    double P2_dtheta = -1.5 * std::sin(2.*th);
    double Ttt = 1.0;
    double Trr = 1.0;

	// compute the ODE:
    // metric components
    double dgtt_dr = 1.;
    double dgrr_dr = 1.;
    // field first derivatives
    double dphi1_dr = Phi1;
    double dphi2_dr = Phi2;
    double dphi3_dr = Phi3;
    // firls second derivatives
    double dPhi1_dr = 1.0;
    double dPhi2_dr = 1.0;
    double dPhi3_dr = 1.0;

	//double da_dr = 0.5* a *      ( (1.-a*a) / r + 8.*M_PI*r*a*a*( etot - 8.*M_PI*s2 ) );
    //double dalpha_dr = 0.5* alpha * ( (a*a-1.) / r + 8.*M_PI*r*a*a*( P - 8.*M_PI*s2 ) );
	//double dP_dr = -(etot + P - 16.*M_PI*s2)/(1. - 8.*M_PI* s2_prime) * dalpha_dr/alpha;

    return vector({dgtt_dr, dgrr_dr, dphi1_dr, dphi2_dr, dphi3_dr, dPhi1_dr, dPhi2_dr, dPhi3_dr});
}

void KFormStar::calculate_star_parameters(const std::vector<integrator::step> &results, const std::vector<integrator::Event> &events) {

    // compute all values of the star, including mass and radius:
    int last_index = results.size() - 1; // last index in the integration
	// we extract the total gravitational mass M_T (asymptotically the solution is just Schwarzschild):
	M_T = results[last_index].first / 2. * (1. - 1./pow(results[last_index].second[0], 2));	// M = R/2 * (1 - 1/ a(R)^2 )
    // ---------------------------------------------------------------------------------
    // update all the global star values:
    this->M_T = M_T;	// total gravitational mass
}

void KFormStar::evaluate_model() {

    std::vector<integrator::step> results;
    this->evaluate_model(results);
}

void KFormStar::evaluate_model(std::vector<integrator::step> &results, std::string filename) {

    // define variables used in the integrator and events during integration:
    integrator::IntegrationOptions intOpts;
    intOpts.save_intermediate = true;
    intOpts.max_stepsize = 1e-3; // smaller stepsize is needed to computethe radius more accurately
    // stop integration if pressure is zero:
    std::vector<integrator::Event> events = {Field_diverging};
    results.clear();

    this->integrate(results, events, this->get_initial_conditions(), intOpts); // integrate the star

    // option to save all the radial profiles into a txt file:
    if (!filename.empty())
    {
		// add two columns for the energy density and restmass density to the results array:
		//for(unsigned i=0; i< results.size(); i++) {
			//results[i].second = vector({results[i].second[0], results[i].second[1], results[i].second[2], this->EOS->get_e_from_P(results[i].second[2]), this->EOS->get_rho_from_P(results[i].second[2])});
		//}
        plotting::save_integration_data(results, {0, 1, 2, 3, 4, 5, 6, 7}, {"gtt", "grr", "phi1", "phi2", "phi3", "Phi1", "Phi2", "Phi3"}, filename);
    }

    this->calculate_star_parameters(results, events);
}

std::ostream& FBS::operator<<(std::ostream &os, const KFormStar &fbs) {

    return os << fbs.M_T << " "					// total gravitational mass in [M_sun]
              << fbs.theta/M_PI << " "			// polar angle/pi
			  << fbs.alpha_Tmunu << " "         // alpha-parameter
              << fbs.beta_Tmunu << " "	        // beta-parameter
              << fbs.phi1_0 << " "				// initial feld component 1
              << fbs.phi2_0 << " "				// initial feld component 2
			  << fbs.phi3_0;                    // initial feld component 3
}

double KFormStar::get_quantity(std::string quantity_label) {
    if (quantity_label == "M_T") { return this->M_T;}
    if (quantity_label == "theta") { return this->theta;}
    if (quantity_label == "alpha_Tmunu") { return this->alpha_Tmunu;}
    if (quantity_label == "beta_Tmunu") { return this->beta_Tmunu;}
    if (quantity_label == "phi1_0") { return this->phi1_0;}
    if (quantity_label == "phi2_0") { return this->phi2_0;}
    if (quantity_label == "phi3_0") { return this->phi3_0;}
    std::cout << "Error: label '" << quantity_label << "' not found. return zero." << std::endl;
    return 0.0;
}

std::vector<std::string> KFormStar::labels() {
	// labels for the Einstein-Cartan NS case:
    return std::vector<std::string>({"M_T", "theta", "alpha_Tmunu", "beta_Tmunu", "phi1_0", "phi2_0", "phi3_0"});
}