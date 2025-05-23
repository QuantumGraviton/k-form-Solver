#include "ns_einstein_cartan_hbar_spin.hpp"

using namespace FBS;

/***********************
 * NSEinsteinCartan *
 ***********************/

// Event to stop the integration when the pressure reaches zero
const integrator::Event NSEinsteinCartanHbarSpin::Pressure_zero = integrator::Event([](const double r, const double dr, const vector &y, const vector &dy, const void *params)
                                                                           { return ((y[2] <= 0.1*P_ns_min) ); }, true);
// Event to stop the integration when the pressure diverges (derivative gets positive)
const integrator::Event NSEinsteinCartanHbarSpin::Pressure_diverging = integrator::Event([](const double r, const double dr, const vector &y, const vector &dy, const void *params)
                                                                           { return ((dy[2] > 110.0) ); }, true);

// initial conditions for two arbitrary fluids
vector NSEinsteinCartanHbarSpin::get_initial_conditions(const double r_init) const {
    //this->prev_step = std::make_pair(r_init, vector({0.0, 0.0, rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(rho_0, 0.) : 0.}));
    return vector({0.0, 0.0, rho_0 > this->EOS->min_rho() ? this->EOS->get_P_from_rho(rho_0, 0.) : 0.});
}

vector NSEinsteinCartanHbarSpin::dy_dr(const double r, const vector &vars, const double r_prev, const vector &vars_prev) const {

	// rename variables for convenience
    const double /*nu = vars[0],*/ lambda = vars[1]; double P = vars[2];
    /*const double prev_nu = vars_prev[0], prev_lambda = vars_prev[1];*/ double P_prev = vars_prev[2];

    // load the EOS for the NS fluid:
    EquationOfState &myEOS = *(this->EOS);

    // call the EOS and compute the wanted values:
    double etot=0.;	// total energy density of the fluid
    double rho=0; // restmass density of the fluid
    double drho_dP=0, dP_drho=0; // derivative restmass w.r.t pressure
    double rho_prev=0; // previous timestep density
    if (P <= 0. || P < myEOS.min_P()) {
        P = 0.;
        etot = 0.;
        rho = 0;
        drho_dP = 0;
        rho_prev=0;
        P_prev=0;
    }
    else {
        etot = myEOS.get_e_from_P(P);
        rho = myEOS.get_rho_from_P(P);
        dP_drho = rho > myEOS.min_rho() ? myEOS.dP_drho(rho,0.) : 0.;
        drho_dP = dP_drho > 0. ? 1./dP_drho : 0.; // valid for barotropic EOS
        // previous timestep:
        rho_prev = myEOS.get_rho_from_P(P_prev);
    }

	// define helper variables:
	double s2 = this->eta_tilde*this->eta_tilde * rho*rho / 2.; // hbar ansatz for the spin fluid: s^2 = 1/2 (eta * hbar/2m_neutron * rho)^2 = 1/2 * eta_tilde^2 * rho^2
    //double s2_prime = this->eta_tilde*this->eta_tilde * rho * drho_dP;
    double dr = r - r_prev;
    double drho = rho - rho_prev;
    if (dr < 1e-14) {dr=1e-14;}
    double s2_dr = this->eta_tilde*this->eta_tilde * rho * (drho/dr); // finite difference derivative
    //std::cout << dr << std::endl;
    //std::cout << "drho: " << drho << std::endl;
    

	// compute the ODE:
    double dnu_dr = (std::exp(lambda) - 1.)/r + 8.*M_PI*r*std::exp(lambda)*( P - 8.*M_PI*s2 );
    double dlambda_dr = (1. - std::exp(lambda))/r + 8.*M_PI*r*std::exp(lambda)*( etot - 8.*M_PI*s2 );
	//double dP_dr = -(etot + P - 16.*M_PI*s2)/(1. - 8.*M_PI* s2_prime) * dnu_dr/2.;
    double dP_dr = -(etot + P - 16.*M_PI*s2) * dnu_dr/2. + 8.*M_PI* s2_dr;
    //std::cout << drho_dP << std::endl;

    //this->prev_step = std::make_pair(r, vector({vars[0], vars[1], vars[2]})); //this->prev_step = std::make_pair(r, vector({vars[0], vars[1], vars[2]}) );
    return vector({dnu_dr, dlambda_dr, dP_dr});
}

void NSEinsteinCartanHbarSpin::calculate_star_parameters(const std::vector<integrator::step> &results, const std::vector<integrator::Event> &events) {

    // compute all values of the star, including mass and radius:
    int last_index = results.size() - 1; // last index in the integration

	// we stop the integration when the pressure reaches zero. This corresponds to the radius of the NS.
	// we extract the total gravitational mass M_T also at this point (outside of the NS is just Schwarzschild):
	//R_NS = results[last_index].first;
    std::cout << "P(R)= " << results[last_index].second[2] << std::endl;
	for (unsigned i = 1; i < results.size(); i++) {
        if (results[i].second[2] < P_ns_min || results[i].second[2] < this->EOS->min_P()) {
            R_NS = results[i].first; // radius of NS
            break;
        }
    }
	M_T = results[last_index].first / 2. * (1. - std::exp(-results[last_index].second[1])); // 1./pow(results[last_index].second[0], 2));	// M = R/2 * (1 - e^-lambda(R))) // (1 - 1/ a(R)^2 )
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
        if (v[2] < P_ns_min || v[2] < this->EOS->min_P()) {rho = 0.;}
        else {this->EOS->callEOS(rho, eps, v[2]);}
        N_F_integrand[i] = 4.*M_PI * std::exp(v[1]/2.) * rho * r[i] * r[i];
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

void NSEinsteinCartanHbarSpin::evaluate_model() {

    std::vector<integrator::step> results;
    this->evaluate_model(results);
}

void NSEinsteinCartanHbarSpin::evaluate_model(std::vector<integrator::step> &results, std::string filename) {

    // define variables used in the integrator and events during integration:
    integrator::IntegrationOptions intOpts;
    intOpts.save_intermediate = true;
    intOpts.max_stepsize = 1e-4; // smaller stepsize is needed to compute the radius more accurately
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
        plotting::save_integration_data(results, {0, 1, 2, 3, 4}, {"nu", "lambda", "P", "e", "rho"}, filename);
    }

    this->calculate_star_parameters(results, events);
}

void NSEinsteinCartanHbarSpin::shooting_constant_Mass(double wanted_mass, std::string quantity_label, double accuracy, int max_steps) {
    // failsafe checks:
    //std::cout << quantity_label << this->get_quantity(quantity_label) << std::endl; std::exit(0);
    std::cout << "Not implemented yet. Quit function" << std::endl; std::exit(0);
    if (!( (quantity_label == "M_T") || (quantity_label == "M_rest") )) { std::cout << "Error: only searches for 'M_T' and 'M_rest' are implemented. quitting function" << std::endl; return;}
    // calc the FBS solution once using an initial rho0 value
    double my_MT = 0.0;
    double rho_0_init = this->rho_0;

    // failsafe check in case that the wanted mass is not attainable with the maximum possible density/pressure, given by 8pi* eta_tilde^2 * rho*drho_dP == 1
    // first compute the star at the highest possible density (rather with density an epsilon smaller to be numerically unproblematic to compute):
    double max_possible_density = 10.;
    if (this->eta_tilde > 0.0) {

        // scan for the critical density:
        max_possible_density = this->EOS->dP_drho(1.,0) / (8.*M_PI * this->eta_tilde*this->eta_tilde) ; // need to do root-finding for this


        this->rho_0 = max_possible_density - 1e-8;
        this->evaluate_model();
        if (this->get_quantity(quantity_label) < wanted_mass) { // in this case, the bisection will automatically not be able to converge
            //std::cout << "failsave activated in NSEinsteinCartan::shooting_constant_Mass" << std::endl;
            this->rho_0 = 0.0;
            this->evaluate_model(); // set to zero
            return;
        }
    }
    

    this->rho_0 = rho_0_init;
    int i = 0;
    while (i<max_steps) {
        i++;
        this->evaluate_model();
        // obtain the current mass
        my_MT = this->get_quantity(quantity_label);
        // check if obtained mass is above the wanted mass
        // if yes, we perform a bisection search in the range [0, rho_0_init]
        // if no, we increase rho_0_init by an amount and perform the above steps again

        if (wanted_mass < my_MT) {
            // the calculated mass is above the wanted mass. We can now perform the bisection search!
            break;
        }
        // the wanted mass is above the calculated mass. Increase the rho_0 for higher mass
        rho_0_init = rho_0_init*1.5;
        if (rho_0_init > max_possible_density) {rho_0_init = max_possible_density- 1e-8;break;}
        this->rho_0 = rho_0_init;
        continue;
    }

    // now perform the bisection until the wanted accuracy is reached:
    // define a few local variables:
    double rho_c_0 = 1e-5;  // this variable might need to be changed for fringe cases
    double rho_c_1 = rho_0_init;
    double rho_c_mid = (rho_c_0 + rho_c_1) / 2.;
    // mass of the lower, mid and upper point in rho0
    double mymass_0;
    double mymass_mid;
    double mymass_1 = my_MT;

    this->rho_0 = rho_c_0;

    this->evaluate_model();
    mymass_0 = this->get_quantity(quantity_label);

    i = 0;
    // continue bisection until the wanted accuracy was reached
    while ( (std::abs(mymass_0 - mymass_1) > accuracy) && (i < max_steps) ) {
        i++;
        rho_c_mid = (rho_c_0 + rho_c_1) / 2.;
        this->rho_0 = rho_c_mid;
        this->evaluate_model();
        // obtain the current mass
        mymass_mid = this->get_quantity(quantity_label);

        if (mymass_mid < wanted_mass) {
            // the mid point is below the wanted ratio and we can adjust the lower bound
            mymass_0 = mymass_mid;
            rho_c_0 = rho_c_mid;
            continue;
        }
        else if (mymass_mid > wanted_mass) {
            // the mid point is above the wanted ratio and we can adjust the upper bound
            mymass_1 = mymass_mid;
            rho_c_1 = rho_c_mid;
            continue;
        }
    }
    // the now obtained rho0 value is now optimized for the wanted gravitational mass and we can quit the function
}

std::ostream& FBS::operator<<(std::ostream &os, const NSEinsteinCartanHbarSpin &fbs) {

    return os << fbs.M_T << " "					// total gravitational mass in [M_sun]
              << fbs.rho_0 << " "				// central density of NS fluid
			  << fbs.R_NS * 1.476625061 << " "  // radius where P(r)=0 in [km]
              << fbs.R_99 * 1.476625061 << " "	// radius (99% matter included) of NS fluid in [km]
              << fbs.M_rest << " "				// total particle number of NS fluid (total restmass) in [M_sun]
              << fbs.C << " "					// Compactness = M_T / R_NS
			  << fbs.eta_tilde;     			// eta tilde for spin fluid model
}

double NSEinsteinCartanHbarSpin::get_quantity(std::string quantity_label) {
    if (quantity_label == "M_T") { return this->M_T;}
    if (quantity_label == "rho_0") { return this->rho_0;}
    if (quantity_label == "R_NS") { return this->R_NS;}
    if (quantity_label == "R_99") { return this->R_99;}
    if (quantity_label == "M_rest") { return this->M_rest;}
    if (quantity_label == "C") { return this->C;}
    if (quantity_label == "eta_tilde") { return this->eta_tilde;}
    std::cout << "Error: label '" << quantity_label << "' not found. return zero." << std::endl;
    return 0.0;
}

std::vector<std::string> NSEinsteinCartanHbarSpin::labels() {
	// labels for the Einstein-Cartan NS case:
    return std::vector<std::string>({"M_T", "rho_0", "R_NS", "R_99", "M_rest", "C", "eta_tilde"});
}