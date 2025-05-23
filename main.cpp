#include <iostream> // for i/o e.g. std::cout etc.
#include <cmath>	// for mathematical functions
#include <vector>	// for std::vector
//#include <iomanip> 	// for std::fixed and std::fixesprecission()
#include <memory>   // shared_ptr
#include <sstream> // string stream

#include "vector.hpp"    // include custom 5-vector class
#include "integrator.hpp"
#include "eos.hpp" // include eos container class
#include "nsmodel.hpp"
#include "fbs_twofluid.hpp"
#include "mr_curves.hpp"
#include "plotting.hpp"     // to use python/matplotlib inside of c++
#include "utilities.hpp"
#include "ns_einstein_cartan.hpp"
#include "ns_einstein_cartan_rotation.hpp"
#include "ns_einstein_cartan_hbar_spin.hpp"
#include "k_form_star.hpp"

// --------------------------------------------------------------------
using namespace FBS;

void k_form_star_single() {
	// define parameters for the 'spin densits EOS': s^2 = beta*P^gamma
	double alpha_Tmunu = 1000.0; //beta_in;
	double theta = 0.6 * M_PI;
	double phi1_0 = -1e10;
    double phi2_0 = 0.1e10;
    double phi3_0 = 1e10;


	// initialize one instance:
	KFormStar KFstar(theta, alpha_Tmunu, 0.0, phi1_0, phi2_0, phi3_0);

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::stringstream stream; std::string tmp;
	std::string plotname = "K_Form_star_profile_alphTmn_"; stream << std::fixed << std::setprecision(5) << alpha_Tmunu;
	stream >> tmp; plotname += (tmp + "_theta_");  stream = std::stringstream(); stream << std::fixed << std::setprecision(5) << theta;
	stream >> tmp; plotname += (tmp + "_initphi_");  stream = std::stringstream(); stream << std::fixed << std::setprecision(5) << phi1_0;
	stream >> tmp; plotname += (tmp + "_"); stream = std::stringstream(); stream << std::fixed << std::setprecision(5) << phi2_0;
	stream >> tmp; plotname += (tmp + "_"); stream = std::stringstream(); stream << std::fixed << std::setprecision(5) << phi3_0;
	stream >> tmp; plotname += tmp;

	// evaluate the model and save the intermediate data into txt file:
	std::vector<integrator::step> results;
    KFstar.evaluate_model(results, "output/" + plotname + ".txt");

	// output global variables:
	std::cout << "calculation complete!" << std::endl;
	std::cout << "global quantities:" << std::endl;
	std::vector<std::string> labels = KFstar.labels();
	std::cout << labels[0] << " " << labels[1] << " " << labels[2] << " " << labels[3] << " " << labels[4] << " " << labels[5] << " " << labels[6] << std::endl;
	std::cout << std::fixed << std::setprecision(10) << KFstar << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}


// computes a single neutron star in Einstein-Cartan gravity and saves the radial profile into a file
// double rho_0 [in saturation density], double beta, double gamma, string EOS_name
void EC_star_single(double rho_0_in, double beta_in, double gamma_in, std::string EOS_in) {

	// define parameters for the 'spin densits EOS': s^2 = beta*P^gamma
	double beta = beta_in;
	double gamma = gamma_in; // setting gamma=0 is somewhat broken, the pressure will not converge to zero...
	// initial density:
	const double sat_to_code = 0.16 * 2.886376934e-6 * 939.565379;	// conversion factor from nuclear saturation density to code units
	double rho_0 = rho_0_in * sat_to_code; // central restmass density of the NS in units of the nucelar saturation density

	// select correct EOS type:
	std::string EOS_filepath = "";
	     if(EOS_in == "EOS_DD2")    { EOS_filepath = "EOS_tables/eos_HS_DD2_with_electrons.beta";}
	else if(EOS_in == "EOS_APR")    { EOS_filepath = "EOS_tables/eos_SRO_APR_SNA_version.beta";}
	else if(EOS_in == "EOS_KDE0v1") { EOS_filepath = "EOS_tables/eos_SRO_KDE0v1_SNA_version.beta";}
	else if(EOS_in == "EOS_LNS")    { EOS_filepath = "EOS_tables/eos_SRO_LNS_SNA_version.beta";}
	else if(EOS_in == "EOS_FSG")    { EOS_filepath = "EOS_tables/eos_HS_FSG_with_electrons.beta";}
	else { std::cout << "Wrong EOS! Supported EOS are: 'EOS_DD2'  'EOS_APR'  'EOS_KDE0v1'  'EOS_LNS'  'EOS_FSG' !" << std::endl; return;}
	auto myEOS = std::make_shared<EoStable>(EOS_filepath); // (load the EOS)

	// initialize one instance:
	NSEinsteinCartan ECstar(myEOS, rho_0, beta, gamma);

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::stringstream stream; std::string tmp;
	std::string plotname = "ECstar_profile_" + EOS_in + "_rho0_"; stream << std::fixed << std::setprecision(10) << rho_0_in; 
	stream >> tmp; plotname += (tmp + "_beta_");  stream = std::stringstream(); stream << std::fixed << std::setprecision(10) << beta;
	stream >> tmp; plotname += (tmp + "_gamma_"); stream = std::stringstream(); stream << std::fixed << std::setprecision(10) << gamma;
	stream >> tmp; plotname += tmp;

	// evaluate the model and save the intermediate data into txt file:
	std::vector<integrator::step> results;
    ECstar.evaluate_model(results, "output/" + plotname + ".txt");

	// output global variables:
	std::cout << "calculation complete!" << std::endl;
	std::cout << "global quantities:" << std::endl;
	std::vector<std::string> labels = ECstar.labels();
	std::cout << labels[0] << " " << labels[1] << " " << labels[2] << " " << labels[3] << " " << labels[4] << " " << labels[5] << " " << labels[6] << " " << labels[7] << " " << std::endl;
	std::cout << std::fixed << std::setprecision(10) << ECstar << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}

// computes multiple neutron stars in Einstein-Cartan gravity and saves the mass and radii into a file
// unsigned Nstars, double rho0_min [in saturation density], double rho0_max [in saturation density], double beta, double gamma, string EOS_name:
void EC_star_curve(unsigned Nstars_in, double rho_0_min_in, double rho_0_max_in, double beta_in, double gamma_in, std::string EOS_in) {
	
	unsigned Nstars = Nstars_in; // number of stars in MR curve
	std::vector<double> rho_c_grid(Nstars, 0.0);

	// define parameters for the 'spin densits EOS': s^2 = beta*P^gamma
	double beta = beta_in;
	double gamma = gamma_in; // setting gamma=0 is somewhat broken, the pressure will not converge to zero...
	// lowest and highest initial densities for the NSs in the MR curve:
	const double sat_to_code = 0.16 * 2.886376934e-6 * 939.565379;	// conversion factor from nuclear saturation density to code units
	double rho_0_min = rho_0_min_in * sat_to_code; // central restmass density of the NSs in units of the nucelar saturation density
	double rho_0_max = rho_0_max_in * sat_to_code;

	// fill array with the initial conditions for every star:
	utilities::fillValuesPowerLaw(rho_0_min, rho_0_max, rho_c_grid, 2);	// power law scaling of 2 or 3 works pretty well

	// select correct EOS type:
	std::string EOS_filepath = "";
	     if(EOS_in == "EOS_DD2")    { EOS_filepath = "EOS_tables/eos_HS_DD2_with_electrons.beta";}
	else if(EOS_in == "EOS_APR")    { EOS_filepath = "EOS_tables/eos_SRO_APR_SNA_version.beta";}
	else if(EOS_in == "EOS_KDE0v1") { EOS_filepath = "EOS_tables/eos_SRO_KDE0v1_SNA_version.beta";}
	else if(EOS_in == "EOS_LNS")    { EOS_filepath = "EOS_tables/eos_SRO_LNS_SNA_version.beta";}
	else if(EOS_in == "EOS_FSG")    { EOS_filepath = "EOS_tables/eos_HS_FSG_with_electrons.beta";}
	else { std::cout << "Wrong EOS! Supported EOS are: 'EOS_DD2'  'EOS_APR'  'EOS_KDE0v1'  'EOS_LNS'  'EOS_FSG' !" << std::endl; return;}
	auto myEOS = std::make_shared<EoStable>(EOS_filepath); // (load the EOS)

	// compute all EC neutron stars:
	std::vector<NSEinsteinCartan> MR_curve;	// holds the stars in the MR curve
	calc_EinsteinCartan_curves(myEOS, rho_c_grid, MR_curve, beta, gamma, 1);	// last argument is verbose

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::stringstream stream; std::string tmp;
	std::string plotname = "ECstar_curve_" + EOS_in + "_beta_"; stream << std::fixed << std::setprecision(10) << beta; 
	stream >> tmp; plotname += (tmp + "_gamma_"); stream = std::stringstream(); stream << std::fixed << std::setprecision(10) << gamma;
	stream >> tmp; plotname += tmp;

	write_MRphi_curve<NSEinsteinCartan>(MR_curve, "output/" + plotname + ".txt");
	std::cout << "calculation complete!" << std::endl;
	std::cout << "parameters: beta gamma" << std::endl;
	std::cout << std::fixed << std::setprecision(10) << beta << " " << gamma << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}

// computes multiple neutron stars in Einstein-Cartan gravity and saves the mass and radii into a file
// unsigned Nstars, double rho0_min [in saturation density], double rho0_max [in saturation density], double beta, double gamma, string EOS_name:
void EC_star_curve_const_mass_with_different_beta(unsigned Nstars_in, double ns_mass_in, std::string mass_quantity_label, double beta_min_in, double beta_max_in, double gamma_in, std::string EOS_in) {
	
	unsigned Nstars = Nstars_in; // number of stars in MR curve
	std::vector<double> beta_grid(Nstars, 0.0);

	// define parameters for the 'spin densits EOS': s^2 = beta*P^gamma
	double gamma = gamma_in; // setting gamma=0 is somewhat broken, the pressure will not converge to zero...
	double beta_min = beta_min_in; // central restmass density of the NSs in units of the nucelar saturation density
	double beta_max = beta_max_in;
	// initial densities for the NSs in the MR curve:
	const double sat_to_code = 0.16 * 2.886376934e-6 * 939.565379;	// conversion factor from nuclear saturation density to code units
	double ns_mass = ns_mass_in;

	// fill array with the initial conditions for every star:
	utilities::fillValuesPowerLaw(beta_min, beta_max, beta_grid, 1); // linear scaling, but other scalings are possible

	// select correct EOS type:
	std::string EOS_filepath = "";
	     if(EOS_in == "EOS_DD2")    { EOS_filepath = "EOS_tables/eos_HS_DD2_with_electrons.beta";}
	else if(EOS_in == "EOS_APR")    { EOS_filepath = "EOS_tables/eos_SRO_APR_SNA_version.beta";}
	else if(EOS_in == "EOS_KDE0v1") { EOS_filepath = "EOS_tables/eos_SRO_KDE0v1_SNA_version.beta";}
	else if(EOS_in == "EOS_LNS")    { EOS_filepath = "EOS_tables/eos_SRO_LNS_SNA_version.beta";}
	else if(EOS_in == "EOS_FSG")    { EOS_filepath = "EOS_tables/eos_HS_FSG_with_electrons.beta";}
	else { std::cout << "Wrong EOS! Supported EOS are: 'EOS_DD2'  'EOS_APR'  'EOS_KDE0v1'  'EOS_LNS'  'EOS_FSG' !" << std::endl; return;}
	auto myEOS = std::make_shared<EoStable>(EOS_filepath); // (load the EOS)

	// compute all EC neutron stars:
	std::vector<NSEinsteinCartan> MR_curve;	// holds the stars in the MR curve
	calc_EinsteinCartan_curves_const_mass(myEOS, 1.0*sat_to_code, MR_curve, beta_grid, gamma, ns_mass, mass_quantity_label, 1);	// last argument is verbose

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::stringstream stream; std::string tmp;
	std::string plotname = "ECstar_curve_const_mass_" + mass_quantity_label + "_with_different_beta_" + EOS_in + "_nsmass_"; stream << std::fixed << std::setprecision(10) << ns_mass; 
	stream >> tmp; plotname += (tmp + "_gamma_"); stream = std::stringstream(); stream << std::fixed << std::setprecision(10) << gamma;
	stream >> tmp; plotname += tmp;

	write_MRphi_curve<NSEinsteinCartan>(MR_curve, "output/" + plotname + ".txt");
	std::cout << "calculation complete!" << std::endl;
	std::cout << "parameters: NSmass gamma" << std::endl;
	std::cout << std::fixed << std::setprecision(10) << ns_mass << " " << gamma << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}

// computes multiple neutron stars in Einstein-Cartan gravity and saves the mass and radii into a file
// used the rotation model as a stand-in to compute the beta value to use with the polytropic model s^2 = beta*P^gamma
// unsigned Nstars, double rho0_min [in saturation density], double rho0_max [in saturation density], double beta, double gamma, string EOS_name:
void EC_star_curve_rotation_rate_model_simple(unsigned Nstars_in, double rho_0_min_in, double rho_0_max_in, double frequency_Hertz_in, std::string EOS_in, double Keplarian = 0.0) {
	
	unsigned Nstars = Nstars_in; // number of stars in MR curve
	std::vector<double> rho_c_grid(Nstars, 0.0);

	// define parameters for the 'spin densits EOS': s^2 = beta*P^gamma
	double Omega_Hertz_to_code_units = 203000.0; // unit conversion: amount of Hz that correspond to Omega=1
	double Omega_rot = 2.*M_PI*frequency_Hertz_in / Omega_Hertz_to_code_units; // angular frequency in code units
	std::vector<double> beta_grid(Nstars, 0.0); // later set this equal to R^4 Omega_rot^2 \Gamma^4, \Gamma = 1/ sqrt(1 - R^2 \Omega^2)

	// lowest and highest initial densities for the NSs in the MR curve:
	const double sat_to_code = 0.16 * 2.886376934e-6 * 939.565379;	// conversion factor from nuclear saturation density to code units
	double rho_0_min = rho_0_min_in * sat_to_code; // central restmass density of the NSs in units of the nucelar saturation density
	double rho_0_max = rho_0_max_in * sat_to_code;

	// fill array with the initial conditions for every star:
	utilities::fillValuesPowerLaw(rho_0_min, rho_0_max, rho_c_grid, 2);	// power law scaling of 2 or 3 works pretty well

	// select correct EOS type:
	std::string EOS_filepath = "";
	     if(EOS_in == "EOS_DD2")    { EOS_filepath = "EOS_tables/eos_HS_DD2_with_electrons.beta";}
	else if(EOS_in == "EOS_APR")    { EOS_filepath = "EOS_tables/eos_SRO_APR_SNA_version.beta";}
	else if(EOS_in == "EOS_KDE0v1") { EOS_filepath = "EOS_tables/eos_SRO_KDE0v1_SNA_version.beta";}
	else if(EOS_in == "EOS_LNS")    { EOS_filepath = "EOS_tables/eos_SRO_LNS_SNA_version.beta";}
	else if(EOS_in == "EOS_FSG")    { EOS_filepath = "EOS_tables/eos_HS_FSG_with_electrons.beta";}
	else { std::cout << "Wrong EOS! Supported EOS are: 'EOS_DD2'  'EOS_APR'  'EOS_KDE0v1'  'EOS_LNS'  'EOS_FSG' !" << std::endl; return;}
	auto myEOS = std::make_shared<EoStable>(EOS_filepath); // (load the EOS)

	std::vector<NSEinsteinCartan> MR_curve;	// holds the stars in the MR curve
	// compute non-rotating neutron stars without torsion:
	calc_EinsteinCartan_curves_beta_grid(myEOS, rho_c_grid, MR_curve, beta_grid, 2.0, 1);	// last argument is verbose
	// use the star results to compute a corresponding beta value:
	for(unsigned int i = 0; i < MR_curve.size(); i++) {
		double Rns = MR_curve[i].R_NS;
		double Mns = MR_curve[i].M_T;
		if(Keplarian > 0.0) { // initialize NS with Keplarian rotation rate
			if(Rns > 0.) {Omega_rot = Keplarian*std::sqrt( Mns/(Rns*Rns*Rns) );} else {Omega_rot = 0.0;}
		} 
        beta_grid[i] = pow(Rns,4) * Omega_rot*Omega_rot / pow(1.- Rns*Rns * Omega_rot*Omega_rot,2);
    }
	calc_EinsteinCartan_curves_beta_grid(myEOS, rho_c_grid, MR_curve, beta_grid, 2.0, 1);	// last argument is verbose

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::string plotname; std::stringstream stream; std::string tmp;
	if(Keplarian > 0.0) {
		plotname = "ECstar_curve_rotation_rate_model_simple_" + EOS_in + "_frequencyKep_"; stream << std::fixed << std::setprecision(10) << Keplarian;
		stream >> tmp; plotname += tmp;
	} else{
		plotname = "ECstar_curve_rotation_rate_model_simple_" + EOS_in + "_frequencyHz_"; stream << std::fixed << std::setprecision(10) << frequency_Hertz_in;
		stream >> tmp; plotname += tmp;
	}

	write_MRphi_curve<NSEinsteinCartan>(MR_curve, "output/" + plotname + ".txt");
	std::cout << "calculation complete!" << std::endl;
	std::cout << "parameters: angular frequency (Hz)" << std::endl;
	std::cout << std::fixed << std::setprecision(10) << frequency_Hertz_in << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}

void EC_star_curve_rotation_rate_model_e_plus_P(unsigned Nstars_in, double rho_0_min_in, double rho_0_max_in, double frequency_Hertz_in, std::string EOS_in, double Keplarian = 0.0) {
	
	unsigned Nstars = Nstars_in; // number of stars in MR curve
	std::vector<double> rho_c_grid(Nstars, 0.0);

	// define parameters for the 'spin densits EOS': s^2 = beta*P^gamma
	double Omega_Hertz_to_code_units = 203000.0; // unit conversion: amount of Hz that correspond to Omega=1
	double Omega_rot = 2.*M_PI*frequency_Hertz_in / Omega_Hertz_to_code_units; // angular frequency in code units
	std::vector<double> beta_grid(Nstars, 0.0); // later set this equal to R^4 Omega_rot^2 \Gamma^4, \Gamma = 1/ sqrt(1 - R^2 \Omega^2)

	// lowest and highest initial densities for the NSs in the MR curve:
	const double sat_to_code = 0.16 * 2.886376934e-6 * 939.565379;	// conversion factor from nuclear saturation density to code units
	double rho_0_min = rho_0_min_in * sat_to_code; // central restmass density of the NSs in units of the nucelar saturation density
	double rho_0_max = rho_0_max_in * sat_to_code;

	// fill array with the initial conditions for every star:
	utilities::fillValuesPowerLaw(rho_0_min, rho_0_max, rho_c_grid, 2);	// power law scaling of 2 or 3 works pretty well

	// select correct EOS type:
	std::string EOS_filepath = "";
	     if(EOS_in == "EOS_DD2")    { EOS_filepath = "EOS_tables/eos_HS_DD2_with_electrons.beta";}
	else if(EOS_in == "EOS_APR")    { EOS_filepath = "EOS_tables/eos_SRO_APR_SNA_version.beta";}
	else if(EOS_in == "EOS_KDE0v1") { EOS_filepath = "EOS_tables/eos_SRO_KDE0v1_SNA_version.beta";}
	else if(EOS_in == "EOS_LNS")    { EOS_filepath = "EOS_tables/eos_SRO_LNS_SNA_version.beta";}
	else if(EOS_in == "EOS_FSG")    { EOS_filepath = "EOS_tables/eos_HS_FSG_with_electrons.beta";}
	else { std::cout << "Wrong EOS! Supported EOS are: 'EOS_DD2'  'EOS_APR'  'EOS_KDE0v1'  'EOS_LNS'  'EOS_FSG' !" << std::endl; return;}
	auto myEOS = std::make_shared<EoStable>(EOS_filepath); // (load the EOS)

	std::vector<NSEinsteinCartanRotation> MR_curve; // holds the stars in the MR curve
	// compute non-rotating neutron stars without torsion:
	calc_EinsteinCartan_curves_rotation_beta_grid(myEOS, rho_c_grid, MR_curve, beta_grid, 1); // last argument is verbose

	// use the star results to compute a corresponding beta value:
	for(unsigned int i = 0; i < MR_curve.size(); i++) {
		double Rns = MR_curve[i].R_NS;
		double Mns = MR_curve[i].M_T;
		if(Keplarian > 0.0) { // initialize NS with Keplarian rotation rate
			if(Rns > 0.) {Omega_rot = Keplarian*std::sqrt( Mns/(Rns*Rns*Rns) );} else {Omega_rot = 0.0;}
		} 
        beta_grid[i] = pow(Rns,4) * Omega_rot*Omega_rot / pow(1.- Rns*Rns * Omega_rot*Omega_rot,2);
    }
	calc_EinsteinCartan_curves_rotation_beta_grid(myEOS, rho_c_grid, MR_curve, beta_grid, 1); // last argument is verbose

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::string plotname; std::stringstream stream; std::string tmp;
	if(Keplarian > 0.0) {
		plotname = "ECstar_curve_rotation_rate_model_e_plus_P_" + EOS_in + "_frequencyKep_"; stream << std::fixed << std::setprecision(10) << Keplarian;
		stream >> tmp; plotname += tmp;
	} else{
		plotname = "ECstar_curve_rotation_rate_model_e_plus_P_" + EOS_in + "_frequencyHz_"; stream << std::fixed << std::setprecision(10) << frequency_Hertz_in;
		stream >> tmp; plotname += tmp;
	}

	write_MRphi_curve<NSEinsteinCartanRotation>(MR_curve, "output/" + plotname + ".txt");
	std::cout << "calculation complete!" << std::endl;
	std::cout << "parameters: angular frequency (Hz)" << std::endl;
	std::cout << std::fixed << std::setprecision(10) << frequency_Hertz_in << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}

void EC_star_curve_rotation_rate_model_e_plus_P_beta_iteration(unsigned Nstars_in, double rho_0_min_in, double rho_0_max_in, double frequency_Hertz_in, std::string EOS_in, double Keplarian = 0.0) {
	
	unsigned Nstars = Nstars_in; // number of stars in MR curve
	std::vector<double> rho_c_grid(Nstars, 0.0);

	// define parameters for the 'spin densits EOS': s^2 = beta*P^gamma
	double Omega_Hertz_to_code_units = 203000.0; // unit conversion: amount of Hz that correspond to Omega=1
	double Omega_rot = 2.*M_PI*frequency_Hertz_in / Omega_Hertz_to_code_units; // angular frequency in code units
	std::vector<double> beta_grid(Nstars, 0.0); // later set this equal to R^4 Omega_rot^2 \Gamma^4, \Gamma = 1/ sqrt(1 - R^2 \Omega^2)

	// lowest and highest initial densities for the NSs in the MR curve:
	const double sat_to_code = 0.16 * 2.886376934e-6 * 939.565379;	// conversion factor from nuclear saturation density to code units
	double rho_0_min = rho_0_min_in * sat_to_code; // central restmass density of the NSs in units of the nucelar saturation density
	double rho_0_max = rho_0_max_in * sat_to_code;

	// fill array with the initial conditions for every star:
	utilities::fillValuesPowerLaw(rho_0_min, rho_0_max, rho_c_grid, 2);	// power law scaling of 2 or 3 works pretty well

	// select correct EOS type:
	std::string EOS_filepath = "";
	     if(EOS_in == "EOS_DD2")    { EOS_filepath = "EOS_tables/eos_HS_DD2_with_electrons.beta";}
	else if(EOS_in == "EOS_APR")    { EOS_filepath = "EOS_tables/eos_SRO_APR_SNA_version.beta";}
	else if(EOS_in == "EOS_KDE0v1") { EOS_filepath = "EOS_tables/eos_SRO_KDE0v1_SNA_version.beta";}
	else if(EOS_in == "EOS_LNS")    { EOS_filepath = "EOS_tables/eos_SRO_LNS_SNA_version.beta";}
	else if(EOS_in == "EOS_FSG")    { EOS_filepath = "EOS_tables/eos_HS_FSG_with_electrons.beta";}
	else { std::cout << "Wrong EOS! Supported EOS are: 'EOS_DD2'  'EOS_APR'  'EOS_KDE0v1'  'EOS_LNS'  'EOS_FSG' !" << std::endl; return;}
	auto myEOS = std::make_shared<EoStable>(EOS_filepath); // (load the EOS)

	std::vector<NSEinsteinCartanRotation> MR_curve; // holds the stars in the MR curve
	// compute non-rotating neutron stars without torsion:
	calc_EinsteinCartan_curves_rotation_beta_grid(myEOS, rho_c_grid, MR_curve, beta_grid, 1); // last argument is verbose

	// use the star results to compute a corresponding beta value:
	unsigned numIterations = 5; // must be an odd number, otherwise the model can return the no-torsion configuration instead
	unsigned j = 0;
	while(j < numIterations) {
		j++;
		for(unsigned int i = 0; i < MR_curve.size(); i++) {
			double Rns = MR_curve[i].R_NS;
			double Mns = MR_curve[i].M_T;
			if(Keplarian > 0.0) { // initialize NS with Keplarian rotation rate
				if(Rns > 0.) {Omega_rot = Keplarian*std::sqrt( Mns/(Rns*Rns*Rns) );} else {Omega_rot = 0.0;}
			} 
        	beta_grid[i] = pow(Rns,4) * Omega_rot*Omega_rot / pow(1.- Rns*Rns * Omega_rot*Omega_rot,2); // update the beta grid
    	}
		// re-compute all configurations with update beta-value
		calc_EinsteinCartan_curves_rotation_beta_grid(myEOS, rho_c_grid, MR_curve, beta_grid, 1); // last argument is verbose
	}
	

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::string plotname; std::stringstream stream; std::string tmp;
	if(Keplarian > 0.0) {
		plotname = "ECstar_curve_rotation_rate_model_e_plus_P_beta_iteration_" + EOS_in + "_frequencyKep_"; stream << std::fixed << std::setprecision(10) << Keplarian;
		stream >> tmp; plotname += tmp;
	} else{
		plotname = "ECstar_curve_rotation_rate_model_e_plus_P_beta_iteration_" + EOS_in + "_frequencyHz_"; stream << std::fixed << std::setprecision(10) << frequency_Hertz_in;
		stream >> tmp; plotname += tmp;
	}

	write_MRphi_curve<NSEinsteinCartanRotation>(MR_curve, "output/" + plotname + ".txt");
	std::cout << "calculation complete!" << std::endl;
	std::cout << "parameters: angular frequency (Hz)" << std::endl;
	std::cout << std::fixed << std::setprecision(10) << frequency_Hertz_in << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}

// NS in EC with hba rspin density prescription:
// double rho_0 [in saturation density], double beta, double gamma, string EOS_name
void EC_star_single_hbar(double rho_0_in, double eta_tilde_in, std::string EOS_in) {

	// define parameters for the 'spin densits EOS': s^2 = beta*P^gamma
	double eta_tilde = eta_tilde_in;
	// initial density:
	const double sat_to_code = 0.16 * 2.886376934e-6 * 939.565379;	// conversion factor from nuclear saturation density to code units
	double rho_0 = rho_0_in * sat_to_code; // central restmass density of the NS in units of the nucelar saturation density

	// select correct EOS type:
	std::string EOS_filepath = "";
	     if(EOS_in == "EOS_DD2")    { EOS_filepath = "EOS_tables/eos_HS_DD2_with_electrons.beta";}
	else if(EOS_in == "EOS_APR")    { EOS_filepath = "EOS_tables/eos_SRO_APR_SNA_version.beta";}
	else if(EOS_in == "EOS_KDE0v1") { EOS_filepath = "EOS_tables/eos_SRO_KDE0v1_SNA_version.beta";}
	else if(EOS_in == "EOS_LNS")    { EOS_filepath = "EOS_tables/eos_SRO_LNS_SNA_version.beta";}
	else if(EOS_in == "EOS_FSG")    { EOS_filepath = "EOS_tables/eos_HS_FSG_with_electrons.beta";}
	else { std::cout << "Wrong EOS! Supported EOS are: 'EOS_DD2'  'EOS_APR'  'EOS_KDE0v1'  'EOS_LNS'  'EOS_FSG' !" << std::endl; return;}
	auto myEOS = std::make_shared<EoStable>(EOS_filepath); // (load the EOS)

	// initialize one instance:
	NSEinsteinCartanHbarSpin ECstar(myEOS, rho_0, eta_tilde);

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::stringstream stream; std::string tmp;
	std::string plotname = "ECstarHbarSpin_profile_fd2_" + EOS_in + "_rho0_"; stream << std::fixed << std::setprecision(10) << rho_0_in; 
	stream >> tmp; plotname += (tmp + "_eta_tilde_drhodr_");  stream = std::stringstream(); stream << std::fixed << std::setprecision(10) << eta_tilde;
	stream >> tmp; plotname += tmp;

	// evaluate the model and save the intermediate data into txt file:
	std::vector<integrator::step> results;
    ECstar.evaluate_model(results, "output/validation/" + plotname + ".txt");

	// output global variables:
	std::cout << "calculation complete!" << std::endl;
	std::cout << "global quantities:" << std::endl;
	std::vector<std::string> labels = ECstar.labels();
	std::cout << labels[0] << " " << labels[1] << " " << labels[2] << " " << labels[3] << " " << labels[4] << " " << labels[5] << " " << labels[6] << " " << std::endl;
	std::cout << std::fixed << std::setprecision(10) << ECstar << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}




int main() {

	k_form_star_single();

    return 0;
}
