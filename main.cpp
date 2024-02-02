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

// --------------------------------------------------------------------
using namespace FBS;

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
void EC_star_curve_rotation_rate(unsigned Nstars_in, double rho_0_min_in, double rho_0_max_in, double frequency_Hertz_in, std::string EOS_in, bool Keplarian = false) {
	
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
		if(Keplarian) { // initialize NS with Keplarian rotation rate
			if(Rns > 0.) {Omega_rot = std::sqrt( Mns/(Rns*Rns*Rns) );} else {Omega_rot = 0.0;}
		} 
        beta_grid[i] = pow(Rns,4) * Omega_rot*Omega_rot / pow(1.- Rns*Rns * Omega_rot*Omega_rot,2);
    }
	calc_EinsteinCartan_curves_beta_grid(myEOS, rho_c_grid, MR_curve, beta_grid, 2.0, 1);	// last argument is verbose

	// name of output textfile. Use stringstream for dynamic naming of output file:
	std::stringstream stream; std::string tmp;
	std::string plotname = "ECstar_curve_rotation_rate_" + EOS_in + "_frequencyHz_"; stream << std::fixed << std::setprecision(10) << frequency_Hertz_in;
	stream >> tmp; plotname += tmp;

	write_MRphi_curve<NSEinsteinCartan>(MR_curve, "output/" + plotname + ".txt");
	std::cout << "calculation complete!" << std::endl;
	std::cout << "parameters: angular frequency (Hz)" << std::endl;
	std::cout << std::fixed << std::setprecision(10) << frequency_Hertz_in << std::endl;
	std::cout << "Name of output file:" << std::endl << "output/" + plotname + ".txt" << std::endl;
}


int main() {

    // integrate a single Einstein-Cartan star:
	// double rho0 [in saturation density], double beta, double gamma, string EOS_name:
    //EC_star_single(1.0, 10.0, 2.0, "EOS_DD2");
	//EC_star_single(4.0, 0.0, 2.0, "EOS_APR");

	// integrate a MR curve of Einstein-Cartan stars:
	// unsigned Nstars, double rho0_min [in saturation density], double rho0_max [in saturation density], double beta, double gamma, string EOS_name:
	/*EC_star_curve(200, 0.2, 12.0, 0.0, 2.0, "EOS_DD2");
	EC_star_curve(200, 0.2, 12.0, 10.0, 2.0, "EOS_DD2");
	EC_star_curve(200, 0.2, 12.0, 20.0, 2.0, "EOS_DD2");
	EC_star_curve(200, 0.2, 12.0, 100.0, 2.0, "EOS_DD2");*/

    // ----------------------------------------------------------------
	// Plots to create data for the figures in the paper:
	// Figure 1:
	//EC_star_single(4.0, 0.0, 2.0, "EOS_DD2");
	//EC_star_single(4.0, 0.0, 2.0, "EOS_APR");
	// Figure 2:
	// ....
	// Figure 3:
	EC_star_curve_const_mass_with_different_beta(200, 0.8, "M_T", 0.0, 10., 1.5, "EOS_DD2");
	EC_star_curve_const_mass_with_different_beta(200, 1.0, "M_T", 0.0, 10., 1.5, "EOS_DD2");
	EC_star_curve_const_mass_with_different_beta(200, 1.4, "M_T", 0.0, 10., 1.5, "EOS_DD2");
	EC_star_curve_const_mass_with_different_beta(200, 2.0, "M_T", 0.0, 10., 1.5, "EOS_DD2");
	// Figure 4:
	EC_star_curve_rotation_rate(200, 0.6, 10.0, 0., "EOS_DD2", false); // non-rotating NS
	EC_star_curve_rotation_rate(200, 0.6, 10.0, 714., "EOS_DD2", false);
	EC_star_curve_rotation_rate(200, 0.6, 10.0, 1., "EOS_DD2", true); // Kelparian rotating NS

	EC_star_curve_rotation_rate(200, 0.7, 10.0, 0., "EOS_APR", false); // non-rotating NS
	EC_star_curve_rotation_rate(200, 0.7, 10.0, 714., "EOS_APR", false);
	EC_star_curve_rotation_rate(200, 0.7, 10.0, 1., "EOS_APR", true); // Kelparian rotating NS

    return 0;
}
