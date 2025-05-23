#pragma once

#include <utility>  // for std::swap

#include "vector.hpp"
#include "eos.hpp"
#include "integrator.hpp"
#include "plotting.hpp"

namespace FBS {

#define R_INIT 1e-10    // the initial integration radius
#define R_MAX 500.      // the general maximum integration radius (might be increased if not sufficient)
#define P_ns_min 1e-15  // the minimum pressure for the "boundary" of the NS

/*  NSmodel
 * this is an abstract class that is supposed to be the backbone for
 *  a physical model of a neutron star
 * This implementation allows the integration of the differential equations
 *  describing the neutron star
 * The r_init and r_end values describe the integration range. This can depend
 *  on the parameters involved, so it is a class member.
 * */

class NSmodel {
protected:
    double r_init, r_end;
public:
    /* A pointer to the Equation of State describing the neutron star matter */
    std::shared_ptr<EquationOfState> EOS;

    /* Constructor */
    NSmodel(std::shared_ptr<EquationOfState> EOS) : r_init(R_INIT), r_end(R_MAX), EOS(EOS) {}
    NSmodel() : r_init(R_INIT), r_end(R_MAX) {}

    /* This function gives the derivatives for the differential equations */
    virtual vector dy_dr(const double r, const vector& vars) const = 0;

    /* This is a static wrapper function, that calls the dy_dr function */
    static vector dy_dr_static(const double r, const vector& y, const void* params);

    /* The initial conditions for a, alpha, phi, Psi, and P */
    virtual vector get_initial_conditions(double r_init=R_INIT) const = 0;

    /* This function calls the integrator and returns the results of the integration */
    int integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, const vector initial_conditions, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=-1., double r_end=-1.) const;

    /* For easy output the class should define an << operator */
    friend std::ostream& operator<<(std::ostream&, const NSmodel&);
    /* The labels for the different parameters that are output by the << operator */
    static std::vector<std::string> labels();
};

/* Same as NSmodel but including previous-step functionality */
class NSmodelv2 {
protected:
    double r_init, r_end;
public:
    /* A pointer to the Equation of State describing the neutron star matter */
    std::shared_ptr<EquationOfState> EOS;

    /* Constructor */
    NSmodelv2(std::shared_ptr<EquationOfState> EOS) : r_init(R_INIT), r_end(R_MAX), EOS(EOS) {}

    /* This function gives the derivatives for the differential equations */
    virtual vector dy_dr(const double r, const vector& vars, const double r_prev, const vector& vars_prev) const = 0;

    /* This is a static wrapper function, that calls the dy_dr function */
    static vector dy_dr_static(const double r, const vector& y, const double r_prev, const vector& y_prev, const void* params);

    /* The initial conditions for a, alpha, phi, Psi, and P */
    virtual vector get_initial_conditions(double r_init=R_INIT) const = 0;

    /* This function calls the integrator and returns the results of the integration */
    int integrate(std::vector<integrator::step>& result, std::vector<integrator::Event>& events, const vector initial_conditions, integrator::IntegrationOptions intOpts = integrator::IntegrationOptions(), double r_init=-1., double r_end=-1.) const;

    /* For easy output the class should define an << operator */
    friend std::ostream& operator<<(std::ostream&, const NSmodel&);
    /* The labels for the different parameters that are output by the << operator */
    static std::vector<std::string> labels();
};

}

