#ifndef _H_LIMBER
#define _H_LIMBER

#include <gsl/gsl_spline.h>
#include <stdbool.h>
#include "interp2d.h"


// The integrated C_ell are in general allowed to be zero or negative if
// they describe cross-correlations. We use these statuses to describe
// errors where a log was requested also.
// LIMBER_STATUS_NEGATIVE is probably always an error
// but LIMBER_STATUS_ZERO is not necessarily
#define LIMBER_STATUS_OK 0
#define LIMBER_STATUS_ZERO 1
#define LIMBER_STATUS_NEGATIVE 2
#define LIMBER_STATUS_ERROR 3

// These are the options you can set for
// the Limber integrator.
typedef struct limber_config{
	bool xlog;  // The output spline will be in terms of log(ell) not ell
	bool ylog;  // The output spline will return log(C_ell) not C_ell
	int n_ell;  // Number of ell values you want in the spline
	double * ell;  // The chosen ell values you want
	double prefactor; //Scaling prefactor
    int status; // did everything go okay?
    double absolute_tolerance;
    double relative_tolerance;
} limber_config;


// Do a flat universe Limber approximation integral, of the form:
// C^{XY}(\ell) 
// = A \int_0^{\chi_{\mathrm{max}}} 
//       \frac{W^X(\chi) W^Y(\chi)}{\chi^2} P(\ell/\chi, \chi)

// The two splines and the matter power passed into the function
// be FUNCTIONS OF CHI NOT Z

// TODO - Fortran version

gsl_spline * limber_integral(limber_config * config, 
	gsl_spline * WX, gsl_spline * WY, Interpolator2D * P);


#endif