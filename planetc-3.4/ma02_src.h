#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_math.h" 
#include "gsl/gsl_sf_ellint.h" //for elliptic integrals
#include "gsl/gsl_sf_gamma.h" //for gamma func, beta func, factorials, Pochhammer symbol

#include "mconf.h"


#define MODE GSL_PREC_DOUBLE
#define TOL 2e-12
#define POW2( x ) ( ( x ) * ( x ) )
#define DOUBLE_EQ( x, v ) ( ( ( v - TOL ) < x ) && ( x < ( v + TOL ) ) )

double F_quad( double NormSep, double RpRs, double gam1, double gam2 ) ;
double F_nonlin( double NormSep, double RpRs, double c1, double c2, double c3, double c4 ) ;

double kappa_0( double NormSep, double RpRs ) ;
double kappa_1( double NormSep, double RpRs ) ;
double lambda_e_pi_funct( double NormSep, double RpRs) ;
double lambda_1( double NormSep, double RpRs, double a, double b, double q, double k ) ;
double lambda_2( double NormSep, double RpRs, double a, double b, double q, double k ) ;
double lambda_3( double RpRs, double k ) ;
double lambda_4( double RpRs, double k ) ;
double lambda_5( double RpRs) ;
double lambda_6( double RpRs) ;
double eta_1( double NormSep, double RpRs, double a, double b ) ;
double eta_2( double NormSep, double RpRs ) ;
double Omega( double * C ) ;
double Nfunc( double NormSep, double RpRs, int n, double a, double b ) ;
double Mfunc( double NormSep, double RpRs, int n, double a, double b ) ;

// Neale's code did some hackwork to define the Gauss hypergeometric and Appell hypergeomtric functions; 
// I'm not quite sure why he didn't use the GSL versions, but he probably had a good reason; I'm going to 
// try them anyway just in case.
double hypappell( double a, double b1, double b2, double c, double x, double y ) ;
void swap_double( double *x1, double *x2 ) ; //swap doubles
#define ETOL_APPELL 1e-8
#define MAX_ITER_APPELL 50
#define ETOL_HYPER 1e-8
#define MAX_ITER_HYPER 50


