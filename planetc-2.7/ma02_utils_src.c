#include "ma02_src.h"

/**************************************************************************************************/
//Functions for flux_quad - elliptic integrals etc

double lambda_e_pi_funct( double NormSep, double RpRs )
{
  double f ;
  f = sqrt( ( 4.*POW2( NormSep ) - POW2( 1. + POW2( NormSep ) - POW2( RpRs ) ) )/4. ) ;
  return ( ( POW2( RpRs ) * kappa_0( NormSep, RpRs ) + kappa_1( NormSep, RpRs ) - f )/M_PI ) ;
}


double kappa_1( double NormSep, double RpRs )
{
  return acos( ( 1. - POW2( RpRs ) + POW2( NormSep ) )/( 2.*NormSep ) ) ;
}


double kappa_0( double NormSep, double RpRs )
{
  return acos( ( POW2( RpRs ) + POW2( NormSep ) - 1. )/( 2.*RpRs*NormSep ) ) ;
}


double lambda_1( double NormSep, double RpRs, double a, double b, double q, double k )
{
  return ( ( 1./( 9.*M_PI*sqrt( RpRs*NormSep ) ) )*( ( ( 1. - b )*( 2.*b + a - 3. ) 
	 - 3.*q*( b-2. ) )*gsl_sf_ellint_Kcomp( k, MODE ) 
         + 4.*RpRs*NormSep*( POW2( NormSep ) 
         + 7.*POW2( RpRs ) - 4. )*gsl_sf_ellint_Ecomp( k, MODE )
	 - 3.*( q/a )*( gsl_sf_ellint_Kcomp( k, MODE )
         - ( ( 1./a - 1. )/3. )*gsl_sf_ellint_RJ( 0., 1. - POW2( k ), 1., 1./a, MODE) ) ) ) ;
}


double lambda_2( double NormSep, double RpRs,double a, double b, double q, double k )
{
  return ( ( 2./( 9.*M_PI*sqrt( 1. - a ) ) )*( ( 1. - 5.*POW2( NormSep )
         + POW2( RpRs ) + POW2( q ) )*gsl_sf_ellint_Kcomp( ( 1./k ), MODE )
	 + ( 1. - a )*( POW2( NormSep ) + 7*POW2( RpRs ) - 4. )
         *gsl_sf_ellint_Ecomp( ( 1./k ), MODE )
	 - 3.*( q/a )*( gsl_sf_ellint_Kcomp( 1/k, MODE ) 
	 - ( ( b/a - 1. )/3. )*gsl_sf_ellint_RJ( 0., 1. - POW2( 1./k ), 1., b/a, MODE) ) ) ) ;
}


double lambda_3( double RpRs, double k )
{
  return ( 1./3. + ( ( 16.*RpRs )/( 9.*M_PI ) )*( 2.*POW2( RpRs ) - 1. )
         *gsl_sf_ellint_Ecomp( ( 1./( 2.*k ) ), MODE )
	 - ( ( 32.*POW2( RpRs ) -20. *RpRs + 3. )/( 9.*M_PI*RpRs ) )
         *gsl_sf_ellint_Kcomp( ( k/2. ), MODE ) ) ;
}


double lambda_4( double RpRs, double k )
{
  return ( 1./3. + ( 2./( 9.*M_PI ) )*( 4.*( 2.*POW2( RpRs ) - 1. )
         *gsl_sf_ellint_Ecomp( 1./k, MODE )
         + ( 1. - 4.*POW2( RpRs ) )*gsl_sf_ellint_Kcomp( 1./k, MODE ) ) ) ;
}


double lambda_5( double RpRs )
{
  return ( ( 2./( 3.*M_PI ) )*( acos( 1. - 2.*RpRs ) ) - ( 4./( 9.*M_PI ) )
         *sqrt( RpRs*( 1. - RpRs ) )*( 3. + 2.*RpRs - 8.*POW2( RpRs ) ) ) ;
}


double lambda_6( double RpRs )
{
  return ( -( 2./3. )*pow( ( 1. - POW2( RpRs ) ), 3./2. ) ) ;
}


double eta_1( double NormSep, double RpRs, double a, double b )
{
  return ( ( 1./( 2.*M_PI ) )*( kappa_1( NormSep, RpRs )
         + 2.*eta_2( NormSep, RpRs )*kappa_0( NormSep, RpRs )
	 - ( 1./4. )*( 1. + 5.*POW2( RpRs ) + POW2( NormSep ) )
         *sqrt( ( 1.-a )*( b-1. ) ) ) ) ;
}


double eta_2( double NormSep, double RpRs )
{
  return ( ( POW2( RpRs )/2. )*( POW2( RpRs ) + 2.*POW2( NormSep ) ) ) ;
}


double Omega( double * C )
{
  double Om = 0. ;
  int n ;
  for ( n=0. ; n<5 ; n++ )
    {
      Om += C[n]/( n + 4. ) ;
    }
  return Om ;
  }


double Nfunc( double NormSep, double RpRs, int n, double a, double b )
{
  // The N function given by Eq 3 in Madel & Agol (2002).

  double N ;
  
  N = pow( 1. - a, ( n + 6. )/4. )/pow( b - a, 0.5 ) ;
  N *= gsl_sf_beta( ( n + 8. )/4. , 0.5 ) ;
  N *= ( POW2( NormSep ) - POW2( RpRs ) )/a
       * hypappell( 0.5, 1., 0.5, ( n + 10. )/4., ( a - 1. )/a, ( 1. - a )/( b - a ) )
       - hyp2f1( 0.5, 0.5, ( n + 10. )/4., ( 1. - a )/( b - a ) ) ;
  
  return N;
}


double Mfunc( double NormSep, double RpRs, int n, double a, double b )
{
  // The M function given by Eq 4 in Mandel & Agol (2002). 

  double M ;

  M = pow( ( 1. - a ), ( n + 4. )/4. ) ;
  M *= ( POW2( NormSep ) - POW2( RpRs ) )/a
       * hypappell( 0.5, -( n + 4. )/4., 1., 1., ( b - a )/( 1. - a), ( a - b )/a )
       - hyp2f1( -( n + 4. )/4., 0.5, 1., ( b - a )/( 1. - a ) ) ;

  return M ;
  }

/* double hyp_appell(double a,double b1,double b2,double c,double x,double y) */
/* { */
/*   // Neale Gibson's custom-written Appell hypergeometric routine. Here are his original comments: */
/*   // */
/*   // This can be written in terms of Gauss hypergeometric functions, which gsl has implementations  */
/*   // of. */
/*   // But gsl version quite slow, and does not work for x>1, therefore I've used the scipy version */
/*   // and adapted the C source  */
/*   // see http://www.eafit.edu.co/revistas/ingenieria-ciencia/Documents/revista10/productInvHyp.pdf  */
/*   // for the equations also looked at mpmath.py module - some useful transformations etc, but can */
/*   // probably still be significantly speeded up with more analytic continuations/transformations  */
/*   // for certain ranges of inputs */
  
/*   double ind_term ; */
/*   double sum ; */
  
/*   // Swap terms so that x<y - the outer loop should  */
/*   // have a smaller value for greater accuracy and  */
/*   // stability: */
/*   if( fabs( x ) > fabs( y ) ) */
/*     { */
/*       swap_double( &x, &y ) ; //ie x smaller */
/*       swap_double( &b1, &b2 ) ; */
/*     } */
  
/*   if( fabs( x ) >= 0.99 ) */
/*     {  */
/*       // Coord transformation if x >= 0.99  */
/*       // (speeds things up?) */
/*       double u = ( x - y )/( x - 1 ) ; */
/*       if( fabs( u ) <= 0.99 ) */
/* 	{ */
/* 	  // But if u > 0.99 ignore this  */
/* 	  // transformation */
/* 	  return pow( 1. - x, -b1 )*pow( 1. - y, c - a - b2 )*hyp_appell( c - a, b1, c - b1 - b2, c, u, y ) ; */
/* 	} */
/*       // Calculate appell hypergeometric func... */
/*       // using eqs frp, http://www.eafit.edu.co/revistas/ingenieria-ciencia/Documents/revista10/productInvHyp.pdf */
/*       // eq A.5 */
/*       double t = 1. ; */
/*       double h ; //first term for r=0 */
/*       int r = 0 ; */
/*       while( 1 ) */
/* 	{ */
/* 	  //sum up the individual terms from r=0->inf */
/* 	  h = hyp2f1( a + r, b2, c + r, y ) ; */
/* 	  ind_term = t*h ; */
/* 	  sum += ind_term ; */
/* 	  //break if max iterations or tolerance reached */
/* 	  if( fabs( ind_term ) < ETOL_APPELL && fabs( h ) > 10.*ETOL_APPELL )  */
/* 	    { */
/* 	      break ; */
/* 	    } */
/* 	  if( r==MAX_ITER_APPELL ) */
/* 	    { */
/* 	      //printf("warning: appellF1 eval - max iterations reached\n"); */
/* 	      break ; */
/* 	    } */
/* 	  //increment r and calculate first term in series */
/* 	  //just need to multiply by these terms to get the new first term (t) */
/* 	  //(ie don't need to recalculate the pochammer vals) */
/* 	  r++ ; */
/* 	  t *= x*( a + r - 1. )*( b1 + r - 1. )/( c + r - 1. )/r ; */
/* 	} */
/*       return sum ; */
/*     } */

/* void swap_double( double *x1, double *x2) */
/* { */
/*     double t = *x1 ; */
/*     *x1 = *x2 ; */
/*     *x2 = t ; */
/* } */
