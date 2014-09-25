#include "ma02_src.h"

double hypappell( double a, double b1, double b2, double c, double x, double y )
{
  //This can be written in terms of Gauss hypergeometric functions, which gsl
  //has implementations of.
  //but gsl version quite slow, and does not work for x>1, therefore I've used the scipy version and adapted the C source
  //see http://www.eafit.edu.co/revistas/ingenieria-ciencia/Documents/revista10/productInvHyp.pdf for the equations
  //also looked at mpmath.py module - some useful transformations etc, but can probably still be significantly
  //speeded up with more analytic continuations/transformations for certain ranges of inputs
  
  double ind_term ;
  double sum = 0. ;
  
  // Swap terms so that x<y - the outer loop should have a smaller value for greater accuracy and stability
  if( fabs( x ) > fabs( y ) ) 
    {
      swap_double( &x, &y ) ; //ie x smaller
      swap_double( &b1, &b2 ) ;
    }
  
  // coord transformation if x >= 0.99 (speeds things up?)
  if( fabs( x ) >= 0.99 )
    { 
      double u = ( x - y )/( x - 1 ) ;
      if( fabs( u ) <= 0.99 ) //if u > 0.99 ignore this transformation
	return pow( 1. - x, -b1 )*pow( 1. - y, c - a - b2 )*hypappell( c - a, b1, c - b1 - b2, c, u, y ) ;
    }
  
  // Calculate appell hypergeometric func...
  // using eqs frp, http://www.eafit.edu.co/revistas/ingenieria-ciencia/Documents/revista10/productInvHyp.pdf
  // eq A.5
  double t = 1. ;
  double h ; //first term for r=0
  int r = 0 ;
  while( 1 )
    {
      // Sum up the individual terms from r=0->inf
      h = hyp2f1( a + r, b2, c + r, y ) ;
      ind_term = t*h ;
      sum += ind_term ;
    
      // Break if max iterations or tolerance reached
      if( fabs( ind_term ) < ETOL_HYPER && fabs( h ) > 10.*ETOL_HYPER )
	{
	  break ;
	}
      if( r == MAX_ITER_HYPER )
	{
	  break ;
	}
      
      //increment r and calculate first term in series
      //just need to multiply by these terms to get the new first term (t)
      //(ie don't need to recalculate the pochammer vals)
      r++ ;
      t *= x*( a + r - 1 )*( b1 + r - 1 )/( c + r - 1 )/r ;
    }
  return sum;
}

void swap_double( double *x1, double *x2 )
{
  double t = *x1 ;
  *x1 = *x2 ;
  *x2 = t ;
}
