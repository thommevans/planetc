#include "keporb_src.h"

// TODO = Explain more fully how the XYZ coordinate
// system is defined. I've made a start here...
// NOTE: To be explicit, in the above coordinate system, 
// the X axis increases along some arbitrary reference 
// direction, defined for convenience such that the 
// longitude of periapse "Omega" (not the argument of 
// periapse "omega") is equal to 180 degrees, or in other
// words, the angle between the ascending node (the point 
// on the planet's orbit at which it crosses the imaginary
// sky plane moving away from the observer) and the
// arbitrary reference direction is 180 degrees:


double NormSep( double MeanAnom, double aRs, double ecc, double omega, double incl )
{
  // This routine returns the sky-projected distance 
  // between the star and planet in units of stellar
  // radii, which can then be used as input to the 
  // Mandel & Agol (2002) equations for calculating 
  // transit lightcurves.
	
  double cos_TrueAnom, sin_TrueAnom, cos_omega, sin_omega, EccAnom, r ;
  double X, Y, D ;
  
  // Ensure the mean anomaly value is between [0,2Pi]:
  while( MeanAnom > 2*M_PI )
    {
      MeanAnom -= 2*M_PI ;
    }
  while( MeanAnom < -2*M_PI )
    {
      MeanAnom += 2*M_PI ;
    }

  if ( ecc == 0. )
    {
      // If the orbit is circular, the eccentirc anomaly
      // is just the same as the mean anomaly:
      EccAnom = MeanAnom ;
    }
  else
    {
      // If the orbit is non-circular, then we must solve
      // the Kepler Equation to get the eccentric anomaly
      // from the true anomaly and eccentricity:
      EccAnom = calc_EccAnom( MeanAnom, ecc ) ;
    }

  cos_TrueAnom = calc_cos_TrueAnom( EccAnom, ecc ) ;
  sin_TrueAnom = calc_sin_TrueAnom( EccAnom, ecc ) ;
  cos_omega = cos( omega ) ;
  sin_omega = sin( omega ) ;
  
  // Calculate the radial distance between the star and
  // planet in units of stellar radii:
  r = aRs * ( 1 - POW2( ecc ) ) / ( 1. + ecc * cos_TrueAnom ) ;
  
  // Calculate the X and Y coordinates of the planet
  // on the sky plane:
  X = - r * ( cos_TrueAnom*cos_omega - sin_TrueAnom*sin_omega ) ;
  Y = - r * ( sin_TrueAnom*cos_omega + cos_TrueAnom*sin_omega ) * cos( incl ) ;
  
  // Finally, we can calculate the sky-plane projected 
  // distance between the star and planet, in units of
  // stellar radii:
  D = sqrt( POW2( X ) + POW2( Y ) ) ;
  
  return D ;
}


double Xcoord( double MeanAnom, double aRs, double ecc, double omega )
{
  // This routine returns the sky-projected X coordinate
  // of the planet. On its own, this is not a particularly
  // useful quantity, unless the aim is to visualise where
  // on the sky the planet would appear. In terms of 
  // computing transit lightcurves, the important quantity
  // is the normalised on-sky separation between the star
  // and planet, which is sqrt( X^2 + Y^2 ), and is calculated
  // by the NormSep() routine above.
  
  double cos_TrueAnom, sin_TrueAnom, cos_omega, sin_omega, EccAnom, r, X ;
  
  // Ensure the mean anomaly value is between [0,2Pi]:
  while( MeanAnom > 2*M_PI)
    {
      MeanAnom -= 2*M_PI;
    }
  while( MeanAnom < -2*M_PI )
    {
      MeanAnom += 2*M_PI;
    }
  
  // Numerically solve the Kepler Equation to get the 
  // eccentric anomaly and use this to calculate the 
  // true anomaly:
  EccAnom = calc_EccAnom( MeanAnom, ecc ) ;
  cos_TrueAnom = calc_cos_TrueAnom( EccAnom, ecc ) ;
  sin_TrueAnom = calc_sin_TrueAnom( EccAnom, ecc ) ;
  cos_omega = cos( omega ) ;
  sin_omega = sin( omega ) ;
  
  // Calculate the radial distance between the star and
  // planet in units of stellar radii:
  r = aRs * ( 1 - POW2( ecc ) ) / ( 1. + ecc * cos_TrueAnom ) ;
  
  // Calculate the X coordinate of the planet on the sky
  // plane:
  X = - r * ( cos_TrueAnom*cos_omega - sin_TrueAnom*sin_omega ) ;

  return X ;
}


double Ycoord( double MeanAnom, double aRs, double ecc, double omega, double incl )
{
  // This routine returns the sky-projected Y coordinate
  // of the planet. On its own, this is not a particularly
  // useful quantity, unless the aim is to visualise where
  // on the sky the planet would appear. In terms of 
  // computing transit lightcurves, the important quantity
  // is the normalised on-sky separation between the star
  // and planet, which is sqrt( X^2 + Y^2 ), and is calculated
  // by the NormSep() routine above.

  double cos_TrueAnom, sin_TrueAnom, cos_omega, sin_omega, EccAnom, r, Y ;
  
  // Ensure the mean anomaly value is between [0,2Pi]:
  while( MeanAnom > 2*M_PI)
    {
      MeanAnom -= 2*M_PI;
    }
  while( MeanAnom < -2*M_PI )
    {
      MeanAnom += 2*M_PI;
    }
  
  // Numerically solve the Kepler Equation to get the 
  // eccentric anomaly and use this to calculate the 
  // true anomaly:
  EccAnom = calc_EccAnom( MeanAnom, ecc ) ;
  cos_TrueAnom = calc_cos_TrueAnom( EccAnom, ecc ) ;
  sin_TrueAnom = calc_sin_TrueAnom( EccAnom, ecc ) ;
  cos_omega = cos( omega ) ;
  sin_omega = sin( omega ) ;

  // Calculate the radial distance between the star and
  // planet in units of stellar radii:
  r = aRs * ( 1 - POW2( ecc ) ) / ( 1. + ecc * cos_TrueAnom ) ;
  
  // Calculate the Y coordinate of the planet on the sky
  // plane:
  Y = - r * ( sin_TrueAnom*cos_omega + cos_TrueAnom*sin_omega ) * cos( incl ) ;

  return Y ;
}


double Zcoord( double MeanAnom, double aRs, double ecc, double omega, double incl )
{
  // This routine returns the Z coordinate, which is perpendicular
  // to the plane of the sky and increases away from the observer.
  // When an on-sky overlap between the stellar disc and planet disc
  // occurs from the point of view of an observer:
  //    i.e. when sqrt( X^2 + Y^2 ) < 1
  // the Z coordinate returned by this function can be used to 
  // distinguish between it being a primary transit (Z<0) or a 
  // secondary eclipse (Z>0).

  double cos_TrueAnom, sin_TrueAnom, cos_omega, sin_omega, EccAnom, r, Z ;
	
  // Ensure the mean anomaly value is between [0,2Pi]:
  while( MeanAnom > 2*M_PI)
    {
      MeanAnom -= 2*M_PI;
    }
  while( MeanAnom < -2*M_PI )
    {
      MeanAnom += 2*M_PI;
    }
  
  // Numerically solve the Kepler Equation to get the 
  // eccentric anomaly and use this to calculate the 
  // true anomaly:
  EccAnom = calc_EccAnom( MeanAnom, ecc ) ;
  cos_TrueAnom = calc_cos_TrueAnom( EccAnom, ecc ) ;
  sin_TrueAnom = calc_sin_TrueAnom( EccAnom, ecc ) ;
  cos_omega = cos( omega ) ;
  sin_omega = sin( omega ) ;
  
  // Calculate the radial distance between the star and
  // planet in units of stellar radii:
  r = aRs * ( 1 - POW2( ecc ) ) / ( 1 + ecc * cos_TrueAnom ) ;
  
  // Calculate the Z coordinate:
  Z = + r * ( sin_TrueAnom*cos_omega + cos_TrueAnom*sin_omega ) * sin( incl ) ;

  return Z ;
}


double calc_EccAnom( double MeanAnom, double ecc )
{
  // Solves the Kepler equation numerically using a Newton-Raphson 
  // iteration method, as outlined in "Solar System Dynamics" by 
  // Murray & Dermott (1999).
  
  double EccAnom, f_EccAnom, fder_EccAnom;
  
  // Initial guess for EccAnom:
  if( sin( MeanAnom ) > 0 )
    {
      EccAnom = MeanAnom + 0.85*ecc ;
    }
  else
    {
    EccAnom = MeanAnom - 0.85*ecc ;
    }

  // Iteratively refine estimate of EccAnom:
  do
    {
      f_EccAnom = EccAnom - ecc*sin( EccAnom ) - MeanAnom ;
      fder_EccAnom = 1 - ecc*cos( EccAnom ) ;
      EccAnom = EccAnom - f_EccAnom/fder_EccAnom ;
    }
  while( f_EccAnom > 0.000001 || f_EccAnom < -0.000001) ;
  
  return EccAnom;
}


double calc_cos_TrueAnom( double EccAnom, double ecc )
{
  // Calculates the cosine of the true anomaly. A derivation of
  // the formula can be found on Wolfram Mathworld in the 
  // "Eccentric Anomaly" article.
  // 
  // The cosine is calculated because the calculation of the
  //  sky-projected normalised separation between the centres 
  // of the planetary and stellar discs only requires the sine
  // and cosine of the true anomaly. If we were to calculate
  // the true anomaly via an inverse trig function, the 
  // computation time is increased a tiny bit, and we also 
  // introduce a tiny bit of additional numerical error.
  return ( cos( EccAnom ) - ecc )/( 1. - ecc*cos( EccAnom ) ) ;
}

double calc_sin_TrueAnom( double EccAnom, double ecc )
{
  // Calculates the sine of the true anomaly. A derivation of
  // the formula can be found on Wolfram Mathworld in the 
  // "Eccentric Anomaly" article.
  // 
  // The cosine is calculated because the calculation of the
  //  sky-projected normalised separation between the centres 
  // of the planetary and stellar discs only requires the sine
  // and cosine of the true anomaly. If we were to calculate
  // the true anomaly via an inverse trig function, the 
  // computation time is increased a tiny bit, and we also 
  // introduce a tiny bit of additional numerical error.
  return sqrt( 1. - POW2( ecc ) )*sin( EccAnom )/( 1. - ecc*cos( EccAnom ) );
}

double calc_TrueAnom( double EccAnom, double ecc )
{
  // Calculates the true anomaly.
  // 
  // Note that we use the sine and cosine of the true anomaly
  // when calculating the normalised on-sky separation of the
  // stellar and planetary discs, and these are calculated by
  // separate routines. However, the present routine could be
  // used if the true anomly itself is desired.

  double xcirc, ycirc ;

  // We use the formula for the true anomaly:
  //
  //   tan( TrueAnom/2 ) = sqrt[ ( 1+ecc )/( 1-ecc ) ]*tan( EccAnom/2 )
  //
  // However, we must allow for the possibility that TrueAnom/2 is
  // outside the range [ -pi/2, +pi/2 ]. This means that we cannot
  // simply use the atan() function, which only maps to this range.
  // Instead, we should remember the definition of sin(theta) as 
  // the y-coordinate of a point on the unit circle with centre at
  // the origin of the xy-plane, and cos(theta) as the x-coordinate,
  // where theta is the angle between the positive x-axis and the
  // vector (x,y), increasing in the counter-clockwise direction. 
  // Then tan(theta) is simply defined as the ratio sin(theta)/cos(theta). 
  // So from above, if we define theta=TrueAnom/2, the corresponding
  // xy coordinates on the unit circle are:
  double sss = sqrt( 1. + ecc ) ;
  ycirc = sss * sin( EccAnom/2. ) ;
  xcirc = sss * cos( EccAnom/2. ) ;
  // where xcirc and ycirc can both obviously be either positive or 
  // negative, depending on the value of TrueAnom/2. We can properly
  // account for this using the atan2() function of the standard
  // maths library:
  
  return 2 * atan2( ycirc, xcirc ) ;
}

