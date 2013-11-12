/*********************************************************************************
This function uses the analytic equations of Mandel and Agol (2002) to compute the
flux of a star during occultation given a normalized separation z, ratio of planet
to star radii, and the two quadratic limb darkening coefficients gam1 and gam2.
*********************************************************************************/

#include "ma02_src.h"

double F_quad( double NormSep, double RpRs, double gam1, double gam2 )
{
	

  // Declare and calculate the variables a, b and q:
  double a = POW2( NormSep - RpRs ) ;
  double b = POW2( NormSep + RpRs ) ;
  double q = ( POW2( RpRs ) - POW2( NormSep ) ) ;
  double k = sqrt( ( 1. - a ) / ( 4.*RpRs*NormSep ) ) ;
  
  // Declare and calculate the limb darkening coeffs:
  double c2 = ( gam1 + 2.*gam2 ) ;
  double c4 = ( -gam2 ) ;
	
  // Declare the main function variables:
  double lambda_e ;
  double lambda_d ;
  double eta_d ;
  double F ;
  
  // Declare and set to -1 a flag that will keep track of
  // whether or not the input RpRs and NormSep values can
  // be matched to one of the cases given in Table 1 of the
  // Mandel & Agol (2002) paper:
  int _case = -1 ;

  // Before going further, check that both the NormSep and
  // RpRs variables have non-negative values:
  if ( RpRs < 0. ) 
    { 
      printf( "RpRs = %+.12f\n", RpRs ) ;
      printf( "RpRs must be >=0 for flux calculation !!\n" ) ;
      exit( -1 ) ;
    }
  if ( NormSep < 0. ) 
    { 
      printf( "NormSep = %+.12f\n", NormSep ) ;
      printf( "NormSep must be >=0 for flux calculation !!\n" ) ;
      exit( -1 ) ;
    }

  // We need to determine which of the 11 configurations 
  // listed in Table 1 of Mandel & Agol (2002) we're dealing
  // with. We can start with the two trivial cases:

  // Case 1  
  // The stellar and planetary discs are not even
  // overlapping or the planetary disc has zero size.
  if ( NormSep > ( 1. + RpRs - TOL ) || DOUBLE_EQ( RpRs, 0. ) )
    {
      _case = 1 ;
      F = 1. ;
      return F ;
    }

  // Case 11 
  // The non-luminous planetary disc is entirely
  // covering the stellar disc.
  else if ( RpRs > 1. - TOL && NormSep < RpRs - 1. + TOL )
    {
      _case = 11 ;
      F = 1. ;
      return F ;
    }

  // Now we proceed to the non-trivial cases, in which the
  // planetary disc is occulting some nonzero fraction of 
  // the limb-darkedned stellar disc. To do this, we need 
  // to compute the relevant lambda_d and eta_d terms in the
  // right-hand column of Table 1 (Mandel & Agol, 2002), as
  // well as the lambda_e value which gives the fraction of
  // the stellar disc that is overlapped by the planet.

  // Case 2
  // The planetary disc overlaps the stellar disc but 
  // without covering the stellar centre.
  else if ( NormSep > 0.5 + fabs( RpRs - 0.5 ) - TOL && 
	    NormSep < 1. + RpRs + TOL )
    {
      _case = 2 ;
      lambda_e = lambda_e_pi_funct( NormSep, RpRs ) ;
      lambda_d = lambda_1( NormSep, RpRs, a, b, q, k ) ;
      eta_d = eta_1( NormSep, RpRs, a, b ) ;
    }
  
  // Case 3
  // The planetary disc is fully overlapping the stellar
  // disc, but the stellar centre is not covered.
  else if ( RpRs < 0.5 + TOL &&
            NormSep > RpRs - TOL &&
            NormSep < 1. - RpRs + TOL )
    {
      _case = 3 ;
      lambda_e = POW2( RpRs ) ;
      lambda_d = lambda_2( NormSep, RpRs, a, b, q, k ) ;
      eta_d = eta_2( NormSep, RpRs ) ;
    }

  // Case 4
  // The same as Case 3 with the added condition that the
  // edges of the planetary disc and stellar disc just touch.
  else if ( RpRs < 0.5 + TOL &&
            DOUBLE_EQ( NormSep, 1. - RpRs ) )
    {
      _case = 4 ;
      lambda_e = POW2( RpRs ) ;
      lambda_d = lambda_5( RpRs ) ;
      eta_d = eta_2( NormSep, RpRs ) ;
    }

  // Case 5
  // Very similar to Case 3 except that the edge of the 
  // planetary disc just touches the stellar centre.
  else if ( RpRs < 0.5 + TOL &&
            DOUBLE_EQ( NormSep, RpRs ) )
    {
      _case = 5 ;
      lambda_e = POW2( RpRs ) ;
      lambda_d = lambda_4( RpRs, k ) ;
      eta_d = eta_2( NormSep, RpRs ) ;      
    }

  // Case 6
  // Similar to cases 3, 4 and 5, but for the even more
  // specific situation where the edges of the planetary
  // disc are simultaneously just touching the centre of
  // the stellar disc and the edge of the stellar disc.
  else if ( DOUBLE_EQ( RpRs, 0.5 ) &&
	    DOUBLE_EQ( NormSep, 0.5 ) )
    {
      _case = 6 ;
      lambda_e = POW2( RpRs ) ;
      lambda_d = 1./3. - 4./( 9.*M_PI ) ;
      eta_d = 3./32. ;
    }

  // Case 7
  // The edge of the planetary disc is simultaneously
  // just touching the centre of the stellar disc and
  // overflowing the edges of the stellar disc.
  else if ( RpRs > 0.5 - TOL &&
            DOUBLE_EQ( NormSep, RpRs ) )
    {
      _case = 7 ;
      lambda_e = lambda_e_pi_funct( NormSep, RpRs ) ;
      lambda_d = lambda_3( RpRs, k ) ;
      eta_d = eta_1( NormSep, RpRs, a, b ) ;
    }

  // Case 8
  // The edge of the planetary disc is simultaneously
  // overlapping the centre of the stellar disc and 
  // overflowing the edges of the stellar disc.
  else if ( RpRs > 0.5 - TOL &&
            NormSep > fabs( 1. - RpRs ) - TOL &&
            NormSep < RpRs + TOL )
    {
      _case = 8 ;
      lambda_e = lambda_e_pi_funct( NormSep, RpRs ) ;
      lambda_d = lambda_1( NormSep, RpRs, a, b, q, k ) ;
      eta_d = eta_1( NormSep, RpRs, a, b ) ;
    }

  // Case 9
  // The planetary disc is fully overlapping the stellar
  // disc, covering the stellar centre, but the midpoints
  // of the two discs do not coincide.
  else if ( RpRs < 1.0 + TOL &&
            NormSep < 0.5 - fabs( RpRs - 0.5 ) + TOL )
    {
      _case = 9 ;
      lambda_e = POW2( RpRs ) ;
      lambda_d = lambda_2( NormSep, RpRs, a, b, q, k ) ;
      eta_d = eta_2( NormSep, RpRs ) ;
    }

  // Case 10
  // The centres of the stellar and planetary discs
  // precisely coincide, and the planetary disc is 
  // smaller than the stellar disc.
  else if ( RpRs < 1.0 + TOL &&
            DOUBLE_EQ( NormSep, 0. ) )
    {
      _case = 10 ;
      lambda_e = POW2( RpRs ) ;
      lambda_d = lambda_6( RpRs ) ;
      eta_d = eta_2( NormSep, RpRs ) ;
    }

  // Confirm a match has been made to one of the valid cases:
  if ( _case == -1 ) 
    { 
      printf( "ERROR: Unable to match RpRs and NormSep values to star-planet configuration case\n" ) ;
      printf( "(probably because it has been inadvertently overlooked in the C function definition)" ) ;
      exit( -1 ) ;
    }

  // The Theta( RpRs - NormSep ) function in the Mandel & Agol (2002) paper
  // is a heaviside step function, which we can implement here as:
  if( NormSep < RpRs ) lambda_d += 2./3. ;
  
  // Calculate the flux now that the lambda and eta functions are evaluated:
  F = 1. - ( 1./( 1. - gam1/3. - gam2/6. ) ) * ( ( 1. - c2 ) * ( lambda_e )
      + c2 * ( lambda_d ) - c4 * ( eta_d ) ) ;
  
  return F ;

}


