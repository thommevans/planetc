/*********************************************************************************
This function uses the analytic equations of Mandel and Agol (2002) to compute the
flux of a star during occultation given a normalized separation NormSep, ratio of planet
to star radii RpRs, and the four nonlinear limb darkening coefficients (c1,c2,c3,c4).
*********************************************************************************/

#include "ma02_src.h"


double F_nonlin( double NormSep, double RpRs, double c1, double c2, double c3, double c4 )
{
	
  //define commonly used parameters
  double F = 0. ;
  double sum = 0. ;
  double L ;
  int n ;
  double c0 = 1. - c1 - c2 - c3 - c4 ;
  double C[] = { c0, c1, c2, c3, c4 } ;

  double a, b ;
  a = POW2( NormSep - RpRs ) ;
  b = POW2( NormSep + RpRs ) ;

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
      for( n=0, sum=0. ; n<5 ; n++ )
	{
	  sum += Nfunc( NormSep, RpRs, n, a, b )*C[n]/( n + 4. ) ;
	}
      F = 1. - pow( 2*M_PI*Omega( C ), -1. )*sum ;
    }
  
  // Case 3
  // The planetary disc is fully overlapping the stellar
  // disc, but the stellar centre is not covered.
  else if ( RpRs < 0.5 + TOL &&
            NormSep > RpRs - TOL &&
            NormSep < 1. - RpRs + TOL )
    {
      _case = 3 ;
      L = POW2( RpRs )*( 1. - POW2( RpRs )/2. - POW2( NormSep ) ) ;
      for( n=1, sum=0. ; n<4 ; n++ )
	{
	  sum += Mfunc( NormSep, RpRs, n, a, b ) * C[n] / ( n + 4. ) ;
	}
      F = 1. - pow( 4.*Omega( C ), -1. )*( C[0]*POW2( RpRs ) + 2.*sum + C[4]*L ) ;
    }

  // Case 4
  // The same as Case 3 with the added condition that the
  // edges of the planetary disc and stellar disc just touch.
  else if ( RpRs < 0.5 + TOL &&
            DOUBLE_EQ( NormSep, 1. - RpRs ) )
    {
      _case = 4 ;
      L = POW2( RpRs )*( 1. - POW2( RpRs )/2. - POW2( NormSep ) ) ;
      for( n=1, sum=0. ; n<4 ; n++ )
	{
	  sum += Mfunc( NormSep, RpRs, n, a, b ) * C[n] / ( n + 4. ) ;
	}
      F = 1. - pow( 4.*Omega( C ), -1. )*( ( C[0]*POW2( RpRs ) ) + 2.*sum + C[4] * L ) ;
    }

  // Case 5
  // Very similar to Case 3 except that the edge of the 
  // planetary disc just touches the stellar centre.
  else if ( RpRs < 0.5 + TOL &&
            DOUBLE_EQ( NormSep, RpRs ) )
    {
      _case = 5 ;
      for( n=0, sum=0. ; n<5 ; n++ )
	{
	  sum += C[n]/( n + 4. )*hyp2f1( 0.5, -( n + 4. )/4., 1., 4.*POW2( RpRs ) ) ;
	}
      F = 0.5 + pow( 2.*Omega( C ), -1. )*sum ;
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
      for( n=0, sum=0. ; n<5 ; n++ )
	{
	  sum += C[n]/( n + 4. )*gsl_sf_gamma( 1.5 + n/4. )/gsl_sf_gamma( 2. + n/4. ) ;
	}
      F = 0.5 + 1./( 2.*M_PI*Omega( C ) )*sum ;
    }

  // Case 7
  // The edge of the planetary disc is simultaneously
  // just touching the centre of the stellar disc and
  // overflowing the edges of the stellar disc.
  else if ( RpRs > 0.5 - TOL &&
            DOUBLE_EQ( NormSep, RpRs ) )
    {
      _case = 7 ;
      for( n=0, sum=0. ; n<5 ; n++ )
	{
	  //sum += C[n]/( n + 4. )*gsl_sf_beta( 0.5, ( n + 8. )/4. )*hyp2f1( 0.5, 0.5, 2.5 + n/4., 1./( 4.*POW2( RpRs ) ) ) ;
	  sum += C[n]/( n + 4. )*gsl_sf_beta( 0.5, ( n + 8. )/4. ) ;
	  sum *= hyp2f1( 0.5, 0.5, 2.5 + n/4., 1./( 4.*POW2( RpRs ) ) ) ;
	}
    
      F = 0.5 + 1./( 4.*RpRs*M_PI*Omega( C ) )*sum ;
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
      for( n=0, sum=0. ; n<5 ; n++ )
	{
	  sum += C[n]*Nfunc( NormSep, RpRs, n, a, b ) * pow( n+4., -1. ) ;
	}
      F = -1./( 2.*M_PI*Omega( C ) )*sum ;
    }

  // Case 9
  // The planetary disc is fully overlapping the stellar
  // disc, covering the stellar centre, but the midpoints
  // of the two discs do not coincide.
  else if ( RpRs < 1.0 + TOL &&
            NormSep < 0.5 - fabs( RpRs - 0.5 ) + TOL )
    {
      _case = 9 ;
      L = POW2( RpRs )*( 1. - POW2( RpRs )/2. - POW2( NormSep ) ) ;
      for( n=1, sum=0. ; n<4 ; n++ )
	{
	  sum += ( C[n]/( n + 4. )*Mfunc( NormSep, RpRs, n, a, b ) ) ;
	}
      F = pow( ( 4.*Omega( C ) ), -1. )*( C[0]*( 1. - POW2( RpRs ) ) + C[4]*( 0.5 - L ) - 2.*sum ) ;
    }

  // Case 10
  // The centres of the stellar and planetary discs
  // precisely coincide, and the planetary disc is 
  // smaller than the stellar disc.
  else if ( RpRs < 1.0 + TOL &&
            DOUBLE_EQ( NormSep, 0. ) )
    {
      _case = 10 ;
      for( n=0, sum=0. ; n<5 ; n++ )
	{
	  sum += C[n]*pow( 1. - POW2( RpRs ), ( n + 4. )/4. )/( n + 4. ) ;
	}
      F = pow( Omega( C ), -1. )*sum ;
    }

  // Confirm a match has been made to one of the valid cases:
  if ( _case == -1 ) 
    { 
      printf( "ERROR: Unable to match RpRs and NormSep values to star-planet configuration case\n" ) ;
      printf( "(probably because it has been inadvertently overlooked in the C function definition)" ) ;
      exit( -1 ) ;
    }
  
  return F ;  
}






















/*   //////////////////////neale's	 */

/*   //Cases 4 and 5 upper and lower limits of 3, respectively */
/*   //Case 4 - planet entirely within star, but does not touch/cross the centre, and touches limb */
/*   //same as Case 3 for nonlin law */
/*   else if( p < 0.5 && DOUBEQ(NormSep,(1.-RpRs)) ) */
/*     { */
/*       Case = 4; */
    
/*       //calculate sum term */
/*       for(n=1,sum=0.;n<4;n++) */
/* 	{ */
/* 	  sum += Mfunc(NormSep, RpRs, n, POW2( NormSep - RpRs ), POW2( NormSep + RpRs ) ) * C[n] / (n+4.); */
/* 	} */
    
/*       L = POW2(RpRs) * (1.-POW2(RpRs)/2.-POW2(NormSep)); */
    
/*       F = 1. - pow(4.*Omega(C),-1.) * ( (C[0] * POW2(RpRs)) + 2.*sum + C[4] * L ); */
/*     } */
  
/*   //Case 5 - planet entirely within star, touches centre of star */
/*   else if( p < 0.5 && DOUBEQ(NormSep,RpRs)) */
/*     { */
/*       Case = 5; */
      
/*       for(n=0,sum=0.;n<5;n++) */
/* 	{ */
/* 	  sum += C[n] / (n+4.) * hyp2f1(0.5,-(n+4.)/4.,1.,4.*POW2(RpRs)); */
/* 	} */
      
/*       F = 0.5 + pow(2.*Omega(C),-1.) * sum; */
/*     } */
  
/*   //Case 6 - planet touches star centre and disk, therefore p = 0.5 */
/*   else if( DOUBEQ(RpRs,0.5) && DOUBEQ(NormSep,0.5) ) */
/*     { */
/*       Case = 6; */
    
/*       for(n=0,sum=0.;n<5;n++) */
/* 	{ */
/* 	  sum += C[n]/(n+4.) * gsl_sf_gamma(1.5+n/4.) / gsl_sf_gamma(2.+n/4.); */
/* 	} */
    
/*       F = 0.5 + 1./(2.*M_PI*Omega(C)) * sum; */
    
/*       flux = 0.; */
/*     } */
  
/*   //Case 7 - planet touches stellar centre, but not entirely within disk */
/*   else if( p > 0.5 && DOUBEQ(NormSep,RpRs) ) */
/*     { */
/*       Case = 7; */
    
/*       for(n=0,sum=0.;n<5;n++) */
/* 	{ */
/* 	  sum += C[n]/(n+4.) * gsl_sf_beta(0.5,(n+8.)/4.) * hyp2f1(0.5,0.5,2.5+n/4.,1./(4.*POW2(RpRs))); */
/* 	} */
    
/*       F = 0.5 + 1./(4.*p*M_PI*Omega(C)) * sum; */
    
/*       flux = 0.; */
/*     } */
  
/*   //Case 2 - planet on limb of star but doesn't touch centre */
/*   else if( NormSep > (0.5+fabs(RpRs-0.5)) && NormSep < (1.+RpRs)) */
/*     { */
/*       Case = 2; */
    
/*       for(n=0,sum=0.;n<5;n++) */
/* 	{ */
/* 	  sum += Nfunc( NormSep, RpRs, n, POW2( NormSep - RpRs ), POW2( NormSep + RpRs ) ) * C[n] / (n+4.); */
/* 	} */
    
/*       F = 1. - pow(2*M_PI*Omega(C),-1.) * sum; */
/*     } */
  
/*   //Case 3 - planet entirely within star, but does not touch/cross the centre, or touch limb */
/*   else if( p < 0.5 && NormSep > p && NormSep < (1.-RpRs)) */
/*     { */
/*       Case = 3; */
      
/*       for(n=1,sum=0.;n<4;n++) */
/* 	{ */
/* 	  sum += Mfunc( NormSep, RpRs, n, POW2( NormSep - RpRs ), POW2( NormSep + RpRs ) ) * C[n] / (n+4.); */
/* 	} */
    
/*       double L = POW2(RpRs) * (1.-POW2(RpRs)/2.-POW2(NormSep)); */
    
/*       F = 1. - pow(4.*Omega(C),-1.) * ( C[0]*POW2(RpRs) + 2.*sum + C[4] * L ); */
/*     } */
  
/*   //Case 8 - planet covers the centre and limb of the star */
/*   else if( p > 0.5 && NormSep >= fabs(1.-RpRs) && NormSep < RpRs) */
/*     { */
/*       Case = 8; */
    
/*       for(n=0,sum=0.;n<5;n++) */
/* 	{ */
/* 	  sum += C[n] * Nfunc( NormSep, RpRs, n, POW2( NormSep - RpRs ), POW2( NormSep + RpRs ) ) * pow(n+4.,-1.); */
/* 	} */
    
/*       F = -1./(2.*M_PI*Omega(C)) * sum; */
/*     } */
  
/*   //Case 9 - planet entirely withing star, and covers the centre */
/*   //if( (RpRs < 1.) && (NormSep > 0.) && (NormSep < (0.5 - abs(RpRs-0.5)))){ */
/*   //had to change this from MA02 paper - it allows situations with p<0.5 and the centre is not covered? */
/*   else if( (RpRs < 1.) && (NormSep > 0.) && (NormSep < RpRs) ) */
/*     { */
/*       Case = 9; */
      
/*       //define L */
/*       double L = POW2(RpRs) * (1.-POW2(RpRs)/2.-POW2(NormSep)); */
    
/*       for(n=1,sum=0.;n<4;n++) */
/* 	{ */
/* 	  sum += (C[n] / (n+4.) * Mfunc( NormSep, RpRs, n, POW2( NormSep - RpRs ), POW2( NormSep + RpRs ) ) ) ; */
/* 	} */
    
/*       F = pow((4.*Omega(C)),-1.) * ( C[0]*(1.-POW2(RpRs)) + C[4]*(0.5-L) -2.*sum ); */
/*     } */
  
/*   //Case 10 - planet centred on star, p<1. */
/*   else if( p < 1. && NormSep == 0.) */
/*     { */
/*       Case = 10; */
    
/*       for(n=0,sum=0.;n<5;n++) */
/* 	{ */
/* 	  sum += C[n] * pow(1.-POW2(RpRs),(n+4.)/4.) / (n+4.); */
/* 	} */
    
/*       F = pow(Omega(C),-1.) * sum; */
/*     } */
  
/*   // Confirm a match has been made to one of the valid cases: */
/*   if ( _case == -1 )  */
/*     {  */
/*       printf( "ERROR: Unable to match RpRs and NormSep values to star-planet configuration case\n" ) ; */
/*       printf( "(probably because it has been inadvertently overlooked in the C function definition)" ) ; */
/*       exit( -1 ) ; */
/*     } */
			
/*   return F; */
/* } */
