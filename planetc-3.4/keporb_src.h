#include <math.h>
#include <stdio.h>

double Xcoord( double MeanAnom, double aRs, double ecc, double omega );
double Ycoord( double MeanAnom, double aRs, double ecc, double omega, double incl );
double Zcoord( double MeanAnom, double aRs, double ecc, double omega, double incl );
double NormSep( double MeanAnom, double aRs, double ecc, double omega, double incl );
double calc_EccAnom( double MeanAnom, double ecc );
double calc_TrueAnom( double EccAnom, double ecc );
double calc_cos_TrueAnom( double EccAnom, double ecc );
double calc_sin_TrueAnom( double EccAnom, double ecc );

#define POW2( x ) ( ( x ) * ( x ) )
