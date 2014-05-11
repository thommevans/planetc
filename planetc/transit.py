import numpy as np
import matplotlib.pyplot as plt
import pdb, sys, os
import keporb, ma02
import phys_consts as consts

# 15Jul2013 TME:
#    Actually, I think I worked out that the 'bug' mentioned
#    in the previous comment was an artefact of the equations,
#    rather than a bug in the code. It might be to do with the
#    fact that the equations implemented in the code are
#    themselves an approximation that is very good for moderate
#    eccentricities, but starts breaking down for higher
#    eccentricity orbits.
# 3Mar2013 TME:
#    There appears to be a bug in the code for calculating
#    eccentric transits that manifests itself as an egress
#    (haven't seen an ingress case yet) that isn't perfectly
#    smooth - it looks basically right but it's definitely
#    a bit funny. Need to fix this bug.
#
# Nov2012 TME:
#    This module provides python wrappers for underlying
#    C routines that compute Mandel & Agol (2002) transit
#    lightcurves.

def ma02_aRs( t, **pars ):
    """

    Uses the Mandel & Agol (2002) analytic equations to compute
    transit lightcurves. This routine takes the semimajor axis
    directly as an input parameter, and therefore does not require
    either the star or planet mass; this differs from the
    ma02_RsMsMp() routine which uses the star and planet masses
    with Kepler's third law to calculate the semimajor axis.

    CALLING
      F = transit.ma02_aRs( t, **pars )

    INPUTS
    ** t - a numpy array containing the times that the lightcurve
    is to be evaluated for.

    KEYWORD INPUTS
    ** pars - a dictionary containing the keyword arguments.
         ['tr_type'] a flag that can be set to either 'primary',
           'secondary' or 'both'; if it is set to 'both', the
           Z-coordinate will be calculated along with the normalised
           separation in order to distinguish between primary transits
           and secondary eclipses, and to scale the alternating flux
           changes accordingly; on the other hand, if it is set to
           either 'primary' or 'secondary', the Z-coordinate will not
           be calculated and all transit events will have the same
           shape and depth depending on which is chosen; the latter
           saves time and should almost always be used when modelling
           on a per-transit basis. 
         ['T0'] time of periapse passage, which this routine will
           force to be the same as the mid-transit time if the orbit
           is circular, by setting the argument of periapse to 3pi/2.
         ['P'] orbital period in same units of time as 'T0'.
         ['aRs'] semimajor axis in units of stellar radii.
         ['RpRs'] planetary radius in units of stellar radii.
         ['incl'] orbital inclination in degrees.
         ['b'] = aRs*cos(i), which can be provided instead of 'incl'
         ['ecc'] orbital eccentricity.
         ['SecDepth'] depth of the secondary eclipse; will be set to
           zero if not explicitly specified.
         ['omega'] argument of periapse in degrees; this must be
           explicitly specified if the eccentricity is nonzero; if
           the orbit is circular, it will be forced to 3pi/2 regardless
           of the input value, to ensure that T0 corresponds to the
           time of mid-transit.
         ['ld'] an optional flag that can be set to either None,
           'quad' or 'nonlin' to specify the type of limb
           darkening law to be used; if it is not set, then no limb
           darkening will be used, and this would also be the case if
           it is set to None.
         ['gam1']+['gam2'] quadratic limb darkening coeffs; required
          if the 'ld' flag is set to 'quad'.
         ['c1']+['c2']+['c3']+['c4'] nonlinear limb darkening coeffs;
           required if the 'ld' flag is set to 'nonlin'.

    OUTPUT
    ** F - a numpy array containing the relative flux values of
    the model transit lightcurve.

    """ 

    # Start unpacking paramteres:
    T0 = pars[ 'T0' ]
    P = pars[ 'P' ]
    aRs = pars[ 'aRs' ]
    RpRs = pars[ 'RpRs' ]
    try:
        incl_rad = np.deg2rad( pars['incl'] )
    except:
        incl_rad = np.arccos( pars['b']/aRs )
        
    ecc = pars[ 'ecc' ]
    try:
        foot = pars[ 'foot' ]
    except:
        foot = 1.
    try:
        grad = pars[ 'grad' ]
    except:
        grad = 0.
    try:
        SecDepth = pars[ 'SecDepth' ]
    except:
        SecDepth = 0.
    try:
        tr_type = pars['tr_type']
    except:
        tr_type = 'both'

    # Following the naming convention of the original
    # Mandel & Agol (2002) paper for the limb darkening
    # coefficients:
    try:
        if pars['ld']=='quad':
            gam1 = pars[ 'gam1' ]
            gam2 = pars[ 'gam2' ]
        elif pars['ld']=='nonlin': 
            c1 = pars[ 'c1' ]
            c2 = pars[ 'c2' ]
            c3 = pars[ 'c3' ]
            c4 = pars[ 'c4' ]
        elif pars['ld']==None:
            pars['ld'] = 'quad'
            pars['gam1'] = 0.
            pars['gam2'] = 0.
    except:
        pars['ld'] = 'quad'
        pars['gam1'] = 0.
        pars['gam2'] = 0.
       
    # Calculate the mean anomaly:
    t = t.flatten()
    MeanAnom = ( 2*np.pi/P )*( t - T0 )

    # Calculate the normalised separation between the
    # planetary and stellar discs:
    if ecc != 0.:
        omega_rad = pars['omega'] * np.pi/180.
        NormSep = keporb.NormSep( MeanAnom, aRs, ecc, omega_rad, incl_rad )
    else:
        omega_rad = 3.*np.pi/2.
        try:
            b = pars['b']
        except:
            b = aRs*np.cos( incl_rad )
        NormSep = np.sqrt( ( ( aRs*np.sin( MeanAnom ) )**2. ) \
                           + ( ( b*np.cos( MeanAnom ) )**2. ) )

    # If we want to model both primary transits and secondary
    # eclipses, we need to compute the Z coordinate to determine
    # when the planet is in front of the star (Z<0) and behind
    # the star (Z>0):
    if tr_type=='both':
        F = np.ones( len( t ) )
        zcoord = keporb.Zcoord( MeanAnom, aRs, ecc, omega_rad, incl_rad )
        ixsf = ( zcoord < 0 )
        if ixsf.max()==True:
            if pars['ld']=='quad':
                F[ixsf] = ma02.F_quad( NormSep[ixsf], RpRs, \
                                       pars['gam1'], pars['gam2'] )
            elif pars['ld']=='nonlin':
                F[ixsf] = ma02.F_nonlin( NormSep[ixsf], RpRs, \
                                         pars['c1'], pars['c2'], \
                                         pars['c3'], pars['c4'] )
            else:
                pdb.set_trace()
        ixsb = ( zcoord >= 0 )
        if ixsb.max()==True:
            temp = ma02.F_quad( NormSep[ixsb], RpRs, 0.0, 0.0 ) - 1.
            F[ixsb] = 1 + temp*SecDepth/( temp.max() - temp.min() ) #/( RpRs**2. )

    # If we're only interested in the primary transits then
    # we must take stellar limb darkening into account
    # while treating the planet as an non-luminous disc:
    elif tr_type=='primary':
        if pars['ld']=='quad':
            F = ma02.F_quad( NormSep, RpRs, \
                             pars['gam1'], pars['gam2'] )
        elif pars['ld']=='nonlin':
            F = ma02.F_nonlin( NormSep, RpRs, \
                               pars['c1'], pars['c2'], \
                               pars['c3'], pars['c4'] )
        else:
            print '\n\n\n{0:s} not recognised as limb darkening type'\
                  .format( pars['ld'] )
            pdb.set_trace() 

    # If we're only interested in the secondary eclipses
    # we treat the planet as a uniform disc with no limb
    # darkening:
    elif tr_type=='secondary':
        temp = ma02.F_quad( NormSep, RpRs, 0.0, 0.0 ) - 1.
        F = 1 + temp*SecDepth/( temp.max() - temp.min() ) #/( RpRs**2. )

    # If requested, re-scale the lightcurve by a linear
    # trend before returning the output:
    if ( grad!=0. )+( foot!=1. ):
        twid = t.max() - t.min()
        tmid = t.min() + 0.5*twid
        F = F * ( foot + grad*( t - tmid ) )
    # NOTE: This will change the absolute value of the
    # eclipse depth, but the fractional value of the
    # eclipse depth will remain the same.

    return F




def ma02_RsMsRpMp( t, **pars ):
    """
    
    Uses the Mandel & Agol (2002) analytic equations to compute
    transit lightcurves. This routine takes the mass and radius
    for both the star and planet as input parameters. The masses
    are used with Kepler's third law to calculate the semimajor
    axis. This differs from the ma02_aRs() routine which takes
    the semimajor axis directly as an input parameter.

    CALLING
      F = transit.ma02_RsMsMp( t, **pars )

    INPUTS
    ** t - a numpy array containing the times that the lightcurve
    is to be evaluated for.

    KEYWORD INPUTS
    ** pars - a dictionary containing the keyword arguments:
         ['tr_type'] a flag that can be set to either 'primary',
           'secondary' or 'both'; if it is set to 'both', the
           Z-coordinate will be calculated along with the normalised
           separation in order to distinguish between primary transits
           and secondary eclipses, and to scale the alternating flux
           changes accordingly; on the other hand, if it is set to
           either 'primary' or 'secondary', the Z-coordinate will not
           be calculated and all transit events will have the same
           shape and depth depending on which is chosen; the latter
           saves time and should almost always be used when modelling
           on a per-transit basis. 
         ['T0'] time of periapse passage, which this routine will
           force to be the same as the mid-transit time if the orbit
           is circular, by setting the argument of periapse to 3pi/2.
         ['P'] orbital period in same units of time as 'T0'.
         ['Rs'] stellar radius in solar radii.
         ['Ms'] stellar mass in solar masses.
         ['Rp'] planetary radius in Jupiter radii.
         ['Mp'] planetary mass in Jupiter masses.
         ['incl'] orbital inclination in degrees.
         ['ecc'] orbital eccentricity.
         ['SecDepth'] depth of the secondary eclipse; will be set to
           zero if not explicitly specified.
         ['omega'] argument of periapse in degrees; this must be
           explicitly specified if the eccentricity is nonzero; if
           the orbit is circular, it will be forced to 3pi/2 regardless
           of the input value, to ensure that T0 corresponds to the
           time of mid-transit.
         ['ld'] an optional flag that can be set to either None,
           'quad' or 'nonlin' to specify the type of limb
           darkening law to be used; if it is not set, then no limb
           darkening will be used, and this would also be the case if
           it is set to None.
         ['gam1']+['gam2'] quadratic limb darkening coeffs; required
          if the 'ld' flag is set to 'quad'.
         ['c1']+['c2']+['c3']+['c4'] nonlinear limb darkening coeffs;
           required if the 'ld' flag is set to 'nonlin'.

    OUTPUT
    ** F - a numpy array containing the relative flux values of
    the model transit lightcurve.

    """

    # Start unpacking paramteres:
    T0 = pars['T0']
    P = pars['P']
    Rs = pars['Rs']
    Ms = pars['Ms']
    Rp = pars['Rp']
    Mp = pars['Mp']
    SecDepth = pars['SecDepth']
    incl_rad = pars['incl'] * np.pi/180.
    ecc = pars['ecc']
    try:
        foot = pars[ 'foot' ]
    except:
        foot = 1.
    try:
        grad = pars[ 'grad' ]
    except:
        grad = 0.
    try:
        SecDepth = pars[ 'SecDepth' ]
    except:
        SecDepth = 0.
    try:
        tr_type = pars['tr_type']
    except:
        tr_type = 'both'

    # Convert some of the units:
    Rs *= consts.RSun
    Ms *= consts.MSun
    Rp *= consts.RJup
    Mp *= consts.MJup
    
    # Assuming a 2-body Keplerian orbit, use Kepler's
    # third law to calculate the semimajor axis:
    a = np.power( ( ( ( ( P*24.*60.*60./( 2*np.pi ) )**2 ) \
                  * consts.G * ( Ms + Mp) ) ) , (1./3.) )
    aRs = a/Rs
    RpRs = Rp/Rs
    
    # Following the naming convention of the original
    # Mandel & Agol (2002) paper for the limb darkening
    # coefficients:
    try:
        if pars['ld']=='quad':
            gam1 = pars[ 'gam1' ]
            gam2 = pars[ 'gam2' ]
        elif pars['ld']=='nonlin': 
            c1 = pars[ 'c1' ]
            c2 = pars[ 'c2' ]
            c3 = pars[ 'c3' ]
            c4 = pars[ 'c4' ]
        elif pars['ld']==None:
            pars['ld'] = 'quad'
            pars['gam1'] = 0.
            pars['gam2'] = 0.
    except:
        pars['ld'] = 'quad'
        pars['gam1'] = 0.
        pars['gam2'] = 0.
        
    # Calculate the mean anomaly:
    MeanAnom = ( 2*np.pi/P )*( t - T0 )

    # Calculate the normalised separation between the
    # planetary and stellar discs:
    if ecc != 0.:
        omega_rad = pars['omega'] * np.pi/180.
        NormSep = keporb.NormSep( MeanAnom, aRs, ecc, omega_rad, incl_rad )
    else:
        omega_rad = 3.*np.pi/2.
        b = aRs*np.cos( incl_rad )
        NormSep = np.sqrt( ( ( aRs*np.sin( MeanAnom ) )**2. ) \
                           + ( ( b*np.cos( MeanAnom ) )**2. ) )

    # If we want to model both primary transits and secondary
    # eclipses, we need to compute the Z coordinate to determine
    # when the planet is in front of the star (Z<0) and behind
    # the star (Z>0):
    if tr_type=='both':
        F = np.ones( len( t ) )
        zcoord = keporb.Zcoord( MeanAnom, aRs, ecc, omega_rad, incl_rad )
        ixsf = ( zcoord < 0 )
        if pars['ld']=='quad':
            F[ixsf] = ma02.F_quad( NormSep[ixsf], RpRs, \
                                   pars['gam1'], pars['gam2'] )
        elif pars['ld']=='nonlin':
            F[ixsf] = ma02.F_nonlin( NormSep[ixsf], RpRs, \
                                     pars['c1'], pars['c2'], \
                                     pars['c3'], pars['c4'] )
        else:
            pdb.set_trace()
        ixsb = ( zcoord >= 0 )
        temp = ma02.F_quad( NormSep[ixsb], RpRs, 0.0, 0.0 ) - 1.
        F[ixsb] = 1 + temp*SecDepth/( temp.max() - temp.min() ) #/( RpRs**2. )

    # If we're only interested in the primary transits then
    # we must take stellar limb darkening into account
    # while treating the planet as an non-luminous disc:
    elif tr_type=='primary':
        if pars['ld']=='quad':
            F = ma02.F_quad( NormSep, RpRs, \
                             pars['gam1'], pars['gam2'] )
        elif pars['ld']=='nonlin':
            F = ma02.F_nonlin( NormSep, RpRs, \
                               pars['c1'], pars['c2'], \
                               pars['c3'], pars['c4'] )
        else:
            pdb.set_trace()

    # If we're only interested in the secondary eclipses
    # we treat the planet as a uniform disc with no limb
    # darkening:
    elif tr_type=='secondary':
        temp = ma02.F_quad( NormSep, RpRs, 0.0, 0.0 ) - 1.
        F = 1 + temp*SecDepth/( temp.max() - temp.min() ) #/( RpRs**2. )

    # If requested, re-scale the lightcurve by a linear
    # trend before returning the output:
    if ( grad!=0. )+( foot!=1. ):
        twid = t.max() - t.min()
        tmid = t.min() + 0.5*twid
        F = F * ( foot + grad*( t - tmid ) )

    return F


def calc_T0( Ttr, P, ecc, omega, transit='primary' ):
    """
    SUMMARY
    Computes time of periapse passage given the mid-time of
    either the primary transit or the secondary eclipse.

    CALLING
      T0 = transit.calc_T0( Ttr, P, ecc, omega, transit='primary' )
    
    INPUTS
    ** Ttr - transit/eclipse mid-time.
    ** P - orbital period.
    ** ecc - orbital eccentricity.
    ** omega - argument of periapse in degrees.

    KEYWORD INPUTS
    ** transit - 'primary' or 'secondary' depending on whether the
      mid-time Ttr is for the primary transit or secondary eclipse;
      default is 'primary'.

    OUTPUT
    ** T0 - the time of periapse passage.

    NOTES:
    ** Ttr and P must have the same units of time, and the output T0
       will also have the same units.
    ** By default, the longitude of the ascending node (big Omega)
       is implicitly taken to be 180 degrees in the calculations.
    """
    
    # Convert omega from degrees to radians:
    omega *= np.pi / 180.
    
    # Make sure omega has been provided between 0-2pi:
    while omega >= 2*np.pi:
        omega -= 2*pi
    while omega < 0:
        omega += 2*np.pi

    # Calculate the true anomaly corresponding to the
    # midpoint of the transit/eclipse, and ensure that
    # the value lies in the 0-2pi range:
    if transit=='primary':
        TrueAnom_tr = 3*np.pi/2. - omega
    elif transit=='secondary':
        TrueAnom_tr = np.pi/2. - omega
    else:
        pdb.set_trace()
    while TrueAnom_tr >= 2*np.pi:
        TrueAnom_tr -= 2*pi
    while TrueAnom_tr < 0:
        TrueAnom_tr += 2*np.pi

    # Calculate the value of the eccentric anomaly at the time
    # of transit/eclipse:
    EccAnom_tr = 2 * np.arctan2( np.sqrt( 1. - ecc )*np.sin( TrueAnom_tr/2. ), \
                                 np.sqrt( 1. + ecc )*np.cos( TrueAnom_tr/2. ) )

    # Use the eccentric anomaly at the time of transit/eclipse
    # to calculate the mean anomaly:
    MeanAnom_tr = EccAnom_tr - ecc*np.sin( EccAnom_tr )

    # Convert the mean anomaly to a time elapsed since the
    # time of periastron passage:
    delt = P * ( MeanAnom_tr/2./np.pi )

    # Calculate the time of periastron passage:
    T0 = Ttr - delt

    return T0


def example():
    """
    Simple routine that demonstrates various ways of
    computing lightcurves, and plots the results.
    """
    
    # Time values:
    t = np.linspace( 0., 7., 10000 )

    # Orbit properties:
    ecc = 0.0
    P = 2.1 # days
    incl = 88.2 # degrees
    omega = 90.0 # degrees
    T_tr = 0.0 # time of transit
    T0 = calc_T0( T_tr, P, ecc, omega, transit='primary' )

    # Star-planet physical:
    Rp = 1.0 # Jupiter radii
    Mp = 1.0 # Jupiter masses
    Rs = 1.0 # solar radii
    Ms = 1.0 # solar masses

    # Lightcurve properties:
    SecDepth = 1e-3
    foot = 1.0
    grad = 0.0

    # Calculate the semimajor axis using Kepler's
    # third law:
    a = np.power( ( ( ( ( P*24.*60.*60./( 2*np.pi ) )**2 ) * consts.G \
                  * ( Ms*consts.MSun + Mp*consts.MJup ) ) ) , (1./3.) )

    RpRs = ( Rp*consts.RJup )/( Rs*consts.RSun )
    aRs = a / ( Rs*consts.RSun )    

    # Nonlinear limb darkening coeffs:
    c1 = -0.1
    c2 = +1.4
    c3 = -1.2
    c4 = +0.5

    # Quadratic limb darkening coeffs:
    gam1 = 0.5
    gam2 = 0.1

    # Quadratic limb darkening + RsMsRpMp parameterisation:
    pars_RsMsRpMp_q = { 'T0':T0, 'P':P, 'Ms':Ms, 'Mp':Mp, 'Rs':Rs, \
                      'Rp':Rp, 'SecDepth':SecDepth, 'incl':incl, \
                      'ecc':ecc, 'omega':omega, 'gam1':gam1, 'gam2':gam2, \
                      'ld':'quad', 'foot':foot, 'grad':grad }
    F_RsMsRpMp_q = ma02_RsMsRpMp( t, **pars_RsMsRpMp_q )
    
    # Nonlinear limb darkening + RsMsRpMp parameterisation:
    pars_RsMsRpMp_nl = { 'T0':T0, 'P':P, 'Ms':Ms, 'Mp':Mp, 'Rs':Rs, \
                      'Rp':Rp, 'SecDepth':SecDepth, 'incl':incl, \
                      'ecc':ecc, 'omega':omega, 'c1':c1, 'c2':c2, \
                       'c3':c3, 'c4':c4, 'ld':'nonlin', \
                       'foot':foot, 'grad':grad }
    F_RsMsRpMp_nl = ma02_RsMsRpMp( t, **pars_RsMsRpMp_nl )

    # No limb darkening + RsMsRpMp parameterisation:
    pars_RsMsRpMp_n = { 'T0':T0, 'P':P, 'Ms':Ms, 'Mp':Mp, 'Rs':Rs, \
                      'Rp':Rp, 'SecDepth':SecDepth, 'incl':incl, \
                      'ecc':ecc, 'omega':omega, 'ld':None, \
                      'foot':foot, 'grad':grad }
    F_RsMsRpMp_n = ma02_RsMsRpMp( t, **pars_RsMsRpMp_n )

    # Quadratic limb darkening + aRs parameterisation:
    pars_aRs_q = { 'T0':T0, 'P':P, 'aRs':aRs, 'RpRs':RpRs, \
                   'SecDepth':SecDepth, 'incl':incl, 'ecc':ecc, \
                   'omega':omega, 'gam1':gam1, 'gam2':gam2,
                   'ld':'quad', 'foot':foot, 'grad':grad }
    F_aRs_q = ma02_aRs( t, **pars_aRs_q )
    
    # Nonlinear limb darkening + aRs parameterisation:
    pars_aRs_nl = { 'T0':T0, 'P':P, 'aRs':aRs, 'RpRs':RpRs, \
                    'SecDepth':SecDepth, 'incl':incl, 'ecc':ecc, \
                    'omega':omega, 'c1':c1, 'c2':c2, 'c3':c3, 'c4':c4,
                    'ld':'nonlin', 'foot':foot, 'grad':grad }
    F_aRs_nl = ma02_aRs( t, **pars_aRs_nl )

    # No limb darkening + aRs parameterisation:
    pars_aRs_n = { 'T0':T0, 'P':P, 'aRs':aRs, 'RpRs':RpRs, \
                   'SecDepth':SecDepth, 'incl':incl, 'ecc':ecc, \
                   'omega':omega, 'ld':None, \
                   'foot':foot, 'grad':grad }
    F_aRs_n = ma02_aRs( t, **pars_aRs_n )

    # Plot the results:
    fig = plt.figure()
    ax1 = fig.add_subplot( 211 )
    ax1.plot( t, F_aRs_n, '--g', lw=1 )    
    ax1.plot( t, F_aRs_nl, '-m', lw=1 )
    ax1.plot( t, F_aRs_q, '-b', lw=1 )
    #ax1.set_ylim( [ 1-1.4*(RpRs**2.), 1+0.2*(RpRs**2.) ] )
    ax1.set_ylim( [ 1. - 1.4*foot*(RpRs**2.), 1. + 0.2*foot*(RpRs**2.) ] )
    ax2 = fig.add_subplot( 212, sharex=ax1 )
    ax2.plot( t, F_RsMsRpMp_n, '--g', lw=1 )    
    ax2.plot( t, F_RsMsRpMp_nl, '-m', lw=1 )
    ax2.plot( t, F_RsMsRpMp_q, '-b', lw=1 )
    #ax2.set_ylim( [ 1-1.4*(RpRs**2.), 1+0.2*(RpRs**2.) ] )    
    ax2.set_ylim( [ 1. - 1.4*foot*(RpRs**2.), 1. + 0.2*foot*(RpRs**2.) ] )
    ax2.set_xlabel( 'Time' )

    # Discrepencies between output from routines
    # using different parameterisations:
    print 'This should be zero --> {0:.10f}'\
          .format( ( F_RsMsRpMp_q - F_aRs_q ).max() )
    print 'This should be zero --> {0:.10f}'\
          .format( ( F_RsMsRpMp_nl - F_aRs_nl ).max() )

    return None
