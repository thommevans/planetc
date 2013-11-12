

def radvel(JD, P, K, T0 = 0, V0 = 0, Ecc = 0, omega = 0):
    """Radial velocity (user-defined semi-amplitude)"""
    Nu = truean(JD, P, T0, Ecc)
    Vr = V0 + K * (scipy.cos(omega + Nu) + Ecc * scipy.cos(omega))
    if (K < 0): Vr[:] = -999
    return Vr

def rv_K(P, f):
    """Calculate RV semi-amplitude for a given mass function f."""
    return (2 * scipy.pi * 6.67e-11 / (P*86400)**2)**(1/3.) * f
