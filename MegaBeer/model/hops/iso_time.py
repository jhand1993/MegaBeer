"""
Mathematics for hop-related calculations for recipes/plotting.
Abbreviations:
    alpha acid: AA
    specific gravity: G
    original gravity: OG
    temperature: T (celsius)
    time: t (minutes)
"""
import numpy as np
from MegaBeer.science import reaction

def maloshell(t, k1, k2):
    """ AA isomizeration rates using the Malowicki & Shellhammer 2005
        model for a fixed gravity. pH fixed to 5.2.

        For a fixed temperature T, the total isomizeration rate is a function of
        time: c_iso(t) = const * f(t)
    Args:
        t (float or numpy.ndarray): Boil time
        k1 (float): Isomerization reaction rate
        k2 (float): Iso-AA degradation rate
    
    Returns:
        float or numpy.ndarray: time component of utilization fraction
    """
    return k1 / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t))


def maloshell_std(t):
    """ Malowicki & Shellhammer 2005 model but with fixed boiling temperature.
        Refer to table 3 in paper: k1 = 0.01141, k2 = 0.00263 for T=100.
    Args:
        t (float or numpy.ndarray): Boil time
    
    Returns:
        float or numpy.ndarray: time component of utilization fraction
    """
    return maloshell(t, 0.01141, 0.00263)

    
def tinseth(t, max_u=0.241, r=0.04):
    """ Temporal component of Tinseth model. Max_u and r taken from Palmer.
        
        Note: 0.241 = 1 / 4.15.
    Args:
        t (float or numpy.ndarray): Boil time.
        max_u (float): Maximum utilization constant.  Default is 0.241.
        r (float): Rate constant of growth.  Default is 0.04.

    Results:
        float or numpy.ndarray: Time component of utilization fraction.
    """
    return max_u * (1. - np.exp(-r * t)) 

