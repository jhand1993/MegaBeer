""" Gravity factor adjustment for total hop utilization.
"""

import numpy as np

def tinseth(G):
    """ Gravity factor component of Tinseth model (1997).
    Args:
        G (float or numpy.ndarray): Standard gravity

    Returns:
        float or numpy.ndarray: Gravity factor.
    """
    return 1.65 * 0.000125**(G - 1.)


def tinseth1050(G):
    """ Same as Tinseth 1997, but normalized to 1 at G=1.050 fo conform 
        with Rager.
    Args:
        G (float or numpy.ndarray): Standard gravity

    Returns:
        float or numpy.ndarray: Gravity factor.
    """
    return 1.5673 * 0.000125 * (G - 1.)

def rager(G):
    """ Gravity factor from Jackie Rager (1990).  Always equal 1.0 if G < 1.050.
    Args:
        G (float or numpy.ndarray): Standard gravity

    Returns:
        float or numpy.ndarray: Gravity factor.
    """
    u = 1. / (1. + 5. * (G - 1.05))

    # Handle array differently
    if type(u) == np.ndarray:
        return np.where(G >= 1.05, u, 1.)
    
    # Return 1.0 for all G < 1.05
    else:
        if G < 1.05:
            return 1.
        
        else:
            return u


def mosher(G):
    """ An order-two polynomial fit result for Randy Mosher's boil gravity factor table 
        fit by Michael L. Hall:
        https://www.homebrewersassociation.org/attachments/0000/2501/IBUs.pdf
    
    Args:
        G (float or numpy.ndarray): Standard gravity

    Returns:
        float or numpy.ndarray: Gravity factor.  
    """
    return 1.0526 * (G - 40. * (G - 1.)**2)