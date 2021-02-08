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
from MegaBeer.science import heat
from scipy.interpolate import RectBivariateSpline as rbs


class MaloShell:
    """ Container class for Malowicki & Shellhammer 2005 results.
    """
    @staticmethod
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

    @staticmethod
    def maloshell_std(t):
        """ Malowicki & Shellhammer 2005 model but with fixed boiling temperature.
            Refer to table 3 in paper: k1 = 0.01141, k2 = 0.00263 for T=100.
        Args:
            t (float or numpy.ndarray): Boil time.
        
        Returns:
            float or numpy.ndarray: time component of utilization fraction
        """
        return MaloShell.maloshell(t, 0.01141, 0.00263)

    @staticmethod
    def maloshell_cooling_curves(t, tau=132.5, T0=21.1):
        """ Approximates isomizeration rate during cooling.  Assumes reactions
            cease at T=80.  Default tau is taken from cooling one gallon of water
            in an 8 quart stainless steel stock pot.
        Args:
            t (float): Time left to cool.
            tau (float): Cooling rate time scale in minutes.  Default is 132.5 min.
            T0 (float): Room temperature water is cooling in.  Default is 21.1 C (70 F).
        """
        # Array of time values ending with t:
        n = 100
        t_arr = np.linspace(0., t, n)

        # Temperature arrays:
        temp_arr = np.array([80., 90., 100., 110.])
        temp_arr_dense = heat.NewtonCooling.T(t, T0, 100., tau)

        # MS utilizations at fixed temperatures:
        t_110 = MaloShell.maloshell(t, 0.03078, 0.0068) # From MS2005 table 3
        t_100 = MaloShell.maloshell_std(t_arr)
        t_90 = MaloShell.maloshell(t, 0.00478, 0.00144) # From MS2005 table 3
        t_80 = np.zeros_like(t) # Assume no isomizeration at 80 deg

        # Stack arrays into matrix
        util_arr = np.concatenate([t_80, t_90, t_100, t_110])

        # LSQ bivariate interpolation to approximate surface
        util_func = rbs(temp_arr, t_arr, util_arr)

        return np.array([util_func(t_arr[i], temp_arr_dense[i]) for i in range(n)], dtype=float)
    
    @staticmethod
    def maloshell_cooling(t, tau=132.5, T0=21.1):
        """ Approximates isomizeration rate during cooling.  Assumes reactions
            cease at T=80.  Default tau is taken from cooling one gallon of water
            in an 8 quart stainless steel stock pot.
        Args:
            t (float): Time left to cool.
            tau (float): Cooling rate time scale in minutes.  Default is 132.5 min.
            T0 (float): Room temperature water is cooling in.  Default is 21.1 C (70 F).
        """
        # Array of time values ending with t:
        n = 100
        t_arr = np.linspace(0., t, n)

        # Temperature arrays:
        temp_arr = np.array([80., 90., 100., 110.])
        temp_final = heat.NewtonCooling.T(t, T0, 100., tau)

        # MS utilizations at fixed temperatures:
        t_110 = MaloShell.maloshell(t, 0.03078, 0.0068) # From MS2005 table 3
        t_100 = MaloShell.maloshell_std(t_arr)
        t_90 = MaloShell.maloshell(t, 0.00478, 0.00144) # From MS2005 table 3
        t_80 = np.zeros_like(t) # Assume no isomizeration at 80 deg

        # Stack arrays into matrix
        util_arr = np.concatenate([t_80, t_90, t_100, t_110])

        # LSQ bivariate interpolation to approximate surface
        util_func = rbs(temp_arr, t_arr, util_arr)

        return util_func(t, temp_final)


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

