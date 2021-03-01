"""
Mathematics for hop-related calculations for recipes/plotting.
Abbreviations:
    alpha acid: AA
    specific gravity: G
    original gravity: OG
    temperature: T
    time: t (minutes)
"""
import numpy as np
from MegaBeer.science import reaction
from MegaBeer.science import heat
from scipy.interpolate import RectBivariateSpline as rbs
from scipy.integrate import quad


class MaloShell:
    """ Container class for Malowicki & Shellhammer 2005 results.
    """
    @staticmethod
    def maloshell(t, k1, k2):
        """ AA isomizeration rates using the Malowicki & Shellhammer 2005
            model for a fixed gravity. pH fixed to 5.2.

            For a fixed temperature T, the total isomizeration rate is a
            function of time: c_iso(t) = const * f(t)
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
    def maloshell_cooling_curves(t, tau, T0=21.1):
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

        return np.array(
            [util_func(t_arr[i], temp_arr_dense[i]) for i in range(n)],
            dtype=float
            )
    
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


class mIBU:
    """ Alchemy Overlords modified Tinseth utilization model
        accounting for cooling of wort after flameout. T(t) has
        been changed Newton's law of cooling for consistency.
    """
    @staticmethod
    def mIBU(
        t, t_boil, t_cool, surface_area, open_area, volume,
        max_u=0.241, r=0.04
    ):
        """ Modified Tinseth from https://alchemyoverlord.wordpress.com/
        Args:
            t (float or numpy.ndarray): Iso time.  Total time hop addition(s)
                is(are) the wort.
            t_boil (float): Boil time.
            t_cool (float): Cooling time.
            surface_area (float): Exposed wort surface area in square 
                centimeters.
            open_area (float): Size of opening of pot in square centimeters.
            volume (float): Volume of wort in liters.
            max_u (float): Maximum utilization constant.  Default is 0.241.
            r (float): Rate constant of growth.  Default is 0.04.

        Results:
            float or numpy.ndarray: Time component of utilization fraction.
        """

        # Convert float to array:
        t_arr = np.asarray(t)

        # Mask that is true is t >= t_cool, false otherwise:
        t_mask = t_arr >= t_cool

        # mIBU model timescale:
        b = mIBU.b(surface_area, open_area, volume)

        # calculate utilization at constant temperature:
        boil_util = np.where(
            t_mask, TinsethTime.tinseth(t_arr, max_u=max_u, r=r), 0.
        )

        # Cooling rate to integrate
        cool_rate = lambda x: TinsethTime.tinseth_rate(x, max_u=max_u, r=r) * \
            mIBU.mIBU_rate_correction(t, b)
        
        # Integrate to calculate cooling utilization.  Note that t_cool - t
        # to t_cool is the total cooling time as t is the total time in the wort:
        cool_util = np.where(
            t_mask,
            quad(cool_rate, 0., t_cool), quad(cool_rate, t_cool - t, t_cool)
        )

        return boil_util + cool_util
    
    @staticmethod
    def b(surface_area, open_area, volume):
        """ Timescale (tau equivalent) from Alchemy Overlord.
        Args:
            surface_area (float): Exposed wort surface area in square 
                centimeters.
            open_area (float): Size of opening of pot in square centimeters.
            volume (float): Volume of wort in liters.
        
        Results:
            float: b
        """
        eff_area = np.sqrt(surface_area * open_area)
        return 2.925e-4 * eff_area / volume + 5.38e-3
    
    @staticmethod
    def mIBU_rate_correction(t, b):
        """ mIBU relative rate differential correction factor.
        Args:
            t (float or numpy.ndarray): time
            b (float): Temperature decay time scale.
        
        Returns:
            float or numpy.ndarray: relative rate correction
        """
        # Constants for mIBU model:
        c1 = 2.39e11
        c2 = 9773.  # Units of E_activation / R
        c3 = 53.7
        c4 = 319.95

        return c1 * np.exp(-c2 / (c3 * np.exp(-b * (t)) + c4))

class TinsethTime():
    """ Container class for Tinseth temporal component calculations
    """
    @staticmethod
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

    @staticmethod
    def tinseth_rate(t, max_u=0.241, r=0.04):
        """ Derivative of the temporal component of Tinseth model.
            Note: 0.241 = 1 / 4.15.
        Args:
            t (float or numpy.ndarray): Boil time.
            max_u (float): Maximum utilization constant.  Default is 0.241.
            r (float): Rate constant of growth.  Default is 0.04.

        Results:
            float or numpy.ndarray: Time component of utilization fraction.
        """
        return max_u * r * np.exp(-r * t)

