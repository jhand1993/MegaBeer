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
from MegaBeer.science.heat import NewtonCooling
# from scipy.interpolate import RectBivariateSpline as rbs
from scipy.integrate import quad, odeint
from scipy.interpolate import interp1d


class MaloShell:
    """ Container class for Malowicki & Shellhammer 2005 function constructors.
        Note: c = [c1, c2, c3] = [AA, iso-AA, degradation products]
    """
    @staticmethod
    def maloshell_constant_temp(k1, k2):
        """ AA isomizeration rates using the Malowicki & Shellhammer 2005
            model for a fixed gravity. pH fixed to 5.2.  Assumes k1 !=k2.

            For a fixed temperature T, the total isomizeration rate is a
            function of time: c_iso(t) = const * f(t)
        Args:
            k1 (float): Isomerization reaction rate
            k2 (float): Iso-AA degradation rate
        
        Returns:
            numpy.ndarray: 3d vector of fractional concentration functions 
        """
        c1 = lambda t: k1 * np.exp(-k1 * t)
        c2 = lambda t: k1 / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t))
        c3 = lambda t: k2 * np.exp(-k2 * t)
        return np.asarray([c1, c2, c3])

    @staticmethod
    def maloshell_boil_temp():
        """ Malowicki & Shellhammer 2005 model but with fixed boiling temperature.
            Refer to table 3 in paper: k1 = 0.01141, k2 = 0.00263 for T=100.
        Args:
            t (float or numpy.ndarray): Boil time.
        
        Returns:
            float or numpy.ndarray: time component of utilization fraction
        """
        k1_boil = 0.01141
        k2_boil = 0.00263
        return MaloShell.maloshell_constant_temp(k1_boil, k2_boil)

    @staticmethod
    def maloshell_cooling(A1, A2, Ea_1, Ea_2, c0, tau=132.5, T_room=21.1):
        """ Approximates isomizeration rate during cooling, and assumes no
            significant isomizeration after one e-fold (tau minutes) have past.
            Default tau is taken from cooling one gallon of water
            in an 8 quart stainless steel stock pot with a steel lid on.
        Args:
            A1 (float): Exponential prefactor for k1.
            A2 (float): Exponential prefactor for k2.
            Ea_1 (float): Activation energy for reaction 1.
            Ea_2 (float): Activation energy for reaction 2.
            c0 (np.ndarray): Initial condition vector.
            tau (float): Cooling rate time scale in minutes.  Default is 132.5 min.
            T_room (float): Room temperature water is cooling in.  Default is 21.1 C (70 F).
        
        Returns:
            function: Utilization function vector for c1, c2, c3
        """
        # Make sure c0 is an array:
        c0 = np.asarray(c0)

        # Create time array with spacings of one minute.
        t_arr = np.linspace(0., tau, np.ceil(tau) + 1)

        # Grab the Newton Cooling function T(t) with T0=T_room and Ti=100 C
        temp_func = NewtonCooling.T(T_room, 100., tau)

        # Get rate functions:
        k1_func = lambda t: reaction.arrhenius(Ea_1, A1)(temp_func(t))
        k2_func = lambda t: reaction.arrhenius(Ea_2, A2)(temp_func(t))

        def dcdt(c, t):
            # Internal function to solve ODE. c = np.array([c1, c2, c3])
            M = np.zeros((3, 3))
            M[0, 0] = -reaction.RateEquations.dkn_dt(1, k1_func(t))(t)
            M[1, 0] = reaction.RateEquations.dkn_dt(1, k1_func(t))(t)
            M[1, 1] = -reaction.RateEquations.dkn_dt(1, k2_func(t))(t)
            M[2, 1] = reaction.RateEquations.dkn_dt(1, k2_func(t))(t)

            # print(M)
            return np.dot(M, c)

        c = odeint(dcdt, c0, t_arr)

        # Extrapolation fill values for each component of vector c.
        ext = [(c[0, i], c[-1, i]) for i in range(3)]

        # Linear interpolation for each vector component.  Assumes no utilization below t=0
        # and fixed utilization above t value given as input.
        return np.asarry(
                [interp1d(t_arr, c[:, i], fill_value=ext[i]) for i in range(3)]
            )


class mIBU:
    """ Alchemy Overlords modified Tinseth utilization model accounting for cooling 
        of wort after flameout.  Note that the piecewise nature of these math functions
        complicated a functional approach, hence the OOP approach here.
    Args:
        max_u (float): Maximum utilization constant.  Default is 0.241.
        r (float): Rate constant of growth.  Default is 0.04.
        surface_area (float): Exposed wort surface area in square 
            centimeters.
        open_area (float): Size of opening of pot in square centimeters.
        volume (float): Volume of wort in liters.
    """
    def __init__(
        self, surface_area, open_area, volume, max_u=0.241, r=0.04
        ):
        self.surface_area = surface_area
        self.open_area = open_area
        self.volume = volume
        self.b = mIBU.calculate_b(surface_area, open_area, volume)
        self.max_u = max_u
        self.r = r

    @staticmethod
    def calculate_b(surface_area, open_area, volume):
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

    def change_b(self, b):
        """ Convenience method to update b attribute
        Args:
            b (float): New b value
        """
        if type(b) is float:
            self.b = b
            return True
        
        else:
            raise ValueError(float)

    def mIBU(self, t, t_boil, t_cool):
        """ Modified Tinseth from https://alchemyoverlord.wordpress.com/
        Args:
            t (float or numpy.ndarray): Iso time.  Total time hop addition(s)
                is(are) the wort.
            t_boil (float): Boil time.
            t_cool (float): Cooling time.

        Results:
            float or numpy.ndarray: Time component of utilization fraction.
        """

        # Convert float to array:
        t_arr = np.asarray(t)

        # Mask that is true is t >= t_cool, false otherwise:
        t_mask = t_arr >= t_cool

        # calculate utilization at constant temperature:
        boil_util = np.where(
            t_mask, TinsethTime.tinseth(max_u=self.max_u, r=self.r)(t_arr), 0.
        )

        # Cooling rate to integrate
        cool_rate = lambda x: TinsethTime.tinseth_rate(max_u=self.max_u, r=self.r)(x) * \
            mIBU.mIBU_rate_correction(t, self.b)
        
        # Integrate to calculate cooling utilization.  Note that t_cool - t
        # to t_cool is the total cooling time as t is the total time in the wort:
        cool_util = np.where(
            t_mask,
            quad(cool_rate, 0., t_cool), quad(cool_rate, t_cool - t, t_cool)
        )

        return boil_util + cool_util
    
    def mIBU_rate_correction(self, t):
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

        return c1 * np.exp(-c2 / (c3 * np.exp(-self.b * (t)) + c4))

class TinsethTime():
    """ Container class for Tinseth temporal component calculations
    """
    @staticmethod
    def tinseth(max_u=0.241, r=0.04):
        """ Temporal component of Tinseth model. Max_u and r taken from Palmer.
            Note: 0.241 = 1 / 4.15.
        Args:
            max_u (float): Maximum utilization constant.  Default is 0.241.
            r (float): Rate constant of growth.  Default is 0.04.

        Results:
            function: Time component of utilization function.
        """
        return lambda t: max_u * (1. - np.exp(-r * t)) 

    @staticmethod
    def tinseth_rate(max_u=0.241, r=0.04):
        """ Derivative of the temporal component of Tinseth model.
            Note: 0.241 = 1 / 4.15.
        Args:
            max_u (float): Maximum utilization constant.  Default is 0.241.
            r (float): Rate constant of growth.  Default is 0.04.

        Results:
            function: Time component of utilization function derivative.
        """
        return lambda t: max_u * r * np.exp(-r * t)

