from numpy import np

class NewtonCooling:
    """ Container class for Newton cooling law methods. Law states that
        cooling of a material over time is T(t) = T0 + DT * exp(-t/tau), 
        where T0 is ambient temperature, DT = T0-Ti is the difference in
        starting temperature and ambient temperature, and tau=C/hA is a
        time scale.  C is the heat capacity, h is the heat transfer
        coefficient, and A is the surface area of the cooling object.
    """
    @staticmethod
    def T(t, T0, Ti, tau):
        """ Base law using tau instead of full C/hA
        Args:
            t (float or numpy.ndarray): Time
            T0 (float): Ambient temperature.
            Ti (float): Initial temperature of material.
            tau (float): Timescale.
        
        Returns:
            float: Temperature
        """
        return T0 + (T0 - Ti) * np.exp(-t / tau)
    
    @staticmethod
    def tau_approximator(
        water_mass, pot_radius, fill_factor=0.5, metric=True
    ):
        """ Approximates tau as:
            1 / (2 * pi) * mass / (radius^2 * [1+fill_factor/2]) * 875.
            875 ~= specific heat / heat transfer coefficient for water
            in a stock pot with a diameter about equal to its height.
            In that case, height = 2 * radius.  This simplifies area
            calculation for a cylinder from 2 * pi * radius^2 (1 + h/r)
            to 2 * pi * radius^2 * (1 + 1/2).  Fill factor considers
            fraction of pot filled with water.
        Args:
            water_mass (float): Mass of water.
            pot_radius (float): Radius of pot.
            fill_factor (float): fraction of pot volume filled with water.
            metric (bool): If True, metric system is used and mass is in kg,
                radius is in meters. Otherwise, lbs and inches are used.
                Default is True.
            
        Returns:
            float: tau
        """
        if metric:
            c_h = 875.
        
        else:
            c_h = 606455.
        
        return c_h * water_mass /\
            (2. * np.pi * pot_radius**2 *(1. + fill_factor / 2.))
