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
        """
        return T0 + (T0 - Ti) * np.exp(-t / tau)