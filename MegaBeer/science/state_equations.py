import numpy as np

# Ideal gas law constant: 
R = 8.314 # J/mol x K

class IdealGasLaw(object):
    """ Ideal gas law EoS container.
    """
    @staticmethod
    def P(n, V, T):
        """ Base pressure function builder
        Args:
            n (float): moles
            V (float): Volume
            T (float): Temperature
        
        Returns:
            Function of n, V, and T
        """
        return lambda n, V, T: n * R * T / V

    @staticmethod
    def V(n, P, T):
        """ Base volume function builder
        Args:
            n (float): moles
            P (float): Pressure
            T (float): Temperature
        
        Returns:
            Function of n, P, and T
        """
        return lambda n, P, T: n * R * T / P

    @staticmethod
    def T(n, V, P):
        """ Base temperature function builder
        Args:
            n (float): moles
            V (float): Volume
            P (float): Pressure
        
        Returns:
            Function of n, V, and P
        """
        return lambda n, V, P: P * V / n / R

