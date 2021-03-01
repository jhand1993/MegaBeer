import numpy as np

class RateEquations(object):
    """ Container class for rate equations. Integer power law
        rate equations orders zero, one, and two are explicitly
        defined.
    """
    @staticmethod
    def order_n(n, k, A0):
        """ n order rate base method 
        Args:
            n (int): Order of reaction.
            k (float): Rate constant.
            A0 (float): Starting amount.
        
        Returns:
            Function of time
        """
        # Handle order one reactions differently:
        if n == 1:
            return lambda t: A0 * np.exp(-k * t)

        g = 1. - float(n) # Convert to float.
        return lambda t: (A0**g - g * k * t)**(1. / g)
    
    @staticmethod
    def order_zero(k, A0):
        """ Order zero rate equation
        Args:
            k (float): Rate constant.
            A0 (float): Starting amount.
        
        Returns:
            Function of time
        """
        return RateEquations.order_n(0, k, A0)

    @staticmethod
    def order_one(k, A0):
        """ Order one rate equation
        Args:
            k (float): Rate constant.
            A0 (float): Starting amount.
        
        Returns:
            Function of time
        """
        return RateEquations.order_n(1, k, A0)

    @staticmethod
    def order_two(k, A0):
        """ Order two rate equation
        Args:
            k (float): Rate constant.
            A0 (float): Starting amount.
        
        Returns:
            Function of time
        """
        return RateEquations.order_n(2, k, A0)


def arrhenius(T, E0, A):
    """ Arrhenius equation modeling temperature dependence on reaction rate.
    Args:
        T (float): Temperature in Celcius.  Converted to Kelvin
        E0 (float): Activation energy in Joules
        A (float): Pre-exponential factor
    
    Returns:
        float: reaction rate k in min^-1
    """
    R = 8.3145
    return A * np.exp( - E0 / (T + 272.15) / R)
    
