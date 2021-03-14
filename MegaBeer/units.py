""" Handle unit system conversion
"""

class Units():
    """ Base class for unit conversion.
    """
    def __init__(
        self, 
        mass, length, volume, temperature,
        mass_back, length_back, volume_back, temperature_back
        ):
        self.mass = mass
        self.length = length
        self.volume = volume
        self.temperature = temperature
        self.mass_back = mass_back
        self.length_back = length_back
        self.volume_back = volume_back
        self.temperature_back = temperature_back
    
    def convert_mass(self, m):
        """ Converts mass units to metric
        """
        return self.mass(m)
    
    def convert_mass_back(self, m):
        """ Converts mass back from metric
        """
        return self.mass_back(m)

    def convert_volume(self, v):
        """ Converts volume units to metric
        """
        return self.volume(v)
    
    def convert_volume_back(self, v):
        """ Converts volume back from metric
        """
        return self.volume_back(v)

    def convert_temperature(self, T):
        """ Converts temperature units to metric
        """
        return self.temperature(T)
    
    def convert_temperature_back(self, T):
        """ Converts temperature back from metric
        """
        return self.temperature_back(T)



class Metric(Units):
    """ Class for metric conversions. Methods are trivial as all
        calculations are performed assuming metric units.
    """
    def __init__(self):
        super().__init__(
            nochange, nochange, nochange, nochange, 
            nochange, nochange, nochange, nochange
            )


class MetricGrams(Units):
    """ Class for metric conversions but with grams instead of kilograms
    """
    def __init__(self):
        super().__init__(
            lambda x: orderchange(x, 3.), nochange, nochange, nochange, 
            lambda x: orderchange(x, -3.), nochange, nochange, nochange
            )




def c_to_f(T):
    """ Temperature in Celsius to Fahrenheit.
    Args:
        T (float or numpy.ndarray): Temperature in Celsius.
    
    Returns:
        float or numpy.ndarray: Temperature in Fahrenheit.
    """
    return T * (9. / 5.) + 32.

def f_to_c(T):
    """ Temperature in Celsius to Fahrenheit.
    Args:
        T (float or numpy.ndarray): Temperature in Fahrenheit.
    
    Returns:
        float or numpy.ndarray: Temperature in Celsius.
    """
    return (T - 32) * 5. / 9.

def nochange(x):
    """ No change mathematical function f(x)=x.
    """
    return x

def orderchange(x, order):
    """ Changes order of magnitude: f(x) = 10**order * x
    """
    return 10.**order * x