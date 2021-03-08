""" IBU calculations
"""

class IBUCalculation():
    """ Container class for IBU calculations
    """
    @staticmethod
    def tinseth_ibu(AA, m, u, c):
        """ IBU calcuation: IBU = c * AA * m * u
        Args:
            AA (float): Alpha acid percentage.
            m (float): Mass of hop addition with alpha acid percentage AA.
            u (float): Hop utilization fraction.
            c (float): Constant to convert mass and per volume to milligrams per liter.  For
                example, if mass is in oz and volume is in gallons, c = 74.89.  If mass
                is in grams and volume in liters, then c = 10. If mass is in grams and volume in
                gallons, c=2.64.
        
        Returns:
            float: IBUs with units milligrams per liter.
        """
        return c * u * m * AA