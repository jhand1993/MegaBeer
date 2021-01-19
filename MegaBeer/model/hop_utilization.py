"""
Mathematics for hop-related calculations for recipes/plotting.
Abbreviations:
    alpha acid: AA
    specific gravity: G
    original gravity: OG
    temperature: T
    time: t
"""
import numpy as np


class MSModel(object):
    """ Class containing methods for AA isomizeration rates using the Malowicki
        & Shellhammer 2005 Model.  pH is fixed for this model.
    """

    @staticmethod
    def ms_basemodel(hop_dict, T, OG):
        pass


class TinsethModel(object):
    """ Class containing methods for AA isomizeration rates using Tinseth 1995
        PhD model.  Temperature is fixed at 100 C (212 F).  This model
        separates utilization rate function into two components f(t, SG) =
        g(t) x h(SG).
    """

    @staticmethod
    def time_component(t, max_u, r):
        """ Temporal component of Tinseth model.
        Args:
            t (float or numpy.ndarray): Boil time.
            max_u (float): Maximum utilization constant.
            r (float): Rate constant of growth.

        Results:
            float or numpy.ndarray: Time component of utilization fraction.
        """
        return (1. - np.exp(-r * t)) / max_u

    @staticmethod
    def G_component(G):
        """ Gravity component of Tinseth model.
        Args:
            G (float or numpy.ndarray): Standard gravity

        Returns:
            float or numpy.ndarray: Gravity component of utilization fraction.
        """
        return 1.65 * 0.000125**(G - 1.)

    @staticmethod
    def utilization(t, G, max_u=4.15, r=0.04):
        """ Utilization fraction for Tenseth model for given time, specific gravity.
        Args:
            t (float or numpy.ndarray): Boil time.
            G (float or numpy.ndarray): Standard gravity.
            max_u (float): Maximum utilization constant.  Default is 4.15.
            r (float): Rate constant of growth.  Default is 0.04.

        Returns:
            float or numpy.ndarray: Gravity component of utilization fraction.
        """
        return TinsethModel.time_component(t, max_u, r) * \
            TinsethModel.G_component(G)

    @staticmethod
    def utilization_grid(t_arr, G_arr, max_u=4.15, r=0.04):
        """ Grid of utilizations for given t_arr, g_arr values.
        Args:
            t_arr (numpy.ndarray): 1-D array of time values.
            G_arr (numpy.ndarray): 1-D array of gravity values.
            max_u (float): Maximum utilization constant.  Default is 4.15.
            r (float): Rate constant of growth.  Default is 0.04.

        Returns:
            numpy.ndarray: 2-D array of utilization values.
        """
        t_grid, G_grid = np.meshgrid(t_arr, G_arr)
        return t_grid, G_grid, TinsethModel.utilization(t_grid, G_grid, max_u=max_u, r=r)
