""" Contains classes for constructing hop utilization models.
"""
import numpy as np

import MegaBeer.calculation.hops.gravity_factor as gravity_factor
import MegaBeer.calculation.hops.ibu_calculations as ibu_calculations
import MegaBeer.calculation.hops.iso_time as iso_time

class Utilization(object):
    """ Base class for all hop utilization classes
    """
    