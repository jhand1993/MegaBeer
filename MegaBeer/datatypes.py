"""
Datatypes for brewing. 
"""

class DBObject:
    """ Datatype which information can be loaded from database
    """
    def __init__(self, ingred_type):
        self.ingred_type = ingred_type

class Hop(DBObject):
    """ Datatype for hops.
    """
    def __init__(self, AA, name=None):
        self.AA = AA
        self.name = name
        super().__init__('hop')
    