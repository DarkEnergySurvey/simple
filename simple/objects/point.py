#!/usr/bin/env python
"""
Classy points
"""
__author__ = "Sid Mau"

# Python libraries

# Ugali libraries
import ugali.utils.healpix
import ugali.utils.projector

# TODO:

########################################################################

class Point:
    """
    Class object for working with data of a point in the sky.
    """
    def __init__(self, ra, dec):
        self.ra  = ra
        self.dec = dec

    @property
    def projector(self):
        """
        Return a projector object for the point's spatial position.
        """
        return ugali.utils.projector.Projector(self.ra, self.dec)

    def healpixel(self, nside):
        """
        Return the healpixel containing the point given an nside.
        """
        return ugali.utils.healpix.angToPix(nside, self.ra, self.dec)

#class HotSpot(Spot):
#    """
#    Class object for working with data of a hotspot from a data search.
#    """
#    def __init__(self, mod=0., sig=0.):
#        self.mod = mod
#        self.sig = sig
