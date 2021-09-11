"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import unicode_literals

"""
NOTE on unit cell constraints determination:
  XDS doesn't handle "real" rhombohedral space group (right?).
  So, No need to support R3 or R32. They are handled as H3 or H32, maybe.
"""

import math
import numpy
import scipy.constants
from cctbx import uctbx

class CellConstraints(object):
    def __init__(self, space_group):
        self.cs = space_group.crystal_system()

        self.free_indices = [0] # index of free parameters
        if not self.is_b_equal_a(): self.free_indices.append(1)
        if not self.is_c_equal_a_b(): self.free_indices.append(2)
        for i, an in enumerate(("alpha", "beta", "gamma")):
            if not self.is_angle_constrained(an):
                self.free_indices.append(i+3)

    # __init__()

    def is_b_equal_a(self): return self.cs in ("Tetragonal", "Hexagonal", "Trigonal", "Cubic")
    def is_c_equal_a_b(self): return self.cs == "Cubic"

    def is_angle_constrained(self, angle):
        assert angle in ("alpha", "beta", "gamma")
        if self.cs == "Triclinic": return False
        if self.cs == "Monoclinic": return angle != "beta"

        return True
    # is_angle_constrained()

    def get_label_for_free_params(self, short=True):
        ret = ["a", "b", "c",
               "al" if short else "alpha",
               "be" if short else "beta",
               "ga" if short else "gamma"]
        ret = [ret[x] for x in self.free_indices]
        return " ".join(ret)
    # get_label_for_free_params()

    def format_free_params(self, uc, lfmt="%6.2f", afmt="%5.1f", sep=" "):
        if hasattr(uc, "parameters"):
            uc = uc.parameters()
            
        return sep.join([(lfmt if x<3 else afmt)%uc[x] for x in self.free_indices])
    # format_free_params()
# class CellConstraints

def v6cell(cell):
    # taken from cctbx/uctbx/determine_unit_cell/target_uc.py
    """ Take a reduced Niggli Cell, and turn it into the G6 representation """
    if hasattr(cell, "parameters"):
        uc = cell.parameters()
    else:
        uc = cell
        
    a = uc[0] ** 2
    b = uc[1] ** 2
    c = uc[2] ** 2
    d = 2 * uc[1] * uc[2] * math.cos(uc[3]/180.*math.pi)
    e = 2 * uc[0] * uc[2] * math.cos(uc[4]/180.*math.pi)
    f = 2 * uc[0] * uc[1] * math.cos(uc[5]/180.*math.pi)
    return [a, b, c, d, e, f]


def is_same_laue_symmetry(sg1, sg2):
    laue = lambda x: x.build_derived_reflection_intensity_group(anomalous_flag=False)
    return laue(sg1) == laue(sg2) # == comparison of space_group object is possible.
# is_same_laue_symmetry()

def is_enantiomorphic_space_group_pair(sg1, sg2):
    if not sg1.type().is_enantiomorphic(): return False
    if not sg2.type().is_enantiomorphic(): return False
    
    return sg1.info().change_hand().group() == sg2
# is_enantiomorphic_space_group_pair()

def is_same_space_group_ignoring_enantiomorph(sg1, sg2):
    if sg1 == sg2: return True
    if sg1.type().is_enantiomorphic() and sg2.type().is_enantiomorphic():
        return sg1.info().change_hand().group() == sg2
    return False
# is_same_space_group_ignoring_enantiomorph()

def abc_convert_real_reciprocal(a, b, c):
    V = numpy.dot(a, numpy.cross(b, c))
    a_ = numpy.cross(b, c) / V
    b_ = numpy.cross(c, a) / V
    c_ = numpy.cross(a, b) / V
    return a_, b_, c_
# abc_convert_real_reciprocal()

def format_unit_cell(uc, lfmt="%6.3f", afmt="%5.2f", sep=" "):
    if hasattr(uc, "parameters"):
        uc = uc.parameters()

    lstr = sep.join([lfmt%x for x in uc[:3]])
    astr = sep.join([afmt%x for x in uc[3:6]])
    return lstr + sep + astr
# format_unit_cell()

def electron_voltage_to_wavelength(voltage):
    h, m0, e, c = scipy.constants.h, scipy.constants.m_e, scipy.constants.e, scipy.constants.c
    wavelen = h/numpy.sqrt(2*m0*e*voltage*(1.+e*voltage/2./m0/c**2)) * 1.e10
    return wavelen
# electron_voltage_to_wavelength()

def shelx_latt(space_group):
    # http://xray.chem.ualberta.ca/xray/shelxl/LATT.htm
    
    latt_type = str(space_group.info())[0]
    latt = dict(P=1, I=2, R=3, F=4, A=5, B=6, C=7).get(latt_type, "") # XXX check :R or :H for R.

    if not space_group.is_centric():
        latt *= -1
    
    return str(latt)
