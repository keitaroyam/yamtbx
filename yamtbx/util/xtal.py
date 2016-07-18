"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

"""
NOTE on unit cell constraints determination:
  XDS doesn't handle "real" rhombohedral space group (right?).
  So, No need to support R3 or R32. They are handled as H3 or H32, maybe.
"""
import math
import numpy
from cctbx import uctbx

class CellConstraints:
    def __init__(self, space_group):
        self.cs = space_group.crystal_system()
    # __init__()

    def is_b_equal_a(self): return self.cs in ("Tetragonal", "Hexagonal", "Trigonal", "Cubic")
    def is_c_equal_a_b(self): return self.cs == "Cubic"

    def is_angle_constrained(self, angle):
        assert angle in ("alpha", "beta", "gamma")
        if self.cs == "Triclinic": return False
        if self.cs == "Monoclinic": return angle != "beta"

        return True
    # is_angle_constrained()
# class CellConstraints

def v6cell(cell):
    # taken from cctbx/uctbx/determine_unit_cell/target_uc.py
    uc = cell.parameters()
    """ Take a reduced Niggli Cell, and turn it into the G6 representation """
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

def abc_convert_real_reciprocal(a, b, c):
    V = numpy.dot(a, numpy.cross(b, c))
    a_ = numpy.cross(b, c) / V
    b_ = numpy.cross(c, a) / V
    c_ = numpy.cross(a, b) / V
    return a_, b_, c_
# abc_convert_real_reciprocal()

def format_unit_cell(uc, lfmt="%6.2f", afmt="%5.1f", sep=" "):
    if hasattr(uc, "parameters"):
        uc = uc.parameters()

    lstr = sep.join(map(lambda x: lfmt%x, uc[:3]))
    astr = sep.join(map(lambda x: afmt%x, uc[3:6]))
    return lstr + sep + astr
# format_unit_cell()
