from __future__ import absolute_import
from __future__ import unicode_literals
from . import geom
from . import hkl
from . import stream
from iotbx import crystal_symmetry_from_any

def get_run_and_tag_from_filename(filename):
    """e.g. SFX_each_193441_79553414.h5"""
    import os
    sp = os.path.splitext(os.path.basename(filename))[0].split("_")
    if len(sp) >= 2:
        run = int(sp[-2])
        tag = int(sp[-1])
        return run, tag
    else:
        return -1, -1

def read_cell_file(filename):
    query = crystal_symmetry_from_any.extract_from(filename)
    if query is not None: return query

    ifs = open(filename)
    fl = ifs.readline().strip()
    assert fl.startswith("CrystFEL unit cell file version ")
    filever = fl.split()[-1]

    cellkeys = ("a", "b", "c", "al", "be", "ga")
    cell = (0, 0, 0, 0, 0, 0)
    latt_type, centering = None, None

    if filever == "1.0":
        for l in ifs:
            sp = [x.strip() for x in l.split("=")]
            if len(sp) != 2: continue
            lhs, rhs = sp
            if lhs in cellkeys:
                cell[cellkeys.index(lhs)] = float(rhs.replace("A","").replace("deg",""))
            elif lhs == "lattice_type":
                latt_type = rhs
            elif lhs == "centering":
                centering = rhs

    return cell, latt_type, centering
# read_cell_file()
