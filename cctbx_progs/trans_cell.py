#!/usr/bin/env cctbx.python

from __future__ import print_function
from __future__ import unicode_literals
from iotbx import crystal_symmetry_from_any
from cctbx import uctbx
from cctbx import sgtbx

if __name__ == "__main__":
    import optparse
    import sys

    parser = optparse.OptionParser(usage="usage: %prog [options] coordinates...")

    parser.add_option("-s", action="store", dest="sym_source", help="Symmetry source (pdb, mtz, etc)")
    parser.add_option("-c", action="store", dest="cell", help="unit cell parameter")

    (opts, args) = parser.parse_args(sys.argv)

    cell = None

    if opts.sym_source is not None:
        symm = crystal_symmetry_from_any.extract_from(opts.sym_source)
        cell = symm.unit_cell()
    else:
        cell = uctbx.unit_cell(opts.cell)

    op = sgtbx.change_of_basis_op(args[1])
    print("   Original cell:", cell)
    print("        Operator:", op)
    print("Transformed cell:", cell.change_basis(op))
