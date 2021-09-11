#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

from yamtbx.dataproc.xds import xscale
from yamtbx.dataproc.xds.command_line import xscale_simple

master_params_str = """\
use_tmpdir_if_available = False
 .type = bool
 .help = "Use RAM disk or temporary directory to run xscale"
cbf_to_dat = True
 .type = bool
 .help = "Convert .cbf files (scale values) to a few .dat files"
aniso_analysis = True
 .type = bool
 .help = "Do anisotropy analysis"
"""

import os
import sys
import iotbx.phil
from yamtbx import util

def inp_file_analysis(inp):
    import numpy
    from yamtbx.util.xtal import format_unit_cell

    cells = []
    for l in open(inp):
        ll = l[:l.index("!")] if "!" in l else l
        if "INPUT_FILE=" in ll: #  and len(l) > 132: # one line is limited to 131 characters!
            filename = ll[ll.index("=")+1:].strip()
            if "*" in filename: filename = filename[filename.index("*")+1:].strip()
            print(ll.strip())
            info = xscale_simple.get_xac_info(filename, get_nframes=False)
            print(" FRIEDEL'S_LAW=", info.get("friedels_law"))
            print(" RESOLUTION_RANGE=", info.get("resol_range"))
            print(" SPACE_GROUP_NUMBER=", info.get("spgr_num"))
            print(" UNIT_CELL_CONSTANTS=", info.get("cell"))
            cells.append(list(map(float, info["cell"].split())))
    
    print()
    print()
    cells = numpy.array(cells)
    
    print("Averaged cell =", format_unit_cell(cells.mean(axis=0)))
    print("    (std dev) =", format_unit_cell(cells.std(axis=0)))
    print("  Median cell =", format_unit_cell(numpy.median(cells, axis=0)))

def run(params, inp):
    try:
        inp_file_analysis(inp)
    finally:
        print()
        print("Starting XSCALE")
        print()
        xscale.run_xscale(inp,
                          cbf_to_dat=params.cbf_to_dat,
                          aniso_analysis=params.aniso_analysis,
                          use_tmpdir_if_available=params.use_tmpdir_if_available)
# run()

def show_help():
    print("""\
This script runs xscale program. No limitation in INPUT_FILE= length.

Usage: yamtbx.run_xscale XSCALE.INP [use_tmpdir_if_available=true]

Parameters:""")
    iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
# show_help()

def run_from_args(args):
    if "-h" in args or "--help" in args:
        show_help()
        return

    cmdline = iotbx.phil.process_command_line(args=args,
                                              master_string=master_params_str)
    params = cmdline.work.extract()

    if len(cmdline.remaining_args) == 0:
        xscale_inp = "XSCALE.INP"
    else:
        xscale_inp = cmdline.remaining_args[0]

    if not os.path.isfile(xscale_inp):
        print("Cannot find %s" % xscale_inp)
        show_help()
        return

    run(params, xscale_inp)
# run_from_args()

if __name__ == "__main__":
    run_from_args(sys.argv[1:])
