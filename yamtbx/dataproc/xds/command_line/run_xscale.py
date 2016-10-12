#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.xds import xscale

master_params_str = """\
use_tmpdir_if_available = False
 .type = bool
 .help = "Use RAM disk or temporary directory to run xscale"
cbf_to_dat = True
 .type = bool
 .help = "Convert .cbf files (scale values) to a few .dat files"
"""

import os
import sys
import iotbx.phil
from yamtbx import util

def run(params, inp):
    xscale.run_xscale(inp,
                      cbf_to_dat=params.cbf_to_dat,
                      use_tmpdir_if_available=params.use_tmpdir_if_available)
# run()

def show_help():
    print """\
This script runs xscale program. No limitation in INPUT_FILE= length.

Usage: yamtbx.run_xscale XSCALE.INP [use_tmpdir_if_available=true]

Parameters:"""
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
        show_help()
        return

    xscale_inp = cmdline.remaining_args[0]

    run(params, xscale_inp)
# run_from_args()

if __name__ == "__main__":
    run_from_args(sys.argv[1:])
