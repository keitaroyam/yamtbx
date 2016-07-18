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
cbf_to_dat = True
 .type = bool
"""

import os
import sys
import iotbx.phil
from yamtbx import util

def run(params, inp):
    xscale.run_xscale(xscale_inp,
                      cbf_to_dat=params.cbf_to_dat,
                      use_tmpdir_if_available=params.use_tmpdir_if_available)
# run()

if __name__ == "__main__":
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    xscale_inp = cmdline.remaining_args[0]

    run(params, xscale_inp)
