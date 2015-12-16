"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

"""
Just do an evaluation part of automatic processing.
"""

import os

from yamtbx.dataproc.xds.command_line import xds_plot_integrate
from yamtbx.dataproc.auto.command_line import run_all_xds_simple
from yamtbx import util

import iotbx.phil
from libtbx import easy_mp

master_params_str = """
topdir = None
 .type = path
nproc = None
 .type = int
 .help = number of processors for single xds job OR number of parallel jobs
"""

def evaluate_run(root):
    integrate_lp = os.path.join(root, "INTEGRATE.LP")
    xds_ascii_hkl = os.path.join(root, "XDS_ASCII.HKL")

    if os.path.isfile(integrate_lp):
        xds_plot_integrate.run(integrate_lp, os.path.join(root, "plot_integrate.log"))

    if os.path.isfile(xds_ascii_hkl):
        ret = run_all_xds_simple.calc_merging_stats(xds_ascii_hkl)
        run_all_xds_simple.run_xdsstat(wdir=root)
# evaluate_run()

def run(params):
    xds_dirs = []
    print "Found xds directories:"
    for root, dirnames, filenames in os.walk(params.topdir, followlinks=True):
        if "XDS.INP" in filenames:
            print "", os.path.relpath(root, params.topdir)
            xds_dirs.append(root)

    print
    print "Start running.."

    npar = util.get_number_of_processors() if params.nproc is None else params.nproc

    fun_local = lambda x: evaluate_run(x)#, params)
    easy_mp.pool_map(fixed_func=fun_local,
                     args=xds_dirs,
                     processes=npar)
# run()


if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    
    for arg in args:
        if os.path.isdir(arg) and params.topdir is None:
            params.topdir = arg

    if params.topdir is None:
        params.topdir = os.getcwd()

    run(params)
