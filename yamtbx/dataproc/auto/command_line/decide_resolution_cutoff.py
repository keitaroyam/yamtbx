from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc.auto.resolution_cutoff import estimate_resolution_based_on_cc_half
from iotbx import reflection_file_reader
import iotbx.phil
import os
import sys

master_params_str = """\
hklin = None
 .type = path
n_bins = 9
 .type = int
cc_one_half_min = 0.5
 .type = float
cc_half_tol = 0.03
 .type = float
anomalous = False
 .type = bool
show_plot = False
 .type = bool
"""

def run(params):
    hkl_in = reflection_file_reader.any_reflection_file(params.hklin)
    miller_arrays = hkl_in.as_miller_arrays(merge_equivalents=False)
    miller_arrays = [x for x in miller_arrays if x.is_xray_intensity_array()]

    print("Reading intensity data from %s" % params.hklin)
    for ma in miller_arrays:
        print(" %s" % ma.info())

    if len(miller_arrays) > 1:
        print("Using %s for analysis" % miller_arrays[0].info())
    elif len(miller_arrays) == 0:
        print("Error: No intensity data!")
        return

    iobs = miller_arrays[0]

    est = estimate_resolution_based_on_cc_half(iobs, params.cc_one_half_min,
                                               params.cc_half_tol, params.n_bins,
                                               params.anomalous, sys.stdout)
    print()
    print("Suggested cutoff= %.2f A (CC1/2= %.4f)" %(est.d_min, est.cc_at_d_min))
    print("# Treated as *%s* data" % ("anomalous" if params.anomalous else "non-anomalous"))

    if params.show_plot:
        est.show_plot(filename="cchalf_plot.pdf")
# run()    

def show_help():
    print("""\
This script determines high resolution cutoff that gives specified CC1/2 value in outer shell.

Usage: kamo.decide_resolution_cutoff xscale.hkl [n_bins=9] [cc_one_half_min=0.5]

Parameters:""")
    iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
# show_help()

def run_from_args(argv):
    if "-h" in argv or "--help" in argv:
        show_help()
        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if not os.path.exists(arg):
            print("Error: Given path does not exist: %s" % arg)
            quit()
        if params.hklin is None:
            params.hklin = arg

    if not params.hklin:
        show_help()
        return

    run(params)
# run_from_args()

if __name__ == "__main__":
    run_from_args(sys.argv[1:])
