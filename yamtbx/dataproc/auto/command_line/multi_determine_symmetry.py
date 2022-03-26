"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

from builtins import range
from yamtbx.dataproc.auto.multi_merging.resolve_reindex import KabschSelectiveBreeding
from yamtbx.util import read_path_list
from libtbx.utils import multi_out
import iotbx.phil
import libtbx.phil
from cctbx import sgtbx
from cctbx import crystal
import os
import sys
import copy
import numpy

master_params_str = """\
lstin = None
 .type = path
 .help = list of XDS_ASCII.HKL
method = pointless *selective_breeding
 .type = choice(multi=False)
 .help = method
from_p1 = False
 .type = bool
 .help = Start from P1 whatever the symmetry of input files.
logfile = "multi_determine_symmetry.log"
 .type = path
 .help = logfile name
nproc = 1
 .type = int
 .help = number of processors

d_min = 3
 .type = float
 .help = high resolution cutoff used in the method
min_ios = None
 .type = float
 .help = minimum I/sigma(I) cutoff used in the method
max_delta = 5
 .type = float
 .help = maximum obliquity used in determining the lattice symmetry, using a modified Le-Page algorithm.

max_cycles = 100
 .type = int(value_min=1)
 .help = Maximum number of cycles for selective_breeding algorithm.

"""

def run(params):
    log_out = multi_out()
    log_out.register("log", open(params.logfile, "w"), atexit_send_to=None)
    log_out.register("stdout", sys.stdout)

    libtbx.phil.parse(master_params_str).format(params).show(out=log_out, prefix=" ")

    xac_files = read_path_list(params.lstin, only_exists=True, err_out=log_out)
    if len(xac_files) == 0:
        print("No (existing) files in the list: %s" % params.lstin, file=log_out)
        return

    if params.method == "selective_breeding":
        rb = KabschSelectiveBreeding(xac_files, max_delta=params.max_delta,
                                     d_min=params.d_min, min_ios=params.min_ios,
                                     nproc=params.nproc, log_out=log_out,
                                     from_p1=params.from_p1)
        xs = rb.representative_crystal_symmetry()

        log_out.write("Starting from:\n")
        xs.show_summary(log_out, "  ")
        log_out.write("\n")

        rb.assign_operators(max_cycle=params.max_cycles)
        rb.show_assign_summary()
        final_cc_means = rb.final_cc_means()
        assert len(final_cc_means) == len(xac_files)
        reidx_ops = rb.reindex_operators()
        sg = copy.copy(xs.space_group())
        unit_cell = xs.unit_cell()

        cc0 = [x[0][1] for x in final_cc_means]
        log_out.write("Analyzing KabschSelectiveBreeding result..\n")

        accepted_ops = []

        for iop in range(1, len(reidx_ops)):
            cci = [x[iop][1] for x in final_cc_means]
            corr = numpy.corrcoef(cc0, cci)[0,1]
            log_out.write("  h,k,l vs %s: corr= %.4f\n" % (reidx_ops[iop].as_hkl(), corr))
            if corr > 0.5:
                accepted_ops.append(reidx_ops[iop])
                sg.expand_smx(reidx_ops[iop].as_hkl())
                unit_cell = unit_cell.change_basis(reidx_ops[iop])
                log_out.write("    this operator accepted. sg= %s\n" % sg.info())

        log_out.write("Solution:\n")
        new_xs = crystal.symmetry(unit_cell, space_group=sg)
        new_xs.show_summary(log_out, "  ")
        log_out.write("As reference setting:\n")
        new_xs.as_reference_setting().show_summary(log_out, "  ")

        log_out.write("Initial:\n")
        xs.show_summary(log_out, "  ")

        log_out.write("""
* Notice *
Here the space group is deduced from the similarity of reflection intensities under the constraint of lattice symmetry.
This could be wrong especially when the crystal is twineed.
Please note that space group is only determined when the structure is solved.
""")

    else:
        raise "invalid method choice (method=%s)" % params.method

# run()

def show_help():
    print("""
""")
    iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
    print() 
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
        if os.path.isfile(arg) and params.lstin is None:
            params.lstin = arg

    if params.lstin is None:
        show_help()
        print("Error: Give .lst of XDS_ASCII files")
        quit()

    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
