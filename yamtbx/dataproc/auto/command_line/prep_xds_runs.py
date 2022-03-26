"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import pickle
from yamtbx.dataproc.bl_logfiles import BssJobLog
from yamtbx.dataproc.dataset import find_existing_files_in_template
from yamtbx.dataproc.xds import xds_inp
from cctbx import sgtbx
from cctbx import uctbx
from cctbx import crystal
import iotbx.phil


master_params_str = """
workdir = None
 .type = path
 .help = working directory for xds runs
pklin = None
 .type = path
dont_overwrite = true
 .type = bool
min_spots_for_ends = None
 .type = int
 .help = minimum number of spots for both ends (to decide DATA_RANGE=)
anomalous = False
 .type = bool
cell = None
 .type = floats
space_group = None
 .type = str
frames_for_index = first *all
 .type = choice(multi=False)
integrate_nimages = None
 .type = int
 .help = the number of frames for 3d profile generation. it affects DELPHI=.
help = False
 .type = bool
"""

def decide_data_range_using_nspots(spotstat, min_spots_for_ends):
    n_spots = [(idx, stat.spots.get_n_spots("hi_pass_resolution_spots"))
               for idx, f, stat in spotstat if stat is not None]

    first, last = None, None

    # Decide first
    for idx, nspots in n_spots:
        if nspots >= min_spots_for_ends:
            first = idx
            break
    
    # Decide last
    for idx, nspots in reversed(n_spots):
        if nspots >= min_spots_for_ends:
            last = idx
            break
    
    return first, last
# decide_data_range_using_nspots()

def make_shikalog(spotstat, out):
    print("   idx filename n_spots", file=out)
    for idx, f, stat in spotstat:
        if stat is None:
            continue
        print("%6d %s %5d" % (idx, f, stat.spots.get_n_spots("hi_pass_resolution_spots")), file=out)
# make_shikalog()     

def create_crystal_symmetry(sg, cell):
    space_group = None if sg is None else sgtbx.space_group_info(sg).group()
    unit_cell = uctbx.infer_unit_cell_from_symmetry(cell, space_group)
    xs = crystal.symmetry(unit_cell, space_group=space_group)
    return xs
# create_crystal_symmetry

def run(params):#data_dir, wdir, use_normalized=False):
    if (params.space_group, params.cell).count(None) == 1:
        print("If you want to specify cell or symmetry, give both.")
        return

    if params.cell is not None and len(params.cell) > 6:
        print("Too many parameters for unit cell:", params.cell)
        return

    xs = None if params.space_group is None else create_crystal_symmetry(params.space_group, params.cell)

    logobjects = pickle.load(open(params.pklin, "rb"))

    # set topdir
    topdir = os.path.dirname(os.path.commonprefix([os.path.dirname(x[0]) for x in logobjects]))

    for log, logobj, spot_stats in logobjects:
        relpath = os.path.relpath(log, topdir)
        if relpath.startswith(os.sep) or relpath.startswith(".."):
            print("Outside:", log)
            continue

        if len(logobj.jobs) > 1:
            print("##############################################################")
            print(relpath)
            print("Why %d jobs!?" % len(logobj.jobs))
            print("May be overwritten? Additional images?? Need to be careful!!") 
            print("Currently, I just pick up the last information.")
            print("##############################################################")

        job, spotstat = logobj.jobs[-1], spot_stats[-1]
        range_inlog = job.get_frame_num_range()
        if None in range_inlog: dirname = "xds_%s" % (job.prefix)
        else: dirname = "xds_%s_%.4d-%.4d" % ((job.prefix,)+job.get_frame_num_range())
        wd = os.path.join(params.workdir, os.path.dirname(relpath), dirname)

        if os.path.isfile(os.path.join(wd, "XDS.INP")) and params.dont_overwrite:
            continue

        print(wd, len(job.images))
        data_range = 1, job.n_images
        if params.min_spots_for_ends is not None:
            data_range = decide_data_range_using_nspots(spotstat, params.min_spots_for_ends)
            if None in data_range or data_range[1]-data_range[0] < 3:
                print("  Oh no!! not useful data!")
                if None in data_range:
                    print("  No images contains more than %d spots." % params.min_spots_for_ends)
                else:
                    print("  Only %d images contains more than %d spots." % (data_range[1]-data_range[0]+1,
                                                                             params.min_spots_for_ends))
                continue
            print(" Data range:", data_range)

        # How about okkake sekibun?
        img_files = find_existing_files_in_template(job.filename, data_range[0], data_range[1],
                                                    datadir=os.path.dirname(log), check_compressed=True)

        os.makedirs(wd)
        xdsinp_str = xds_inp.generate_xds_inp(img_files=img_files,
                                              inp_dir=os.path.abspath(wd),
                                              reverse_phi=True, anomalous=params.anomalous,
                                              spot_range=params.frames_for_index, minimum=False,
                                              crystal_symmetry=xs,
                                              integrate_nimages=params.integrate_nimages)
        ofs = open(os.path.join(wd, "XDS.INP"), "w")
        ofs.write(xdsinp_str)

        if spotstat != []:
            make_shikalog(spotstat, open(os.path.join(wd, "shika.log"), "w"))
# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    
    if "-h" in args: params.help = True

    if params.help:
        print("Parameter syntax:\n")
        iotbx.phil.parse(master_params_str).show(prefix="  ")
        print()
        print("Usage: prep_xds_runs.py logobjects.pkl [min_spots_for_ends=25] [workdir=run1]")
        quit()

    for arg in args:
        if os.path.isdir(arg) and params.workdir is None:
            params.workdir = arg
        if params.pklin is None and os.path.isfile(arg) and arg.endswith(".pkl"):
            params.pklin = arg

    if None in (params.pklin,):
        print("Provide single .pkl file")
        quit()

    if params.workdir is None:
        params.workdir = os.getcwd()

    run(params)

    print()
    print("Want to start all XDS runs?")
    print("run_all_xds_simple.py %s" % params.workdir)
