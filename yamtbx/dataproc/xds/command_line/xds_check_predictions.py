#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
Compare located spots and predicted positions.

TODO: for fine-slicing, frame-width parameter should be prepared.
"""
from __future__ import print_function
from __future__ import unicode_literals

from yamtbx.dataproc.xds.command_line import xds_predict_mitai
from yamtbx.dataproc.xds.colspot import get_colspot_result
from yamtbx.dataproc.xds import integrate_hkl_as_flex
from yamtbx.dataproc.xds import files as xds_files
from yamtbx.dataproc.xds import make_backup, revert_files, modify_xdsinp
import iotbx.phil

import scipy.spatial
import numpy
import os
import shutil
import subprocess

master_params_str = """
xdsdir = None
 .type = path
spotfinder = *xds shika
 .type = choice(multi=False)
frames = []
 .type = ints
 .help = frame number to check predictions
d_min = 5
 .type = float
 .help = High resolution limit for check
distance_limit_in_px = 5
 .type = float
 .help = maximum distance to match spot locations
sigmar = 0.7
 .type = float
sigmab = 0.03
 .type = float
"""

def calc_matches(numer, denom, dist_limit, out):
    kdtree = scipy.spatial.cKDTree(numpy.array(numer))
    nmatch = 0
    for pred in denom:
        dist, idx = kdtree.query(pred, k=1, p=2)
        if idx==kdtree.n:
            continue
        if dist < dist_limit:
            print("%.2f %.2f "%tuple(numer[idx]), file=out)
            nmatch += 1
    return nmatch
# calc_matches()

def run(params, out):
    print("Frames:", params.frames, file=out)
    backup_needed = xds_files.generated_by_DEFPIX + ("XDS.INP","BKGINIT.cbf",)
    bk_prefix = make_backup(backup_needed, wdir=params.xdsdir)

    ret = {} # {frame: [matches, spots, predicted]}

    try:
        # run DEFPIX to limit resolution.
        modify_xdsinp(os.path.join(params.xdsdir, "XDS.INP"),
                      [("JOB", "DEFPIX"), ("INCLUDE_RESOLUTION_RANGE", "50 %.2f"%params.d_min)])

        p = subprocess.Popen("xds", cwd=params.xdsdir)
        p.wait()

        # copy BKGPIX.cbf -> BKGINIT.cbf (for COLSPOT)
        shutil.copyfile(os.path.join(params.xdsdir, "BKGPIX.cbf"),
                        os.path.join(params.xdsdir, "BKGINIT.cbf"))
        
        for frame in params.frames:
            print("Frame %d" % frame, file=out)
            print("====================\n", file=out)
            # search spots
            if params.spotfinder == "xds":
                spotxds = get_colspot_result(frame_ranges=[[frame, frame],], wdir=params.xdsdir)
                spots = [x[:2] for x in spotxds.collected_spots(with_resolution=False)]
            else:
                raise "Sorry!"

            # run INTEGRATE to collect predicted coords
            integrate_results = xds_predict_mitai.run(param_source=os.path.join(params.xdsdir, "XPARM.XDS"), 
                                                      frame_num=frame,
                                                      wdir=params.xdsdir, need_adx=False,
                                                      sigmar=params.sigmar, sigmab=params.sigmab)

            # read predicted coords
            tmp = [x for x in integrate_results if x.endswith(".HKL")]
            if len(tmp) == 0:
                print("Integration failed!", file=out)
                ret[frame] = (0, len(spots), 0)
                continue

            integrate_hkl = tmp[0]
            cols = integrate_hkl_as_flex.reader(integrate_hkl, [], False).get_column_names()
            i_xcal, i_ycal = cols.index("XCAL"), cols.index("YCAL")
            predicted = []

            for l in open(integrate_hkl):
                if l.startswith("!"): continue
                sp = l.split()
                predicted.append(list(map(float, (sp[i_xcal], sp[i_ycal]))))

            # compare them
            nmatch = calc_matches(spots, predicted, params.distance_limit_in_px,
                                  open(os.path.join(params.xdsdir, "matched_predicted_%.4d.adx"%frame), "w"))
            #nmatch = calc_matches(predicted, spots, params.distance_limit_in_px,
            #                      open(os.path.join(params.xdsdir, "matched_located_%.4d.adx"%frame), "w"))

            ret[frame] = (nmatch, len(spots), len(predicted))

    finally:
        revert_files(backup_needed, bk_prefix, wdir=params.xdsdir)

    print(file=out)
    for frame in sorted(ret):
        nmatch, nspots, npredicted = ret[frame]
        print("Frame %4d Located/Predicted: %d/%d= %.2f%%" % (frame, nmatch, npredicted,
                                                                     100.*float(nmatch)/npredicted if npredicted>0 else float("nan")), file=out)
        print("Frame %4d Predicted/Located: %d/%d= %.2f%%" % (frame, nmatch, nspots,
                                                                     100.*float(nmatch)/nspots if nspots>0 else float("nan")), file=out)
        print(file=out)


    return ret
# run()    

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if os.path.isdir(arg) and params.xdsdir is None:
            params.xdsdir = arg

    if params.xdsdir is None:
        params.xdsdir = "."

    if len(params.frames) == 0:
        cmdline.work.format(python_object=params).show(out=sys.stdout, prefix=" ")
        print()
        print("Give frames=.")
        quit()

    run(params, sys.stdout)
