"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os
import subprocess

from yamtbx.dataproc.xds.command_line.xds_predict_mitai import check_needed_files
from yamtbx.dataproc.xds import files as xds_files
from yamtbx.dataproc.xds.idxreflp import SpotXds
from yamtbx.dataproc.xds import modify_xdsinp, make_backup, revert_files

def get_colspot_result(frame_ranges, wdir):
    assert all(map(lambda x: len(x)==2 and x[0]<=x[1], frame_ranges))

    # Check all needed files exist
    if not check_needed_files(xds_files.needed_by_COLSPOT, wdir):
        return

    xdsinp = os.path.join(wdir, "XDS.INP")

    backup_needed = xds_files.generated_by_COLSPOT + ("XDS.INP",)

    # 1. Backup XDS.INP, etc. (Make copies; not renaming)
    bk_prefix = make_backup(backup_needed, wdir=wdir)

    try:
        # 2. Modify XDS.INP
        spot_ranges = map(lambda x: ("SPOT_RANGE", "%d %d"%tuple(x)), frame_ranges)
        data_range = "%d %d" % (min(map(lambda x: min(x), frame_ranges)),
                                max(map(lambda x: max(x), frame_ranges)))
        modify_xdsinp(xdsinp, [("JOB", "COLSPOT"), ("DATA_RANGE", data_range)]+spot_ranges)

        # 3. Run xds
        p = subprocess.Popen("xds", cwd=wdir)
        p.wait()

        spotxds = SpotXds(os.path.join(wdir, "SPOT.XDS"))
            
    finally:
        # 6. Revert XDS.INP, etc.
        revert_files(backup_needed, bk_prefix, wdir=wdir)

    return spotxds
# get_colspot_result()
