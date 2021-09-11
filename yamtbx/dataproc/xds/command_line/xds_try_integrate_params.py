#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
Try several parameters in INTEGRATE step of XDS, and optionally run CORRECT and make mtz files which are ready for refinement.

Usage:
  PHENIX_TRUST_OTHER_ENV=1 phenix.python xds_try_integrate_params.py delphi=5,10,15,...
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals


master_params_str = """\
delphi = []
    .type = floats
all_possible_delphi = False
    .type = bool

all_refine_combinations = False
    .type = bool

run_correct = True
    .type = bool
    .help = Run CORRECT (and make mtz) or not.
mtz_free = None
    .type = path
"""

import os
import shutil
import itertools
import traceback
import iotbx.phil

from yamtbx.dataproc.xds import *
from yamtbx.util import call
from yamtbx.dataproc.xds.command_line import xds2mtz
from yamtbx.dataproc.xds import files
from .xds_try_scale_params import get_digit, copy_testflag

def make_all_refine_combinations():
    all_keywords = ["DISTANCE","BEAM","AXIS","ORIENTATION","CELL"]
    result = []
    for i in range(len(all_keywords)+1):
        for c in itertools.combinations(all_keywords, i):
            result.append(c)

    return result
# make_all_refine_combinations()

def get_delphi_defined(xdsinp):
    org_kwd = [x for x in get_xdsinp_keyword(xdsinp) if x[0].strip().startswith("DELPHI")]
    if len(org_kwd) == 0:
        return 5.
    else:
        return float(org_kwd[-1][1])
# get_delphi_defined()

def get_refine_defined(xdsinp):
    org_kwd = [x for x in get_xdsinp_keyword(xdsinp) if x[0].strip().startswith("REFINE(INTEGRATE)")]
    if len(org_kwd) == 0:
        return [("DISTANCE", "BEAM", "ORIENTATION", "CELL")]
    else:
        return [org_kwd[-1][1].strip().split()]
# get_refine_defined(xdsinp)

def run(params):
    xdsinp = "XDS.INP"
    workdir = os.getcwd()

    # concerning DELPHI=
    if len(params.delphi) == 0:
        params.delphi = [get_delphi_defined(xdsinp)]

    params.delphi.sort()
    digit_delphi = str(get_digit(params.delphi))

    # concerning REFINE(INTEGRATE)=
    if params.all_refine_combinations:
        refine_params = make_all_refine_combinations()
    else:
        refine_params = get_refine_defined(xdsinp)

    backup_needed = files.generated_by_INTEGRATE + ("XDS.INP",)
    if params.run_correct:
        backup_needed += files.generated_by_CORRECT

    bk_prefix = None

    try:
        for rp in refine_params:
            rp_str = "refine_none" if len(rp) == 0 else "refine_" + "+".join(rp)

            for delphi in params.delphi:
                delphi_str = ("delphi_%."+digit_delphi+"f") % delphi

                bk_prefix = make_backup(backup_needed) # Backup existing XDS.INP and others

                work_name = "INTEGRATE_%s_%s" % (rp_str, delphi_str)

                inp_params = [("JOB","INTEGRATE"),
                              ("DELPHI", delphi),
                              ("REFINE(INTEGRATE)", " ".join(rp))
                              ]
                if params.run_correct:
                    inp_params[0] = ("JOB","INTEGRATE CORRECT")

                modify_xdsinp(xdsinp, inp_params=inp_params)
                call("xds_par", stdout=sys.stdout)

                # Move files into a new directory.
                os.mkdir(work_name)
                for f in backup_needed:
                    shutil.move(f, work_name)

                if params.run_correct:
                    # make .mtz file
                    try:
                        call("xdsstat", stdin="\n", stdout=open(os.path.join(work_name, "XDSSTAT.LP"),"w"), wdir=work_name)

                        xds2mtz.xds2mtz(os.path.abspath(os.path.join(work_name, "XDS_ASCII.HKL")),
                                        dir_name=os.path.join(workdir, work_name, "ccp4"),
                                        run_ctruncate=False, run_xtriage=True
                                        )

                        if params.mtz_free is not None:
                            copy_testflag(mtzfree=params.mtz_free,
                                          mtzin=os.path.join(workdir, work_name, "ccp4", "XDS_ASCII.mtz"))
                    except:
                        print(traceback.format_exc())
                        print("Ignoring xds2mtz error..")
                        print()

                revert_files(backup_needed, bk_prefix) # Revert XDS.INP and others

    finally:
        # 6. Revert XDS.INP, etc.
        revert_files(backup_needed, bk_prefix)

# run()


if __name__ == "__main__":

    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()

    if params.mtz_free is not None:
        params.mtz_free = os.path.abspath(params.mtz_free)

    print("Parameters:")
    cmdline.work.format(python_object=params).show(out=sys.stdout, prefix=" ")

    #if len(params.delphi) == 1:
    #    print "Specify params!"
    #    sys.exit(1)

    run(params)
