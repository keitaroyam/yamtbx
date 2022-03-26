#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
"""
Helper script for setting up mtz files (?) for comparing results with scala

Usage:
  PHENIX_TRUST_OTHER_ENV=1 xds_setup_scala_exp.py

What we want to compare:
 1. XDS CORRECT
   a. all defaults
   b. change WFAC1?
   c.
 2. Scala after INTEGRATE
 3. Scala after CORRECT
   a. all defaults
   b. CORRECTIONS=            NBATCH= 1 MINIMUM_I/SIGMA= 50 WFAC1= 2
   c. CORRECTIONS= MODULATION NBATCH= 1
  How about REFINE(CORRECT)=???

CORRECT related parameters
  OVERLOAD=, SILICON=, SENSOR_THICKNESS=, MINIMUM_ZETA=, FRACTION_OF_POLARIZATION=, POLARIZATION_PLANE_NORMAL=
  AIR=, REFINE(CORRECT)=, MINPK=, WFAC1=, STRICT_ABSORPTION_CORRECTION=, PATCH_SHUTTER_PROBLEM=, NBATCH=, REJECT_ALIEN=

This script MUST be used after geometric parameters are converged. Do recycling GXPARM->XPARM; INTEGRATE->CORRECT !
  Scala's schemes:

Scheme.
 1. Backup XDS.INP, CORRECT.LP, XDS_ASCII.HKL (and all other CORRECT outputs?)
 2. Make INTEGRATE_forscale.mtz with pointless.
 3. Run CORRECT with all default parameters and make xds_CORRECT.mtz (with xdsconv), CORRECT_a_forscale.mtz (with pointless)
 4. Run CORRECT with setting b and make CORRECT_b_forscale.mtz (with pointless)
 5. Run CORRECT with setting c and make CORRECT_c_forscale.mtz (with pointless)
 6. Run Scala or Aimless for
   - INTEGRATE_forscale.mtz => aimless_INTEGRATE.mtz
   - CORRECT_a_forscale.mtz => aimless_CORRECT_a.mtz
   - CORRECT_b_forscale.mtz => aimless_CORRECT_b.mtz
   - CORRECT_c_forscale.mtz => aimless_CORRECT_c.mtz
"""

import sys, os, optparse, shutil, re, subprocess, random, glob, string
import iotbx.mtz
from cctbx import r_free_utils
from cctbx import miller
from iotbx.reflection_file_utils import get_r_free_flags_scores

from yamtbx.dataproc.xds import *
from yamtbx.util import call
from yamtbx.dataproc.xds.command_line import xds2mtz

backup_needed = ("XDS.INP", "XDS_ASCII.HKL", "CORRECT.LP")
needed_files = ("BKGPIX.cbf", "BLANK.cbf", "GAIN.cbf", "X-CORRECTIONS.cbf", "Y-CORRECTIONS.cbf",
                "XDS.INP", "XPARM.XDS")

def check_needed_files():
    not_exists = [ f for f in needed_files if not os.path.isfile(f) ]

    if not_exists != []:
        print("We need these files!")
        print("  " + ",".join(not_exists))
        print()

    return len(not_exists) == 0
# check_needed_files()

def is_anomalous(xdsinp):
    kwds = get_xdsinp_keyword(xdsinp)
    for kwd, val in kwds:
        if kwd == "FRIEDEL'S_LAW" and val.strip() == "FALSE":
            return True

    return False
# is_anomalous()

def copy_testflag(mtzfree, mtzin, mtzout):
    call(cmd="copy_free_R_flag.py -r %s -o %s %s" % (mtzfree, mtzout, mtzin))
# copy_testflag()

def run_aimless(prefix, wdir):
    script_body = """\
pref=%s
mtzin=../${pref}_forscale.mtz
scaled=${pref}_aimless.mtz
truncated=${pref}_aimless_truncate.mtz

#lowres=$2
#highres=$3

aimless \
 HKLIN $mtzin \
 HKLOUT $scaled \
 SCALES ${pref}.scala \
 ROGUES ${pref}_scala_rogues.log \
 NORMPLOT ${pref}_scala_normplot.xmgr \
 ANOMPLOT ${pref}_scala_anomplot.xmgr \
 PLOT ${pref}_scala_surface_plot.plt \
 CORRELPLOT ${pref}_scala_correlplot.xmgr \
 ROGUEPLOT ${pref}_scala_rogueplot.xmgr \
<<eoi
!resolution low $lowres high $highres
!exclude batch ** to **
bins 20
!run 1 batch 1 to 750
!scales rotation spacing 5.000000 bfactor on tails
!cycles 5
""" % (prefix)

    if anomalous:
        script_body += """\
anomalous on
eoi

ctruncate -hklin $scaled -hklout $truncated -colin "/*/*/[IMEAN,SIGIMEAN]" -colano "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]"
"""
    else:
        script_body += """\
anomalous off
eoi

ctruncate -hklin $scaled -hklout $truncated -colin "/*/*/[IMEAN,SIGIMEAN]"
"""
    open(os.path.join(wdir, "run_aimless.sh"),"w").write(script_body)

    call("sh run_aimless.sh",
         stdout=open(os.path.join(wdir, "aimless_truncate.log"), "w"),
         wdir=wdir)


if __name__ == "__main__":

    parser = optparse.OptionParser(usage="usage: %prog [options] pdbfile")
    #parser.add_option("--wdir", action="store", dest="wdir", default="see_predict", type=str,
    #                  help="Working directory")
    parser.add_option("--mtz-free","-r", action="store", dest="mtz_free", type=str, help="MTZ file for test flag")
    parser.add_option("--free-fraction","-f", action="store", dest="free_fraction", type=float, help="Test set fraction if newly set")


    (opts, args) = parser.parse_args(sys.argv)

    if opts.mtz_free is not None and opts.free_fraction is not None:
        print("mtz-free and free-fraction cannot be specified at the same time.")
        quit()
    if 0:
        parser.print_help()
        quit()

    # Check all needed files exist
    if not check_needed_files():
        sys.exit(1)

    xdsinp = "XDS.INP"

    workdir = "scale_exp"
    os.mkdir(workdir)
    os.mkdir(os.path.join(workdir, "scaled"))


    anomalous = is_anomalous(xdsinp)

    # *. INTEGRATE_forscale.mtz
    call("pointless -copy xdsin INTEGRATE.HKL hklout %s/INTEGRATE_forscale.mtz" % workdir)

    # 1. Backup XDS.INP, etc. (Make copies; not renaming)
    bk_prefix = make_backup(backup_needed)

    try:
        # CORRECT with a
        modify_xdsinp(xdsinp, inp_params=[("JOB","CORRECT"),
                                          ("CORRECTIONS", None),
                                          ("NBATCH", None),
                                          ("MINIMUM_I/SIGMA", None),
                                          ("WFAC1", None)
                                          ])

        call("xds", stdout=open("CORRECT_a_xds.log", "w"))

        # make xds_CORRECT.mtz
        xds2mtz.xds2mtz(os.path.abspath("XDS_ASCII.HKL"), dir_name=os.path.join(workdir, "xdsconv_CORRECT"),
                        hklout="xds_CORRECT.mtz")

        # make CORRECT_a_forscale.mtz
        call("pointless -copy xdsin XDS_ASCII.HKL hklout %s/CORRECT_a_forscale.mtz" % workdir)
        make_backup(backup_needed, "CORRECT_a_")

        # CORRECT with b
        #modify_xdsinp(xdsinp, mode="b")
        modify_xdsinp(xdsinp, inp_params=[("JOB","CORRECT"),
                                          ("CORRECTIONS", ""),
                                          ("NBATCH", "1"),
                                          ("MINIMUM_I/SIGMA", "50"),
                                          ("WFAC1", "2")
                                          ])

        call("xds", stdout=open("CORRECT_b_xds.log", "w"))
        call("pointless -copy xdsin XDS_ASCII.HKL hklout %s/CORRECT_b_forscale.mtz" % workdir)
        make_backup(backup_needed, "CORRECT_b_")
        modify_xdsinp(xdsinp, inp_params=[("JOB","CORRECT"),
                                          ("CORRECTIONS", "MODULATION"),
                                          ("NBATCH", "1"),
                                          ("MINIMUM_I/SIGMA", None),
                                          ("WFAC1", None)
                                          ])
        # CORRECT with c
        #modify_xdsinp(xdsinp, mode="c")

        call("xds", stdout=open("CORRECT_c_xds.log", "w"))
        call("pointless -copy xdsin XDS_ASCII.HKL hklout %s/CORRECT_c_forscale.mtz" % workdir)
        make_backup(backup_needed, "CORRECT_c_")

    finally:
        # 6. Revert XDS.INP, etc.
        revert_files(backup_needed, bk_prefix)


    # Run aimless
    print("Running aimless")
    for prefix in ("INTEGRATE", "CORRECT_a", "CORRECT_b", "CORRECT_c"):
        print("running aimless for", prefix)
        wdir = os.path.join(workdir, "aimless_%s"%prefix)
        os.mkdir(wdir)
        run_aimless(prefix, wdir)

        truncate_mtz = os.path.join(wdir, prefix+"_aimless_truncate.mtz")
        if not os.path.exists(truncate_mtz):
            print("WARNING: mtz file was not created:", truncate_mtz)
            continue

        if opts.mtz_free:
            copy_testflag(mtzfree=opts.mtz_free,
                          mtzin=truncate_mtz,
                          mtzout=os.path.join(workdir, "scaled", prefix+"_aimless_truncate.mtz"))
        else:
                shutil.copy(truncate_mtz,
                            os.path.join(workdir, "scaled"))



    if opts.mtz_free:
        copy_testflag(mtzfree=opts.mtz_free,
                      mtzin=os.path.join(workdir, "xdsconv_CORRECT", "xds_CORRECT.mtz"),
                      mtzout=os.path.join(workdir, "scaled", "xds_CORRECT.mtz"))
    else:
        shutil.copy(os.path.join(os.path.join(workdir, "xdsconv_CORRECT", "xds_CORRECT.mtz")),
                    os.path.join(workdir, "scaled"))
