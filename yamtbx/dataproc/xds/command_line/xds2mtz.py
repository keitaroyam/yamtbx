#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

"""
xds2mtz.py : convert XDS_ASCII.HKL or xscale output to mtz file. MTZ file will include all possible columns user would need.

When FRIDEL'S_LAW= FALSE,
 MTZ columns will be F(+), F(-), I(+), I(-), IMEAN, FP, DANO, ISYM

When FRIDEL'S_LAW= TRUE,
 MTZ columns will be IMEAN, FP

With -x option, phenix.xtriage is executed automatically.
With -t option, ctruncate is used for converting I to F.
"""

import sys, os, optparse, subprocess, re
import traceback
from yamtbx.util import mtzutil
from yamtbx.util import call
from yamtbx.dataproc.xds import re_xds_kwd
from yamtbx.dataproc.command_line import copy_free_R_flag
from yamtbx.dataproc.command_line import create_free_R_flag
import iotbx.mtz
from cctbx import sgtbx

def run_xtriage_in_module_if_possible(args,  wdir):
    try:
        from mmtbx.scaling import xtriage
        cwd_org = os.getcwd()
        try:
            os.chdir(wdir)
            xtriage.run(args, command_name="mmtbx.xtriage")
        except:
            print(traceback.format_exc())
        finally:
            os.chdir(cwd_org)
    except ImportError:
        call("phenix.xtriage", arg=" ".join(args),
             stdin=None, stdout=sys.stdout, wdir=wdir)
# run_xtriage()

def unique(mtzin, mtzout, sg, wdir, logout):
    ##
    # unique -> cad -> mtzutil (EXCLUDE FUNI SIGFUNI)
    #
    m = mtzutil.MtzFile(os.path.join(wdir,mtzin))
    cell = m.get_cell_str()
    sg_org = m.get_spacegroup()[1].decode()
    resol = min(m.get_resolution())

    call(cmd="unique",
         arg="hklout unique.mtz",
         wdir=wdir,
         stdin="CELL %s\nSYMMETRY '%s'\nLABOUT F=FUNI SIGF=SIGFUNI\nRESOLUTION %f" % (cell, sg_org, resol),
         expects_in=[],
         expects_out=["unique.mtz"],
         stdout=logout,
         )

    call(cmd="cad",
         arg="hklin1 %s hklin2 %s hklout %s" %(mtzin, "unique.mtz", "unique_cad.mtz"),
         wdir=wdir,
         stdin="labin file 1 all\nlabin file 2 all\nend\n",
         expects_in=[mtzin, "unique.mtz"],
         expects_out=["unique_cad.mtz"],
         stdout=logout,
         )

    call(cmd="mtzutils",
         arg="hklin %s hklout %s" % ("unique_cad.mtz", mtzout),
         wdir=wdir,
         stdin="EXCLUDE FUNI SIGFUNI\nSYMMETRY %s\nRUN\n" % sg,
         expects_in=["unique_cad.mtz"],
         expects_out=[mtzout],
         stdout=logout,
         )

    os.remove(os.path.join(wdir, "unique.mtz"))
    os.remove(os.path.join(wdir, "unique_cad.mtz"))

# unique()

def prep_xdsconv_inp(wdir, hklin, out_type, anom_flag, dmin, dmax):
    ofs = open(os.path.join(wdir, "XDSCONV.INP"), "w")
    ofs.write("OUTPUT_FILE=tmp.hkl %s\n" % out_type)
    ofs.write("INPUT_FILE=%s\n" % hklin)
    ofs.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.0\n")
    ofs.write("WILSON_STATISTICS= TRUE\n")
    ofs.write("FRIEDEL'S_LAW= %s\n" % ("FALSE" if anom_flag else "TRUE"))

    if None not in (dmin, dmax):
        ofs.write("INCLUDE_RESOLUTION_RANGE= %s %s\n" % (dmax, dmin))

    ofs.close()
# prep_xdsconv_inp()

def xds2mtz_work(refl, mtzout, sg, wavelen, logout, anom_flag, use_ctruncate=False, dmin=None, dmax=None):
    wdir = os.path.dirname(mtzout)

    if not os.path.exists(os.path.join(wdir, "original")):
        os.symlink(refl, os.path.join(wdir, "original"))

    ##
    # prepare XDSCONV.INP and run
    #

    # for I
    logout.write("generating MTZ for %s\n" % ("I(+),I(-),SIGI(+),SIGI(-)" if anom_flag else "IMEAN,SIGIMEAN"))
    logout.flush()

    prep_xdsconv_inp(wdir, "original", "CCP4_I", anom_flag, dmin, dmax)

    call(cmd="xdsconv",
         wdir=wdir,
         expects_in=["original"],
         expects_out=["F2MTZ.INP", "tmp.hkl"],
         stdout=logout
         )

    call(cmd="f2mtz",
         arg="hklout CCP4_I.mtz",
         stdin=open(os.path.join(wdir, "F2MTZ.INP")).read(),
         wdir=wdir,
         expects_in=["tmp.hkl"],
         expects_out=["CCP4_I.mtz"],
         stdout=logout
         )

    # for F
    logout.write("generating MTZ for %s\n" % ("F(+),F(-),SIGF(+),SIGF(-)" if anom_flag else "F,SIGF"))
    logout.flush()
    ctruncate_ok = False
    if use_ctruncate:
        cmd_add = " -colano '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]'" if anom_flag else ""
        call(cmd="ctruncate -hklin CCP4_I.mtz -hklout ctruncate.mtz -colin '/*/*/[IMEAN,SIGIMEAN]'"+cmd_add,
             wdir=wdir,
             expects_in=["CCP4_I.mtz"],
             stdout=open(os.path.join(wdir, "ctruncate.log"), "w")
             )

        ctruncate_ok = os.path.isfile(os.path.join(wdir, "ctruncate.mtz"))

        if ctruncate_ok:
            call(cmd="cad",
                 arg="hklin1 ctruncate.mtz hklout CCP4_FI.mtz",
                 stdin="""\
labin file 1 all
xname file 1 ALL=XDS
dname file 1 ALL=XDS
dwavelength file 1 XDS XDS %s
end
    """%(wavelen),
                 wdir=wdir,
                 expects_in=["ctruncate.mtz"],
                 expects_out=["CCP4_FI.mtz"],
                 stdout=logout
                 )
        else:
            logout.write("Ctruncate failed. Use xdsconv instead.\n")
            logout.flush()

    if not ctruncate_ok:
        prep_xdsconv_inp(wdir, "original", "CCP4_F", anom_flag, dmin, dmax)

        call(cmd="xdsconv",
             wdir=wdir,
             expects_in=["original"],
             expects_out=["F2MTZ.INP", "tmp.hkl"],
             stdout=logout
             )

        call(cmd="f2mtz",
             arg="hklout CCP4_F.mtz",
             stdin=open(os.path.join(wdir, "F2MTZ.INP")).read().replace("FP","F"), # for ctruncate-comatibility
             wdir=wdir,
             expects_in=["tmp.hkl"],
             expects_out=["CCP4_F.mtz"],
             stdout=logout
             )

        if anom_flag:
            # for DANO, ISYM
            logout.write("generating MTZ for DANO, ISYM\n")
            logout.flush()

            
            prep_xdsconv_inp(wdir, "original", "CCP4", anom_flag, dmin, dmax)

            call(cmd="xdsconv",
                 wdir=wdir,
                 expects_in=["original"],
                 expects_out=["F2MTZ.INP", "tmp.hkl"],
                 stdout=logout
                 )

            call(cmd="f2mtz",
                 arg="hklout CCP4.mtz",
                 stdin=open(os.path.join(wdir, "F2MTZ.INP")).read(),
                 wdir=wdir,
                 expects_in=["tmp.hkl"],
                 expects_out=["CCP4.mtz"],
                 stdout=logout
                 )
            
            # CAD all mtz files
            logout.write("concatenating MTZ files\n")
            logout.flush()
        
            call(cmd="cad",
                 arg="hklin1 CCP4_I.mtz hklin2 CCP4_F.mtz hklin3 CCP4.mtz hklout CCP4_FI.mtz",
                 stdin="""\
labin file 1 all
labin file 2 all
labin file 3 E1=DANO E2=SIGDANO E3=ISYM
xname file 1 ALL=XDS
xname file 2 ALL=XDS
xname file 3 ALL=XDS
dname file 1 ALL=XDS
dname file 2 ALL=XDS
dname file 3 ALL=XDS
dwavelength file 1 XDS XDS %s
end
"""%(wavelen),
                 wdir=wdir,
                 expects_in=["CCP4_I.mtz", "CCP4_F.mtz", "CCP4.mtz"],
                 expects_out=["CCP4_FI.mtz"],
                 stdout=logout
                 )
        else:
            ##
            # CAD all mtz files
            logout.write("concatenating MTZ files\n")
            logout.flush()
            
            call(cmd="cad",
                 arg="hklin1 CCP4_I.mtz hklin2 CCP4_F.mtz hklout CCP4_FI.mtz",
                 stdin="""\
labin file 1 all
labin file 2 all
xname file 1 ALL=XDS
xname file 2 ALL=XDS
dname file 1 ALL=XDS
dname file 2 ALL=XDS
dwavelength file 1 XDS XDS %s
end
"""%(wavelen),
                 wdir=wdir,
                 expects_in=["CCP4_I.mtz", "CCP4_F.mtz"],
                 expects_out=["CCP4_FI.mtz"],
                 stdout=logout
            )

    ##
    # Generate all unique reflections
    print("Genrating all unique reflections")
    unique(mtzin="CCP4_FI.mtz", mtzout=os.path.basename(mtzout), sg=sg, wdir=wdir, logout=logout)

    # remove files
    os.remove(os.path.join(wdir, "CCP4_I.mtz"))
    os.remove(os.path.join(wdir, "CCP4_FI.mtz"))
    os.remove(os.path.join(wdir, "tmp.hkl"))
    os.remove(os.path.join(wdir, "XDSCONV.INP"))
    os.remove(os.path.join(wdir, "XDSCONV.LP"))
    os.remove(os.path.join(wdir, "F2MTZ.INP"))
    os.remove(os.path.join(wdir, "original"))

    if ctruncate_ok:
        os.remove(os.path.join(wdir, "ctruncate.mtz"))
    else:
        os.remove(os.path.join(wdir, "CCP4_F.mtz"))
        if anom_flag: os.remove(os.path.join(wdir, "CCP4.mtz"))

# xds2mtz_work()


def add_multi(xds_file, workmtz, dmin=None, dmax=None, force_anomalous=False):
    from yamtbx.dataproc.xds import xds_ascii

    print("Adding multiplicity for each reflection")

    xac = xds_ascii.XDS_ASCII(xds_file)
    iobs = xac.i_obs(anomalous_flag=True if force_anomalous else None)
    iobs = iobs.resolution_filter(d_min=float(dmin) if dmin is not None else None,
                                  d_max=float(dmax) if dmax is not None else None)
    iobs = iobs.select(iobs.sigmas() > 0).map_to_asu()
    
    merge = iobs.merge_equivalents(use_internal_variance=False)
    array_merged = merge.array()
    reject_sel = (array_merged.data() < -3*array_merged.sigmas())
    #print " rejecting %d reflections (<<I>/sd(<I>)> < -3)" % reject_sel.count(True)
    iobs = iobs.delete_indices(other=array_merged.select(reject_sel))
    
    # merge again after rejection
    merge = iobs.merge_equivalents(use_internal_variance=False)

    mtz_object = iotbx.mtz.object(workmtz)
    crystals = mtz_object.crystals()
    crystals[-1].datasets()[-1].add_miller_array(miller_array=merge.redundancies(), column_root_label="MULTIPLICITY")
    mtz_object.write(file_name=workmtz)
# generate_multi()

def copy_test_flag(workmtz, flag_source, flag_name=None, flag_value=None, log_out=sys.stdout):
    import iotbx.file_reader

    f = iotbx.file_reader.any_file(flag_source, force_type="hkl", raise_sorry_if_errors=True)
    flag_array, flag_name = copy_free_R_flag.get_flag_array(f.file_server.miller_arrays, flag_name, log_out=log_out)
    print("Copying test flag from", file=log_out)
    flag_array.show_summary(log_out, " ")

    if flag_value is None:
        flag_scores = copy_free_R_flag.get_r_free_flags_scores(miller_arrays=[flag_array],
                                                               test_flag_value=flag_value)
        flag_value = flag_scores.test_flag_values[0]
        print(" Guessing flag number:", flag_value, file=log_out)

    copy_free_R_flag.copy_flag_to_mtz(flag_array, flag_name, flag_value,
                                      workmtz, workmtz, log_out)
# copy_test_flag()

def add_test_flag(workmtz, preferred_fraction=0.05, max_fraction=0.1, max_free=2000,
                  ccp4_style=True, use_lattice_symmetry=True, log_out=sys.stdout):
    print("Adding test flag to %s" % workmtz, file=log_out)
    mtz_object = iotbx.mtz.object(workmtz)

    nref = mtz_object.n_reflections()
    fraction = preferred_fraction
    print(" preferred fraction= %.4f" % preferred_fraction, file=log_out)
    print(" total number of reflections= %d" % nref, file=log_out)

    if preferred_fraction*nref < 500:
        fraction = min(500./nref, max_fraction)
        print("As test reflections are too few (%d < 500), fraction changed to %.4f" % (preferred_fraction*nref, fraction), file=log_out)
    elif preferred_fraction*nref > max_free:
        fraction = max_free/float(nref)
        print("As test reflections are too many (%d > 2000), fraction changed to %.4f" % (preferred_fraction*nref, fraction), file=log_out)

    print(file=log_out)
    create_free_R_flag.run(workmtz, workmtz, fraction, "FreeR_flag",
                           ccp4=ccp4_style, use_lattice_symmetry=use_lattice_symmetry, log_out=log_out)
# add_test_flag()

def xds2mtz(xds_file, dir_name, hklout=None, run_xtriage=False, run_ctruncate=False, dmin=None, dmax=None, force_anomalous=False, with_multiplicity=False, space_group=None, flag_source=None, add_flag=False):
    if hklout is None:
        hklout = os.path.splitext(os.path.basename(xds_file))[0] + ".mtz"

    # if output file already exists, exit.
    if os.path.isfile(os.path.join(dir_name, hklout)):
        raise Exception(os.path.join(dir_name, hklout), "already exists.")

    # read header
    header = {}
    for l in open(xds_file):
        if l.startswith("!END_OF_HEADER"):
            break

        headers = re_xds_kwd.findall(l[l.index("!")+1:])
        is_iset = headers and headers[0][0] == "ISET"
        for k, v in headers:
            if not is_iset:
                if k == "FRIEDEL'S_LAW":
                    header["FRIEDEL'S_LAW"] = v.strip()
                if k == "SPACE_GROUP_NUMBER":
                    header["SPACE_GROUP_NUMBER"] = v.strip()
            if k == "X-RAY_WAVELENGTH":
                header["X-RAY_WAVELENGTH"] = v.strip() # XXX could be wrong if XSCALE result

    anom_flag = header["FRIEDEL'S_LAW"] == "FALSE"
    if force_anomalous: anom_flag = True


    sginfo_org = sgtbx.space_group_info(header["SPACE_GROUP_NUMBER"])
    
    if space_group:
        sginfo = sgtbx.space_group_info(space_group)
    else:
        sginfo = sginfo_org
        
    # make output directory
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)

    print("Header information read from", xds_file)
    for k in header:
        print(k, "=", header[k])
    print()
    
    logout = open(os.path.join(dir_name, "xds2mtz.log"), "w")
    print("xds2mtz.py running in %s" % os.getcwd(), file=logout)
    print("output directory: %s" % dir_name, file=logout)
    print("original file: %s" % xds_file, file=logout)
    print("flag_source: %s" % flag_source, file=logout)
    print("space group: %s (SPACE_GROUP_NUMBER=%s; requested space_group=%s)" % (sginfo, header.get("SPACE_GROUP_NUMBER",""), space_group), file=logout)
    if sginfo_org.group().build_derived_reflection_intensity_group(False) != sginfo.group().build_derived_reflection_intensity_group(False):
        print("  WARNING!! specified space group is incompatible with original file (%s)." % sginfo_org, file=logout)
    print("anomalous: %s (FRIEDEL'S_LAW=%s force_anomalous=%s)" % (anom_flag, header["FRIEDEL'S_LAW"], force_anomalous), file=logout)
    print("", file=logout)
    logout.flush()

    ##
    # convert to MTZ
    xds2mtz_work(xds_file,
                 mtzout=os.path.join(dir_name, hklout),
                 sg=str(sginfo).replace(" ",""),
                 wavelen=header.get("X-RAY_WAVELENGTH","0"),
                 logout=logout,
                 anom_flag=anom_flag,
                 use_ctruncate=run_ctruncate,
                 dmin=dmin, dmax=dmax)

    if run_xtriage:
        print("Running xtriage..")
        args = [hklout]
        if anom_flag: args.append('input.xray_data.obs_labels="I(+),SIGI(+),I(-),SIGI(-)"')
        run_xtriage_in_module_if_possible(args=args, wdir=dir_name)

    if with_multiplicity:
        add_multi(xds_file, os.path.join(dir_name, hklout),
                  dmin=dmin, dmax=dmax, force_anomalous=anom_flag)

    if flag_source is not None:
        copy_test_flag(os.path.join(dir_name, hklout), flag_source, log_out=logout)
    elif add_flag:
        add_test_flag(os.path.join(dir_name, hklout), log_out=logout)

# xds2mtz()

if __name__ == "__main__":
    parser = optparse.OptionParser(usage="usage: %prog [options] [XDS_ASCII.HKL]")

    parser.add_option("--dir","-d", action="store", type=str, dest="dir", default="ccp4",
                      help="output directory")
    parser.add_option("--xtriage","-x", action="store_true", dest="run_xtriage", help="run phenix.xtriage")
    parser.add_option("--truncate","-t", action="store_true", dest="run_ctruncate", help="use ctruncate to estimate F")
    parser.add_option("--multiplicity","-m", action="store_true", dest="make_mtzmulti", help="Add multiplicity info")
    parser.add_option("--anomalous","-a", action="store_true", dest="anomalous", help="force anomalous")
    parser.add_option("--dmin", action="store", dest="dmin", help="high resolution cutoff") # as str
    parser.add_option("--dmax", action="store", dest="dmax", help="low resolution cutoff") # as str
    parser.add_option("--copy-test-flag","-r", action="store", dest="flag_source", help="")
    parser.add_option("--add-test-flag", action="store_true", dest="add_flag", help="")
    parser.add_option("--space-group", action="store", type=str, dest="sg", help="Space group number or name")

    (opts, args) = parser.parse_args(sys.argv)

    ##
    # the default input file name is "XDS_ASCII.HKL"
    #

    if len(args) < 2:
        xds_file = os.path.abspath("XDS_ASCII.HKL")
    else:
        xds_file = os.path.abspath(args[1])

    if not os.path.isfile(xds_file):
        print("Cannot open", xds_file)
        sys.exit(1)

    xds2mtz(xds_file, dir_name=opts.dir,
            run_xtriage=opts.run_xtriage, run_ctruncate=opts.run_ctruncate,
            dmin=opts.dmin, dmax=opts.dmax, force_anomalous=opts.anomalous,
            with_multiplicity=opts.make_mtzmulti,
            flag_source=opts.flag_source, add_flag=opts.add_flag,
            space_group=opts.sg)

