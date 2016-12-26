#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

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

def run_xtriage_in_module_if_possible(args,  wdir):
    try:
        from mmtbx.scaling import xtriage
        cwd_org = os.getcwd()
        try:
            os.chdir(wdir)
            xtriage.run(args, command_name="mmtbx.xtriage")
        except:
            print traceback.format_exc()
        finally:
            os.chdir(cwd_org)
    except ImportError:
        call("phenix.xtriage", arg=" ".join(args),
             stdin=None, stdout=sys.stdout, wdir=wdir)
# run_xtriage()

def unique(mtzin, mtzout, wdir, logout):
    ##
    # unique -> cad -> mtzutil (EXCLUDE FUNI SIGFUNI)
    #
    m = mtzutil.MtzFile(os.path.join(wdir,mtzin))
    cell = m.get_cell_str()
    sg = m.get_spacegroup()[1]
    resol = min(m.get_resolution())

    call(cmd="unique",
         arg="hklout unique.mtz",
         wdir=wdir,
         stdin="CELL %s\nSYMMETRY '%s'\nLABOUT F=FUNI SIGF=SIGFUNI\nRESOLUTION %f" % (cell, sg, resol),
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
         stdin="EXCLUDE FUNI SIGFUNI\nRUN\n",
         expects_in=["unique_cad.mtz"],
         expects_out=[mtzout],
         stdout=logout,
         )

    os.remove(os.path.join(wdir, "unique.mtz"))
    os.remove(os.path.join(wdir, "unique_cad.mtz"))

# unique()

def xds2mtz_normal(refl, mtzout, sg, wavelen, logout, use_ctruncate=False, dmin=None, dmax=None):
    wdir = os.path.dirname(mtzout)

    if not os.path.exists(os.path.join(wdir, "original")):
        os.symlink(refl, os.path.join(wdir, "original"))


    ##
    # prepare XDSCONV.INP and run
    #

    # for I
    print "generating MTZ for IMEAN,SIGIMEAN"

    ofs = open(os.path.join(wdir, "XDSCONV.INP"), "w")
    ofs.write("OUTPUT_FILE=tmp.hkl CCP4_I\n")
    ofs.write("INPUT_FILE=original\n")
    ofs.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.0\n")
    ofs.write("WILSON_STATISTICS= TRUE\n")
    if None not in (dmin, dmax):
        ofs.write("INCLUDE_RESOLUTION_RANGE= %s %s\n" % (dmax, dmin))

    ofs.close()

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
    print "generating MTZ for FP,SIGFP"
    ctruncate_ok = False
    if use_ctruncate:
        call(cmd="ctruncate -hklin CCP4_I.mtz -hklout ctruncate.mtz -colin '/*/*/[IMEAN,SIGIMEAN]'",
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
    symmetry %s
    end
    """%(wavelen,sg),
                 wdir=wdir,
                 expects_in=["ctruncate.mtz"],
                 expects_out=["CCP4_FI.mtz"],
                 stdout=logout
                 )
        else:
            print "Ctruncate failed. Use xdsconv instead."

    if not ctruncate_ok:
        ofs = open(os.path.join(wdir, "XDSCONV.INP"), "w")
        ofs.write("OUTPUT_FILE=tmp.hkl CCP4_F\n")
        ofs.write("INPUT_FILE=original\n")
        ofs.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.0\n")
        ofs.write("WILSON_STATISTICS= TRUE\n")
        if None not in (dmin, dmax):
            ofs.write("INCLUDE_RESOLUTION_RANGE= %s %s\n" % (dmax, dmin))

        ofs.close()

        call(cmd="xdsconv",
             wdir=wdir,
             expects_in=["original"],
             expects_out=["F2MTZ.INP", "tmp.hkl"],
             stdout=logout
             )

        call(cmd="f2mtz",
             arg="hklout CCP4_F.mtz",
             stdin=open(os.path.join(wdir, "F2MTZ.INP")).read(),
             wdir=wdir,
             expects_in=["tmp.hkl"],
             expects_out=["CCP4_F.mtz"],
             stdout=logout
             )

        ##
        # CAD all mtz files
        print "concatenating MTZ files"

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
symmetry %s
end
"""%(wavelen,sg),
             wdir=wdir,
             expects_in=["CCP4_I.mtz", "CCP4_F.mtz"],
             expects_out=["CCP4_FI.mtz"],
             stdout=logout
             )

    ##
    # Generate all unique reflections
    print "Genrating all unique reflections"
    unique(mtzin="CCP4_FI.mtz", mtzout=os.path.basename(mtzout), wdir=wdir, logout=logout)


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


# xds2mtz_anom()

def xds2mtz_anom(refl, mtzout, sg, wavelen, logout, use_ctruncate=False, dmin=None, dmax=None):
    wdir = os.path.dirname(mtzout)

    if not os.path.exists(os.path.join(wdir, "original")):
        os.symlink(refl, os.path.join(wdir, "original"))


    ##
    # prepare XDSCONV.INP and run

    # for I(+), I(-), SIGI(+), SIGI(-)
    print "generating MTZ for I(+), I(-), SIGI(+), SIGI(-)"

    ofs = open(os.path.join(wdir, "XDSCONV.INP"), "w")
    ofs.write("OUTPUT_FILE=tmp.hkl CCP4_I\n")
    ofs.write("INPUT_FILE=original\n")
    ofs.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.0\n")
    ofs.write("WILSON_STATISTICS= TRUE\n")
    ofs.write("FRIEDEL'S_LAW= FALSE\n")
    if None not in (dmin, dmax):
        ofs.write("INCLUDE_RESOLUTION_RANGE= %s %s\n" % (dmax, dmin))

    ofs.close()

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


    # for F(+), F(-), SIGF(+), SIGF(-)
    print "generating MTZ for F(+), F(-), SIGF(+), SIGF(-)"
    ctruncate_ok = False
    if use_ctruncate:
        call(cmd="ctruncate -hklin CCP4_I.mtz -hklout ctruncate.mtz -colin '/*/*/[IMEAN,SIGIMEAN]' -colano '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]'",
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
    symmetry %s
    end
    """%(wavelen,sg),
                 wdir=wdir,
                 expects_in=["ctruncate.mtz"],
                 expects_out=["CCP4_FI.mtz"],
                 stdout=logout
                 )
        else:
            print "Ctruncate failed. Use xdsconv instead."

    if not ctruncate_ok:
        ofs = open(os.path.join(wdir, "XDSCONV.INP"), "w")
        ofs.write("OUTPUT_FILE=tmp.hkl CCP4_F\n")
        ofs.write("INPUT_FILE=original\n")
        ofs.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.0\n")
        ofs.write("WILSON_STATISTICS= TRUE\n")
        ofs.write("FRIEDEL'S_LAW= FALSE\n")
        if None not in (dmin, dmax):
            ofs.write("INCLUDE_RESOLUTION_RANGE= %s %s\n" % (dmax, dmin))

        ofs.close()

        call(cmd="xdsconv",
             wdir=wdir,
             expects_in=["original"],
             expects_out=["F2MTZ.INP", "tmp.hkl"],
             stdout=logout
             )

        call(cmd="f2mtz",
             arg="hklout CCP4_F.mtz",
             stdin=open(os.path.join(wdir, "F2MTZ.INP")).read(),
             wdir=wdir,
             expects_in=["tmp.hkl"],
             expects_out=["CCP4_F.mtz"],
             stdout=logout
             )


        # for DANO, ISYM
        print "generating MTZ for DANO, ISYM"

        ofs = open(os.path.join(wdir, "XDSCONV.INP"), "w")
        ofs.write("OUTPUT_FILE=tmp.hkl CCP4\n")
        ofs.write("INPUT_FILE=original\n")
        ofs.write("GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.0\n")
        ofs.write("WILSON_STATISTICS= TRUE\n")
        ofs.write("FRIEDEL'S_LAW= FALSE\n")
        ofs.close()

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

        ##
        # CAD all mtz files
        print "concatenating MTZ files"

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
symmetry %s
end
"""%(wavelen,sg),
             wdir=wdir,
             expects_in=["CCP4_I.mtz", "CCP4_F.mtz", "CCP4.mtz"],
             expects_out=["CCP4_FI.mtz"],
             stdout=logout
             )

    ##
    # Generate all unique reflections
    print "Genrating all unique reflections"
    unique(mtzin="CCP4_FI.mtz", mtzout=os.path.basename(mtzout), wdir=wdir, logout=logout)

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
        os.remove(os.path.join(wdir, "CCP4.mtz"))


# xds2mtz_anom()

def add_multi(xds_file, workmtz, dmin=None, dmax=None, force_anomalous=False):
    from yamtbx.dataproc.xds import xds_ascii
    import iotbx.mtz

    print "Adding multiplicity for each reflection"

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

def xds2mtz(xds_file, dir_name, hklout=None, run_xtriage=False, run_ctruncate=False, dmin=None, dmax=None, force_anomalous=False, with_multiplicity=False, sgnum=None):
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
        for k, v in headers:
            if k == "FRIEDEL'S_LAW":
                header["FRIEDEL'S_LAW"] = v.strip()
            if k == "SPACE_GROUP_NUMBER":
                header["SPACE_GROUP_NUMBER"] = v.strip()
            if k == "X-RAY_WAVELENGTH":
                header["X-RAY_WAVELENGTH"] = v.strip() # XXX could be wrong if XSCALE result

    if force_anomalous:
        header["FRIEDEL'S_LAW"] = "FALSE"

    if sgnum is not None:
        header["SPACE_GROUP_NUMBER"] = str(sgnum)

    # make output directory
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)

    logout = open(os.path.join(dir_name, "xds2mtz.log"), "w")
    print >>logout, "xds2mtz.py running in %s" % os.getcwd()
    print >>logout, "output directory: %s" % dir_name
    print >>logout, "original file: %s" % xds_file

    print "Header information read from", xds_file

    for k in header:
        print k, "=", header[k]

    print

    ##
    # convert to MTZ
    print >>logout, "FRIEDEL'S_LAW= %s" % header["FRIEDEL'S_LAW"]
    print >>logout, ""
    logout.flush()

    if header["FRIEDEL'S_LAW"] == "TRUE":
        print xds_file, "is not anomalous dataset."
        print

        xds2mtz_normal(xds_file,
                       mtzout=os.path.join(dir_name, hklout),
                       sg=header["SPACE_GROUP_NUMBER"],
                       wavelen=header.get("X-RAY_WAVELENGTH","0"),
                       logout=logout,
                       use_ctruncate=run_ctruncate,
                       dmin=dmin, dmax=dmax)
        if run_xtriage:
            print "Running xtriage.."
            run_xtriage_in_module_if_possible(args=[hklout], wdir=dir_name)

        if with_multiplicity:
            add_multi(xds_file, os.path.join(dir_name, hklout),
                      dmin=dmin, dmax=dmax, force_anomalous=False)

    else:
        print xds_file, "is anomalous dataset."
        print
        xds2mtz_anom(xds_file,
                     mtzout=os.path.join(dir_name, hklout),
                     sg=header["SPACE_GROUP_NUMBER"],
                     wavelen=header.get("X-RAY_WAVELENGTH","0"),
                     logout=logout,
                     use_ctruncate=run_ctruncate,
                     dmin=dmin, dmax=dmax)
        if run_xtriage:
            print "Running xtriage.."
            run_xtriage_in_module_if_possible(args=[hklout, 'input.xray_data.obs_labels="I(+),SIGI(+),I(-),SIGI(-)"'],
                        wdir=dir_name)

        if with_multiplicity:
            add_multi(xds_file, os.path.join(dir_name, hklout),
                      dmin=dmin, dmax=dmax, force_anomalous=True)

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

    (opts, args) = parser.parse_args(sys.argv)

    ##
    # the default input file name is "XDS_ASCII.HKL"
    #

    if len(args) < 2:
        xds_file = os.path.abspath("XDS_ASCII.HKL")
    else:
        xds_file = os.path.abspath(args[1])

    if not os.path.isfile(xds_file):
        print "Cannot open", xds_file
        sys.exit(1)

    xds2mtz(xds_file, dir_name=opts.dir,
            run_xtriage=opts.run_xtriage, run_ctruncate=opts.run_ctruncate,
            dmin=opts.dmin, dmax=opts.dmax, force_anomalous=opts.anomalous,
            with_multiplicity=opts.make_mtzmulti)

