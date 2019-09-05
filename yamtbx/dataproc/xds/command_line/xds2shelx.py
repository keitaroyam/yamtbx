import os
from yamtbx.util import call
from yamtbx.util import xtal
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from cctbx import sgtbx

def xds2shelx(xds_file, dir_name, prefix=None, dmin=None, dmax=None, force_anomalous=False, space_group=None, flag_source=None, add_flag=False):
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(xds_file))[0]

    hklout = prefix + ".hkl"
    
    # if output file already exists, exit.
    if os.path.isfile(os.path.join(dir_name, hklout)):
        raise Exception(os.path.join(dir_name, hklout), "already exists.")

    # read header
    xac = XDS_ASCII(xds_file, read_data=False)

    wavelength = xac.wavelength
    if xac.wavelength is None and xac.input_files:
        wavelength = float(xac.input_files.values()[0][1])
    else:
        wavelength = 1.0

    anom_flag = xac.anomalous
    if force_anomalous: anom_flag = True

    sginfo_org = xac.symm.space_group_info()
    
    if space_group:
        sginfo = sgtbx.space_group_info(space_group)
    else:
        sginfo = sginfo_org
        
    # make output directory
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)

    logout = open(os.path.join(dir_name, "xds2shelx.log"), "w")
    print >>logout, "xds2shelx.py running in %s" % os.getcwd()
    print >>logout, "output directory: %s" % dir_name
    print >>logout, "original file: %s" % xds_file
    print >>logout, "flag_source: %s" % flag_source
    print >>logout, "space group: %s (original=%s, requested space_group=%s)" % (sginfo, sginfo_org, space_group)
    if sginfo_org.group().build_derived_reflection_intensity_group(False) != sginfo.group().build_derived_reflection_intensity_group(False):
        print >>logout, "  WARNING!! specified space group is incompatible with original file (%s)." % sginfo_org
    print >>logout, "anomalous: %s (original=%s force_anomalous=%s)" % (anom_flag, xac.anomalous, force_anomalous)
    print >>logout, ""
    logout.flush()

    ##
    if not os.path.exists(os.path.join(dir_name, "original")):
        os.symlink(xds_file, os.path.join(dir_name, "original"))

    ##
    # prepare XDSCONV.INP and run
    #

    with open(os.path.join(dir_name, "XDSCONV.INP"), "w") as ofs:
        ofs.write("OUTPUT_FILE=%s SHELX\n" % hklout)
        ofs.write("INPUT_FILE=original\n")
        ofs.write("MERGE= FALSE\n")
        ofs.write("FRIEDEL'S_LAW= %s\n" % ("FALSE" if anom_flag else "TRUE"))

        if None not in (dmin, dmax):
            ofs.write("INCLUDE_RESOLUTION_RANGE= %s %s\n" % (dmax, dmin))


    call(cmd="xdsconv",
         wdir=dir_name,
         expects_in=["original"],
         expects_out=[hklout],
         stdout=logout
         )

    acgrp = sginfo.group().build_derived_acentric_group()
    
    cell_str = xtal.format_unit_cell(xac.symm.unit_cell(), lfmt="%8.4f", afmt="%7.3f")
    with open(os.path.join(dir_name, "%s.ins"%prefix), "w") as ofs:
        ofs.write("CELL %.4f %s\n" % (wavelength, cell_str))
        ofs.write("ZERR 1 0 0 0 0 0 0\n")
        ofs.write("LATT %s\n" % xtal.shelx_latt(sginfo.group()))
        for op in acgrp.all_ops():
            if op.is_unit_mx(): continue
            ofs.write("SYMM %s\n" % op.as_xyz(decimal=True, t_first=True, symbol_letters="XYZ"))
        ofs.write("SFAC C N O S\n")
        ofs.write("UNIT 6 6 6 6\n")
        ofs.write("FIND 10\n") # TODO more intelligent
        ofs.write("NTRY 1000\n")
        ofs.write("HKLF 4\n")
        ofs.write("END\n")

    #if flag_source is not None:
    #    copy_test_flag(os.path.join(dir_name, hklout), flag_source, log_out=logout)
    #elif add_flag:
    #    add_test_flag(os.path.join(dir_name, hklout), log_out=logout)


if __name__ == "__main__":
    import sys
    import optparse
    parser = optparse.OptionParser(usage="usage: %prog [options] [XDS_ASCII.HKL]")

    parser.add_option("--dir","-d", action="store", type=str, dest="dir", default="shelx",
                      help="output directory")
    parser.add_option("--anomalous","-a", action="store_true", dest="anomalous", help="force anomalous")
    parser.add_option("--dmin", action="store", dest="dmin", help="high resolution cutoff") # as str
    parser.add_option("--dmax", action="store", dest="dmax", help="low resolution cutoff") # as str
    #parser.add_option("--copy-test-flag","-r", action="store", dest="flag_source", help="")
    #parser.add_option("--add-test-flag", action="store_true", dest="add_flag", help="")
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
        print "Cannot open", xds_file
        sys.exit(1)

    xds2shelx(xds_file, dir_name=opts.dir,
              dmin=opts.dmin, dmax=opts.dmax, force_anomalous=opts.anomalous,
              #flag_source=opts.flag_source, add_flag=opts.add_flag,
              space_group=opts.sg)

