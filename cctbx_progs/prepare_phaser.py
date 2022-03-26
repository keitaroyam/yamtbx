#!/usr/bin/env phenix.python

from __future__ import print_function
from __future__ import unicode_literals
import sys, os, optparse
import iotbx.mtz
from iotbx.reflection_file_editor import is_rfree_array
from iotbx.reflection_file_utils import get_r_free_flags_scores

if __name__ == "__main__":
    parser = optparse.OptionParser(usage="usage: %prog [options] mtzfile pdbfiles [seqfile]",
                                   description="Prepare script for phaser.")

    (opts, args) = parser.parse_args(sys.argv)

    pdbfiles, mtzfile = [], None

    for arg in args[1:]:
        if arg.endswith(".pdb"):
            pdbfiles.append(arg)
        elif arg.endswith(".mtz"):
            mtzfile = arg

    if mtzfile is None or len(pdbfiles)==0:
        print("Provide pdb and mtz")
        quit()

    arrays = iotbx.mtz.object(mtzfile).as_miller_arrays()

    col_target = []

    for array in arrays:
        labs = tuple(array.info().labels)

        if array.is_xray_intensity_array():
            """
            if "(+)" in array.info().label_string(): #array.anomalous_flag():
                col_target.append("I+=%s SIGI+=%s I-=%s SIGI-=%s" % tuple(array.info().labels))
            else:
                if len(labs) < 2:
                    col_target.append("IP=%s" % labs)
                else:
                    col_target.append("IP=%s SIGIP=%s" % labs)
                    """
            pass
        elif array.is_xray_amplitude_array():
            if array.sigmas() is None:
                print("Not having sigma:", array.info().labels)
                continue
            if array.anomalous_flag():
                #if "+" not in array.info().label_string():
                #    print "non-supported anomalous array:", array.info().labels
                #    continue
                #col_target.append("F+=%s SIGF+=%s F-=%s SIGF-=%s" % tuple(array.info().labels))
                pass
            else:
                col_target.append("F=%s SIGF=%s" % tuple(array.info().labels))

        else:
            print("Unrecognized:", array.info().labels)

    if os.path.exists("run_phaser.sh"):
        print()
        print("Error: Script file already exists. aborting.")
        sys.exit(1)

    sh_str = """\
pref=mr_1
phenix.phaser <<+ > phaser_${pref}.log
MODE MR_AUTO
HKLIN %(mtzfile)s
""" % dict(mtzfile=mtzfile)

    for trg in col_target:
        sh_str += "labin " + trg + "\n"

    str_search = ""
    for i, pdb in enumerate(pdbfiles):
        sh_str += "ENSEMBLE ens%d PDBFILE %s IDENTITY 0.3\n" % (i+1, pdb)
        str_search += "SEARCH ENSEMBLE ens%d NUM 1\n" % (i+1)

    sh_str += """\
!COMPOSITION PROTEIN MW 58000 NUM 4
%s
ROOT $pref
+
""" % (str_search)

    open("run_phaser.sh", "w").write(sh_str)

