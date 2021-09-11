#!/usr/bin/env phenix.python

from __future__ import print_function
from __future__ import unicode_literals
import sys, os, optparse
import iotbx.mtz
from iotbx.reflection_file_editor import is_rfree_array
from iotbx.reflection_file_utils import get_r_free_flags_scores

if __name__ == "__main__":
    parser = optparse.OptionParser(usage="usage: %prog [options] mtzfile pdbfile",
                                   description="Prepare script for refmac5.")

    (opts, args) = parser.parse_args(sys.argv)

    pdbfile, mtzfile = None, None

    for arg in args[1:]:
        if arg.endswith(".pdb"):
            pdbfile = arg
        elif arg.endswith(".mtz"):
            mtzfile = arg

    if None in (pdbfile, mtzfile):
        print("Provide pdb and mtz")
        quit()



    arrays = iotbx.mtz.object(mtzfile).as_miller_arrays()

    col_free = []
    col_target = []
    col_ph = []
    wavelengths = []
    for array in arrays:
        wavelengths.append(array.info().wavelength)
        labs = tuple(array.info().labels)

        if is_rfree_array(array, array.info()):
            flag_scores = get_r_free_flags_scores(miller_arrays=[array], test_flag_value=None)
            flag_value = flag_scores.test_flag_values[0]
            col_free.append(("FREE="+array.info().label_string(), flag_value))

        elif array.is_xray_intensity_array():
            if "(+)" in array.info().label_string(): #array.anomalous_flag():
                col_target.append("I+=%s SIGI+=%s I-=%s SIGI-=%s" % tuple(array.info().labels))
            else:
                if len(labs) < 2:
                    col_target.append("IP=%s" % labs)
                else:
                    col_target.append("IP=%s SIGIP=%s" % labs)

        elif array.is_xray_amplitude_array():
            if array.sigmas() is None:
                print("Not having sigma:", array.info().labels)
                continue
            if array.anomalous_flag():
                if "+" not in array.info().label_string():
                    print("non-supported anomalous array:", array.info().labels)
                    col_target.append("FP=%s SIGFP=%s" % tuple(array.info().labels[:2]))
                else:
                    col_target.append("F+=%s SIGF+=%s F-=%s SIGF-=%s" % tuple(array.info().labels))
            else:
                col_target.append("FP=%s SIGFP=%s" % tuple(array.info().labels))

        elif array.is_hendrickson_lattman_array():
            col_ph.append("HLA=%s HLB=%s HLC=%s HLD=%s" % tuple(array.info().labels))

        else:
            print("Unrecognized:", array.info().labels)

    if os.path.exists("go_refmac_1.sh") or os.path.exists("go_refmac_2.sh"):
        print()
        print("Error: Script file(s) already exists. aborting.")
        sys.exit(1)

    sh_1_str = """\
pref=refmac_1
refmac5 \\
 hklin %(mtzfile)s \\
 xyzin %(pdbfile)s \\
 hklout $pref.mtz \\
 xyzout $pref.pdb \\
<<+ > $pref.log
""" % dict(mtzfile=mtzfile, pdbfile=pdbfile)

    wavelengths = list(set([x for x in wavelengths if x>0]))
    wavelengths.sort()
    wavelengths_str = " ! ".join(["%.6f"%w for w in wavelengths])

    sh_common_str = ""
    for free, val in col_free:
        for trg in col_target:
            if len(col_ph) > 0 and "+=" not in trg:
                for ph in col_ph:
                    sh_common_str += "labin " + trg + " " + ph + " " + free + "\n"
            else:
                sh_common_str += "labin " + trg + " " + free + "\n"
            sh_common_str += "free %d\n" % val

    sh_common_str += """\

MAKE HYDR Y
ridg dist sigm 0.01
ncsr local
NCYC 30
! For twinned data:
! twin
! For SAD function:
! anomalous wavelength %s
+
""" % (wavelengths_str)

    open("go_refmac_1.sh", "w").write(sh_1_str + sh_common_str)

    sh_2_str = """\
for cyc in `seq 2 10`
do
prev=refmac_$((cyc-1))
pref=refmac_${cyc}
refmac5 \
 hklin %(mtzfile)s \\
 xyzin $prev.pdb \\
 hklout $pref.mtz \\
 xyzout $pref.pdb \\
<<+ > $pref.log
""" % dict(mtzfile=mtzfile)

    open("go_refmac_2.sh", "w").write(sh_2_str + sh_common_str + "done\n")
