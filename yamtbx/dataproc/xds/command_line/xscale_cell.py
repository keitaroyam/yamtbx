"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from cctbx import sgtbx
import numpy
import shutil
import os

def run(xscale_inp):
    inp_dir = os.path.dirname(xscale_inp)

    files = [y[1].replace("*","") for y in [x for x in get_xdsinp_keyword(xscale_inp) if x[0]=="INPUT_FILE"]]
    files = [os.path.join(inp_dir, x) if not os.path.isabs(x) else x for x in files]
    symms = [XDS_ASCII(x,read_data=False).symm for x in files]
    cells = numpy.array([x.unit_cell().parameters() for x in symms])
    sgs = [str(x.space_group_info()) for x in symms]
    laues = [str(x.space_group().build_derived_reflection_intensity_group(False).info()) for x in symms]

    median_cell = [numpy.median(cells[:,i]) for i in range(6)]
    mean_cell = [cells[:,i].mean() for i in range(6)]
    cell_sd = [numpy.std(cells[:,i]) for i in range(6)]

    print("%4d files loaded" % len(files))
    print("Space groups:", ", ".join(["%s (%d files)"%(x,sgs.count(x)) for x in set(sgs)]))
    print(" Laue groups:", ", ".join(["%s (%d files)"%(x,laues.count(x)) for x in set(laues)]))
    print(" Median cell:", " ".join(["%7.3f"%x for x in median_cell]))
    print("   Mean cell:", " ".join(["%7.3f"%x for x in mean_cell]))
    print("          SD:", " ".join(["%7.1e"%x for x in cell_sd]))

    # for BLEND $CCP4/share/blend/R/blend0.R
    # names(macropar) <- c("cn","a","b","c","alpha","beta","gamma","mosa","ctoddist","wlength")
    ofs = open("forR_macropar.dat", "w")
    for i, cell in enumerate(cells):
        print("%4d" % (i+1), end=' ', file=ofs)
        print(" ".join(["%7.3f"%x for x in cell]), end=' ', file=ofs)
        print(" 0 0 0", file=ofs)
    ofs.close()

    shutil.copyfile("forR_macropar.dat", "forR_macropar.dat.bak")
    print()
    print("Run BLEND?")
    print("Rscript $CCP4/share/blend/R/blend0.R")
# run()

if __name__ == "__main__":
    import sys
    xscale_inp = sys.argv[1]
    run(xscale_inp)
