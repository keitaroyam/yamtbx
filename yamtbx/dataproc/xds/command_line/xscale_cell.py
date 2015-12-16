"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from cctbx import sgtbx
import numpy
import shutil
import os

def run(xscale_inp):
    inp_dir = os.path.dirname(xscale_inp)

    files = map(lambda y: y[1].replace("*",""), filter(lambda x: x[0]=="INPUT_FILE", get_xdsinp_keyword(xscale_inp)))
    files = map(lambda x: os.path.join(inp_dir, x) if not os.path.isabs(x) else x, files)
    symms = map(lambda x: XDS_ASCII(x,read_data=False).symm, files)
    cells = numpy.array(map(lambda x: x.unit_cell().parameters(), symms))
    sgs = map(lambda x: str(x.space_group_info()), symms)
    laues = map(lambda x: str(x.space_group().build_derived_reflection_intensity_group(False).info()), symms)

    median_cell = map(lambda i: numpy.median(cells[:,i]), xrange(6))
    mean_cell = map(lambda i: cells[:,i].mean(), xrange(6))
    cell_sd = map(lambda i: numpy.std(cells[:,i]), xrange(6))

    print "%4d files loaded" % len(files)
    print "Space groups:", ", ".join(map(lambda x: "%s (%d files)"%(x,sgs.count(x)), set(sgs)))
    print " Laue groups:", ", ".join(map(lambda x: "%s (%d files)"%(x,laues.count(x)), set(laues)))
    print " Median cell:", " ".join(map(lambda x: "%7.3f"%x, median_cell))
    print "   Mean cell:", " ".join(map(lambda x: "%7.3f"%x, mean_cell))
    print "          SD:", " ".join(map(lambda x: "%7.1e"%x, cell_sd))

    # for BLEND $CCP4/share/blend/R/blend0.R
    # names(macropar) <- c("cn","a","b","c","alpha","beta","gamma","mosa","ctoddist","wlength")
    ofs = open("forR_macropar.dat", "w")
    for i, cell in enumerate(cells):
        print >>ofs, "%4d" % (i+1),
        print >>ofs, " ".join(map(lambda x: "%7.3f"%x, cell)),
        print >>ofs, " 0 0 0"
    ofs.close()

    shutil.copyfile("forR_macropar.dat", "forR_macropar.dat.bak")
    print
    print "Run BLEND?"
    print "Rscript $CCP4/share/blend/R/blend0.R"
# run()

if __name__ == "__main__":
    import sys
    xscale_inp = sys.argv[1]
    run(xscale_inp)
