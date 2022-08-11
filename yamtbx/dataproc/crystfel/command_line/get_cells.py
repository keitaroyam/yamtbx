"""
Usage: python ~/yamashita/work_scripts/get_cells.py all.stream  > cells.dat

R
library(reshape2)
library(ggplot2)

d <- read.table("cells.dat",h=T)
md <- melt(d, id=c("filename","tag","index.method"))
ggplot(md, aes(x=tag, y=value)) + geom_point() + facet_grid(variable~., scales="free") + geom_smooth(colour="darkgoldenrod1", size=1.5, method="loess", degree=0, span=0.1, se=FALSE)

"""
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import re
import shlex

# Cell parameters 4.89018 7.91817 8.63370 nm, 90.00000 90.00002 90.00002 deg
re_cell = re.compile("([0-9\.]+) ([0-9\.]+) ([0-9\.]+) nm, ([0-9\.]+) ([0-9\.]+) ([0-9\.]+) deg")

def parse_table(tablein):
    lines = open(tablein).readlines()
    keys = lines[0].strip().split()
    i_tag = keys.index("tag")
    data = {}
    for l in lines[1:]:
        sp = shlex.split(l)
        data[int(sp[i_tag])] = sp

    return keys, data

def run(streamin, tablein):

    keys, table = parse_table(tablein) if tablein is not None else (None, None)
    i_timestamp, i_photonE, i_pulse = None, None, None
    if table is not None:
        i_timestamp = keys.index("timestamp")
        i_photonE = keys.index("oh2_photonE")
        i_pulseE = keys.index("eh2_bm1_pulseE")

    imgf, tag, index_meth, pr = None, None, None, "nan"
    cells = []

    if table is None:
        print("filename tag index.method a b c alpha beta gamma pr dlim")
    else:
        print("filename tag index.method a b c alpha beta gamma pr dlim timestamp photonE pulseE")

    for l in open(streamin):
        if l.startswith("----- Begin chunk -----"):
            imgf, tag, index_meth, pr = None, None, None, "nan"
        elif l.startswith("Image filename:"):
            imgf = l[l.index(":")+1:].strip()
            try: tag = int(os.path.splitext(os.path.basename(imgf))[0].split("-")[-1])
            except: tag = -1
        elif l.startswith("Event: tag-"):
            tag = int(l[l.index("-")+1:l.index("/")])
        elif l.startswith("indexed_by ="):
            index_meth = l.split("=")[1].strip()
        elif l.startswith("profile_radius ="):
            pr, unit = l.split("=")[1].split()
            assert unit == "nm^-1"
        elif l.startswith("Cell parameters"):
            cell = [float(x) for x in re_cell.search(l).groups()]
            cells.append(cell)
        elif l.startswith("diffraction_resolution_limit = "):# diffraction_resolution_limit = 6.01 nm^-1 or 1.66 A
            dlim = l.split()[-2]
        elif l.startswith("Reflections measured after indexing"):
            ex_str = ""
            if table is not None:
                d = table[tag]
                ex_str = '"%s" %s %s' % (d[i_timestamp], d[i_photonE], d[i_pulseE])

            print("%s %d %s"%(imgf, tag, index_meth), "%.2f %.2f %.2f" % tuple([x*10 for x in cell[:3]]), "%.3f %.3f %.3f" % tuple(cell[3:]), pr, dlim, ex_str)

        #if l.startswith("----- End chunk -----"):

    # write median_cell.pdb if requested



if __name__ == "__main__":
    streamin = sys.argv[1]
    tablein = sys.argv[2] if len(sys.argv) > 2 else None
    run(streamin, tablein)
