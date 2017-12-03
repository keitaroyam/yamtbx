#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
TODO: plot differences in direct beam and rotation axis
"""

import sys
from yamtbx.dataproc.xds import integratelp
from yamtbx.util.xtal import CellConstraints

def make_plot(lp, log_out):
    ofs = open(log_out, "w")

    ofs.write("$TABLE: Parameters estimated for each frame:\n")
    ofs.write("$GRAPHS\n")
    ofs.write(":scales")
    ofs.write(":A:1,2:\n")
    ofs.write(":number of overloaded/strong/unexpected reflections")
    ofs.write(":A:1,3,4,5:\n")
    ofs.write(":SIGMAB (beam divergence e.s.d.)")
    ofs.write(":A:1,6:\n")
    ofs.write(":SIGMAR (reflecting range e.s.d.)")
    ofs.write(":A:1,7:\n")
    ofs.write("$$\n")
    ofs.write("Frame scale overlods nstrong nrej sigmaD sigmaM $$\n$$\n")
    for f, scale, novl, nstr, nrej, sd, sm in zip(lp.frames, lp.scales, lp.overloads, lp.strongs, lp.rejecteds, lp.sigmads, lp.sigmars):
        ofs.write("%5d %s %d %d %d %s %s\n" % (f, scale, novl, nstr, nrej, sd, sm))

    ofs.write("$$\n")
    ofs.write("\n\n\n")

    ofs.write("$TABLE: Parameters estimated for each block:\n")
    ofs.write("$GRAPHS\n")
    ofs.write(":unit cell length a")
    ofs.write(":A:1,2:\n")

    cellconstr = CellConstraints(lp.space_group)
    if not cellconstr.is_b_equal_a():
        ofs.write(":unit cell length b")
        ofs.write(":A:1,3:\n")
    if not cellconstr.is_c_equal_a_b():
        ofs.write(":unit cell length c")
        ofs.write(":A:1,4:\n")
    if not cellconstr.is_angle_constrained("alpha"):
        ofs.write(":unit cell angle alpha")
        ofs.write(":A:1,5:\n")
    if not cellconstr.is_angle_constrained("beta"):
        ofs.write(":unit cell angle beta")
        ofs.write(":A:1,6:\n")
    if not cellconstr.is_angle_constrained("gamma"):
        ofs.write(":unit cell angle gamma")
        ofs.write(":A:1,7:\n")
    ofs.write(":rotations off from initial orientation")
    ofs.write(":A:1,8,9,10:\n")
    ofs.write(":distance")
    ofs.write(":A:1,11:\n")
    ofs.write(":deviations from predicted positions")
    ofs.write(":A:1,12,13:\n")
    ofs.write(":beam center")
    ofs.write(":A:1,14,15:\n")
    ofs.write(":missetting angles")
    ofs.write(":A:1,16,17,18:\n")
    ofs.write("$$\n")
    ofs.write("#image a b c alpha beta gamma rotx roty rotz dist spot spindle orgx orgy phix phiy phiz$$\n$$\n")
    for images, param in sorted(lp.blockparams.items()):
        for i in images:
            print >>ofs, "%4d " % i, "  ".join(param.get("cell", ["D"]*6)), " ".join(param.get("rotation", ["D"]*3)), param.get("dist","D"), param.get("spot","D"), param.get("spindle","D"), " ".join(param.get("orig",["D"]*2)), " ".join(param.get("misset",["D"]*3))

    ofs.write("$$\n")
    ofs.write("\n\n\n")


    ofs.write("$TABLE: sigmaB and sigmaR on 9 areas for each block:\n")
    ofs.write("$GRAPHS\n")
    ofs.write(":SIGMAB")
    ofs.write(":A:1,2,3,4,5,6,7,8,9,10:\n")
    ofs.write(":SIGMAR")
    ofs.write(":A:1,11,12,13,14,15,16,17,18,19:\n")
    ofs.write("$$\n")
    ofs.write("#image %s %s$$\n$$\n" % (" ".join(["sigmab%d"%x for x in range(1,10)]), " ".join(["sigmar%d"%x for x in range(1,10)])))
    for images, param in sorted(lp.blockparams.items()):
        for i in images:
            print >>ofs, "%4d " % i, " ".join(param["sigmab9"]), " ".join(param["sigmar9"])

    ofs.write("$$\n")
    ofs.write("\n\n\n")
# make_plot()

def run(int_lp, log_out="plot_integrate.log"):
    lpobj = integratelp.IntegrateLp(int_lp)
    make_plot(lpobj, log_out)
# run()

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        int_lp = sys.argv[1]
    else:
        int_lp = "INTEGRATE.LP"

    log_out = "plot_integrate.log"

    run(int_lp, log_out)

    print
    print "Run:"
    print "loggraph", log_out
