#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

import os
from yamtbx.dataproc.xds import xscalelp

def run(wdir):
    xscale_lp = os.path.join(wdir, "XSCALE.LP")

    corrs = xscalelp.get_pairwise_correlations(xscale_lp)
    
    ofs_corrs = open("xscale_correlations.dat", "w")
    print("  i   j n.common  cc ratio.common      b", file=ofs_corrs)
    for i, j, ncommon, corr, ratio, bfac in corrs:
        print("%3d %3d %5d % .3f %9.4f % 9.4f" % (i, j, ncommon, corr, ratio, bfac), file=ofs_corrs)

    print("""
# Want to plot correlations?
R
library(ggplot2)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
d <- read.table("xscale_correlations.dat", h=T)
ggplot(d, aes(x=i, y=j, colour=cc, shape=cc<0.8)) + geom_point(size=4) + scale_color_gradientn(colours=jet.colors(10))
##

""")

if __name__ == "__main__":
    run(".")
