"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc.xds import xds_ascii
from yamtbx.dataproc.xds.command_line.eval_resolution_hkl_files import eval_resolution
from cctbx.array_family import flex


def run(hklin):
    xscaled = xds_ascii.XDS_ASCII(hklin) # Must be XSCALE output
    merged_iobs = xscaled.i_obs().merge_equivalents(use_internal_variance=False).array()
    binner = merged_iobs.setup_binner(n_bins=100)

    isets = set(xscaled.iset)
    for_plot = {}

    cut_ios = (2, 1, 0.5, 0)
    print "iset file",
    for cut in cut_ios: print "cut_ios_%.2f" % cut,
    print

    for iset in isets:
        sel = (xscaled.iset == iset)
        data_i = xscaled.i_obs().select(sel).merge_equivalents(use_internal_variance=False).array()
        cutoffs = eval_resolution(data_i, 100, cut_ios)
        print "%3d %s %s" % (iset, xscaled.input_files[iset][0], " ".join(map(lambda x: "%.2f"%x, cutoffs)))

        name = "data%.2d" % iset
        for i_bin in binner.range_used():
            dmax, dmin = binner.bin_d_range(i_bin)
            Isel = data_i.resolution_filter(d_max=dmax, d_min=dmin)
            for_plot.setdefault(name, []).append(flex.mean(Isel.data()) if Isel.size()>0 else 0)

    import matplotlib
    matplotlib.use('Agg') # Allow to work without X
    import pylab
    import math
    from matplotlib.ticker import FuncFormatter
    s2_formatter = lambda x,pos: "inf" if x == 0 else "%.2f" % (1./math.sqrt(x))

    fig, ax1 = pylab.plt.subplots()
    plot_x = [binner.bin_d_range(i)[1]**(-2) for i in binner.range_used()]

    plots = {}
    for name, vals in for_plot.items():
        name = name[-30:].lstrip("_")
        plots[name] = pylab.plot(plot_x, vals, label=name)

    pylab.legend()
    pylab.xlabel('resolution (d^-3)')
    pylab.ylabel('<I>')
    pylab.setp(pylab.gca().get_legend().get_texts(), fontsize="small")
    plot_title = ""
    pylab.title(plot_title)

    pylab.gca().xaxis.set_major_formatter(FuncFormatter(s2_formatter))
    pylab.savefig("test.pdf")
    pylab.show()

if __name__ == "__main__":
    import sys

    hklin = sys.argv[1]
    run(hklin)
    
