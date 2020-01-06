"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os
import math
import numpy
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
    print "iset file cmpl redun",
    for cut in cut_ios: print "cut_ios_%.2f" % cut,
    print

    for iset in isets:
        sel = (xscaled.iset == iset)
        merge_i = xscaled.i_obs().select(sel).merge_equivalents(use_internal_variance=False)
        redun_i = flex.mean(merge_i.redundancies().data().as_double())
        data_i = merge_i.array()
        cmpl_i = data_i.completeness()
        cutoffs = eval_resolution(data_i, 100, cut_ios)
        print "%3d %s %5.1f%% %.1f %s" % (iset, xscaled.input_files[iset][0], cmpl_i*100, redun_i, " ".join(map(lambda x: "%.2f"%x, cutoffs)))

        for i_bin in binner.range_used():
            dmax, dmin = binner.bin_d_range(i_bin)
            Isel = data_i.resolution_filter(d_max=dmax, d_min=dmin)
            for_plot.setdefault(iset, []).append(flex.mean(Isel.data()) if Isel.size()>0 else 0)

    import matplotlib
    matplotlib.use('Agg') # Allow to work without X
    import pylab
    import math
    from matplotlib.ticker import FuncFormatter
    s2_formatter = lambda x,pos: "inf" if x == 0 else "%.2f" % (1./numpy.sqrt(x))
    exp_formatter = lambda x,pos: "%.1e" % x

    plot_x = [binner.bin_d_range(i)[1]**(-2) for i in binner.range_used()]

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages("test.pdf")

    keys = sorted(for_plot)
    for names in (keys[i:i+100] for i in xrange(0, len(keys), 100)):
        ncols = 5
        nrows = int(math.ceil(len(names)/float(ncols)))
        
        fig, axes = pylab.plt.subplots(ncols=ncols, nrows=nrows, figsize=(5*ncols, 5*nrows),
                                       sharex=False, sharey=False)
        axes = axes.flatten()
        for name, ax in zip(names, axes):
            ax.plot(plot_x, for_plot[name], linewidth=1, )
            ax.axhline(y=0, color="red", linestyle="-")
            ax.set_xlabel('(d^-2)')
            ax.set_ylabel('<I>')
            ax.xaxis.set_major_formatter(FuncFormatter(s2_formatter))
            ax.yaxis.set_major_formatter(FuncFormatter(exp_formatter))
            ax.set_title("data%.4d"%name)
            ax.grid(True)

        pylab.plt.tight_layout()
        plot_title = ""
        pylab.title(plot_title)
        pp.savefig()
        #pylab.savefig("test_%.3d.png"%i)
        #pylab.show()

    pp.close()

if __name__ == "__main__":
    import sys

    hklin = sys.argv[1]
    run(hklin)
    
