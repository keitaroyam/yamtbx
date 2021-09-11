#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

import sys, os, optparse, math
from collections import OrderedDict

from yamtbx.dataproc.scale_data import kBdecider

import iotbx.mtz
import iotbx.phil
from cctbx.array_family import flex
from cctbx import miller

master_params_str="""\
logscale = False
  .type = bool
  .help = log scale
dmin = None
  .type = float
dmax = None
  .type = float
nbins = 120
  .type = int
show = true
  .type = bool
  .help = show graphics
output = "plot.pdf"
  .type = path
take_anom_diff = False
  .type = bool
  .help = Use anomalous difference, I(+) - I(-)
over_sigma = False
  .type = bool
  .help = plot I/sigma rather than I
noscale = False
  .type = bool
  .help = Do not apply linear scale
extra = *cc rfactor rsplit no
  .type = choice(multi=False)
  .help = Plot CC or R-factor when two data given.
force_same_cell = False
  .type = bool
  .help = Use the same cell for all data
scale {
  dmin = None
    .type = float
  dmax = None
    .type = float
  bscale = False
    .type = bool
    .help = B-factor scaling
}
"""


def commonalize(Is_in):
    Is = [s[1] for s in Is_in]

    new_Is = []
    Is0 = Is[0]
    for I in Is[1:]:
        Is0, I = Is0.common_sets(I, assert_is_similar_symmetry=False)
        new_Is.append(I)

    Is = []

    for I in new_Is:
        I = I.common_set(Is0, assert_is_similar_symmetry=False)
        Is.append(I)

    res = [Is0,] + Is

    return [[o[0], r, o[2]] for o, r in zip(Is_in, res)] # putting back original data
# commonalize()

def decide_label(labels):
    labs = [x for x in labels if not x.upper().startswith(("SIG","PHI"))]
    return ",".join(labs)

if __name__ == "__main__":
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    if len(args) == 0:
        print("Usage: %s mtz1 lab1 mtz2 lab2 [mtz3 lab3...] param=value" % sys.argv[0])
        print()
        print("Defaultparamters:")
        print(cmdline)
        cmdline.work.format(python_object=params).show(out=sys.stdout, prefix=" ", attributes_level=1)
        quit()

    print("Paramters:")
    cmdline.work.format(python_object=params).show(out=sys.stdout, prefix=" ")
    print()

    if params.over_sigma:
        assert params.noscale


    Is = [] # [[name, miller_array, scale], ..]

    for mtzfile, label in ((args[2*i],args[2*i+1]) for i in range((len(args))//2)):
        mtzobj = iotbx.mtz.object(file_name=mtzfile)
        arrays = [s for s in mtzobj.as_miller_arrays() if label in s.info().labels]

        if len(arrays) == 0:
            print("ERROR! %s does not have column %s"%(mtzfile, label))
            print("Candidates:", [x.info().labels for x in mtzobj.as_miller_arrays()])
            quit()

        labels = arrays[0].info().labels
        if params.take_anom_diff:
            assert arrays[0].anomalous_flag()
            data = arrays[0].as_intensity_array().anomalous_differences()
        elif arrays[0].is_complex_array() or arrays[0].is_xray_reconstructed_amplitude_array() or arrays[0].is_xray_amplitude_array():
            print("Warning - amplitude or complex array %s" % labels)
            data = arrays[0].as_intensity_array().as_non_anomalous_array().merge_equivalents(use_internal_variance=False).array()
        elif arrays[0].is_integer_array() or arrays[0].is_real_array():
            print("Warning - no experimental array %s" % labels)
            data = arrays[0].as_non_anomalous_array().as_double()
        else:
            raise "Can't plot %s" % labels

        data = data.resolution_filter(d_max=params.dmax, d_min=params.dmin)
        Is.append([mtzfile+":"+decide_label(labels), data, (1.0,0.0)])
        print("loaded:", mtzfile, labels)

        #if "hkl" in mtzfile:
        #    Is[-1][1] = miller.array(miller_set=Is[-1][1], data= Is[-1][1].data() * flex.exp(4.8*Is[-1][1].d_star_sq().data()))

    if params.force_same_cell:
        for x in Is[1:]:
            x[1] = x[1].customized_copy(crystal_symmetry=Is[0][1])
        
    # Take common sets
    Is = commonalize(Is) ####

    # Decide scale
    if not params.noscale:
        for i in range(1, len(Is)):
            I = Is[i][1].resolution_filter(d_max=params.scale.dmax, d_min=params.scale.dmin)
            I0 = Is[0][1].resolution_filter(d_max=params.scale.dmax, d_min=params.scale.dmin)
            I, I0 = I.common_sets(I0, assert_is_similar_symmetry=False)

            if params.scale.bscale:
                Is[i][2] = kBdecider(I0, I).run()
                print("Scale for", Is[i][0], "is", Is[i][2])
            else:
                scale = flex.sum(I0.data()*I.data()) / flex.sum(flex.pow2(I.data()))
                Is[i][2] = scale, 0
                print("Scale for", Is[i][0], "is", scale)

    print(Is[0][1].data().size())
    # Prepare plot data
    for_plot = OrderedDict() # {name: [mean, ...], ..}
    binner = Is[0][1].setup_binner(n_bins=params.nbins)#reflections_per_bin=50)
    for i_bin in binner.range_used():
        for name, I, (scale, b) in Is:
            dmax, dmin = binner.bin_d_range(i_bin)
            Isel = I.resolution_filter(d_max=dmax, d_min=dmin)

            #Isel = I.select(binner.bin_indices() == i_bin) # crash if not common
            if params.over_sigma:
                data = Isel.data() / Isel.sigmas()
            else:
                bfac = flex.exp(-b * Isel.d_star_sq().data()) if b != 0 else 1.
                data = Isel.data() *scale*bfac
            if len(data)==0:
                print("WARNING: ", name, "No data in %f .. %f" % binner.bin_d_range(i_bin))
                for_plot.setdefault(name, []).append(float("nan"))
            elif params.logscale:
                for_plot.setdefault(name, []).append(math.log(flex.mean(data))) # taking log<I>
            else:
                for_plot.setdefault(name, []).append(flex.mean(data))

    # If only two data in, calc CC.
    extra = []
    if len(Is) == 2 and params.extra.lower() != "no":
        for i_bin in binner.range_used():
            dmax, dmin = binner.bin_d_range(i_bin)
            #Isel0, Isel1 = map(lambda x:x[1].select(binner.bin_indices() == i_bin), Is)
            Isel0, Isel1 = [x[1].resolution_filter(d_max=dmax, d_min=dmin) for x in Is]
            Isel0, Isel1 = Isel0.common_sets(Isel1, assert_is_similar_symmetry=False)
            scale, b = Is[1][2][0], Is[1][2][1]
            #bfac = flex.exp(-b * Isel.d_star_sq().data()) if b != 0 else 1.
            #data = Isel.data() * scale*bfac
            bfac = flex.exp(-b * Isel1.d_star_sq().data()) if b != 0 else 1.

            if params.extra.lower() == "cc":
                corr = flex.linear_correlation(Isel0.data(), Isel1.data()*bfac)
                if corr.is_well_defined():
                    extra.append(corr.coefficient())
                else:
                    extra.append(float("nan"))
            elif params.extra.lower() == "rfactor":
                # This R-factor is like model R-factor.
                # but maybe this should be sum(|I1-I2|)/2sum(I1+I2) or something like that..?
                denom = flex.sum(Isel0.data())
                numer = flex.sum(flex.abs(Isel0.data() - Isel1.data()*scale*bfac))
                if denom != 0:
                    extra.append(numer/denom)
                else:
                    extra.append(float("nan"))
            elif params.extra.lower() == "rsplit":
                denom = flex.sum((Isel0.data()+Isel1.data()*scale*bfac)/2.)
                numer = flex.sum(flex.abs(Isel0.data() - Isel1.data()*scale*bfac))
                if denom != 0:
                    extra.append(1./math.sqrt(2.)*numer/denom)
                else:
                    extra.append(float("nan"))
            else:
                raise "Never reaches here"

        # Calc overall CC
        Isel0, Isel1 = Is[0][1].common_sets(Is[1][1], assert_is_similar_symmetry=False)
        scale, b = Is[1][2][0], Is[1][2][1]
        bfac = flex.exp(-b * Isel1.d_star_sq().data()) if b != 0 else 1.
        corr = flex.linear_correlation(Isel0.data(), Isel1.data()*bfac)
        if corr.is_well_defined():
            print("Overall CC=", corr.coefficient(), "with %d reflections" % len(Isel0.data()))
        # Calc overall R
        denom = flex.sum(Isel0.data())
        numer = flex.sum(flex.abs(Isel0.data() - Isel1.data()*scale*bfac))
        print("Overall R=", numer/denom)

    plot_x = [binner.bin_d_range(i)[1]**(-2) for i in binner.range_used()]

    print("       d %s" % " ".join(list(for_plot.keys())))
    f = ["%"+"%d"%len(x)+".2f" for x in list(for_plot.keys())]
    for i, x in enumerate(plot_x):
        line = "%8.3f " % (1./math.sqrt(x))
        line += " ".join([fj%for_plot[k][i] for fj, k in zip(f, for_plot)])
        if extra != []:
            line += " %.4f" % extra[i]
        print(line)

    # Plot
    import matplotlib
    matplotlib.use('Agg') # Allow to work without X
    from pylab import *
    from matplotlib.ticker import FuncFormatter
    s2_formatter = lambda x,pos: "inf" if x<1e-10 else "%.2f" % (1./math.sqrt(x))

    #from matplotlib.backends.backend_pdf import PdfPages

    fig, ax1 = plt.subplots()

    plots = {}
    for name, vals in list(for_plot.items()):
        name = name[-30:].lstrip("_")
        plots[name] = plot(plot_x, vals, label=name)

    legend()
    xlabel('resolution (d^-3)')
    if params.logscale:
        ylabel('log<I>')
    else:
        ylabel('<I/sigma>' if params.over_sigma else '<I>')
    setp(gca().get_legend().get_texts(), fontsize="small")
    plot_title = ""
    if params.scale.bscale:
        plot_title += ' Scaled with B-factors'
        plot_title += ' (' + ", ".join(["%.1f"%x[2][1] for x in Is]) + ')'

    if params.take_anom_diff:
        plot_title += ' With anom diff'

    title(plot_title)

    gca().xaxis.set_major_formatter(FuncFormatter(s2_formatter))

    if extra != []:
        ax2 = ax1.twinx()
        ax2.plot(plot_x, extra, 'r')
        ax2.set_ylabel('CC' if params.extra.lower()=="cc" else "R-factor", color='r')

    savefig(params.output)
    if params.show:
        show()
