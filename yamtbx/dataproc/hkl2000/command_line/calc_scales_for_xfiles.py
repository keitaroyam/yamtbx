"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.hkl2000.xfile import DenzoXfile
from yamtbx.dataproc.scale_data import kBdecider, kBdecider3

import iotbx.phil
import iotbx.file_reader
from cctbx.array_family import flex

import traceback

master_params_str = """
lstin = None
 .type = path

d_min = None
 .type = float
d_max = None
 .type = float
sigma_cutoff = None
 .type = float

reference {
 file = None
  .type = path
 label = None
  .type = str
 d_min = None
  .type = float
 d_max = None
  .type = float
}

show_plot = False
 .type = bool
help = False
 .type = bool
"""


def run(params, xfiles):
    # read reference
    arrays = iotbx.file_reader.any_file(params.reference.file).file_server.miller_arrays
    arrays = filter(lambda ar: ar.is_xray_data_array(), arrays)
    if params.reference.label is not None:
        arrays = filter(lambda ar: ar.info().label_string() == params.reference.label, arrays)

    if len(arrays) != 1:
        print "Can't decide data to use in reference file:", params.reference.file
        print "Choose label"
        for ar in arrays: print ar.info().label_string()
        return

    refdata = arrays[0].as_intensity_array()
    refdata = refdata.resolution_filter(d_max=params.reference.d_max, d_min=params.reference.d_min)

    print "file n.common k b cc.org cc.mean cc.scaled a b c al be ga"
    
    for xf in xfiles:
        print "# Reading", xf
        try:
            xfile = DenzoXfile(xf)
        except:
            traceback.print_exc()
            continue
        a = xfile.miller_array(anomalous_flag=refdata.anomalous_flag())
        a = a.select(a.sigmas() > 0)
        a = a.resolution_filter(d_min=params.d_min, d_max=params.d_max)
        if params.sigma_cutoff is not None:
            a = a.select(a.data()/a.sigmas() >= params.sigma_cutoff)

        a = a.merge_equivalents(use_internal_variance=False).array()

        tmp, a = refdata.common_sets(a, assert_is_similar_symmetry=False)
        n_common = tmp.size()

        if n_common == 0:
            print "# No useful reflection in this file. skip."
            continue

        corr = flex.linear_correlation(tmp.data(), a.data())
        cc_org = corr.coefficient() if corr.is_well_defined() else float("nan")

        # Calc CC in resolution bin and average
        tmp.setup_binner(auto_binning=True)
        cc_bins = []
        for i_bin in tmp.binner().range_used():
            sel = tmp.binner().selection(i_bin)
            corr = flex.linear_correlation(tmp.select(sel).data(), a.select(sel).data())
            if not corr.is_well_defined(): continue
            cc_bins.append(corr.coefficient())

        cc_mean = sum(cc_bins) / float(len(cc_bins)) if len(cc_bins) > 0 else float("nan")
            
        # Determine scale and B
        k, b = kBdecider(tmp, a).run()

        bfac = flex.exp(-b * a.d_star_sq().data()) if b != 0 else 1.
        corr = flex.linear_correlation(tmp.data(), a.data() * k*bfac)
        cc_scaled = corr.coefficient() if corr.is_well_defined() else float("nan")

        print "%s %5d %.3e %.3e %.4f %.4f %.4f" % (xf, n_common, k, b, cc_org, cc_mean, cc_scaled),
        print ("%.3f "*6)%a.unit_cell().parameters()

        if params.show_plot:
            import pylab
            from matplotlib.ticker import FuncFormatter
            s3_formatter = lambda x,pos: "inf" if x == 0 else "%.2f" % (x**(-1/3))

            fig, ax1 = pylab.plt.subplots()

            plot_x = map(lambda i: tmp.binner().bin_d_range(i)[1]**(-3), tmp.binner().range_used())

            #for name, ar in (("reference", tmp), ("data", a)):
            vals = map(lambda i: flex.mean(tmp.data().select(tmp.binner().selection(i))), tmp.binner().range_used())
            pylab.plot(plot_x, vals, label="reference")

            scale = flex.sum(tmp.data()*a.data()) / flex.sum(flex.pow2(a.data()))
            print "Linear-scale=", scale
            vals = map(lambda i: scale*flex.mean(a.data().select(tmp.binner().selection(i))), tmp.binner().range_used())
            pylab.plot(plot_x, vals, label="data")
            vals = map(lambda i: flex.mean((a.data()*k*bfac).select(tmp.binner().selection(i))), tmp.binner().range_used())
            pylab.plot(plot_x, vals, label="data_scaled")

            """
            from mmtbx.scaling import absolute_scaling, relative_scaling
            ls_scaling = relative_scaling.ls_rel_scale_driver(tmp, tmp.customized_copy(data=a.data(),sigmas=a.sigmas()), use_intensities=True, scale_weight=True, use_weights=True)
            ls_scaling.show()
            vals = map(lambda i: flex.mean(ls_scaling.derivative.resolution_filter(*tmp.binner().bin_d_range(i)).data()), tmp.binner().range_used())
            pylab.plot(plot_x, vals, label="data_scaled2")
            """
            
            pylab.legend()
            pylab.xlabel('resolution (d^-3)')
            pylab.ylabel('<I>')
            pylab.setp(pylab.gca().get_legend().get_texts(), fontsize="small")
            pylab.title('Scaled with B-factors (%.2f)' % b)

            pylab.gca().xaxis.set_major_formatter(FuncFormatter(s3_formatter))

            ax2 = ax1.twinx()
            ax2.plot(plot_x, cc_bins, "black")
            ax2.set_ylabel('CC')
            pylab.show()
            


if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    if "-h" in args: params.help = True

    if params.help:
        print "Parameter syntax:\n"
        iotbx.phil.parse(master_params_str).show(prefix="  ")
        print
        print "Usage: calc_scales_for_xfiles.py [.x files] [lstin=xfiles.lst] reference.sca"
        quit()

    reffile = filter(lambda x: x.endswith((".sca",".mtz")), args)
    assert len(reffile) == 1
    params.reference.file = reffile[0]
        
    xfiles = filter(lambda x: x.endswith(".x"), args)
    if params.lstin is not None:
        for l in open(params.lstin):
            if "#" in l: l = l[:l.index("#")]
            l = l.strip()
            if l != "": xfiles.append(l)
    
    run(params, xfiles)
