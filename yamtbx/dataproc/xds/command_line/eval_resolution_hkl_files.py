"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.xds import integrate_hkl_as_flex
from yamtbx.dataproc.xds import xds_ascii

import iotbx.phil
from cctbx.array_family import flex
import os

master_params_str = """\
d_min = None
 .type = float
d_max = None
 .type = float
n_bins = 100
 .type = int

cut_ios = 2, 1, 0.5, 0
 .type = floats
 .help = I/sigma cutoff points

variance_model = 4, 0.00010
 .type = floats(size=2)
fix_variance_model = false
 .type = bool
"""

def eval_resolution(i_obs, n_bins, cut_ios):
    binner = i_obs.setup_binner(n_bins=n_bins)
    cutoffs = map(lambda x: float("nan"), xrange(len(cut_ios)))
    d_min_last = i_obs.d_max_min()[0]
    for i_bin in binner.range_used():
        sel = binner.selection(i_bin)
        iobs = i_obs.select(selection=sel)
        if iobs.size() == 0: continue
        d_max, d_min = iobs.d_max_min()
        ios = flex.mean(iobs.data()/iobs.sigmas())
        for i, cut in enumerate(cut_ios):
            if ios < cut and cutoffs[i]!=cutoffs[i]: # NaN test
                cutoffs[i] = d_min_last

        d_min_last = d_min
    return cutoffs
# eval_resolution()

def run(files, params):
    print "filename",
    for cut in params.cut_ios: print "cut_ios_%.2f" % cut,
    print

    for f in files:
        is_xac = xds_ascii.is_xds_ascii(f)
        i_obs = None

        if is_xac:
            xac = xds_ascii.XDS_ASCII(f, read_data=True, i_only=True)
            xac.remove_rejected()
            i_obs = xac.i_obs().resolution_filter(d_min=params.d_min, d_max=params.d_max)

            if params.fix_variance_model:
                ao, bo = xac.variance_model
                an, bn = params.variance_model
                i_obs = i_obs.customized_copy(sigmas = flex.sqrt(flex.abs(an * (i_obs.sigmas()**2/ao + (bn-bo)*flex.pow2(i_obs.data())))))
        else:
            ihkl = integrate_hkl_as_flex.reader(f, read_columns=("IOBS","SIGMA"))
            i_obs = ihkl.i_obs().resolution_filter(d_min=params.d_min, d_max=params.d_max)

            if params.fix_variance_model:
                a, b = params.variance_model
                i_obs = i_obs.customized_copy(sigmas = flex.sqrt(a * (i_obs.sigmas()**2 + b*flex.pow2(i_obs.data()))))

        cutoffs = eval_resolution(i_obs, params.n_bins, params.cut_ios)

        print "%s %s" % (f, " ".join(map(lambda x: "%.2f"%x, cutoffs)))
        

if __name__ == "__main__":
    import sys
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    files = filter(lambda x: integrate_hkl_as_flex.is_integrate_hkl(x) or xds_ascii.is_xds_ascii(x), args)
    for arg in args:
        if arg.endswith(".lst"):
            for l in open(arg):
                if "#" in l: l = l[:l.index("#")]
                files.append(l.strip())

    run(files, params)
