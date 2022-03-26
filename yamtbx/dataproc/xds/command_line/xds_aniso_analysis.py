"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import re
import os
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.anisotropic_analysis import calc_principal_vectors, calc_stats_along_axes
from yamtbx.dataproc.auto import resolution_cutoff
import iotbx.file_reader
import iotbx.phil
from libtbx.utils import null_out
import numpy

master_params_str = """\
hklin = None
 .type = path
 .help = "XDS unmerged file"
hklin_merged = None
 .type = path
 .help = "Merged file if exists"
anomalous = None
 .type = bool
cone_angle = 20.
 .type = float
 .help = "Cone angle"
n_bins = 10
 .type = int
 .help = "Number of resolution bin"
fit_curve = True
 .type = bool
 .help = "Fit function for CC1/2 curve"
log_out = None
 .type = path
 .help = "Output filename"
"""

def parse_logfile(logfile):
    re_reso_cut = re.compile("CC1/2=([^ ]+) resolution for ([^ ]+): ([^ ]+) Angstrom")
    ret = {}
    if not os.path.isfile(logfile): return ret

    flag_eigen = False
    ret["aniso_cutoffs"] = []
    ret["eigen_values"] = []
    ret["has_anisotropy"] = True
    
    for l in open(logfile):
        if "No anisotropy in this symmetry" in l:
            ret["has_anisotropy"] = False
        elif "Eigenvalues/vectors:" in l:
            flag_eigen = True
        elif flag_eigen and not l.strip():
            flag_eigen = False
        elif flag_eigen:
            val = float(l.split("(")[0])
            label = l.split(")")[-1].strip()
            ret["eigen_values"].append((val, label))
        else:
            r = re_reso_cut.search(l)
            if r:
                cchalf_cutoff, label, resol = float(r.group(1)), r.group(2), float(r.group(3))
                ret["aniso_cutoffs"].append([cchalf_cutoff, label, resol, None]) # last None will be eigen value

    if ret["aniso_cutoffs"]:
        ret["aniso_cutoffs"].sort(key=lambda x:x[2])
        ret["d_min_best"] = ret["aniso_cutoffs"][0][2]
        ret["d_min_worst"] = ret["aniso_cutoffs"][-1][2]

        # Fill eigen values
        eigen_vals = dict([x[::-1] for x in ret["eigen_values"]])
        for x in ret["aniso_cutoffs"]:
            x[-1] = eigen_vals.get(x[1])
            if x[-1] is None:
                tmp = [k for k in list(eigen_vals.keys()) if k!="c*"] # XXX does this always work??
                if tmp: x[-1] = eigen_vals[tmp[0]]
    else:
        ret["d_min_best"] = ret["d_min_worst"] = float("nan")

    return ret
# parse_logfile()

def make_aniso_stats_table(i_obs, array_merged, cone_angle, n_bins, do_fit, log_out):
    """
    Rejected data must be removed from i_obs
    array_merged must have positive sigma values only.
    """
    i_obs = i_obs.map_to_asu()

    v_stars = calc_principal_vectors(array_merged, log_out)
    print("", file=log_out)

    if not v_stars:
        return

    binner = array_merged.setup_binner(n_bins=n_bins)
    
    cchalf = calc_stats_along_axes(i_obs, v_stars, binner, angle=cone_angle, kind="cchalf")
    ios = calc_stats_along_axes(array_merged, v_stars, binner, angle=cone_angle, kind="ioversigma")

    ret = dict(binner=binner, cchalf=cchalf, ios=ios)
    
    print("Anisotropic stats with cone angle of %.1f deg:" % cone_angle, file=log_out)
    for vc, vi in zip(cchalf, ios):
        label = vc[1]
        cc_ovl = vc[3]
        ios_ovl = vi[3]
        print("%20s Overall CC1/2=%.4f Mn(I/sigma)=%.4f" %(label, cc_ovl, ios_ovl), file=log_out)

    s2_max_min = i_obs.min_max_d_star_sq()
    print("""
$TABLE: Anisotropic statistics:
$GRAPHS: weighted CC1/2 v resolution:%(xmin).4f|%(xmax).4fx0|1:1,%(ccnumbers)s:
: Mn(I/sigma) v resolution:%(xmin).4f|%(xmax).4fx%(iosmin).4f|%(iosmax).4f:1,%(isnumbers)s:
: Number of unique reflections used in analysis v resolution:%(xmin).4f|%(xmax).4fx0|%(nrmax)d:1,%(nrnumbers)s:
$$ 1/d^2 %(cclabels)s %(islabels)s %(nrlabels)s$$
$$""" % dict(ccnumbers=",".join([str(x+2) for x in range(len(cchalf))]),
             isnumbers=",".join([str(x+2+len(cchalf)) for x in range(len(ios))]),
             nrnumbers=",".join([str(x+2+len(cchalf)+len(ios)) for x in range(len(ios))]),
             cclabels=" ".join(["CC1/2:"+x[1] for x in cchalf]),
             islabels=" ".join(["I/sd:"+x[1] for x in ios]),
             nrlabels=" ".join(["N:"+x[1] for x in ios]),
             iosmin=min(0,min([min(x[-1]) for x in ios])),
             iosmax=max([max(x[-1]) for x in ios]),
             nrmax=max([max(x[-2]) for x in ios]),
             xmin=s2_max_min[0], xmax=s2_max_min[1]), file=log_out)

    for i, i_bin in enumerate(binner.range_used()):
        d_max, d_min = binner.bin_d_range(i_bin)
        log_out.write("%.4f " % ((1./d_min**2+1./d_max**2)/2))
        for vstar, label, pn, cc_ovl, binnref, binstats in cchalf: log_out.write("% .4f "%binstats[i])
        for vstar, label, pn, cc_ovl, binnref, binstats in ios: log_out.write("% 6.2f "%binstats[i])
        for vstar, label, pn, cc_ovl, binnref, binstats in ios: log_out.write("%6d "%binnref[i])
        log_out.write("\n")
    print("$$", file=log_out)
    

    if do_fit:
        log_out.write("\n\nFitting function for CC1/2 curves\n\n")
        fitted_vals = {}
        ret["aniso_d_min"] = []
        for vstar, label, pn, cc_ovl, binnref, binstats in cchalf:
            #s2_list = map(lambda x: numpy.mean(1/numpy.array(binner.bin_d_range(x))**2), binner.range_used())
            s2_list = numpy.array([1/binner.bin_d_min(x)**2 for x in binner.range_used()])
            d0, r = resolution_cutoff.fit_curve_for_cchalf(s2_list, binstats, log_out, verbose=False)
            d_min_est = resolution_cutoff.resolution_fitted(d0, r, cchalf_min=0.5)
            ret["aniso_d_min"].append((label, d_min_est))
            fitted_vals[label] = (d_min_est, resolution_cutoff.fun_ed_aimless(s2_list, d0, r))
            log_out.write("CC1/2=0.5 resolution for %s: %.4f Angstrom\n" % (label, d_min_est))
        log_out.write("\n")

        print("""
$TABLE: Anisotropic statistics - CC1/2 fitted:
$GRAPHS: fitted curve for weighted CC1/2 v resolution:%(xmin).4f|%(xmax).4fx0|1:1,%(ccnumbers)s:
$$ 1/d^2 %(cclabels)s $$
$$""" % dict(ccnumbers=",".join([str(x+2) for x in range(len(cchalf))]),
             cclabels=" ".join(["CC1/2:%s:d=%.2fA"%(x[1],fitted_vals[x[1]][0]) for x in cchalf]),
             xmin=s2_max_min[0], xmax=s2_max_min[1]), file=log_out)

        for i, i_bin in enumerate(binner.range_used()):
            d_max, d_min = binner.bin_d_range(i_bin)
            log_out.write("%.4f " % ((1./d_min**2+1./d_max**2)/2))
            for x in cchalf: log_out.write("% .4f "%fitted_vals[x[1]][1][i])
            log_out.write("\n")
        print("$$", file=log_out)

    return ret
# make_aniso_stats_table()

def run(hklin, hklin_merged=None, cone_angle=20., n_bins=10, anomalous=None, do_fit=True, log_out=null_out()):
    if 1:
        xac = XDS_ASCII(hklin, i_only=True)
        xac.remove_rejected()
        i_obs = xac.i_obs()
    #else:
    #    import iotbx.mtz
    #    i_obs = filter(lambda x: "SIGI" in x.info().label_string(), iotbx.mtz.object(hklin).as_miller_arrays(merge_equivalents=False))[0]

    print("Unmerged intensity read from", hklin, file=log_out)
    i_obs.show_summary(log_out, prefix=" ")
    print("", file=log_out)

    if anomalous is not None and i_obs.anomalous_flag() != anomalous:
        print("Changing anomalous flag based on user's input", file=log_out)
        i_obs = i_obs.customized_copy(anomalous_flag=anomalous)

    if hklin_merged is not None:
        f = iotbx.file_reader.any_file(hklin)
        array_merged = f.file_server.get_xray_data(file_name=None,
                                                   labels=None,
                                                   ignore_all_zeros=True,
                                                   parameter_scope="",
                                                   prefer_anomalous=False,
                                                   prefer_amplitudes=False)
        print("Merged intensity read from", hklin_merged, file=log_out)
        array_merged.show_summary(log_out, prefix=" ")
    else:
        array_merged = i_obs.merge_equivalents(use_internal_variance=False).array()
        print("Merged intensity calculated", file=log_out)

    print("", file=log_out)

    bad_data = array_merged.select(array_merged.data() < -3*array_merged.sigmas()) # FIXME What if already omitted..
    i_obs = i_obs.delete_indices(other=bad_data)

    array_merged = array_merged.select(array_merged.sigmas()>0)

    if anomalous is not None and not anomalous and array_merged.anomalous_flag():
        print("Converting to non-anomalous data..\n", file=log_out)
        array_merged = array_merged.average_bijvoet_mates()


    return make_aniso_stats_table(i_obs, array_merged, cone_angle, n_bins, do_fit, log_out)
# run()

def run_from_args(args):
    import sys
    import libtbx.phil

    if not args or "-h" in args or "--help" in args:
        print("""\
Perform anisotropy analysis for XDS unmerged data.
Parameters:""")
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        return

    cmdline = iotbx.phil.process_command_line(args=args,
                                              master_string=master_params_str)
    params = cmdline.work.extract()

    if params.hklin is None and cmdline.remaining_args:
        params.hklin = cmdline.remaining_args.pop(0)

    if params.hklin_merged is None and cmdline.remaining_args:
        params.hklin_merged = cmdline.remaining_args.pop(0)

    log_out = sys.stdout if params.log_out is None else open(params.log_out, "w")

    libtbx.phil.parse(master_params_str).format(params).show(log_out, prefix=" ")
    log_out.write("\n")

    run(params.hklin, params.hklin_merged,
        cone_angle=params.cone_angle, n_bins=params.n_bins, anomalous=params.anomalous,
        do_fit=params.fit_curve, log_out=log_out)

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])

