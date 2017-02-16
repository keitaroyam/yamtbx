"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.anisotropic_analysis import calc_principal_vectors, calc_stats_along_axes
import iotbx.file_reader
import iotbx.phil
from libtbx.utils import null_out

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
log_out = None
 .type = path
 .help = "Output filename"
"""

def make_aniso_stats_table(i_obs, array_merged, cone_angle, n_bins, log_out):
    """
    Rejected data must be removed from i_obs
    array_merged must have positive sigma values only.
    """
    i_obs = i_obs.map_to_asu()

    v_stars = calc_principal_vectors(array_merged, log_out)
    print >>log_out, ""

    if not v_stars:
        return

    binner = array_merged.setup_binner(n_bins=n_bins)
    
    cchalf = calc_stats_along_axes(i_obs, v_stars, binner, angle=cone_angle, kind="cchalf")
    ios = calc_stats_along_axes(array_merged, v_stars, binner, angle=cone_angle, kind="ioversigma")
    
    print >>log_out, "Anisotropic stats with cone angle of %.1f deg:" % cone_angle
    for vc, vi in zip(cchalf, ios):
        label = vc[1]
        cc_ovl = vc[3]
        ios_ovl = vi[3]
        print >>log_out, "%20s Overall CC1/2=%.4f Mn(I/sigma)=%.4f" %(label, cc_ovl, ios_ovl)

    s2_max_min = i_obs.min_max_d_star_sq()
    print >>log_out, """
$TABLE: Anisotropic statistics:
$GRAPHS: weighted CC1/2 v resolution:%(xmin).4f|%(xmax).4fx0|1:1,%(ccnumbers)s:
: Mn(I/sigma) v resolution:%(xmin).4f|%(xmax).4fx%(iosmin).4f|%(iosmax).4f:1,%(isnumbers)s:
: Number of unique reflections used in analysis v resolution:%(xmin).4f|%(xmax).4fx0|%(nrmax)d:1,%(nrnumbers)s:
$$ 1/d^2 %(cclabels)s %(islabels)s %(nrlabels)s$$
$$""" % dict(ccnumbers=",".join(map(lambda x: str(x+2), xrange(len(cchalf)))),
             isnumbers=",".join(map(lambda x: str(x+2+len(cchalf)), xrange(len(ios)))),
             nrnumbers=",".join(map(lambda x: str(x+2+len(cchalf)+len(ios)), xrange(len(ios)))),
             cclabels=" ".join(map(lambda x: "CC1/2:"+x[1], cchalf)),
             islabels=" ".join(map(lambda x: "I/sd:"+x[1], ios)),
             nrlabels=" ".join(map(lambda x: "N:"+x[1], ios)),
             iosmin=min(0,min(map(lambda x:min(x[-1]), ios))),
             iosmax=max(map(lambda x:max(x[-1]), ios)),
             nrmax=max(map(lambda x:max(x[-2]), ios)),
             xmin=s2_max_min[0], xmax=s2_max_min[1])

    for i, i_bin in enumerate(binner.range_used()):
        d_max, d_min = binner.bin_d_range(i_bin)
        log_out.write("%.4f " % ((1./d_min**2+1./d_max**2)/2))
        for vstar, label, pn, cc_ovl, binnref, binstats in cchalf: log_out.write("% .4f "%binstats[i])
        for vstar, label, pn, cc_ovl, binnref, binstats in ios: log_out.write("% 6.2f "%binstats[i])
        for vstar, label, pn, cc_ovl, binnref, binstats in ios: log_out.write("%6d "%binnref[i])
        log_out.write("\n")

    print >>log_out, "$$"
# make_aniso_stats_table()

def run(hklin, hklin_merged=None, cone_angle=20., n_bins=10, anomalous=None, log_out=null_out()):
    if 1:
        xac = XDS_ASCII(hklin, i_only=True)
        xac.remove_rejected()
        i_obs = xac.i_obs()
    #else:
    #    import iotbx.mtz
    #    i_obs = filter(lambda x: "SIGI" in x.info().label_string(), iotbx.mtz.object(hklin).as_miller_arrays(merge_equivalents=False))[0]

    print >>log_out, "Unmerged intensity read from", hklin
    i_obs.show_summary(log_out, prefix=" ")
    print >>log_out, ""

    if anomalous is not None and i_obs.anomalous_flag() != anomalous:
        print >>log_out, "Changing anomalous flag based on user's input"
        i_obs = i_obs.customized_copy(anomalous_flag=anomalous)

    if hklin_merged is not None:
        f = iotbx.file_reader.any_file(hklin)
        array_merged = f.file_server.get_xray_data(file_name=None,
                                                   labels=None,
                                                   ignore_all_zeros=True,
                                                   parameter_scope="",
                                                   prefer_anomalous=False,
                                                   prefer_amplitudes=False)
        print >>log_out, "Merged intensity read from", hklin_merged
        array_merged.show_summary(log_out, prefix=" ")
    else:
        array_merged = i_obs.merge_equivalents(use_internal_variance=False).array()
        print >>log_out, "Merged intensity calculated"

    print >>log_out, ""

    bad_data = array_merged.select(array_merged.data() < -3*array_merged.sigmas()) # FIXME What if already omitted..
    i_obs = i_obs.delete_indices(other=bad_data)

    array_merged = array_merged.select(array_merged.sigmas()>0)

    if anomalous is not None and not anomalous and array_merged.anomalous_flag():
        print >>log_out, "Converting to non-anomalous data..\n"
        array_merged = array_merged.average_bijvoet_mates()


    make_aniso_stats_table(i_obs, array_merged, cone_angle, n_bins, log_out)
# run()

def run_from_args(args):
    import sys
    import libtbx.phil

    if not args or "-h" in args or "--help" in args:
        print """\
Perform anisotropy analysis for XDS unmerged data.
Parameters:"""
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
        log_out=log_out)

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])

