"""
Apply sigma cutoff for merged intensity.

Usage:

PHENIX_TRUST_OTHER_ENV=1 phenix.python apply_sigma_cutoff_for_merged.py hoge.mtz cutoff=-2
"""

master_params_str = """\
hklin = None
    .type=str
    .short_caption = file name of mtz including I
hklref = None
    .type=str
    .short_caption = mtz file after phenix.refine
    .help = If specified, calculate R-cryst factor between rejected F-obs and provided R-model.
xds_ascii = None
    .type=str
    .short_caption = XDS_ASCII.HKL
    .help = If specified, discarded reflections will be mapped on detector surface.
hklout = None
    .type=str
intensity_labels = None
    .type = strings
amplitude_labels = None
    .type = strings
    .help = If specified, remove amplitudes of reflections whose I/sigI are less than specified.
cutoff = None
    .type = float
d_min = None
    .type = float
    .help = d_min for applying cutoff (reflections with smaller d will be kept)
"""

import iotbx.phil
import iotbx.file_reader
import iotbx.mtz
import mmtbx.utils
import os

if __name__ == "__main__":
    import sys

    parsed = iotbx.phil.parse(master_params_str, process_includes=True)

    processed_args = mmtbx.utils.process_command_line_args(args = sys.argv[1:],
                                                           log = sys.stdout,
                                                           master_params = parsed)

    working_phil = processed_args.params
    params = working_phil.extract()

    if params.hklin is None:
        mtz_files = filter(lambda x:x.endswith(".mtz"), processed_args.reflection_file_names)
        if len(mtz_files) != 1:
            print "Exactly one mtz file must be given."
            sys.exit(1)
        params.hklin = mtz_files[0]

    if params.xds_ascii is None:
        xds_ascii_files = filter(lambda x:"XDS_ASCII" in x and not x.endswith(".mtz"), processed_args.reflection_file_names)
        if len(xds_ascii_files) > 1:
            print "Exactly one XDS_ASCII file must be given."
            sys.exit(1)
        elif len(xds_ascii_files) == 1:
            params.xds_ascii = xds_ascii_files[0]

    if params.cutoff is None:
        print "Specify cutoff value."
        sys.exit(1)

    if params.hklout is None:
        params.hklout = os.path.splitext(os.path.basename(params.hklin))[0] + "_cut%.2f.mtz" % params.cutoff

    working_phil = parsed.format(python_object=params)
    print "Parameters:"
    working_phil.show(out = sys.stdout, prefix=" ")


    fs = iotbx.file_reader.any_file(params.hklin).file_server

    i_obs = fs.get_xray_data(file_name = params.hklin,
                             labels = params.intensity_labels,
                             ignore_all_zeros = True,
                             parameter_scope = '',
                             parameter_name = 'intensity_labels')

    f_obs = None
    if params.amplitude_labels is not None:
        f_obs = fs.get_xray_data(file_name = params.hklin,
                                 labels = params.amplitude_labels,
                                 ignore_all_zeros = True,
                                 parameter_scope = '',
                                 parameter_name = 'amplitude_labels',
                                 prefer_amplitudes = True)
        i_obs, f_obs = i_obs.common_sets(f_obs)

    print "Using:", i_obs.info()

    remove_sel = i_obs.data() / i_obs.sigmas() <= params.cutoff

    if params.d_min is not None:
        print "Applying resolution cutoff for remove selection", params.d_min
        remove_sel &= i_obs.d_spacings().data() > params.d_min

    # Show by bins
    binner = i_obs.setup_binner(n_bins=20)
    for i_bin in binner.range_used():
        sel = remove_sel.select(binner.bin_indices() == i_bin)
        print "%6.3f - %6.3f" % binner.bin_d_range(i_bin),
        print "", sum(sel), "removed"

    print "Totally", sum(remove_sel), "reflections are removed."

    mtz_dataset = iotbx.mtz.object(params.hklin).\
        add_crystal("crystal", "project", i_obs.unit_cell()).\
        add_dataset(name="dataset", wavelength=0)

    mtz_dataset.add_miller_array(miller_array=i_obs.select(~remove_sel), column_root_label="ICUT")
    mtz_dataset.add_miller_array(miller_array=i_obs.select(remove_sel), column_root_label="IREMOVED")

    if f_obs is not None:
        mtz_dataset.add_miller_array(miller_array=f_obs.select(~remove_sel), column_root_label="FCUT")
        mtz_dataset.add_miller_array(miller_array=f_obs.select(remove_sel), column_root_label="FREMOVED")

    mtz_dataset.mtz_object().write(file_name=params.hklout)

    if params.xds_ascii is not None:
        # XXX Need to check unit cell compatiblity
        from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
        from cctbx import miller
        xa = XDS_ASCII(params.xds_ascii, sys.stdout)
        miller.map_to_asu(xa.symm.space_group_info().type(), False, xa.indices)
        removed_indices = i_obs.indices().select(remove_sel)

        out = open("removed_positions.dat", "w")
        for hkl, x, y, z, i, sigi in zip(xa.indices, xa.xd, xa.yd, xa.zd, xa.iobs, xa.sigma_iobs):
            if sigi <= 0:
                print "sigi<=0", x, y, z, i, sigi
                continue
            if hkl in removed_indices:
                print >>out, x, y, z, i, sigi

    if params.hklref is not None:
        #from eval_Rfree_factors_with_common_reflections import get_flag
        from cctbx.array_family import flex
        calc_r = lambda f_obs, f_model: flex.sum(flex.abs(f_obs.data() - f_model.data())) / flex.sum(f_obs.data())
        hklref_arrays = iotbx.mtz.object(params.hklref).as_miller_arrays()
        f_model = filter(lambda x:x.info().labels[0]=="F-model", hklref_arrays)[0].as_amplitude_array()
        #test_flag = get_flag(hkref_arrays)
        for i_bin in binner.range_used():
            f_obs_sel = f_obs.select((binner.bin_indices() == i_bin) & remove_sel)
            f_obs_sel, f_model_sel = f_obs_sel.common_sets(f_model)
            print "%6.3f - %6.3f" % binner.bin_d_range(i_bin),
            if f_obs_sel.data().size() > 0:
                print " %4d Rcryst=%.4f" % (f_obs_sel.data().size(), calc_r(f_obs_sel, f_model_sel))
            else:
                print " %4d Rcryst=-1" % f_obs_sel.data().size()
