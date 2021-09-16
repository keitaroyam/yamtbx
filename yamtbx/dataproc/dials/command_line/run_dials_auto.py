from __future__ import print_function
from yamtbx import util
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.auto import resolution_cutoff
from yamtbx.dataproc.dataset import find_existing_files_in_template
from yamtbx.dataproc import pointless
import pickle 
import json
import os
from cctbx import sgtbx
from cctbx import crystal
from libtbx.utils import Sorry
import dxtbx.model.crystal

def get_most_possible_symmetry(workdir):
    try:
        pointless_log = os.path.join(workdir, "pointless.log")
        xs = pointless.parse_pointless_output_for_symm(open(pointless_log).read()).get("symm")
        if xs is not None: return xs        
    except:
        pass

    try:
        xs = get_crystal_symmetry_from_json(os.path.join(workdir, "integrated_experiments.json"))
        if xs is not None: return xs
    except:
        pass

    try:
        xac = XDS_ASCII(os.path.join(workdir, "DIALS.HKL"),read_data=False)
        return xac.symm
    except:
        pass

# get_most_possible_symmetry()

def get_crystal_symmetry_from_json(jsonin):
    x = json.load(open(jsonin))
    if "crystal" not in x: return None
    assert len(x["crystal"]) == 1

    cf = dxtbx.model.crystal.CrystalFactory()
    c = cf.from_dict(x["crystal"][0])
    
    return crystal.symmetry(c.get_unit_cell(), space_group=c.get_space_group())
# get_crystal_symmetry_from_json()

def calc_merging_stats(hklin, cut_resolution=True):
    import iotbx.merging_statistics
    from cctbx.array_family import flex

    wdir = os.path.dirname(hklin)
    #pklout = os.path.join(wdir, "merging_stats.pkl")
    logout = open(os.path.join(wdir, "merging_stats.log"), "w")

    print(hklin, file=logout)
    print("", file=logout)
    print("Estimate cutoff", file=logout)
    print("================", file=logout)

    if not os.path.isfile(hklin):
        print("Error: does not exist: %s" % hklin, file=logout)
        return {}

    i_obs = iotbx.merging_statistics.select_data(hklin, data_labels=None)
    d_min = None
    if i_obs.size() < 10: return {}

    """
    # Test I/sigI in 100 shells
    cutoff_isigi = None
    try:
        i_obs_merged = i_obs.select(i_obs.sigmas()>0).merge_equivalents(use_internal_variance=False).array()
        binner = i_obs_merged.setup_binner(n_bins=100)
        flag_last = False
        min_ios = 1
        for i_bin in binner.range_used():
            d_min_bin = f1.binner().bin_d_range(i_bin)[0]
            sel = i_obs_merged.select(binner.selection(i_bin))
            if sel.size() == 0: continue
            ios = flex.mean(sel.data()/sel.sigmas())
            if flag_last:
                return d_min_bin
            if ios < min_ios:
                flag_last = True
                continue

    except Exception, e:
        print >>logout, e.message
        """

    try:
        cutoffs = resolution_cutoff.estimate_crude_resolution_cutoffs(i_obs=i_obs, i_over_sigma_min=1)
        cutoffs.show(out=logout)

        #if cutoffs.cc_one_half_cut != float("inf") and cut_resolution:
        #    d_min = cutoffs.cc_one_half_cut
        if cut_resolution and cutoffs.i_over_sigma_cut is not None:
            d_min = cutoffs.i_over_sigma_cut
    except Sorry as e:
        print(str(e), file=logout)

    print("", file=logout)
    print("Merging statistics", file=logout)
    print("===================", file=logout)

    try:
        stats = iotbx.merging_statistics.dataset_statistics(i_obs=i_obs,
                                                            d_min=d_min,
                                                            d_max=None,
                                                            n_bins=10,
                                                            anomalous=False, # OK?
                                                            log=logout)
        stats.show(out=logout)
    except (Sorry, RuntimeError) as e:
        print(str(e), file=logout)
        return {}

    #ret = dict(cutoff=d_min, cutoffs=cutoffs, stats=stats)
    #pickle.dump(ret, open(pklout, "w"))
    return dict(d_min=d_min, cutoffs=cutoffs, stats=stats)
# calc_merging_stats()

def run_dials_sequence(filename_template, prefix, nr_range, wdir, known_xs, overrides, scan_varying, nproc):
    log_out = open(os.path.join(wdir, "dials_sequence.log"), "w")
    pointless_log = os.path.join(wdir, "pointless.log")

    # Prepare
    img_files = find_existing_files_in_template(filename_template, nr_range[0], nr_range[1],
                                                datadir=os.path.dirname(prefix), check_compressed=True)
    if len(img_files) == 0:
        mylog.error("No files found for %s %s" % (filename_template, nr_range))
        return

    nproc_str = "nproc=%d"%nproc

    log_out.write("Importing %s range=%s\n" % (img_files, nr_range))
    log_out.write(" Overrides: %s\n" % overrides)
    log_out.flush()

    override_str = "" # TODO support other stuff.. (wavelength, distance, osc_range, rotation_axis,..)
    if "orgx" in overrides and "orgy" in overrides:
        override_str += "slow_fast_beam_centre=%.2f,%.2f " % (overrides["orgy"], overrides["orgx"])

    if len(img_files) == 1 and img_files[0].endswith(".h5"):
        util.call('dials.import "%s" %s image_range=%d,%d output.experiments=imported.expt' % (img_files[0], override_str,
                                                                                               nr_range[0], nr_range[1]),
                  wdir=wdir, stdout=log_out,
                  expects_out=[os.path.join(wdir, "imported.expt")])
    else:
        util.call('dials.import %s template="%s" image_range=%d,%d output.experiments=imported.expt' % (override_str,
                                                                                                        filename_template.replace("?","#"),
                                                                                                        nr_range[0], nr_range[1]),
                  wdir=wdir, stdout=log_out,
                  expects_out=[os.path.join(wdir, "imported.expt")])

    util.call("dials.find_spots imported.expt filter.d_max=30 %s output.reflections=strong.refl" % nproc_str, # global_threshold=200
              wdir=wdir, stdout=log_out,
              expects_out=[os.path.join(wdir, "strong.refl")])

    util.call("dials.export strong.refl format=xds xds.directory=.",
              wdir=wdir, stdout=log_out)

    index_ok = False
    for index_meth in ("fft3d", "fft1d", "real_space_grid_search"):
        for index_assi in ("local", "simple"):
            if index_ok: break
            cmd = "dials.index imported.expt strong.refl output.experiments=indexed.expt output.reflections=indexed.refl "
            cmd += "indexing.method=%s index_assignment.method=%s " % (index_meth, index_assi)
            if known_xs is not None:# not in (known.space_group, known.unit_cell):
                cmd += "unit_cell=%s space_group=%d " % (",".join(["%.3f"%x for x in known_xs.unit_cell().parameters()]),
                                                        known_xs.space_group().type().number())
            elif index_meth == "real_space_grid_search":
                continue

            log_out.write("Trying indexing.method=%s index_assignment.method=%s\n" % (index_meth, index_assi))
            log_out.flush()
            util.call(cmd, wdir=wdir, stdout=log_out)
            if os.path.isfile(os.path.join(wdir, "indexed.expt")):
                index_ok = True
            else:
                for f in ("dials.index.log", "dials.index.debug.log"):
                    util.rotate_file(os.path.join(wdir, f))

    if not index_ok:
        return

    files_for_integration = "indexed.expt indexed.refl"

    if scan_varying:
        util.call("dials.refine indexed.expt indexed.refl scan_varying=true output.experiments=refined.expt output.reflections=refined.refl",
                  wdir=wdir, stdout=log_out)
        if os.path.isfile(os.path.join(wdir, "refined.refl")):
            files_for_integration = "refined.expt refined.refl"
        else:
            log_out.write("dials.refine failed. using intedexed results.\n")

    util.call("dials.integrate %s min_spots.per_degree=10 output.experiments=integrated.expt output.reflections=integrated.refl %s" % (files_for_integration, nproc_str),
              wdir=wdir, stdout=log_out)
    util.call("dials.export integrated.refl integrated.expt mtz.hklout=integrated.mtz",
              wdir=wdir, stdout=log_out)
    util.call("pointless integrated.mtz hklout pointless.mtz",
              wdir=wdir, stdin="SETTING SYMMETRY-BASED\ntolerance 10\n", stdout=open(pointless_log, "w"))
    util.call("dials.export integrated.expt integrated.refl format=xds_ascii xds_ascii.hklout=DIALS.HKL",
              wdir=wdir, stdout=log_out)
    util.call("aimless hklin pointless.mtz hklout aimless.mtz",
              wdir=wdir, stdin="output UNMERGED\n", stdout=open(os.path.join(wdir, "aimless.log"), "w"))
              
    #job_str += "touch dials_job_finished\n"

    ret = calc_merging_stats(os.path.join(wdir, "aimless_unmerged.mtz"))
    ret["symm"] = get_most_possible_symmetry(wdir)

    pickle.dump(ret, open(os.path.join(wdir, "kamo_dials.pkl"), "wb"), -1)

    # TODO config.params.xds.exclude_resolution_range config.params.reverse_phi

# run_dials()
