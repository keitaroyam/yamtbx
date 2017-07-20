"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.xds import xds_ascii
from yamtbx.dataproc.crystfel import hkl as crystfel_hkl
from yamtbx.util import read_path_list
import yamtbx_utils_ext
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx import easy_mp
import iotbx.scalepack.merge
import iotbx.phil
import libtbx.phil
import os
import sys
import numpy
import random
import multiprocessing
import threading
import cPickle as pickle

master_params_str = """\
lstin = None
 .type = path
dials_phil = None
 .type = path
dials_phil_selection_lst = None
 .type = path
 .help = can be used when subset of dials_phil= is merged. Give a list of reflection files
start_after = None
 .type = int
stop_after = None
 .type = int
method = *mc xscale
 .type = choice(multi=False)
 .help = mc: Monte-Carlo (same as CrystFEL) or xscale
nproc = 1
 .type = int(value_min=1)

usecell = *given mean median
 .type = choice(multi=False)
unit_cell = None
 .type = floats(size=6)
space_group = None
 .type = str
anomalous_flag = true
 .type = bool
sigma_calculation = *population experimental
  .type = choice(multi=False)

min_peak = None
 .type = float
min_peak_percentile = None
 .type = float
dials_data = *prf sum
 .type = choice(multi=False)
skip_rejected = true
 .type = bool
 .help = skip data with negative sigmas
correct_peak = false
 .type = bool
 .help = If true, intensity will be multiplied by PEAK
cancel_rlp = false
 .type = bool
 .help = divide intensities by RLP

calc_split = false
 .type = bool
split_method = *alternate random
 .type = choice(multi=False)
random_seed = 1234
 .type = int

output = *mtz *sca *crystfel
 .type = choice(multi=True)
prefix = merged
 .type = str
 .help = output prefix
dmin = None
 .type = float
dmax = None
 .type = float

scaling {
 parameter = k kb
  .type = choice(multi=False)
 bscale_ref = min mean max *asis
  .type = choice(multi=False)
  .help = when min, all data will be scaled down.
 reference = None
  .type = path
  .help = scaling reference. if None, two-pass scaling.
 reference_label = None
  .type = str
 calc_cc = false
  .type = bool
}

polarization {
 correct = false
  .type = bool
 plane_normal = 0,1,0
  .type = floats(size=3)
 incident_beam_direction = 0,0,1
  .type = floats(size=3)
 fraction = 0.99
  .type = float
}
"""

def read_list(lstin):
    ret = []
    for l in open(lstin):
        if "#" in l: l = l[:l.index("#")]
        l = l.strip()
        if l == "": continue
        if not os.path.isfile(l):
            print "WARNING:: not exists: %s" % l
            continue
        ret.append(l)
    return ret
# read_list()

def read_dials_phil(phil_in, dials_phil_selection_lst):
    master_str = """\
input {
  experiments = None
   .type = path
   .multiple = true
  reflections = None
   .type = path
   .multiple = true
}
"""
    params = iotbx.phil.process_command_line(args=[phil_in], master_string=master_str).work.extract()

    assert len(params.input.experiments) == len(params.input.reflections)
    assert all(map(lambda x: os.path.isfile(x), params.input.experiments))
    assert all(map(lambda x: os.path.isfile(x), params.input.reflections))

    ret = zip(params.input.experiments, params.input.reflections)

    if dials_phil_selection_lst:
        selected_files = set(read_path_list(dials_phil_selection_lst))
        ret = filter(lambda x: x[1] in selected_files, ret)

    return ret
# read_dials_phil()

#@profile
def get_data_from_xac(params, xac):
    if xac.endswith(".pkl"):
        tmp = pickle.load(open(xac))
    else:
        tmp = xds_ascii.XDS_ASCII(xac)
        
    sel_remove = flex.bool(tmp.iobs.size(), False)

    if params.min_peak is not None:
        sel = tmp.peak < params.min_peak
        sel_remove |= sel
    elif params.min_peak_percentile is not None:
        q = numpy.percentile(tmp.peak, params.min_peak_percentile)
        print "percentile %.2f %s" % (q, xac)
        sel = tmp.peak < q
        sel_remove |= sel

    if params.skip_rejected: sel_remove |= (tmp.sigma_iobs <= 0)
    if params.dmin is not None:
        sel_remove |= ~tmp.as_miller_set().resolution_filter_selection(d_min=params.dmin)

    if params.correct_peak:
        sel_remove |= (tmp.peak < 1) # remove PEAK==0

    # Remove selected
    print "DEBUG:: removing %d reflections" % sel_remove.count(True) #sum(sel_remove)#
    tmp.remove_selection(sel_remove)

    if not params.skip_rejected: tmp.sigma_iobs = flex.abs(tmp.sigma_iobs)

    # Correct I,sigI if needed
    if params.correct_peak:
        tmp.iobs *= tmp.peak * .01
        tmp.sigma_iobs *= tmp.peak * .01

    if params.cancel_rlp:
        tmp.iobs /= tmp.rlp
        tmp.sigma_iobs /= tmp.rlp

        if params.polarization.correct:
            # Only works with single-panel detector!!
            # Assumes detector fast = (1,0,0), slow = (0,1,0)
            sin_sq_2theta  = tmp.symm.unit_cell().sin_sq_two_theta(tmp.indices, tmp.wavelength)
            cos_sq_2theta = 1. - sin_sq_2theta
            sin_theta = tmp.wavelength / tmp.symm.unit_cell().d(tmp.indices) / 2.
            Eppi = numpy.cross(params.polarization.plane_normal, params.polarization.incident_beam_direction)
            Eppi /= numpy.linalg.norm(Eppi)
            S = flex.vec3_double(tmp.xd - tmp.orgx, tmp.yd - tmp.orgy,
                                 flex.double(tmp.xd.size(), tmp.distance/tmp.qx))
            S /= S.norms()

            zp = S.dot(Eppi.tolist()) * 2. * sin_theta
            cosrho = zp / flex.sqrt(sin_sq_2theta)
            P0 = 0.5 * (1. + cos_sq_2theta)
            PP = (params.polarization.fraction - 0.5) * (2.*cosrho**2 - 1.) * sin_sq_2theta
            P = P0 - PP

            # Apply correction
            tmp.iobs /= P
            tmp.sigma_iobs /= P

            if 0: # debug
                for x, y, p in zip(tmp.xd, tmp.yd, P):
                    print "pdebug:: %.2f %.2f %.4e" % (x, y, p)

    return tmp.i_obs().customized_copy(anomalous_flag=params.anomalous_flag,
                                       space_group_info=sgtbx.space_group_info(params.space_group))
# get_data_from_xac()

def get_data_from_dials(params, files):
    from dials.util.options import Importer, flatten_reflections, flatten_experiments
    
    importer = Importer(files, read_experiments=True, read_reflections=True)
    reflections = flatten_reflections(importer.reflections)
    experiments = flatten_experiments(importer.experiments)
    assert len(reflections) == len(experiments) == 1

    xs = crystal.symmetry(experiments[0].crystal.get_unit_cell(), space_group=experiments[0].crystal.get_space_group())

    tmp = reflections[0].select(reflections[0]["id"] >= 0)
    assert max(tmp["id"]) == 0

    if params.dials_data == "sum":
        assert "intensity.sum.value" in tmp
        assert "intensity.sum.variance" in tmp
    else:
        assert "intensity.prf.value" in tmp
        assert "intensity.prf.variance" in tmp

    intensity_key = "intensity.sum.value" if params.dials_data == "sum" else "intensity.prf.value"
    variance_key = "intensity.sum.variance" if params.dials_data == "sum" else "intensity.prf.variance"

    sel_remove = flex.bool(tmp.size(), False)
    
    if params.min_peak is not None:
        sel = tmp["partiality"] < params.min_peak/100.
        sel_remove |= sel
    elif params.min_peak_percentile is not None:
        q = numpy.percentile(tmp["partiality"], params.min_peak_percentile)
        print "percentile %.2f %s" % (q*100., xac)
        sel = tmp["partiality"] < q
        sel_remove |= sel

    if params.skip_rejected: sel_remove |= tmp[variance_key] <= 0
    if params.dmin is not None:
        sel_remove |= xs.unit_cell().d(tmp["miller_index"]) < params.dmin

    if params.correct_peak:
        sel_remove |= (tmp["partiality"] < .01) # remove PEAK==0

    # Remove selected
    print "DEBUG:: removing %d reflections" % sel_remove.count(True) #sum(sel_remove)#
    tmp = tmp.select(~sel_remove)

    ret = miller.array(miller.set(xs, tmp["miller_index"], params.anomalous_flag),
                       data=tmp[intensity_key],
                       sigmas=flex.sqrt(flex.abs(tmp[variance_key])))

    scale = flex.double(ret.size(), 1.)
    if not params.cancel_rlp and "lp" in tmp: scale *= tmp["lp"]
    if "dqe" in tmp: scale /= tmp["dqe"]
    if params.correct_peak: scale *= tmp["partiality"]

    ret = ret.apply_scaling(factor=scale)

    if params.cancel_rlp and params.polarization.correct:
        raise "Not implemented"

    return ret
# get_data_from_dials()

#@profile
def scale_data(indices, iobs, scale_ref, parameter, calc_cc):
    k, b, cc = 1, float("nan"), float("nan")

    sortp = yamtbx_utils_ext.sort_permutation_fast_less(indices)
    indices = indices.select(sortp)
    iobs = iobs.select(sortp)

    sel0, sel1 = yamtbx_utils_ext.my_common_indices(scale_ref.indices(), indices)
    #indices = indices.select(sel1)
    iobs_c = iobs.select(sel1)
    ref_c = scale_ref.data().select(sel0)

    if iobs_c.size() < 10 and ref_c.size() < 10:
        return k, b, cc

    if parameter == "k":
        k = flex.sum(ref_c*iobs_c) / flex.sum(flex.pow2(iobs_c))
    elif parameter == "kb":
        from yamtbx.dataproc.scale_data import kBdecider
        kbd = kBdecider(scale_ref,
                        miller.array(scale_ref.customized_copy(indices=indices),data=iobs))
        k, b = kbd.run()
    else:
        raise "Never reaches here"
    
    if calc_cc:
        corr = flex.linear_correlation(ref_c, iobs_c)
        if corr.is_well_defined(): cc = corr.coefficient()

    return k, b, cc
# scale_data()

#@profile
def mc_integration(params, input_files, input_type, scale_ref=None, split_idxes=None):
    """
    split_idxes: [i]->{0,1} if index given, returns 0 or 1, the group indicator.
    """
    import yamtbx_dataproc_crystfel_ext

    assert input_type in ("xds_ascii", "dials")

    sgtype = sgtbx.space_group_info(params.space_group).type()

    manager = multiprocessing.Manager()
    scale_strs = manager.list()

    def load_xds_or_dials(i, file_in, scout):
        print "Reading %5d %s" %(i, file_in if input_type=="xds_ascii" else file_in[1])
        if input_type=="xds_ascii":
            tmp = get_data_from_xac(params, file_in)
        else:
            tmp = get_data_from_dials(params, file_in)

        miller.map_to_asu(sgtype,
                          params.anomalous_flag,
                          tmp.indices())

        k, b = None, None
        if scale_ref is not None:
            k, b, cc = scale_data(tmp.indices(), tmp.data(), scale_ref, params.scaling.parameter, params.scaling.calc_cc)
            scout.append("%s %.4e %.4e %.4f\n" % (file_in if input_type=="xds_ascii" else file_in[1], k, b, cc))
            #tmp.iobs *= k
            #tmp.sigma_iobs *= k

            #if b == b: # nan if not calculated
            #    d_star_sq = tmp.symm.unit_cell().d_star_sq(tmp.indices)
            #    tmp.iobs *= flex.exp(-b*d_star_sq)
            #    tmp.sigma_iobs *= flex.exp(-b*d_star_sq)

        return i, tmp, k, b
    # load_xds_or_dials()

    if params.nproc > 1:
        # It may cost much memory..
        xds_data = easy_mp.pool_map(fixed_func=lambda x: load_xds_or_dials(x[0],x[1], scale_strs),
                                    args=map(lambda x:x, enumerate(input_files)),
                                    processes=params.nproc)
    else:
        # create generator
        xds_data = (load_xds_or_dials(i, x, scale_strs) for i, x in enumerate(input_files))

    merger = yamtbx_dataproc_crystfel_ext.merge_equivalents_crystfel()
    merger_split = None
    if split_idxes is not None:
        merger_split= map(lambda x: yamtbx_dataproc_crystfel_ext.merge_equivalents_crystfel(), xrange(2))

    print "Start merging"
    cells = []
    bs = [] # b-factor list
    bs_split = [[], []] # b-factor list for split data
    for i, x, k, b in xds_data:
        sys.stdout.write("Merging %7d\r" % (i+1))
        sys.stdout.flush()

        data, sigmas = x.data(), x.sigmas()

        if None not in (k,b):
            data *= k
            sigmas *= k

            if b == b: # nan if not calculated
                d_star_sq = x.unit_cell().d_star_sq(x.indices())
                data *= flex.exp(-b*d_star_sq)
                sigmas *= flex.exp(-b*d_star_sq)
                bs.append(b)

        if params.sigma_calculation == "population":
            merger.add_observations(x.indices(), data)
        else: # experimental
            merger.add_observations(x.indices(), data, sigmas)

        cells.append(x.unit_cell().parameters())
        if split_idxes is not None:
            if b is not None and b==b: bs_split[split_idxes[i]].append(b)
            if params.sigma_calculation == "population":
                merger_split[split_idxes[i]].add_observations(x.indices(), data)
            else: # experimental
                merger_split[split_idxes[i]].add_observations(x.indices(), data, sigmas)

    print "\nDone."

    if scale_ref is not None:
        scales_out = open(params.prefix+"_scales.dat", "w")
        print >>scales_out, "file k b cc"
        for s in scale_strs: scales_out.write(s)

    # Merge
    if params.sigma_calculation == "population":
        merger.merge()
    else: # experimental
        merger.merge(sigma="experimental sigma")

    # Construct arrays
    cells = numpy.array(cells)
    if params.usecell=="mean":
        cell = map(lambda i: cells[:,i].mean(), xrange(6))
    elif params.usecell=="median":
        cell = map(lambda i: numpy.median(cells[:,i]), xrange(6))
    else:
        cell = params.unit_cell

    miller_set = miller.set(crystal_symmetry=crystal.symmetry(cell, params.space_group),
                            indices=merger.indices,
                            anomalous_flag=params.anomalous_flag)
    iobs = miller.array(miller_set, data=merger.data, sigmas=merger.sigmas).set_observation_type_xray_intensity()
    reds = miller.array(miller_set, data=merger.redundancies)

    if len(bs) > 0 and params.scaling.bscale_ref != "asis":
        if params.scaling.bscale_ref == "mean": bref = sum(bs)/float(len(bs))
        elif params.scaling.bscale_ref == "min": bref = min(bs)
        elif params.scaling.bscale_ref == "max": bref = max(bs)
        else: raise "Never reaches here"

        print "Changing B-factor reference to %s (apply %.3f to all)" % (params.scaling.bscale_ref, -bref)
        scale = flex.exp(bref*iobs.unit_cell().d_star_sq(iobs.indices()))
        iobs = iobs.customized_copy(data=iobs.data()*scale,
                                    sigmas=iobs.sigmas()*scale)
        
    extra = []
    if split_idxes is not None:
        for m, bs in zip(merger_split, bs_split):
            if params.sigma_calculation == "population":
                m.merge()
            else: # experimental
                m.merge(sigma="experimental sigma")

            miller_set = miller.set(crystal_symmetry=crystal.symmetry(cell, params.space_group),
                                    indices=m.indices,
                                    anomalous_flag=params.anomalous_flag)
            a = miller.array(miller_set, data=m.data, sigmas=m.sigmas).set_observation_type_xray_intensity()
            r = miller.array(miller_set, data=m.redundancies)

            if len(bs) > 0 and params.scaling.bscale_ref != "asis":
                if params.scaling.bscale_ref == "mean": bref = sum(bs)/float(len(bs))
                elif params.scaling.bscale_ref == "min": bref = min(bs)
                elif params.scaling.bscale_ref == "max": bref = max(bs)
                else: raise "Never reaches here"

                print "Changing B-factor reference to %s (apply %.3f to all) for split data" % (params.scaling.bscale_ref, -bref)

                scale = flex.exp(bref*a.unit_cell().d_star_sq(a.indices()))
                a = a.customized_copy(data=a.data()*scale,
                                      sigmas=a.sigmas()*scale)

            extra.append((a,r))

    return iobs, reds, extra
# mc_integration()

def write_out(params, array, redundancies, suffix=None):
    prefix = params.prefix
    if suffix is not None: prefix += suffix

    if "sca" in params.output:
        iotbx.scalepack.merge.write(file_name="%s.sca" % prefix,
                                    miller_array=array,
                                    scale_intensities_for_scalepack_merge=True)

    if "mtz" in params.output and array.size() > 0:
        mtz_dataset = array.as_mtz_dataset(column_root_label="I")
        mtz_dataset.add_miller_array(miller_array=redundancies, column_root_label="MULT")
        mtz_object = mtz_dataset.mtz_object()
        mtz_object.write(file_name="%s.mtz" % prefix)

    if "crystfel" in params.output:
        crystfel_hkl.write_file(filename="%s.hkl" % prefix,
                                i_obs=array, nmeas=redundancies)
# write_out()

#@profile
def run(params):
    libtbx.phil.parse(master_params_str).format(params).show(prefix=" ")
    print

    if params.space_group is None:
        print "Give space_group!"
        quit()

    if params.usecell == "given" and params.unit_cell is None:
        print "Give unit_cell if usecell=given! Otherwise give usecell=mean"
        quit()

    if params.lstin:
        input_files = read_list(params.lstin)
        input_type = "xds_ascii"
    elif params.dials_phil:
        input_files = read_dials_phil(params.dials_phil, params.dials_phil_selection_lst)
        input_type = "dials"

    if params.start_after is not None:
        input_files = input_files[params.start_after:]
    if params.stop_after is not None:
        input_files = input_files[:params.stop_after]

    print "%3d %s files will be merged" % (len(input_files), input_type)

    scale_ref = None

    # Limitations
    if params.scaling.reference is not None:
        scale_ref = None
        raise "Sorry, not supported."

    random.seed(params.random_seed)

    # For split stats
    split_idxes = None
    if params.calc_split:
        split_idxes = ([0,1]*(len(input_files)//2+1))[:len(input_files)]
        if params.split_method == "alternate":
            pass
        elif params.split_method == "random":
            random.shuffle(split_idxes)
        else:
            raise "Never reaches here"

        ofs = open(params.prefix+"_split_idxes.dat", "w")
        for i, f in zip(split_idxes, input_files): ofs.write("%d %s\n"%(i,f if input_type=="xds_ascii" else f[1]))
        ofs.close()

    if params.method == "mc":
        iobs, reds, sp = mc_integration(params, input_files, input_type, scale_ref, split_idxes)
        if params.scaling.parameter is not None and params.scaling.reference is None:
            write_out(params, iobs, reds, suffix="_beforescale")
            if params.calc_split:
                for i, s in enumerate(sp):
                    write_out(params, s[0], s[1], suffix="_beforescale_sp%d"%(i+1))

            print "Second pass with scaling"

            sortp = yamtbx_utils_ext.sort_permutation_fast_less(iobs.indices())
            iobs = iobs.select(sortp)

            iobs, reds, sp = mc_integration(params, input_files, input_type, iobs, split_idxes)
    elif params.method == "xscale":
        raise "Not supported yet"
    else:
        raise "Never reaches here"

    write_out(params, iobs, reds)

    if "crystfel" in params.output:
        import iotbx.pdb
        s = iotbx.pdb.format_cryst1_record(iobs)
        open(params.prefix+"_cell.pdb", "w").write(s)

    if params.calc_split:
        for i, s in enumerate(sp):
            write_out(params, s[0], s[1], suffix="_sp%d"%(i+1))
# run()

def run_from_args(argv):
    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()

    run(params)

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
