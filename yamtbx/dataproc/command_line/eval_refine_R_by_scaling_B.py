"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
Assume that we have more than one mtz file given by phenix.refine.
Then we want to evaluate R-work/free factors using the common reflections.

Usage:
 phenix.python eval_Rfree_factors_with_common_reflections.py refine_001.mtz refine2_001.mtz ...
"""

import iotbx.mtz
from cctbx.array_family import flex
from cctbx import miller
from iotbx.reflection_file_editor import is_rfree_array, get_best_resolution, get_original_array_types
from iotbx.reflection_file_utils import get_r_free_flags_scores
from mmtbx.scaling import absolute_scaling
from yamtbx.dataproc.scale_data import kBdecider
import StringIO
import numpy

master_params_str="""\
dmin = None
  .type = float
dmax = None
  .type = float
take_common = False
  .type = bool
reference = first *bmin bmax bmed
  .type = choice(multi=False)
"""

def get_flag(arrays, flag_name=None, flag_value=None):
    flag_array = None

    # Get flag array
    if flag_name is None:
        flags = filter(lambda x: is_rfree_array(x, x.info()), arrays)
        if len(flags) == 0:
            print " No R free flags like column found."
            return None
        elif len(flags) > 1:
            print " More than one column which looks like R free flag:"
            for f in flags:
                print " ", f.info().label_string()
            return None
        else:
            flag_name = flags[0].info().label_string()
            flag_array = flags[0]
            print "# Guessing R free flag:", flag_name
    else:
        flags = filter(lambda x: flag_name==x.info().label_string(), arrays)
        if len(flags) == 0:
            print " Specified flag name not found:", flag
            quit()
        else:
            print "# Use specified flag:", flag
            flag_array = flags[0]

    # Get flag number
    if flag_value is None:
        flag_scores = get_r_free_flags_scores(miller_arrays=[flag_array], test_flag_value=flag_value)
        flag_value = flag_scores.test_flag_values[0]
        print "# Guessing flag number:", flag_value
    else:
        print "# Specified flag number:", flag_value

    return flag_array.customized_copy(data=flag_array.data() == flag_value)
# get_flag()

def commonalize(arrays):
    new_arrays = []
    a0 = arrays[0]
    for f, f_obs, f_model, flag in arrays[1:]:
        pairs = miller.match_indices(a0[1].indices(), f_obs.indices()).pairs()
        a0[1] = a0[1].select(pairs.column(0))
        a0[2] = a0[2].select(pairs.column(0))
        a0[3] = a0[3].select(pairs.column(0))
        f_obs = f_obs.select(pairs.column(1))
        f_model = f_model.select(pairs.column(1))
        flag = flag.select(pairs.column(1))
        new_arrays.append([f, f_obs, f_model, flag])

    new_arrays2 = []

    for f, f_obs, f_model, flag in new_arrays:
        pairs = miller.match_indices(a0[1].indices(), f_obs.indices()).pairs()
        f_obs = f_obs.select(pairs.column(1))
        f_model = f_model.select(pairs.column(1))
        flag = flag.select(pairs.column(1))
        new_arrays2.append([f, f_obs, f_model, flag])

    return [a0] + new_arrays2
# commonalize()

def get_arrays(mtzfiles, d_min=None, d_max=None):
    ret = []
    for f in mtzfiles:
        arrays = iotbx.mtz.object(f).as_miller_arrays()
        #print "labs:", map(lambda x:x.info().labels, arrays)
        f_obs = filter(lambda x:x.info().labels[0]=="F-obs-filtered", arrays)[0]
        f_model = filter(lambda x:x.info().labels[0]=="F-model", arrays)[0].as_amplitude_array()
        flag = get_flag(arrays)

        f_obs, f_model, flag = map(lambda x: x.resolution_filter(d_min=d_min, d_max=d_max),
                                   (f_obs, f_model, flag))

        print "#", f, "includes %d reflections. (d= %.2f - %.2f)" % ((f_obs.data().size(),)+ f_obs.d_max_min())

        f_obs, flag = f_obs.common_sets(flag)
        f_model, flag = f_model.common_sets(flag)
        f_model, flag = f_model.common_sets(flag)
        # Assuming f_obs and f_model is an already common set.
        ret.append([f, f_obs, f_model, flag])

    return ret
# get_arrays()

def calc_r(f_obs, f_model, scale):
    return flex.sum(scale*flex.abs(f_obs.data() - f_model.data())) / flex.sum(scale*f_obs.data())
# calc_r()

def calc_cc(obs, model, as_intensity=True):
    if as_intensity:
        obs, model = map(lambda x:x.as_intensity_array(), (obs, model))

    corr = flex.linear_correlation(obs.data(), model.data())
    if corr.is_well_defined(): return corr.coefficient()
    else: return float("nan")
# calc_cc()

def calc_sigmaa(obs, model, flag):
    from mmtbx.scaling.sigmaa_estimation import sigmaa_estimator
    se = sigmaa_estimator(obs, model, flag, kernel_width_free_reflections=100)
    return se
# calc_sigmaa()    

def run(params, mtzfiles):
    arrays = get_arrays(mtzfiles, d_min=params.dmin, d_max=params.dmax)

    if params.take_common:
        arrays = commonalize(arrays)

    maxlen_f = max(map(lambda x: len(x[0]), arrays))

    ref_f_obs = arrays[0][1]

    scales = []
    for f, f_obs, f_model, flag in arrays:
        if ref_f_obs == f_obs: k, B = 1., 0
        else: k, B = kBdecider(ref_f_obs, f_obs).run()

        scales.append((k, B))

    if params.reference != "first":
        if params.reference == "bmin": # scale to strongest
            kref, bref = max(scales, key=lambda x:x[1])
        elif params.reference == "bmax": # scale to most weak
            kref, bref = min(scales, key=lambda x:x[1])
        elif params.reference == "bmed": # scale to most weak
            perm = range(len(scales))
            perm.sort(key=lambda i:scales[i][1])
            kref, bref = scales[perm[len(perm)//2]]
        else:
            raise "Never reaches here"

        print "# Set K=%.2f B=%.2f as reference" % (kref,bref)
        scales = map(lambda x: (x[0]/kref, x[1]-bref), scales) # not bref-x[1], because negated later

    print ("%"+str(maxlen_f)+"s r_work r_free cc_work.E cc_free.E sigmaa fom k B") % "filename"
    for (f, f_obs, f_model, flag), (k, B) in zip(arrays, scales):
        d_star_sq = f_obs.d_star_sq().data()
        scale = k * flex.exp(-B*d_star_sq)
        
        # Normalized
        #f_obs.setup_binner(auto_binning=True)
        #f_model.setup_binner(auto_binning=True)
        #e_obs, e_model = map(lambda x:x.quasi_normalize_structure_factors(), (f_obs, f_model))
        e_obs = absolute_scaling.kernel_normalisation(f_obs.customized_copy(data=f_obs.data()*scale, sigmas=None), auto_kernel=True)
        e_obs = e_obs.normalised_miller_dev_eps.f_sq_as_f()
        e_model = absolute_scaling.kernel_normalisation(f_model.customized_copy(data=f_model.data()*scale, sigmas=None), auto_kernel=True)
        e_model = e_model.normalised_miller_dev_eps.f_sq_as_f()

        f_obs_w, f_obs_t = f_obs.select(~flag.data()), f_obs.select(flag.data())
        f_model_w, f_model_t = f_model.select(~flag.data()), f_model.select(flag.data())

        e_obs_w, e_obs_t = e_obs.select(~flag.data()), e_obs.select(flag.data())
        e_model_w, e_model_t = e_model.select(~flag.data()), e_model.select(flag.data())

        r_work = calc_r(f_obs_w, f_model_w, scale.select(~flag.data()))
        r_free = calc_r(f_obs_t, f_model_t, scale.select(flag.data()))

        cc_work_E = calc_cc(e_obs_w, e_model_w, False)
        cc_free_E = calc_cc(e_obs_t, e_model_t, False)
        #cc_work_E2 = calc_cc(e_obs_w, e_model_w, True)
        #cc_free_E2 = calc_cc(e_obs_t, e_model_t, True)

        se = calc_sigmaa(f_obs, f_model, flag)
        sigmaa = flex.mean(se.sigmaa().data())
        fom = flex.mean(se.fom().data())

        print ("%"+str(maxlen_f)+"s %.4f %.4f % 7.4f % 7.4f %.4e %.4e %.3e %.3e") % (f, r_work, r_free, cc_work_E, cc_free_E, sigmaa, fom, k, B)
# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    run(params, args)
