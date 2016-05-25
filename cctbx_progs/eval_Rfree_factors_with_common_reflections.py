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
            print " Guessing R free flag:", flag_name
    else:
        flags = filter(lambda x: flag_name==x.info().label_string(), arrays)
        if len(flags) == 0:
            print " Specified flag name not found:", flag
            quit()
        else:
            print " Use specified flag:", flag
            flag_array = flags[0]

    # Get flag number
    if flag_value is None:
        flag_scores = get_r_free_flags_scores(miller_arrays=[flag_array], test_flag_value=flag_value)
        flag_value = flag_scores.test_flag_values[0]
        print " Guessing flag number:", flag_value
    else:
        print " Specified flag number:", flag_value

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

def get_arrays(mtzfiles):
    ret = []
    for f in mtzfiles:
        arrays = iotbx.mtz.object(f).as_miller_arrays()
        f_obs = filter(lambda x:x.info().labels[0].startswith("F-obs-filtered"), arrays)[0]
        f_model = filter(lambda x:x.info().labels[0].startswith("F-model"), arrays)[0].as_amplitude_array()
        flag = get_flag(arrays)

        print f, "includes %d reflections. (d= %.2f - %.2f)" % ((f_obs.data().size(),)+ f_obs.d_max_min())

        # XXX better treatment!!
        if f_obs.anomalous_flag(): f_obs = f_obs.average_bijvoet_mates()
        if f_model.anomalous_flag(): f_model = f_model.average_bijvoet_mates()

        f_obs, flag = f_obs.common_sets(flag)
        f_model, flag = f_model.common_sets(flag)
        f_model, flag = f_model.common_sets(flag)
        # Assuming f_obs and f_model is an already common set.
        ret.append([f, f_obs, f_model, flag])

    return ret
# get_arrays()

def calc_r(f_obs, f_model):
    return flex.sum(flex.abs(f_obs.data() - f_model.data())) / flex.sum(f_obs.data())
# calc_r()

if __name__ == "__main__":
    import sys

    mtzfiles = sys.argv[1:]
    arrays = get_arrays(mtzfiles)
    arrays = commonalize(arrays)

    maxlen_f = max([len(f) for f, fobs, fmodel, flag in arrays])

    print
    print "number of common reflections:", arrays[0][1].data().size()
    print

    # self-test
    for f, f_obs, f_model, flag in arrays[1:]:
        if not f_obs.is_similar_symmetry(arrays[0][1]):
            print "WARNING: unsimilar symmetry"
            print f_obs.crystal_symmetry().unit_cell(),
            print arrays[0][1].crystal_symmetry().unit_cell()
        assert (f_obs.indices() == arrays[0][1].indices()).count(False) == 0
        print f_model.indices().size()
        print arrays[0][2].indices().size()
        assert (f_model.indices() == arrays[0][2].indices()).count(False) == 0
        assert (flag.indices() == arrays[0][3].indices()).count(False) == 0
        assert (f_obs.indices() == f_model.indices()).count(False) == 0
        assert (f_obs.indices() == flag.indices()).count(False) == 0


    print ("%"+str(maxlen_f)+"s r_work r_free") % "filename"
    for f, f_obs, f_model, flag in arrays:
        f_obs_w, f_obs_t = f_obs.select(~flag.data()), f_obs.select(flag.data())
        f_model_w, f_model_t = f_model.select(~flag.data()), f_model.select(flag.data())

        r_work, r_free = calc_r(f_obs_w, f_model_w), calc_r(f_obs_t, f_model_t)

        print ("%"+str(maxlen_f)+"s %.4f %.4f") % (f, r_work, r_free)

    print
    print "=== By-shell"
    print 

    # By shell

    print ("%"+str(maxlen_f)+"s d_max d_min r_work r_free") % "filename"
    for f, f_obs, f_model, flag in arrays:
        binner = f_obs.setup_binner(n_bins=20)
        for i_bin in binner.range_used():
            sel = binner.selection(i_bin)
            d_max, d_min = binner.bin_d_range(i_bin)

            f_obs_w, f_obs_t = f_obs.select(~flag.data()&sel), f_obs.select(flag.data()&sel)
            f_model_w, f_model_t = f_model.select(~flag.data()&sel), f_model.select(flag.data()&sel)
            
            r_work, r_free = calc_r(f_obs_w, f_model_w), calc_r(f_obs_t, f_model_t)

            print ("%"+str(maxlen_f)+"s %6.2f %6.2f %.4f %.4f") % (f, d_max, d_min, r_work, r_free)
