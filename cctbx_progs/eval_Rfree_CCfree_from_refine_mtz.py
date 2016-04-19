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

def get_arrays(f):
    arrays = iotbx.mtz.object(f).as_miller_arrays()
    f_obs = filter(lambda x:x.info().labels[0]=="F-obs-filtered", arrays)[0]
    f_model = filter(lambda x:x.info().labels[0]=="F-model", arrays)[0].as_amplitude_array()
    flag = get_flag(arrays)

    print "#",f, "includes %d reflections. (d= %.2f - %.2f)" % ((f_obs.data().size(),)+ f_obs.d_max_min())

    f_obs, flag = f_obs.common_sets(flag)
    f_model, flag = f_model.common_sets(flag)
    f_model, flag = f_model.common_sets(flag)
    # Assuming f_obs and f_model is an already common set.

    return f_obs, f_model, flag
# get_arrays()

def calc_r(f_obs, f_model):
    return flex.sum(flex.abs(f_obs.data() - f_model.data())) / flex.sum(f_obs.data())
# calc_r()

def calc_cc(obs, model, as_intensity=True):
    if as_intensity:
        obs, model = map(lambda x:x.as_intensity_array(), (obs, model))

    corr = flex.linear_correlation(obs.data(), model.data())
    if corr.is_well_defined(): return corr.coefficient()
    else: return float("nan")
# calc_cc()

def run(mtzfile):
    f_obs, f_model, flag = get_arrays(mtzfile)

    f_obs.setup_binner(auto_binning=True)
    f_model.setup_binner(auto_binning=True)
    e_obs, e_model = map(lambda x:x.quasi_normalize_structure_factors(), (f_obs, f_model))

    binner = f_obs.setup_binner(reflections_per_bin=400)

    for i_bin in binner.range_used():
        dmax, dmin = binner.bin_d_range(i_bin)
        sel = binner.selection(i_bin)
        sel_work = sel & ~flag.data()
        sel_free = sel & flag.data()
        #f_obs_w, f_obs_t = f_obs.select(~flag.data()), f_obs.select(flag.data())
        #f_model_w, f_model_t = f_model.select(~flag.data()), f_model.select(flag.data())

        f_obs_w_b, f_model_w_b = f_obs.select(sel_work), f_model.select(sel_work)
        f_obs_t_b, f_model_t_b = f_obs.select(sel_free), f_model.select(sel_free)
        
        e_obs_w_b, e_model_w_b = e_obs.select(sel_work), e_model.select(sel_work)
        e_obs_t_b, e_model_t_b = e_obs.select(sel_free), e_model.select(sel_free)

        r_work, r_free = calc_r(f_obs_w_b, f_model_w_b), calc_r(f_obs_t_b, f_model_t_b)
        cc_work_E, cc_free_E = calc_cc(e_obs_w_b, e_model_w_b, False), calc_cc(e_obs_t_b, e_model_t_b, False)
        cc_work, cc_free = calc_cc(f_obs_w_b, f_model_w_b, True), calc_cc(f_obs_t_b, f_model_t_b, True)

        tmp = dmax, dmin, r_free, r_work, cc_free, cc_work, cc_free_E, cc_work_E, mtzfile
        print "%6.3f %6.3f %.4f %.4f % 7.4f % 7.4f % 7.4f % 7.4f %s" % tmp
# run()

if __name__ == "__main__":
    import sys

    print "  dmax   dmin  rfree  rwork  ccfree  ccwork ccfree.E ccwork.E file"

    for f in sys.argv[1:]:
        run(f)
