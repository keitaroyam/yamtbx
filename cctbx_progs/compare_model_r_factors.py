import math
import iotbx.mtz
import iotbx.phil
from cctbx.array_family import flex
from cctbx import miller
from iotbx.reflection_file_editor import is_rfree_array
from iotbx.reflection_file_utils import get_r_free_flags_scores

from get_ccano_mtz_mtz import remove_phase
from eval_Rfree_factors_with_common_reflections import get_flag

from pylab import *
from matplotlib.ticker import FuncFormatter

master_params_str="""\
dmin = None
  .type = float
dmax = None
  .type = float
nbins = 20
  .type = int
flag_name = None
  .type = str
flag_value = None
  .type = int
plot = *rfree *rwork
  .type = choice(multi=True)
"""

def calc_r(fo_sel, fm_sel):
    return flex.sum(flex.abs(fo_sel.data()-fm_sel.data())) / flex.sum(fo_sel.data())
# calc_r()

def run(mtz_files, params):
    fig, ax1 = plt.subplots()

    for mtzfile in mtz_files:
        arrays = iotbx.mtz.object(file_name=mtzfile).as_miller_arrays()
        remove_phase(arrays)
        f_obs = filter(lambda s: "F-obs-filtered" in s.info().label_string(), arrays)
        f_model = filter(lambda s: "F-model" in s.info().label_string(), arrays)

        if not (len(f_obs) == len(f_model) == 1):
            print "File %s does not contain single F-obs-filtered and F-model" % mtzfile
            continue

        #flag = get_flag(arrays, params.flag_name, params.flag_value)
        flag = get_flag(arrays, params.flag_name, params.flag_value)
        if flag is None:
            print "File %s does not contain test flag" % mtzfile
            continue

        assert f_obs[0].anomalous_flag() == f_model[0].anomalous_flag()
        if f_obs[0].anomalous_flag():
            flag = flag.generate_bijvoet_mates()

        flag = flag.resolution_filter(d_max=params.dmax, d_min=params.dmin)
        f_obs, f_model = f_obs[0], f_model[0]
        #f_obs, f_model = map(lambda d: d.resolution_filter(d_max=params.dmax, d_min=params.dmin),
        #                     (f_obs, f_model))
        f_model, flag = f_model.common_sets(flag)
        f_obs, flag = f_obs.common_sets(flag)
        f_model, flag = f_model.common_sets(flag)
        #f_model = f_model.common_set(flag)
        #f_obs = f_obs.common_set(flag)

        binner = f_obs.setup_binner(n_bins=params.nbins)
        plot_data = [[], [], []]
        for i_bin in binner.range_used():
            dmax, dmin = binner.bin_d_range(i_bin)
            flag_sel = flag.resolution_filter(d_max=dmax, d_min=dmin)
            fo_sel = f_obs.common_set(flag_sel)
            fm_sel = f_model.common_set(flag_sel)
            print fo_sel.size(), fm_sel.size(), flag_sel.size()
            r_free = calc_r(fo_sel.select(flag_sel.data()), fm_sel.select(flag_sel.data()))
            r_work = calc_r(fo_sel.select(~flag_sel.data()), fm_sel.select(~flag_sel.data()))
            plot_data[0].append(1./dmin**2)
            plot_data[1].append(r_free)
            plot_data[2].append(r_work)

        if "rfree" in params.plot:
            plot(plot_data[0], plot_data[1], label=mtzfile[-30:]+":free")
        if "rwork" in params.plot:
            plot(plot_data[0], plot_data[2], label=mtzfile[-30:]+":work")

    legend(loc="upper left")
    xlabel('d^-2')
    ylabel("R-factors")
    s2_formatter = lambda x,pos: "inf" if x == 0 else "%.2f" % (1./math.sqrt(x))
    gca().xaxis.set_major_formatter(FuncFormatter(s2_formatter))

    show()

# run()

if __name__ == "__main__":
    import sys
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    mtz_files = filter(lambda s:s.endswith(".mtz"), args)

    run(mtz_files, params)
