import os
import cPickle as pickle
import numpy
from yamtbx.dataproc import crystfel
from yamtbx.util import read_path_list
from yamtbx.dataproc.xds import xds_ascii
from yamtbx.util.xtal import v6cell
import iotbx.phil
from cctbx.array_family import flex
from cctbx import uctbx
from mmtbx.scaling.absolute_scaling import ml_iso_absolute_scaling
from libtbx import easy_mp


master_params_str = """\
lstin = None
 .type = path
 .help = List of XDS_ASCII.HKL(.pkl)
datout = None
 .type = path
 .help = Output dat file, which is useful for plotting
nproc = 1
 .type = int

min_peak = None
 .type = float
min_peak_percentile = None
 .type = float
correct_peak = false
 .type = bool
 .help = If true, intensity will be multiplied by PEAK

ref_cell = 48.07, 77.45, 84.77, 90., 90., 90.
 .type = floats(size=6)
 .help = used for comparing unit cells
n_residues = 308
 .type = int
 .help = used for Wilson B calculation

stats = *reslimit *ioversigma *resnatsnr1 *pr *wilsonb *abdist
 .type = choice(multi=True)
 .help = "statistics. reslimit: diffraction_resolution_limit reported by CrystFEL; ioversigma: <I/sigma(I)> in each pattern; resnatsnr1: resolution where <I/sigma(I)>_bin drops below 1; pr: profile_radius reported by CrystFEL; wilsonb: ML-estimate of Wilson B value; abdist: Andrews-Bernstein distance (unit cell dissimilarity); ccref: CC with reference data"
"""

def calc_stats(xac_file, stat_choice, n_residues=None, ref_v6cell=None,
               min_peak=None, min_peak_percentile=None, correct_peak=None):
    # Open XDS_ASCII
    if xac_file.endswith(".pkl"): xac = pickle.load(open(xac_file))
    else: xac = xds_ascii.XDS_ASCII(xac_file)
    
    sel_remove = flex.bool(xac.iobs.size(), False)
    if min_peak is not None:
        sel = xac.peak < min_peak
        sel_remove |= sel
    elif min_peak_percentile is not None:
        q = numpy.percentile(xac.peak, min_peak_percentile)
        print "percentile %.2f %s" % (q, xac)
        sel = xac.peak < q
        sel_remove |= sel

    if correct_peak: sel_remove |= (xac.peak < 1) # remove PEAK==0

    xac.remove_selection(sel_remove)

    if params.correct_peak:
        xac.iobs *= xac.peak * .01
        xac.sigma_iobs *= xac.peak * .01

    iobs = xac.i_obs(anomalous_flag=False)
    iobs = iobs.select(iobs.sigmas()>0).merge_equivalents(use_internal_variance=False).array()
    
    stats = dict(filename=xac_file, cell=iobs.unit_cell().parameters())

    if iobs.size() == 0:
        return stats

    if "ioversigma" in stat_choice or "resnatsnr1" in stat_choice:
        binner = iobs.setup_binner(auto_binning=True)

        if "ioversigma" in stat_choice: stats["ioversigma"] = flex.mean(iobs.data()/iobs.sigmas())

        if "resnatsnr1" in stat_choice:
            res = float("nan")
            for i_bin in binner.range_used():
                sel = binner.selection(i_bin)
                tmp = iobs.select(sel)
                if tmp.size() == 0: continue
                sn = flex.mean(tmp.data()/tmp.sigmas())
                if sn <= 1:
                    res = binner.bin_d_range(i_bin)[1]
                    break

            stats["resnatsnr1"] = res

    if "abdist" in stat_choice:
        from cctbx.uctbx.determine_unit_cell import NCDist
        G6a, G6b = ref_v6cell, v6cell(iobs.unit_cell().niggli_cell())
        abdist = NCDist(G6a, G6b)
        stats["abdist"] = abdist

    if "wilsonb" in stat_choice:
        iso_scale_and_b = ml_iso_absolute_scaling(iobs, n_residues, 0)
        stats["wilsonb"] = iso_scale_and_b.b_wilson

    print stats
    return stats
# calc_stats()

def run(params):
    if params.datout is None: params.datout = os.path.basename(params.lstin)+".dat"

    xac_files = read_path_list(params.lstin)
    ofs_dat = open(params.datout, "w")

    ref_v6cell = None
    if params.ref_cell is not None:
        ref_v6cell = v6cell(uctbx.unit_cell(params.ref_cell).niggli_cell())
        ofs_dat.write("# ref_cell= %s\n" % params.ref_cell)

    if params.n_residues is not None: ofs_dat.write("# n_residues= %d\n" % params.n_residues)

    ofs_dat.write("file ioversigma resnatsnr1 wilsonb abdist a b c al be ga\n")

    ret = easy_mp.pool_map(fixed_func=lambda x: calc_stats(x, params.stats, params.n_residues, ref_v6cell,
                                                           params.min_peak, params.min_peak_percentile,
                                                           params.correct_peak),
                           args=xac_files,
                           processes=params.nproc)
    
    for stat in ret:
        getornan = lambda x: stat.get(x, float("nan")) # get or nan
        ofs_dat.write("%s %.3f %.3f %.3f %.3e"%(stat["filename"],
                                                getornan("ioversigma"), getornan("resnatsnr1"),
                                                getornan("wilsonb"), getornan("abdist")))
        ofs_dat.write(" %.3f %.3f %.3f %.2f %.2f %.2f\n" % stat["cell"])

    ofs_dat.close()

# run()

if __name__ == "__main__":
    import sys

    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print "All parameters:\n"
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        quit()

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if os.path.isfile(arg) and params.lstin is None:
            params.lstin = arg

    if params.lstin is None:
        print "Give .lst file"
        quit()

    run(params)
