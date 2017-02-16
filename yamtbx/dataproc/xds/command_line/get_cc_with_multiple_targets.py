import iotbx.phil
import iotbx.file_reader
import iotbx.reflection_file_utils
from cctbx.array_family import flex
from mmtbx.scaling.absolute_scaling import kernel_normalisation

from yamtbx import util
from yamtbx.util.xtal import CellConstraints
from yamtbx.dataproc.xds import xds_ascii

import collections

master_params_str = """\
lstin = None
 .type = path
 .help = "List of XDS_ASCII.HKL"
dat_out = "cc_with_targets.dat"
 .type = path
 .help = "Output file"
normalization = *no E rel_B wilson_B
 .type = choice(multi=False)
 .help = "Normalization of each intensities"
d_min = None
 .type = float
d_max = None
 .type = float
"""

def calc_cc(ari, arj):
  ari, arj = ari.common_sets(arj, assert_is_similar_symmetry=False)
  corr = flex.linear_correlation(ari.data(), arj.data())
  if corr.is_well_defined():
      return corr.coefficient(), ari.size()
  else:
      return float("nan"), ari.size()
# calc_cc()

def read_target_files(target_files, d_min, d_max, normalization, log_out):
    ret = collections.OrderedDict()
    for i, f in enumerate(target_files):
        f = iotbx.file_reader.any_file(f, force_type="hkl", raise_sorry_if_errors=True)
        arrays = f.file_server.get_miller_arrays(None)
        scores = iotbx.reflection_file_utils.get_xray_data_scores(arrays, ignore_all_zeros=True,
                                                                  prefer_anomalous=False, prefer_amplitudes=False)
        array = arrays[scores.index(max(scores))]

        log_out.write("# target%.3d = %s %s\n" % (i, array.info(), array.d_max_min()))
        
        if array.anomalous_flag(): array = array.average_bijvoet_mates()
        array = array.as_intensity_array().resolution_filter(d_max=d_max, d_min=d_min)
        
        if normalization == "E":
            normaliser = kernel_normalisation(array, auto_kernel=True)
            ret[f] = array.customized_copy(data=array.data()/normaliser.normalizer_for_miller_array,
                                           sigmas=array.sigmas()/normaliser.normalizer_for_miller_array if array.sigmas() else None)
        else:
            ret[f] = array

    return ret
# read_target_files()

def run(params, target_files):
    assert params.normalization in ("no", "E")
    ofs = open(params.dat_out, "w")

    xac_files = util.read_path_list(params.lstin)
    targets = read_target_files(target_files, params.d_min, params.d_max, params.normalization, ofs)

    cellcon = CellConstraints(targets.values()[0].space_group())
    
    #for i, t in enumerate(targets): ofs.write("# target%.3d = %s\n" % (i,t))
    ofs.write("# normalization = %s\n" % params.normalization)
    ofs.write("# d_min, d_max = %s, %s\n" % (params.d_min, params.d_max))
    ofs.write("file %s " % cellcon.get_label_for_free_params())
    ofs.write(" ".join(map(lambda x: "cc.%.3d nref.%.3d"%(x,x), xrange(len(targets)))))
    ofs.write("\n")
    
    for xac_file in xac_files:
        print "reading", xac_file
        xac = xds_ascii.XDS_ASCII(xac_file)
        xac.remove_rejected()
        iobs = xac.i_obs(anomalous_flag=False).merge_equivalents(use_internal_variance=False).array()
        ofs.write("%s %s" % (xac_file, cellcon.format_free_params(iobs.unit_cell())))
        fail_flag = False
        if params.normalization == "E":
            try:
                normaliser = kernel_normalisation(iobs, auto_kernel=True)
                iobs = iobs.customized_copy(data=iobs.data()/normaliser.normalizer_for_miller_array,
                                            sigmas=iobs.sigmas()/normaliser.normalizer_for_miller_array)
            except:
                fail_flag = True

        for i, ta in enumerate(targets.values()):
            if fail_flag:
                ofs.write(" % .4f %4d" % cc_num)
            else:
                cc_num = calc_cc(iobs, ta)
                ofs.write(" % .4f %4d" % cc_num)
            
        ofs.write("\n")
# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    targets = cmdline.remaining_args
    
    run(params, targets)
