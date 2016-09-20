"""
Usage:
yamtbx.python test_indexing_mode_with_reference.py \
  all-amb.stream \
  refine_001.mtz \
  "F-model,PHIF-model" \
  "h,k,l" "k,h,-l" \
| tee test_index.dat

"""

from cctbx import sgtbx
from cctbx.array_family import flex
import iotbx.file_reader
import re
from yamtbx.dataproc.crystfel.command_line.stats_stream_savememory import read_stream

def get_reference_data(refdata, reflabel):
    arrays = iotbx.file_reader.any_file(refdata).file_server.miller_arrays
    array = filter(lambda x: x.info().label_string()==reflabel, arrays)[0]
    print "# Read", array.info()
    array = array.resolution_filter(d_min=2)
    if array.anomalous_flag(): array = array.average_bijvoet_mates()
    return array.as_intensity_array()
# get_reference_data()

def calc_cc(a1, a2):
    a1, a2 = a1.common_sets(a2, assert_is_similar_symmetry=False)
    corr = flex.linear_correlation(a1.data(), a2.data())

    if corr.is_well_defined():# and a1.size() > 20:
        return corr.coefficient(), a1.size()
    else:
        return float("nan"), a1.size()
# calc_cc()

def run(stream_in, refmtz_in, refmtz_lab, ops):
    refarr = get_reference_data(refmtz_in, refmtz_lab)

    print "# Stream_in=", stream_in

    print "file event dmin",
    for op in ops: print "%s n_%s" % (op,op),
    print

    ops = map(lambda x: sgtbx.change_of_basis_op(x), ops)

    refarrs = map(lambda op: refarr.customized_copy(indices=op.apply(refarr.indices())).map_to_asu(), ops)

    for chunk in read_stream(stream_in):
        datarr = chunk.data_array(refarr.space_group(), False).merge_equivalents(use_internal_variance=False).array()
        #datarr = datarr.resolution_filter(d_min=chunk.res_lim)
        ccs = []
        for r in refarrs:
            #if op.is_identity_op(): tmp = datarr
            #else: tmp = datarr.customized_copy(indices=op.apply(datarr.indices())).map_to_asu()
            ccs.append(calc_cc(datarr, r))

        print "%s %s %.2f %s" % (chunk.filename, chunk.event, chunk.res_lim, " ".join(map(lambda x: "% .4f %4d" % x, ccs)))
# run()        

if __name__ == "__main__":
    import sys
    stream_in = sys.argv[1]
    refmtz_in = sys.argv[2]
    refmtz_lab = sys.argv[3]
    ops = sys.argv[4:]
    run(stream_in, refmtz_in, refmtz_lab, ops)
