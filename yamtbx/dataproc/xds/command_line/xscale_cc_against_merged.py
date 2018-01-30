"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os
import collections
import copy
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc.xds import xds_ascii
from cctbx.array_family import flex
from libtbx import easy_mp

def eval_cc_with_original_file(f, merged_iobs):
    print "reading",f

    xac = xds_ascii.XDS_ASCII(f)
    iobs = xac.i_obs(anomalous_flag=merged_iobs.anomalous_flag()).merge_equivalents(use_internal_variance=False).array()

    n_all = iobs.size()

    m, i = merged_iobs.common_sets(iobs, assert_is_similar_symmetry=False)
    n_common = m.size()

    corr = flex.linear_correlation(m.data(), i.data())
    cc = corr.coefficient() if corr.is_well_defined() else float("nan")

    ret1 = (n_all, n_common, cc)
    ret2 = []

    for frame in xrange(min(xac.iframe), max(xac.iframe)+1):
        sel = xac.iframe == frame
        iobs = xac.i_obs(anomalous_flag=merged_iobs.anomalous_flag()).select(sel)
        iobs = iobs.merge_equivalents(use_internal_variance=False).array()
        n_all = iobs.size()

        m, i = merged_iobs.common_sets(iobs, assert_is_similar_symmetry=False)
        n_common = m.size()

        corr = flex.linear_correlation(m.data(), i.data())
        cc = corr.coefficient() if corr.is_well_defined() else float("nan")

        ret2.append([frame, n_all, n_common, cc])
    return ret1, ret2
# eval_cc_with_original_file()

def eval_cc_internal(merged_iobs, merged, setno):
    """
    merged_iobs: i_obs array where all data in xscale hkl were merged
    merged: xscale hkl object
    setno: iset number in interest
    """
    merged_sel = copy.copy(merged)
    merged_sel.remove_selection(merged.iset != setno)
    
    iobs_set = merged_sel.i_obs(anomalous_flag=merged_iobs.anomalous_flag())
    iobs_set_merged = iobs_set.merge_equivalents(use_internal_variance=False).array()

    n_all = iobs_set_merged.size()

    m, i = merged_iobs.common_sets(iobs_set_merged, assert_is_similar_symmetry=False)
    n_common = m.size()

    corr = flex.linear_correlation(m.data(), i.data())
    cc = corr.coefficient() if corr.is_well_defined() else float("nan")

    ret1 = (n_all, n_common, cc)
    ret2 = []

    for frame in xrange(min(merged_sel.iframe), max(merged_sel.iframe)+1):
        iobs = iobs_set.select(merged_sel.iframe == frame)
        iobs = iobs.merge_equivalents(use_internal_variance=False).array()
        n_all = iobs.size()

        m, i = merged_iobs.common_sets(iobs, assert_is_similar_symmetry=False)
        n_common = m.size()

        corr = flex.linear_correlation(m.data(), i.data())
        cc = corr.coefficient() if corr.is_well_defined() else float("nan")

        ret2.append([frame, n_all, n_common, cc])
    return ret1, ret2
# eval_cc_internal()


def run(hklin, output_dir=None, eval_internal=True):
    if output_dir is None: output_dir = os.getcwd()

    merged = xds_ascii.XDS_ASCII(hklin)
    merged_iobs = merged.i_obs().merge_equivalents(use_internal_variance=False).array()

    fwidth = max(map(lambda x: len(x[0]), merged.input_files.values()))
    formatf = "%"+str(fwidth)+"s"

    out_files = open(os.path.join(output_dir, "cc_files.dat"), "w")
    out_frames = open(os.path.join(output_dir, "cc_frames.dat"), "w")

    print >>out_files, "file name n.all n.common cc"
    print >>out_frames, "file name frame n.all n.common cc"

    cutforname1 = len(os.path.commonprefix(map(lambda x: x[0], merged.input_files.values())))
    cutforname2 = len(os.path.commonprefix(map(lambda x: x[0][::-1], merged.input_files.values())))
    formatn = "%"+str(fwidth-cutforname1-cutforname2)+"s"

    if eval_internal:
        results = map(lambda x: eval_cc_internal(merged_iobs, merged, x), sorted(set(merged.iset)))
    else:
        results = map(lambda x: eval_cc_with_original_file(x, merged_iobs),
                      map(lambda x: x[0] if os.path.isabs(x[0]) else os.path.join(os.path.dirname(hklin), x[0]),
                          merged.input_files.values()))

    ret = collections.OrderedDict()

    for (f, wavelen, cell), (ret1, ret2) in zip(merged.input_files.values(), results):
        name = f[cutforname1+1:-cutforname2]
        n_all, n_common, cc = ret1
        print >>out_files, formatf%f, formatn%name, "%5d %5d %.4f" % (n_all, n_common, cc)
        ret[f] = []

        for frame, n_all, n_common, cc in ret2:
            print >>out_frames, formatf%f, formatn%name, "%6d %5d %5d %.4f" % (frame, n_all, n_common, cc)
            ret[f].append([frame, n_all, n_common, cc])

    return ret
# run()

if __name__ == "__main__":
    import sys

    hklin = sys.argv[1]
    run(hklin)
    
