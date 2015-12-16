"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.auto.command_line import multi_check_cell_consistency
from yamtbx.dataproc.auto.command_line.run_all_xds_simple import run_xds
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.xds import files as xds_files
from yamtbx.dataproc.xds import correctlp
from yamtbx.dataproc.xds import modify_xdsinp, make_backup, revert_files
from yamtbx.dataproc.xds.xparm import XPARM
from yamtbx import util

import iotbx.phil
import libtbx.phil
from libtbx.utils import multi_out
from cctbx.crystal import reindex
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex

import os
import sys
import shutil
import numpy

master_params_str = """
topdir = None
 .type = path
 .help = top directory containing xds results
lstout = "formerging.lst"
 .type = path
 .help = data list for merging
space_group = None
 .type = str
 .help = Choose the space group (name or number).
group_choice = None
 .type = int
 .help = Choose the group (run once and choose).
reference_for_reindex = None
 .type = path
 .help = Reference reflection data for resolving indexing ambiguity
"""

def rescale_with_specified_symm(topdir, dirs, symms, out, sgnum=None, reference_symm=None):
    assert (sgnum, reference_symm).count(None) == 1

    if sgnum is not None:
        sgnum_laue = sgtbx.space_group_info(sgnum).group().build_derived_reflection_intensity_group(False).type().number()

        matches = filter(lambda x:x.reflection_intensity_symmetry(False).space_group_info().type().number()==sgnum_laue, symms)
        matched_cells = numpy.array(map(lambda x: x.unit_cell().parameters(), matches))
        median_cell = map(lambda x: numpy.median(matched_cells[:,x]), xrange(6))

        reference_symm = crystal.symmetry(median_cell, sgnum)
    else:
        sgnum = reference_symm.space_group_info().type().number()
        sgnum_laue = reference_symm.space_group().build_derived_reflection_intensity_group(False).type().number()
    
    print >>out
    print >>out,  "Re-scaling with specified symmetry:", reference_symm.space_group_info().symbol_and_number()
    print >>out,  " reference cell:", reference_symm.unit_cell()
    print >>out
    print >>out

    cells = {} # cell and file
    for sym, wd in zip(symms, dirs):
        print >>out,  os.path.relpath(wd, topdir),

        # Find appropriate data
        xac_file = util.return_first_found_file(("XDS_ASCII.HKL_noscale.org", "XDS_ASCII.HKL_noscale", 
                                                 "XDS_ASCII_fullres.HKL.org", "XDS_ASCII_fullres.HKL",
                                                 "XDS_ASCII.HKL.org", "XDS_ASCII.HKL"),
                                                wd=wd)
        if xac_file is None:
            print >>out, "Can't find XDS_ASCII file in %s" % wd
            continue

        xac = XDS_ASCII(xac_file, read_data=False)
        print >>out, "%s %s (%s)" % (os.path.basename(xac_file), xac.symm.space_group_info(),
                                     ",".join(map(lambda x: "%.2f"%x, xac.symm.unit_cell().parameters())))

        if xac.symm.reflection_intensity_symmetry(False).space_group_info().type().number() == sgnum_laue:
            if xac.symm.unit_cell().is_similar_to(reference_symm.unit_cell(), 0.1, 10):
                print >>out,  "  Already scaled with specified symmetry"
                cells[wd] = (numpy.array(xac.symm.unit_cell().parameters()), xac_file)
                continue

        xdsinp = os.path.join(wd, "XDS.INP")
        cosets = reindex.reindexing_operators(reference_symm, xac.symm, 0.2, 20)

        if len(cosets.combined_cb_ops())==0:
            print >>out, "Can't find operator:"
            sym.show_summary(out, " ")
            reference_symm.show_summary(out, " ")
            continue

        newcell = reference_symm.space_group().average_unit_cell(xac.symm.change_basis(cosets.combined_cb_ops()[0]).unit_cell())
        newcell = " ".join(map(lambda x: "%.3f"%x, newcell.parameters()))
        print >>out,  "Scaling with transformed cell:", newcell

        #for f in xds_files.generated_by_CORRECT:
        #    util.rotate_file(os.path.join(wd, f))
        bk_prefix = make_backup(xds_files.generated_by_CORRECT, wdir=wd, quiet=True)

        modify_xdsinp(xdsinp, inp_params=[("JOB", "CORRECT"),
                                          ("SPACE_GROUP_NUMBER", "%d"%sgnum),
                                          ("UNIT_CELL_CONSTANTS", newcell),
                                          ("INCLUDE_RESOLUTION_RANGE", "50 0"),
                                          ("CORRECTIONS", ""),
                                          ("NBATCH", "1"),
                                          ("MINIMUM_I/SIGMA", "50"),
                                          ("REFINE(CORRECT)", "")])
        run_xds(wd)
        for f in ("CORRECT.LP", "XDS_ASCII.HKL", "GXPARM.XDS"):
            if os.path.exists(os.path.join(wd, f)):
                shutil.copyfile(os.path.join(wd, f), os.path.join(wd, f+"_rescale"))

        revert_files(xds_files.generated_by_CORRECT, bk_prefix, wdir=wd, quiet=True)

        new_xac = os.path.join(wd, "XDS_ASCII.HKL_rescale")
        new_gxparm = os.path.join(wd, "GXPARM.XDS_rescale")
        if os.path.isfile(new_xac) and os.path.isfile(new_gxparm):
            cells[wd] = (XPARM(new_gxparm).unit_cell, new_xac)
            print "OK:", cells[wd][0]
        else:
            print >>out, "Error: rescaling failed (Can't find XDS_ASCII.HKL)"
            continue

    return cells, reference_symm
# rescale_with_specified_symm()

def reindex_with_specified_symm(topdir, reference_symm, dirs, out):
    print >>out
    print >>out,  "Re-index to specified symmetry:"
    reference_symm.show_summary(out, "  ")
    print >>out
    print >>out

    cells = {} # cell and file

    sgnum_laue = reference_symm.space_group().build_derived_reflection_intensity_group(False).type().number()

    for wd in dirs:
        print >>out, "%s:" % os.path.relpath(wd, topdir),

        # Find appropriate data
        xac_file = util.return_first_found_file(("XDS_ASCII.HKL_noscale.org", "XDS_ASCII.HKL_noscale", 
                                                 "XDS_ASCII_fullres.HKL.org", "XDS_ASCII_fullres.HKL",
                                                 "XDS_ASCII.HKL.org", "XDS_ASCII.HKL"),
                                                wd=wd)
        if xac_file is None:
            print >>out, "Can't find XDS_ASCII file in %s" % wd
            continue

        if xac_file.endswith(".org"): xac_file_org, xac_file = xac_file, xac_file[:-4]
        else: xac_file_org = xac_file+".org"

        if not os.path.isfile(xac_file_org):
            os.rename(xac_file, xac_file_org)

        xac = XDS_ASCII(xac_file_org, read_data=False)
        print >>out, "%s %s (%s)" % (os.path.basename(xac_file), xac.symm.space_group_info(),
                                   ",".join(map(lambda x: "%.2f"%x, xac.symm.unit_cell().parameters())))

        if xac.symm.reflection_intensity_symmetry(False).space_group_info().type().number() == sgnum_laue:
            if xac.symm.unit_cell().is_similar_to(reference_symm.unit_cell(), 0.1, 10):
                print >>out,  "  Already scaled with specified symmetry"
                os.rename(xac_file_org, xac_file) # rename back
                cells[wd] = (numpy.array(xac.symm.unit_cell().parameters()), xac_file)
                continue

        cosets = reindex.reindexing_operators(reference_symm, xac.symm, 0.2, 20)

        if len(cosets.combined_cb_ops())==0:
            print >>out, "Can't find operator:"
            xac.symm.show_summary(out, " ")
            reference_symm.show_summary(out, " ")
            continue

        newcell = xac.write_reindexed(op=cosets.combined_cb_ops()[0],
                                      space_group=reference_symm.space_group(),
                                      hklout=xac_file)
        cells[wd] = (numpy.array(newcell.parameters()), xac_file)

        newcell = " ".join(map(lambda x: "%.3f"%x, newcell.parameters()))
        print >>out,  "  Reindexed to transformed cell: %s with %s" % (newcell, cosets.combined_cb_ops()[0].as_hkl())

    return cells
# reindex_with_specified_symm()

def read_strong_i_from_xds_ascii(xds_ascii_in):
    tmp = XDS_ASCII(xds_ascii_in, i_only=True).i_obs(anomalous_flag=False)
    sel = tmp.sigmas() > 0
    sel &= tmp.data()/tmp.sigmas() > 2
    sel &= tmp.d_spacings() > 3
    if sel.count(True) < 10:
        return None
    tmp = tmp.select(sel)
    merge = tmp.merge_equivalents(use_internal_variance=False)
    return merge.array()
# read_strong_i_from_xds_ascii()

def get_cc(lhs, rhs):
    di, dj = lhs.common_sets(rhs, assert_is_similar_symmetry=False)
    corr = flex.linear_correlation(di.data(), dj.data())
    if corr.is_well_defined():
        return corr.coefficient()
    return None
# get_cc()

def resolve_indexing_ambiguity(dirs, reidx_ops):
    """
    This very simple implementation would be buggy and gangerous.
    Don't use this now!!
    """

    data = {}

    # Read all (strong) data
    for wd in dirs:
        print "reading", wd
        tmp = read_strong_i_from_xds_ascii(os.path.join(wd, "XDS_ASCII.HKL"))
        if tmp is None:
            continue
        data[wd] = tmp

    # Pairwise analysis
    best_op = map(lambda x: (0, 0), xrange(len(dirs)))
    for i in xrange(len(dirs)-1):
        for j in xrange(i+1, len(dirs)):
            wd_i, wd_j = dirs[i], dirs[j]
            if wd_i not in data or wd_j not in data:
                continue
            data_i, data_j = data[wd_i], data[wd_j]
            cc_list = []
            idx_min, idx_max = 0, 0
            for k, rop in enumerate(reidx_ops):
                if k == 0:
                    data_j_ = data_j
                else:
                    data_j_ = data_j.customized_copy(indices=rop.apply(data_j.indices()))
                cc_list.append(get_cc(data_i, data_j_))
                if k > 0 and cc_list[k] is not None:
                    if cc_list[k] < cc_list[idx_min]: idx_min = k
                    if cc_list[k] > cc_list[idx_max]: idx_max = k

            if len(cc_list) - cc_list.count(None) < 2:
                print i,j,"skip"
                continue
            cc_min, cc_max = cc_list[idx_min], cc_list[idx_max]
            score = (cc_max+1)/(cc_min+1) if cc_min+1 != 0 else float("nan")
            print i, j, 
            for cc in cc_list: print cc, 
            print score

            if score > best_op[j][1]:
                best_op[j] = (idx_max, score)
            if score > best_op[i][1]:
                best_op[i] = (0, score)

    for i in xrange(len(dirs)):
        print "%.3d op= %d score= %.2f" % (i, best_op[i][0], best_op[i][1])

    # Apply reindexing

# resolve_indexing_ambiguity()

def resolve_indexing_ambiguity_using_reference(dirs, reidx_ops, reference_file):
    from iotbx import reflection_file_reader

    hkl_ref = filter(lambda x: x.is_xray_intensity_array(),
                     reflection_file_reader.any_reflection_file(reference_file).as_miller_arrays(merge_equivalents=False))
    if len(hkl_ref) == 0:
        raise Exception("No intensity data in %s"%reference_file)
    
    ref_data = hkl_ref[0].merge_equivalents(use_internal_variance=False).array()
    data = {}
    
    # Read all (strong) data
    for wd in dirs:
        print "reading", wd
        tmp = read_strong_i_from_xds_ascii(os.path.join(wd, "XDS_ASCII.HKL"))
        if tmp is None:
            continue
        data[wd] = tmp

    # Compare with reference
    best_op = map(lambda x: 0, xrange(len(dirs)))
    for i in xrange(len(dirs)):
        if dirs[i] not in data:
            continue

        data_i = data[dirs[i]]
        cc_list = []
        idx_max = 0
        for k, rop in enumerate(reidx_ops):
            if k == 0:
                data_i_ = data_i
            else:
                data_i_ = data_i.customized_copy(indices=rop.apply(data_i.indices()))
            cc_list.append(get_cc(data_i_, ref_data))
            if cc_list[k] is not None:
                if cc_list[k] > cc_list[idx_max]: idx_max = k

        if len(cc_list) - cc_list.count(None) < 2:
            print i,"skip"
            continue
        cc_max = cc_list[idx_max]
        print i,
        for cc in cc_list: print cc, 
        print
        best_op[i] = idx_max
    
    print "\n Best reindexing operators"
    print "=================================\n"
    for i in xrange(len(dirs)):
        print "%.3d op= %d %s" % (i, best_op[i], dirs[i])

    print "\nApplying.."
    for i, opi in enumerate(best_op):
        if opi == 0: continue
        # XXX treat!!
    print "Done."

# resolve_indexing_ambiguity_using_reference()

def run(params):
    out = multi_out()
    out.register("log", open(os.path.join(os.path.dirname(params.lstout), "multi_prep_merging.log"), "w"), atexit_send_to=None)
    out.register("stdout", sys.stdout)

    cell_params = libtbx.phil.parse(input_string=multi_check_cell_consistency.master_params_str).extract()
    cell_params.topdir = params.topdir
    cm = multi_check_cell_consistency.run(cell_params, out=out)

    if len(cm.groups) == 0:
        print "Oh, no. No data."
        return

    if params.space_group is not None and params.group_choice is None:
        params.group_choice = 1 # maybe the largest group.

    if params.group_choice < 1 or len(cm.groups) < params.group_choice:
        print "Invalid group_choice=. Should be in 1..%d" % len(cm.groups)
        return

    possible_sgs = set(map(lambda x: cm.symms[x].space_group(), cm.groups[params.group_choice-1]))
    if params.space_group is None:
        print "Please specify space_group=. The possible choice is:"
        print "\n".join(map(lambda x: x.info().symbol_and_number(), possible_sgs))
        return
    try:
        if sgtbx.space_group_info(params.space_group).group() not in possible_sgs:
            print "Invalid space group choice. The possible choice is:"
            print "\n".join(map(lambda x: x.info().symbol_and_number(), possible_sgs))
            return
    except RuntimeError:
        print "Invalid space group name or number (%s)" % params.space_group
        return

    symms = map(lambda i: cm.symms[i], cm.groups[params.group_choice-1])
    dirs = map(lambda i: cm.dirs[i], cm.groups[params.group_choice-1])

    sgnum = sgtbx.space_group_info(params.space_group).group().type().number()

    # 1. Scale with specified symmetry
    rescale_with_specified_symm(params.topdir, dirs, symms, out, sgnum=sgnum)

    # 2. Resolve reindexing problem
    # TODO: reconstruct unit cell by averaging
    reference_symm = filter(lambda x:x.reflection_intensity_symmetry(False).space_group_info().type().number()==sgnum_laue, symms)[0]

    cosets = reindex.reindexing_operators(reference_symm, reference_symm)
    reidx_ops = cosets.combined_cb_ops()
    print " Reference symmetry:", reference_symm.unit_cell(), reference_symm.space_group_info().symbol_and_number(), 
    print " Possible reindex operators:", map(lambda x: str(x.as_hkl()), reidx_ops)
    if len(reidx_ops) == 1:
        print " No reindexing problem exists."
    elif params.reference_for_reindex is not None:
        resolve_indexing_ambiguity_using_reference(dirs, reidx_ops, params.reference_for_reindex)
    else:
        resolve_indexing_ambiguity(dirs, reidx_ops)

    ofs = open(params.lstout, "w")
    for wd in dirs:
        xas_full = os.path.join(wd, "XDS_ASCII_fullres.HKL")
        if os.path.isfile(xas_full): # if _fullres.HKL exists, use it.
            ofs.write("%s\n" % xas_full)
        else:
            ofs.write("%s\n" % os.path.join(wd, "XDS_ASCII.HKL"))
# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    
    for arg in args:
        if os.path.isdir(arg) and params.topdir is None:
            params.topdir = arg

    if params.topdir is None:
        params.topdir = os.getcwd()

    run(params)
