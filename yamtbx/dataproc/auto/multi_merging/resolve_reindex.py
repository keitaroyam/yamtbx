"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.util.xtal import format_unit_cell
from cctbx import crystal
from cctbx.crystal import reindex
from cctbx.array_family import flex
from cctbx import sgtbx
from libtbx.utils import null_out
from libtbx import easy_mp
from libtbx import adopt_init_args
from cctbx.merging import brehm_diederichs

import os
import copy
import multiprocessing
import time
import numpy

def calc_cc(a1, a2):
    a1, a2 = a1.common_sets(a2, assert_is_similar_symmetry=False)
    corr = flex.linear_correlation(a1.data(), a2.data())

    if corr.is_well_defined():# and a1.size() > 20:
        return corr.coefficient()
    else:
        return float("nan")
# calc_cc()

class ReindexResolver:
    def __init__(self, xac_files, d_min=3, min_ios=3, nproc=1, max_delta=5, log_out=null_out()):
        adopt_init_args(self, locals())
        self.arrays = []
        self.best_operators = None
        self._representative_xs = None
        self.bad_files = []
    # __init__()

    def representative_crystal_symmetry(self): return self._representative_xs

    def read_xac_files(self, from_p1=False):
        op_to_p1 = None
        if from_p1:
            """
            This option is currently for multi_determine_symmetry.
            Do not use this for ambiguity resolution! op_to_p1 is not considered when writing new HKL files.
            """
            self.log_out.write("\nAveraging symmetry of all inputs..\n")
            cells = []
            sgs = []
            for f in self.xac_files:
                xac =  XDS_ASCII(f, read_data=False)
                cells.append(xac.symm.unit_cell().parameters())
                sgs.append(xac.symm.space_group())
            assert len(set(sgs)) < 2
            avg_symm = crystal.symmetry(list(numpy.median(cells, axis=0)), space_group=sgs[0])
            op_to_p1 = avg_symm.change_of_basis_op_to_niggli_cell()
            self.log_out.write("  Averaged symmetry: %s (%s)\n" % (format_unit_cell(avg_symm.unit_cell()), sgs[0].info()))
            self.log_out.write("  Operator to Niggli cell: %s\n" % op_to_p1.as_hkl())
            self.log_out.write("        Niggli cell: %s\n" % format_unit_cell(avg_symm.unit_cell().change_basis(op_to_p1)))

        print >>self.log_out, "\nReading"
        cells = []
        bad_files, good_files = [], []
        for i, f in enumerate(self.xac_files):
            print >>self.log_out, "%4d %s" % (i, f)
            xac = XDS_ASCII(f, i_only=True)
            self.log_out.write("     d_range: %6.2f - %5.2f" % xac.i_obs().resolution_range())
            self.log_out.write(" n_ref=%6d" % xac.i_obs().size())
            xac.remove_rejected()
            a = xac.i_obs().resolution_filter(d_min=self.d_min)
            if self.min_ios is not None: a = a.select(a.data()/a.sigmas()>=self.min_ios)
            self.log_out.write(" n_ref_filtered=%6d" % a.size())
            if from_p1:
                a = a.change_basis(op_to_p1).customized_copy(space_group_info=sgtbx.space_group_info("P1"))
            a = a.as_non_anomalous_array().merge_equivalents(use_internal_variance=False).array()
            self.log_out.write(" n_ref_merged=%6d\n" % a.size())
            if a.size() < 2:
                self.log_out.write("     !! WARNING !! number of reflections is dangerously small!!\n")
                bad_files.append(f)
            else:
                self.arrays.append(a)
                cells.append(a.unit_cell().parameters())
                good_files.append(f)

        if bad_files:
            self.xac_files = good_files
            self.bad_files = bad_files

        assert len(self.xac_files) == len(self.arrays) == len(cells)
            
        print >>self.log_out, ""

        self._representative_xs = crystal.symmetry(list(numpy.median(cells, axis=0)),
                                                   space_group_info=self.arrays[0].space_group_info())
    # read_xac_files()

    def show_assign_summary(self, log_out=None):
        if not log_out: log_out = self.log_out
        if not self.best_operators:
            log_out.write("ERROR: Operators not assigned.\n")
            return

        unique_ops = set(self.best_operators)
        op_count = map(lambda x: (x, self.best_operators.count(x)), unique_ops)
        op_count.sort(key=lambda x: x[1])
        log_out.write("Assigned operators:\n")
        for op, num in reversed(op_count):
            log_out.write("  %10s: %4d\n" % (op.as_hkl(), num))
        log_out.write("\n")
    # show_assign_summary()

    def find_reindex_ops(self):
        symm = self.representative_crystal_symmetry()
        cosets = reindex.reindexing_operators(symm, symm, max_delta=self.max_delta)
        reidx_ops = cosets.combined_cb_ops()
        return reidx_ops
    # find_reindex_ops()

    def modify_xds_ascii_files(self, suffix="_reidx", cells_dat_out=None):
        #ofs_lst = open("for_merge_new.lst", "w")
        if cells_dat_out: cells_dat_out.write("file a b c al be ga\n")

        new_files = []
        print >>self.log_out, "Writing reindexed files.."
        assert len(self.xac_files) == len(self.best_operators)
        for i, (f, op) in enumerate(zip(self.xac_files, self.best_operators)):
            xac = XDS_ASCII(f, read_data=False)
            if op.is_identity_op():
                new_files.append(f)
                if cells_dat_out:
                    cell = xac.symm.unit_cell().parameters()
                    cells_dat_out.write(f+" "+" ".join(map(lambda x:"%7.3f"%x, cell))+"\n")

                continue

            newf = f.replace(".HKL", suffix+".HKL") if ".HKL" in f else os.path.splitext(f)[0]+suffix+".HKL"
            print >>self.log_out, "%4d %s" % (i, newf)

            cell_tr = xac.write_reindexed(op, newf, space_group=self.arrays[0].crystal_symmetry().space_group())
            #ofs_lst.write(newf+"\n")
            new_files.append(newf)

            if cells_dat_out:
                cells_dat_out.write(newf+" "+" ".join(map(lambda x:"%7.3f"%x, cell_tr.parameters()))+"\n")

        return new_files
    # modify_xds_ascii_files()

    def debug_write_mtz(self):
        arrays = self.arrays
        
        merge_org = arrays[0].deep_copy()
        for a in arrays[1:]: merge_org = merge_org.concatenate(a, assert_is_similar_symmetry=False)
        merge_org = merge_org.merge_equivalents(use_internal_variance=False).array()
        merge_org.as_mtz_dataset(column_root_label="I").mtz_object().write(file_name="noreindex.mtz")

        merge_new = None
        for a, op in zip(arrays, self.best_operators):
            if not op.is_identity_op(): a = a.customized_copy(indices=op.apply(a.indices()))

            if merge_new is None: merge_new = a
            else: merge_new = merge_new.concatenate(a, assert_is_similar_symmetry=False)

        merge_new = merge_new.merge_equivalents(use_internal_variance=False).array()
        merge_new.as_mtz_dataset(column_root_label="I").mtz_object().write(file_name="reindexed.mtz")
    # debug_write_mtz()
# class ReindexResolver

def kabsch_selective_breeding_worker(args):
    ref, i, tmp, new_ops = args
    global kabsch_selective_breeding_worker_dict
    nproc = kabsch_selective_breeding_worker_dict["nproc"]
    arrays = kabsch_selective_breeding_worker_dict["arrays"]
    reindexed_arrays = kabsch_selective_breeding_worker_dict["reindexed_arrays"]
    reidx_ops = kabsch_selective_breeding_worker_dict["reidx_ops"]

    if ref==i: return None
    
    if nproc > 1:
        tmp2 = reindexed_arrays[new_ops[ref]][ref]
    else:                                
        if reidx_ops[new_ops[ref]].is_identity_op(): tmp2 = arrays[ref]
        else: tmp2 = arrays[ref].customized_copy(indices=reidx_ops[new_ops[ref]].apply(arrays[ref].indices())).map_to_asu()

    cc = calc_cc(tmp, tmp2)
    if cc==cc: return cc
    return None
# work_local()

class KabschSelectiveBreeding(ReindexResolver):
    """
    Reference: W. Kabsch "Processing of X-ray snapshots from crystals in random orientations" Acta Cryst. (2014). D70, 2204-2216
    http://dx.doi.org/10.1107/S1399004714013534
    If I understand correctly...
    """

    def __init__(self, xac_files, d_min=3, min_ios=3, nproc=1, max_delta=5, from_p1=False, log_out=null_out()):
        ReindexResolver.__init__(self, xac_files, d_min, min_ios, nproc, max_delta, log_out)
        self._final_cc_means = [] # list of [(op_index, cc_mean), ...]
        self._reidx_ops = []
        self.read_xac_files(from_p1=from_p1)
    # __init__()

    def final_cc_means(self): return self._final_cc_means
    def reindex_operators(self): return self._reidx_ops

    def assign_operators(self, reidx_ops=None, max_cycle=100):
        arrays = self.arrays
        self.best_operators = None

        if reidx_ops is None: reidx_ops = self.find_reindex_ops()

        print >>self.log_out, "Reindex operators:"
        for i, op in enumerate(reidx_ops): print >>self.log_out, " %2d: %s" % (i, op.as_hkl())
        print >>self.log_out, ""

        reidx_ops.sort(key=lambda x: not x.is_identity_op()) # identity op to first
        self._reidx_ops = reidx_ops

        if self.nproc > 1:
            # consumes much memory.. (but no much benefits)
            reindexed_arrays = [arrays]
            for op in reidx_ops[1:]:
                reindexed_arrays.append(map(lambda x: x.customized_copy(indices=op.apply(x.indices())).map_to_asu(), arrays))
        else:
            reindexed_arrays = None

        old_ops = map(lambda x:0, xrange(len(arrays)))
        new_ops = map(lambda x:0, xrange(len(arrays)))

        global kabsch_selective_breeding_worker_dict
        kabsch_selective_breeding_worker_dict = dict(nproc=self.nproc, arrays=arrays,
                                                     reindexed_arrays=reindexed_arrays,
                                                     reidx_ops=reidx_ops) # constants during cycles
        pool = multiprocessing.Pool(self.nproc)
        for ncycle in xrange(max_cycle):
            #new_ops = copy.copy(old_ops) # doesn't matter
            self._final_cc_means = []

            for i in xrange(len(arrays)):
                cc_means = []
                a = arrays[i]
                #ttt=time.time()
                for j, op in enumerate(reidx_ops):
                    cc_list = []
                    if self.nproc > 1:
                        tmp = reindexed_arrays[j][i]
                    else:
                        if op.is_identity_op(): tmp = a
                        else: tmp = a.customized_copy(indices=op.apply(a.indices())).map_to_asu()
                    
                    """
                    def work_local(ref): # XXX This function is very slow when nproc>1...
                        if ref==i: return None

                        if self.nproc > 1:
                            tmp2 = reindexed_arrays[new_ops[ref]][ref]
                        else:                                
                            if reidx_ops[new_ops[ref]].is_identity_op(): tmp2 = arrays[ref]
                            else: tmp2 = arrays[ref].customized_copy(indices=reidx_ops[new_ops[ref]].apply(arrays[ref].indices())).map_to_asu()

                        cc = calc_cc(tmp, tmp2)
                        if cc==cc: return cc
                        return None
                    # work_local()

                    cc_list = easy_mp.pool_map(fixed_func=work_local,
                                               args=range(len(arrays)),
                                               processes=self.nproc)
                    """
                    cc_list = pool.map(kabsch_selective_breeding_worker,
                                       ((k, i, tmp, new_ops) for k in range(len(arrays))))
                    cc_list = filter(lambda x: x is not None, cc_list)

                    if len(cc_list) > 0:
                        cc_means.append((j, sum(cc_list)/len(cc_list)))
                        #print  >>self.log_out, "DEBUG:", i, j, cc_list, cc_means[-1]

                if cc_means:
                    max_el = max(cc_means, key=lambda x:x[1])
                    print >>self.log_out, "%3d %s" % (i, " ".join(map(lambda x: "%s%d:% .4f" % ("*" if x[0]==max_el[0] else " ", x[0], x[1]), cc_means)))
                    self._final_cc_means.append(cc_means)
                    #print "%.3f sec" % (time.time()-ttt)
                    new_ops[i] = max_el[0]
                else:
                    print >>self.log_out, "%3d %s Error! cannot calculate CC" % (i, " ".join(map(lambda x: " %d:    nan" % x, xrange(len(reidx_ops)))))
                    # XXX append something to self._final_cc_means?

            print >>self.log_out, "In %4d cycle" % (ncycle+1)
            print >>self.log_out, "  old",old_ops
            print >>self.log_out, "  new",new_ops
            print >>self.log_out, "  number of different assignments:", len(filter(lambda x:x[0]!=x[1], zip(old_ops,new_ops)))
            print >>self.log_out, ""
            if old_ops==new_ops:
                self.best_operators = map(lambda x: reidx_ops[x], new_ops)
                print >>self.log_out, "Selective breeding is finished in %d cycles" % (ncycle+1)
                return

            old_ops = copy.copy(new_ops)

        print >>self.log_out, "WARNING:: Selective breeding is not finished. max cycles reached."
        self.best_operators = map(lambda x: reidx_ops[x], new_ops) # better than nothing..
    # assign_operators()
# class KabschSelectiveBreeding

class ReferenceBased(ReindexResolver):
    def __init__(self, xac_files, ref_array, d_min=3, min_ios=3,  nproc=1, max_delta=5, log_out=null_out()):
        ReindexResolver.__init__(self, xac_files, d_min, min_ios, nproc, max_delta, log_out)
        self.read_xac_files()
        self.ref_array = ref_array.resolution_filter(d_min=d_min).as_non_anomalous_array().merge_equivalents(use_internal_variance=False).array()
        
    # __init__()

    def assign_operators(self, reidx_ops=None):
        arrays = self.arrays
        self.best_operators = None

        if reidx_ops is None: reidx_ops = self.find_reindex_ops()

        print >>self.log_out, "Reindex operators:", map(lambda x: str(x.as_hkl()), reidx_ops)
        print >>self.log_out, ""

        reidx_ops.sort(key=lambda x: not x.is_identity_op()) # identity op to first

        new_ops = map(lambda x:0, xrange(len(arrays)))

        for i, a in enumerate(arrays):
            cc_list = []
            a = arrays[i]

            for j, op in enumerate(reidx_ops):
                if op.is_identity_op(): tmp = a
                else: tmp = a.customized_copy(indices=op.apply(a.indices())).map_to_asu()

                cc = calc_cc(tmp, self.ref_array)
                if cc==cc: cc_list.append((j,cc))

            max_el = max(cc_list, key=lambda x:x[1])
            print >>self.log_out, "%3d %s" % (i, " ".join(map(lambda x: "%s%d:% .4f" % ("*" if x[0]==max_el[0] else " ", x[0], x[1]), cc_list)))
            new_ops[i] = max_el[0]

        print >>self.log_out, "  operator:", new_ops
        print >>self.log_out, "  number of different assignments:", len(filter(lambda x:x!=0, new_ops))
        print >>self.log_out, ""

        self.best_operators = map(lambda x: reidx_ops[x], new_ops)
    # assign_operators()
# class ReferenceBased

class BrehmDiederichs(ReindexResolver):
    def __init__(self, xac_files, d_min=3, min_ios=3, nproc=1, max_delta=5, log_out=null_out()):
        ReindexResolver.__init__(self, xac_files, d_min, min_ios, nproc, max_delta, log_out)
        self.read_xac_files()
    # __init__()

    def assign_operators(self, reidx_ops=None):
        arrays = self.arrays
        self.best_operators = None

        if reidx_ops is None: reidx_ops = self.find_reindex_ops()

        print >>self.log_out, "Reindex operators:", map(lambda x: str(x.as_hkl()), reidx_ops)
        print >>self.log_out, ""

        reidx_ops.sort(key=lambda x: not x.is_identity_op()) # identity op to first

        data = None
        latt_id = flex.int([])

        for i, a in enumerate(arrays):
            if data is None: data = a
            else: data = data.concatenate(a, assert_is_similar_symmetry=False)

            latt_id.extend(flex.int(a.size(), i))

        latt_id = data.customized_copy(data=latt_id.as_double())
        result = brehm_diederichs.run(L=[data, latt_id], nproc=self.nproc, verbose=True)

        self.best_operators = map(lambda x: None, xrange(len(arrays)))
        for op in result:
            idxes = map(int, result[op])
            print >>self.log_out, " %s num=%3d idxes= %s" %(op, len(result[op]), idxes)
            for idx in idxes:
                self.best_operators[idx] = sgtbx.change_of_basis_op(op)
    # assign_operators()

# class BrehmDiederichs

if __name__ == "__main__":
    import sys
    lst = sys.argv[1]
    xac_files = map(lambda x:x.strip(), open(lst))

    ksb = KabschSelectiveBreeding(xac_files, log_out=sys.stdout)

    
    if 1: # debug code
        from cctbx import sgtbx
        import random
        debug_op = sgtbx.change_of_basis_op("k,h,l")
        idxes = range(len(ksb.arrays))
        random.shuffle(idxes)
        for i in idxes[:len(ksb.arrays)//2]:
            ksb.arrays[i] = ksb.arrays[i].customized_copy(indices=debug_op.apply(ksb.arrays[i].indices()))

        print "altered:", idxes

    ksb.assign_operators([debug_op, sgtbx.change_of_basis_op("h,k,l")])
    print "right?:", [i for i, x in enumerate(ksb.best_operators) if not x.is_identity_op()]
    #ksb.debug_write_mtz()
    #ksb.modify_xds_ascii_files()

    quit()

    arrays = []
    for f in xac_files:
        print "Reading", f
        xac = XDS_ASCII(f, i_only=True)
        xac.remove_rejected()
        a = xac.i_obs().resolution_filter(d_min=3)
        a = a.merge_equivalents(use_internal_variance=False).array()
        arrays.append(a)

    symm = arrays[0].crystal_symmetry()
    cosets = reindex.reindexing_operators(symm, symm)
    reidx_ops = cosets.combined_cb_ops()
    reidx_ops.sort(key=lambda x: not x.is_identity_op())
    print " Possible reindex operators:", map(lambda x: str(x.as_hkl()), reidx_ops)

    determined = set([0,])
    old_ops = map(lambda x:0, xrange(len(arrays)))

    for ncycle in xrange(100):  # max cycle
        new_ops = map(lambda x:0, xrange(len(arrays)))
        for i in xrange(len(arrays)):
            cc_list = []
            a = arrays[i]
            for j, op in enumerate(reidx_ops):
                tmp = a.customized_copy(indices=op.apply(a.indices())).map_to_asu()
                for ref in determined:
                    if ref==i: continue
                    tmp2 = arrays[ref].customized_copy(indices=reidx_ops[new_ops[ref]].apply(arrays[ref].indices())).map_to_asu()
                    cc = calc_cc(tmp, tmp2)
                    #print "%d reindex= %10s cc=.%4f" % (i, op.as_hkl(), cc)
                    if cc==cc:
                        #cc_list.setdefault(op, []).append(cc)
                        cc_list.append((j,cc))
            if len(cc_list) == 0: continue
            max_el = max(cc_list, key=lambda x:x[1])
            print i, max_el, sum(map(lambda x:x[1], cc_list))/len(cc_list)
            new_ops[i] = max_el[0]
            #arrays[i] = a.customized_copy(indices=reidx_ops[max_el[0]].apply(a.indices())).map_to_asu()
            determined.add(i)

        print "In %4d cycle" % ncycle
        print "old",old_ops
        print "new",new_ops
        print "eq?", old_ops==new_ops
        print 
        if old_ops==new_ops: break
        old_ops = new_ops

    # Junk
    merge_org = arrays[0].deep_copy()
    for a in arrays[1:]: merge_org = merge_org.concatenate(a, assert_is_similar_symmetry=False)
    merge_org = merge_org.merge_equivalents(use_internal_variance=False).array()
    merge_org.as_mtz_dataset(column_root_label="I").mtz_object().write(file_name="noreindex.mtz")

    merge_new = None
    for a, opi in zip(arrays, new_ops):
        a = a.customized_copy(indices=reidx_ops[opi].apply(a.indices()))
        if merge_new is None: merge_new = a
        else: merge_new = merge_new.concatenate(a, assert_is_similar_symmetry=False)

    merge_new = merge_new.merge_equivalents(use_internal_variance=False).array()
    merge_new.as_mtz_dataset(column_root_label="I").mtz_object().write(file_name="reindexed.mtz")
