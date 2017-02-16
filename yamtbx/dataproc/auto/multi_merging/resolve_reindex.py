"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII

from cctbx.crystal import reindex
from cctbx.array_family import flex
from cctbx import sgtbx
from libtbx.utils import null_out
from libtbx import easy_mp
from cctbx.merging import brehm_diederichs

import os
import copy

def calc_cc(a1, a2):
    a1, a2 = a1.common_sets(a2, assert_is_similar_symmetry=False)
    corr = flex.linear_correlation(a1.data(), a2.data())

    if corr.is_well_defined():# and a1.size() > 20:
        return corr.coefficient()
    else:
        return float("nan")
# calc_cc()

class ReindexResolver:
    def __init__(self, xac_files, d_min=3, min_ios=3, nproc=1, max_delta=3, log_out=null_out()):
        self.xac_files = xac_files
        self.log_out = log_out
        self.nproc = nproc
        self.arrays = []
        self.max_delta = max_delta
        self.best_operators = None

        print >>self.log_out, "Reading"
        for i, f in enumerate(self.xac_files):
            print >>self.log_out, "%4d %s" % (i, f)
            xac = XDS_ASCII(f, i_only=True)
            xac.remove_rejected()
            a = xac.i_obs().resolution_filter(d_min=d_min)
            if min_ios is not None: a = a.select(a.data()/a.sigmas()>=min_ios)
            a = a.as_non_anomalous_array().merge_equivalents(use_internal_variance=False).array()
            self.arrays.append(a)

        print >>self.log_out, ""
    # __init__()

    def find_reindex_ops(self):
        symm = self.arrays[0].crystal_symmetry()
        cosets = reindex.reindexing_operators(symm, symm, max_delta=self.max_delta)
        reidx_ops = cosets.combined_cb_ops()
        return reidx_ops
    # find_reindex_ops()

    def modify_xds_ascii_files(self, suffix="_reidx", cells_dat_out=None):
        #ofs_lst = open("for_merge_new.lst", "w")
        if cells_dat_out: cells_dat_out.write("file a b c al be ga\n")

        new_files = []
        print >>self.log_out, "Writing reindexed files.."
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

class KabschSelectiveBreeding(ReindexResolver):
    """
    Reference: W. Kabsch "Processing of X-ray snapshots from crystals in random orientations" Acta Cryst. (2014). D70, 2204-2216
    http://dx.doi.org/10.1107/S1399004714013534
    If I understand correctly...
    """

    def __init__(self, xac_files, d_min=3, min_ios=3, nproc=1, max_delta=3, log_out=null_out()):
        ReindexResolver.__init__(self, xac_files, d_min, min_ios, nproc, max_delta, log_out)
    # __init__()
    
    def assign_operators(self, reidx_ops=None, max_cycle=100):
        arrays = self.arrays
        self.best_operators = None

        if reidx_ops is None: reidx_ops = self.find_reindex_ops()

        print >>self.log_out, "Reindex operators:"
        for i, op in enumerate(reidx_ops): print >>self.log_out, " %2d: %s" % (i, op.as_hkl())
        print >>self.log_out, ""

        reidx_ops.sort(key=lambda x: not x.is_identity_op()) # identity op to first

        if self.nproc > 1:
            # consumes much memory..
            reindexed_arrays = [arrays]
            for op in reidx_ops[1:]:
                reindexed_arrays.append(map(lambda x: x.customized_copy(indices=op.apply(x.indices())).map_to_asu(), arrays))
        else:
            reindexed_arrays = None

        old_ops = map(lambda x:0, xrange(len(arrays)))
        new_ops = map(lambda x:0, xrange(len(arrays)))

        for ncycle in xrange(max_cycle):
            #new_ops = copy.copy(old_ops) # doesn't matter
            for i in xrange(len(arrays)):
                cc_means = []
                a = arrays[i]

                for j, op in enumerate(reidx_ops):
                    cc_list = []
                    if self.nproc > 1:
                        tmp = reindexed_arrays[j][i]
                    else:
                        if op.is_identity_op(): tmp = a
                        else: tmp = a.customized_copy(indices=op.apply(a.indices())).map_to_asu()
                    
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
                                               processes=1)#self.nproc)
                    cc_list = filter(lambda x: x is not None, cc_list)

                    if len(cc_list) > 0:
                        cc_means.append((j, sum(cc_list)/len(cc_list)))
                        #print  >>self.log_out, "DEBUG:", i, j, cc_list, cc_means[-1]

                max_el = max(cc_means, key=lambda x:x[1])
                print >>self.log_out, "%3d %s" % (i, " ".join(map(lambda x: "%s%d:% .4f" % ("*" if x[0]==max_el[0] else " ", x[0], x[1]), cc_means)))
                new_ops[i] = max_el[0]

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
    def __init__(self, xac_files, ref_array, d_min=3, min_ios=3,  nproc=1, max_delta=3, log_out=null_out()):
        ReindexResolver.__init__(self, xac_files, d_min, min_ios, nproc, max_delta, log_out)
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

            cc_list.sort(key=lambda x:-x[1])
            max_el = cc_list[0]
            print >>self.log_out, "%4d"%i, " ".join(map(lambda x:"% .4f"%x[1],cc_list))
            new_ops[i] = max_el[0]

        print >>self.log_out, "  operator:", new_ops
        print >>self.log_out, "  number of different assignments:", len(filter(lambda x:x!=0, new_ops))
        print >>self.log_out, ""

        self.best_operators = map(lambda x: reidx_ops[x], new_ops)
    # assign_operators()
# class ReferenceBased

class BrehmDiederichs(ReindexResolver):
    def __init__(self, xac_files, d_min=3, min_ios=3, nproc=1, max_delta=3, log_out=null_out()):
        ReindexResolver.__init__(self, xac_files, d_min, min_ios, nproc, max_delta, log_out)
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
