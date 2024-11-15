"""
(c) RIKEN 2015-2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from cctbx.array_family import flex
from cctbx import miller
from libtbx import easy_mp
from libtbx.utils import null_out
from libtbx.utils import Sorry
from yamtbx.util import call
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.auto.blend import load_xds_data_only_indices
import os
import numpy
import collections
import scipy.cluster
import scipy.spatial
import json
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

def calc_cc(ari, arj):
  ari, arj = ari.common_sets(arj, assert_is_similar_symmetry=False)
  corr = flex.linear_correlation(ari.data(), arj.data())
  if corr.is_well_defined():
      return corr.coefficient(), ari.size()
  else:
      return float("nan"), ari.size()
# calc_cc()

def read_xac_files(xac_files, d_min=None, d_max=None, min_ios=None):
    arrays = collections.OrderedDict()

    for f in xac_files:
        xac = XDS_ASCII(f, i_only=True)
        xac.remove_rejected()
        a = xac.i_obs().resolution_filter(d_min=d_min, d_max=d_max)
        a = a.as_non_anomalous_array().merge_equivalents(use_internal_variance=False).array()
        if min_ios is not None: a = a.select(a.data()/a.sigmas()>=min_ios)
        arrays[f] = a

    return arrays
# read_xac_files()

class CCClustering(object):
    def __init__(self, wdir, xac_files, d_min=None, d_max=None, min_ios=None):
        self.arrays = read_xac_files(xac_files, d_min=d_min, d_max=d_max, min_ios=min_ios)
        self.wdir = wdir
        self.clusters = {}
        self.all_cc = {} # {(i,j):cc, ...}
        
        if not os.path.exists(self.wdir): os.makedirs(self.wdir)

        open(os.path.join(self.wdir, "filenames.lst"), "w").write("\n".join(xac_files))
    # __init__()

    def do_clustering(self, nproc=1, b_scale=False, use_normalized=False, cluster_method="ward", distance_eqn="sqrt(1-cc)", min_common_refs=3, html_maker=None):
        """
        Using correlation as distance metric (for hierarchical clustering)
        https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering

        Correlation "Distances" and Hierarchical Clustering
        http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm
        """

        self.clusters = {}
        prefix = os.path.join(self.wdir, "cctable")
        assert (b_scale, use_normalized).count(True) <= 1

        distance_eqns = {"sqrt(1-cc)": lambda x: numpy.sqrt(1.-x),
                         "1-cc": lambda x: 1.-x,
                         "sqrt(1-cc^2)": lambda x: numpy.sqrt(1.-x**2),
        }
        cc_to_distance = distance_eqns[distance_eqn] # Fail when unknown options
        assert cluster_method in ("single", "complete", "average", "weighted", "centroid", "median", "ward") # available methods in scipy
        
        
        if len(self.arrays) < 2:
            print("WARNING: less than two data! can't do cc-based clustering")
            self.clusters[1] = [float("nan"), [0]]
            return

        # Absolute scaling using Wilson-B factor 
        if b_scale:
            from mmtbx.scaling.matthews import p_vm_calculator
            from mmtbx.scaling.absolute_scaling import ml_iso_absolute_scaling
            
            ofs_wilson = open("%s_wilson_scales.dat"%prefix, "w")
            n_residues = p_vm_calculator(list(self.arrays.values())[0], 1, 0).best_guess
            ofs_wilson.write("# guessed n_residues= %d\n" % n_residues)
            ofs_wilson.write("file wilsonB\n")
            for f in self.arrays:
                arr = self.arrays[f]
                iso_scale_and_b = ml_iso_absolute_scaling(arr, n_residues, 0)
                wilson_b = iso_scale_and_b.b_wilson
                ofs_wilson.write("%s %.3f\n" % (f, wilson_b))
                if wilson_b > 0: # Ignoring data with B<0? is a bad idea.. but how..?
                    tmp = flex.exp(-2. * wilson_b * arr.unit_cell().d_star_sq(arr.indices())/4.)
                    self.arrays[f] = arr.customized_copy(data=arr.data()*tmp,
                                                         sigmas=arr.sigmas()*tmp)
            ofs_wilson.close()

        elif use_normalized:
            from mmtbx.scaling.absolute_scaling import kernel_normalisation
            failed = {}
            for f in self.arrays:
                arr = self.arrays[f]
                try:
                    normaliser = kernel_normalisation(arr, auto_kernel=True)
                    self.arrays[f] = arr.customized_copy(data=arr.data()/normaliser.normalizer_for_miller_array,
                                                         sigmas=arr.sigmas()/normaliser.normalizer_for_miller_array)
                except Exception as e:
                    failed.setdefault(str(e), []).append(f)

            if failed:
                msg = ""
                for r in failed: msg += " %s\n%s\n" % (r, "\n".join(["  %s"%x for x in failed[r]]))
                raise Sorry("intensity normalization failed by following reason(s):\n%s"%msg)
                    
        # Prep 
        args = []
        for i in range(len(self.arrays)-1):
            for j in range(i+1, len(self.arrays)):
                args.append((i,j))
           
        # Calc all CC
        worker = lambda x: calc_cc(list(self.arrays.values())[x[0]], list(self.arrays.values())[x[1]])
        results = easy_mp.pool_map(fixed_func=worker,
                                   args=args,
                                   processes=nproc)

        # Check NaN and decide which data to remove
        idx_bad = {}
        nans = []
        cc_data_for_html = []
        for (i,j), (cc,nref) in zip(args, results):
            cc_data_for_html.append((i,j,cc,nref))
            if cc==cc and nref>=min_common_refs: continue
            idx_bad[i] = idx_bad.get(i, 0) + 1
            idx_bad[j] = idx_bad.get(j, 0) + 1
            nans.append([i,j])

        if html_maker is not None:
            html_maker.add_cc_clustering_details(cc_data_for_html)

        idx_bad = list(idx_bad.items())
        idx_bad.sort(key=lambda x:x[1])
        remove_idxes = set()
        
        for idx, badcount in reversed(idx_bad):
            if len([x for x in nans if idx in x]) == 0: continue
            remove_idxes.add(idx)
            nans = [x for x in nans if idx not in x]
            if len(nans) == 0: break

        use_idxes = [x for x in range(len(self.arrays)) if x not in remove_idxes]
        N = len(use_idxes)

        # Make table: original index (in file list) -> new index (in matrix)
        count = 0
        org2now = collections.OrderedDict()
        for i in range(len(self.arrays)):
            if i in remove_idxes: continue
            org2now[i] = count
            count += 1

        if len(remove_idxes) > 0:
            open("%s_notused.lst"%prefix, "w").write("\n".join([list(self.arrays.keys())[x] for x in remove_idxes]))

        # Make matrix
        mat = numpy.zeros(shape=(N, N))
        self.all_cc = {}
        ofs = open("%s.dat"%prefix, "w")
        ofs.write("   i    j      cc  nref\n")
        for (i,j), (cc,nref) in zip(args, results):
            ofs.write("%4d %4d %13.11f %4d\n" % (i,j,cc,nref))
            self.all_cc[(i,j)] = cc
            if i not in remove_idxes and j not in remove_idxes:
                mat[org2now[j], org2now[i]] = cc_to_distance(min(cc, 1.)) #numpy.sqrt(1.-min(cc, 1.)) # safety guard (once I encounterd..
        ofs.close()

        # Perform cluster analysis
        D = scipy.spatial.distance.squareform(mat+mat.T) # convert to reduced form (first symmetrize)
        Z = scipy.cluster.hierarchy.linkage(D, cluster_method) # doesn't work with SciPy 0.17 (works with 0.18 or newer?)
        pyplot.figure(figsize=(max(5, min(N,200)/10.),)*2, dpi=100)
        pyplot.title("%s, %s" % (distance_eqn, cluster_method))
        hclabels = [x+1 for x in list(org2now.keys())]
        scipy.cluster.hierarchy.dendrogram(Z, orientation="left", labels=hclabels)
        pyplot.savefig(os.path.join(self.wdir, "tree.png"))
        pyplot.savefig(os.path.join(self.wdir, "tree.pdf"))

        def traverse(node, results, dendro):
            # Cluster id starts with the number of data files
            leaves = [hclabels[x] for x in sorted(node.pre_order())] # file numbers
            
            if not node.right and not node.left: # this must be leaf
                dendro["children"].append(dict(name=str(hclabels[node.id]))) # file number
            else:
                dendro["children"].append(dict(name=str(node.id-N+1), children=[])) # cluster number
                results.append((node.id-N+1, node.dist, node.count, leaves))

            if node.right: traverse(node.right, results, dendro["children"][-1])
            if node.left: traverse(node.left, results, dendro["children"][-1])
        # traverse()

        results = []
        dendro = dict(name="root", children=[])
        node = scipy.cluster.hierarchy.to_tree(Z)
        traverse(node, results, dendro)
        dendro = dendro["children"][0]
        results.sort(key=lambda x: x[0])
        json.dump(dendro, open(os.path.join(self.wdir, "dendro.json"), "w"))

        # Save CLUSTERS.txt and set self.clusters
        ofs = open(os.path.join(self.wdir, "CLUSTERS.txt"), "w")
        ofs.write("ClNumber             Nds         Clheight   IDs\n")
        for clid, dist, ncount, leaves in results:
            leavestr = " ".join(map(str, leaves))
            ofs.write("%04d %4d %7.3f %s\n" % (clid, ncount, dist, leavestr))
            self.clusters[int(clid)] = [float(dist), leaves]
        ofs.close()
            
    # do_clustering()

    def cluster_completeness(self, clno, anomalous_flag, d_min, calc_redundancy=True):
        if clno not in self.clusters:
            print("Cluster No. %d not found" % clno)
            return

        cls = self.clusters[clno][-1]
        msets = [self.miller_sets[list(self.arrays.keys())[x-1]] for x in cls]
        #msets = map(lambda x: self.arrays.values()[x-1], cls)
        num_idx = sum([x.size() for x in msets])
        all_idx = flex.miller_index()
        all_idx.reserve(num_idx)
        for mset in msets: all_idx.extend(mset.indices())

        # Calc median cell
        cells = numpy.array([x.unit_cell().parameters() for x in msets])
        median_cell = [numpy.median(cells[:,i]) for i in range(6)]
        symm = msets[0].customized_copy(unit_cell=median_cell)

        assert anomalous_flag is not None
        # XXX all must belong to the same Laue group and appropriately reindexed..
        all_set = miller.set(indices=all_idx,
                             crystal_symmetry=symm, anomalous_flag=anomalous_flag)
        all_set = all_set.resolution_filter(d_min=d_min)

        # dummy for redundancy calculation. dirty way..
        if calc_redundancy:
            dummy_array = miller.array(miller_set=all_set, data=flex.int(all_set.size()))
            merge = dummy_array.merge_equivalents()
            cmpl = merge.array().completeness()
            redun = merge.redundancies().as_double().mean()
            return cmpl, redun
        else:
            cmpl = all_set.unique_under_symmetry().completeness()
            return cmpl
    # cluster_completeness()

    def get_all_cc_in_cluster(self, clno):
      IDs = self.clusters[clno][1]
      ret = []
      
      for i in range(len(IDs)-1):
          for j in range(i+1, len(IDs)):
            ids = IDs[i]-1, IDs[j]-1
            ret.append(self.all_cc[(min(ids), max(ids))])
      return ret
    # get_all_cc_in_cluster()
    
    def show_cluster_summary(self, d_min, out=null_out()):
        tmp = []
        self.miller_sets = load_xds_data_only_indices(xac_files=list(self.arrays.keys()), d_min=d_min)

        for clno in self.clusters:
            cluster_height, IDs = self.clusters[clno]
            cmpl, redun = self.cluster_completeness(clno, anomalous_flag=False, d_min=d_min)
            acmpl, aredun = self.cluster_completeness(clno, anomalous_flag=True, d_min=d_min)
            all_cc = self.get_all_cc_in_cluster(clno)
            ccmean, ccmin = numpy.mean(all_cc), min(all_cc)
            tmp.append((clno, IDs, cluster_height, cmpl*100., redun, acmpl*100., aredun, ccmean, ccmin))

        self.miller_sets = None # clear memory

        tmp.sort(key=lambda x: (-x[4], -x[3])) # redundancy & completeness
        out.write("# d_min= %.3f\n" % (d_min))
        out.write("# Sorted by redundancy & completeness\n")
        out.write("Cluster Number   CLh   Cmpl Redun  ACmpl ARedun CCmean CCmin\n")
        for clno, IDs, clh, cmpl, redun, acmpl, aredun, ccmean, ccmin in tmp:
            out.write("%7d %6d %5.1f %6.2f %5.1f %6.2f %5.1f %.4f %.4f\n" % (clno, len(IDs), clh, cmpl, redun,
                                                                             acmpl, aredun, ccmean, ccmin))
        return tmp
    # show_cluster_summary()

# class CCClustering

