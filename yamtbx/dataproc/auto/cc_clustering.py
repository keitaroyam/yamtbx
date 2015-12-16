"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from cctbx.array_family import flex
from cctbx import miller
from libtbx import easy_mp
from libtbx.utils import null_out
from yamtbx.util import call
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
from yamtbx.dataproc.auto.blend import load_xds_data_only_indices
import os
import numpy
import collections

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

class CCClustering:
    def __init__(self, wdir, xac_files, d_min=None, d_max=None, min_ios=None):
        self.arrays = read_xac_files(xac_files, d_min=d_min, d_max=d_max, min_ios=min_ios)
        self.wdir = wdir
        self.clusters = {}
        
        if not os.path.exists(self.wdir): os.makedirs(self.wdir)

        open(os.path.join(self.wdir, "filenames.lst"), "w").write("\n".join(xac_files))
    # __init__()

    def do_clustering(self, nproc=1, b_scale=False, use_normalized=False):
        self.clusters = {}
        prefix = os.path.join(self.wdir, "cctable")
        assert (b_scale, use_normalized).count(True) <= 1

        if len(self.arrays) < 2:
            print "WARNING: less than two data! can't do cc-based clustering"
            self.clusters[1] = [float("nan"), [0]]
            return

        # Absolute scaling using Wilson-B factor 
        if b_scale:
            from mmtbx.scaling.matthews import p_vm_calculator
            from mmtbx.scaling.absolute_scaling import ml_iso_absolute_scaling
            
            ofs_wilson = open("%s_wilson_scales.dat"%prefix, "w")
            n_residues = p_vm_calculator(self.arrays.values()[0], 1, 0).best_guess
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
            for f in self.arrays:
                arr = self.arrays[f]
                normaliser = kernel_normalisation(arr, auto_kernel=True)
                self.arrays[f] = arr.customized_copy(data=arr.data()/normaliser.normalizer_for_miller_array,
                                                     sigmas=arr.sigmas()/normaliser.normalizer_for_miller_array)
        # Prep 
        args = []
        for i in xrange(len(self.arrays)-1):
            for j in xrange(i+1, len(self.arrays)):
                args.append((i,j))
           
        # Calc all CC
        worker = lambda x: calc_cc(self.arrays.values()[x[0]], self.arrays.values()[x[1]])
        results = easy_mp.pool_map(fixed_func=worker,
                                   args=args,
                                   processes=nproc)

        # Check NaN and decide which data to remove
        idx_bad = {}
        nans = []
        for (i,j), (cc,nref) in zip(args, results):
            if cc==cc: continue
            idx_bad[i] = idx_bad.get(i, 0) + 1
            idx_bad[j] = idx_bad.get(j, 0) + 1
            nans.append([i,j])

        idx_bad = idx_bad.items()
        idx_bad.sort(key=lambda x:x[1])
        remove_idxes = set()
        
        for idx, badcount in reversed(idx_bad):
            if len(filter(lambda x: idx in x, nans)) == 0: continue
            remove_idxes.add(idx)
            nans = filter(lambda x: idx not in x, nans)
            if len(nans) == 0: break

        use_idxes = filter(lambda x: x not in remove_idxes, xrange(len(self.arrays)))

        # Make table: original index (in file list) -> new index (in matrix)
        count = 0
        org2now = collections.OrderedDict()
        for i in xrange(len(self.arrays)):
            if i in remove_idxes: continue
            org2now[i] = count
            count += 1

        if len(remove_idxes) > 0:
            open("%s_notused.lst"%prefix, "w").write("\n".join(map(lambda x: self.arrays.keys()[x], remove_idxes)))

        # Make matrix
        mat = numpy.zeros(shape=(len(use_idxes), len(use_idxes)))
        for (i,j), (cc,nref) in zip(args, results):
            if i in remove_idxes or j in remove_idxes: continue
            mat[org2now[j], org2now[i]] = cc
            
        open("%s.matrix"%prefix, "w").write(" ".join(map(lambda x:"%.4f"%x, mat.flatten())))

        ofs = open("%s.dat"%prefix, "w")
        ofs.write("   i    j     cc  nref\n")
        for (i,j), (cc,nref) in zip(args, results):
            ofs.write("%4d %4d %.4f %4d\n" % (i,j,cc,nref))

        open("%s_ana.R"%prefix, "w").write("""\
treeToList2 <- function(htree)
{  # stolen from $CCP4/share/blend/R/blend0.R
 groups <- list()
 itree <- dim(htree$merge)[1]
 for (i in 1:itree)
 { 
  il <- htree$merge[i,1]
  ir <- htree$merge[i,2]
  if (il < 0) lab1 <- htree$labels[-il]
  if (ir < 0) lab2 <- htree$labels[-ir]
  if (il > 0) lab1 <- groups[[il]]
  if (ir > 0) lab2 <- groups[[ir]]
  lab <- c(lab1,lab2)
  lab <- as.integer(lab)
  groups <- c(groups,list(lab))
 }
 return(groups)
}

cc<-scan("%(prefix)s.matrix")
md<-matrix(1-cc, ncol=%(ncol)d, byrow=TRUE)
hc <- hclust(as.dist(md),method="ward")
pdf("tree.pdf")
plot(hc)
dev.off()
png("tree.png",height=1000,width=1000)
plot(hc)
dev.off()

hc$labels <- c(%(hclabels)s)
groups <- treeToList2(hc)
cat("ClNumber             Nds         Clheight   IDs\\n",file="./CLUSTERS.txt")
for (i in 1:length(groups))
{
 sorted_groups <- sort(groups[[i]])
 linea <- sprintf("%%04d %%4d %%7.3f %%s\\n",
                  i,length(groups[[i]]),hc$height[i], paste(sorted_groups,collapse=" "))
 cat(linea, file="./CLUSTERS.txt", append=TRUE)
}

# reference: http://www.coppelia.io/2014/07/converting-an-r-hclust-object-into-a-d3-js-dendrogram/
library(rjson)
HCtoJSON<-function(hc){
  labels<-hc$labels
  merge<-data.frame(hc$merge)
  for (i in (1:nrow(merge))) {
    if (merge[i,1]<0 & merge[i,2]<0) {eval(parse(text=paste0("node", i, "<-list(name=\\"", i, "\\", children=list(list(name=labels[-merge[i,1]]),list(name=labels[-merge[i,2]])))")))}
    else if (merge[i,1]>0 & merge[i,2]<0) {eval(parse(text=paste0("node", i, "<-list(name=\\"", i, "\\", children=list(node", merge[i,1], ", list(name=labels[-merge[i,2]])))")))}
    else if (merge[i,1]<0 & merge[i,2]>0) {eval(parse(text=paste0("node", i, "<-list(name=\\"", i, "\\", children=list(list(name=labels[-merge[i,1]]), node", merge[i,2],"))")))}
    else if (merge[i,1]>0 & merge[i,2]>0) {eval(parse(text=paste0("node", i, "<-list(name=\\"", i, "\\", children=list(node",merge[i,1] , ", node" , merge[i,2]," ))")))}
  }
  eval(parse(text=paste0("JSON<-toJSON(node",nrow(merge), ")")))
  return(JSON)
}

JSON<-HCtoJSON(hc)
cat(JSON, file="dendro.json")

q(save="yes")
""" % dict(prefix=os.path.basename(prefix),
           ncol=len(self.arrays),
           hclabels=",".join(map(lambda x: "%d"%(x+1), org2now.keys()))))

        call(cmd="Rscript", arg="%s_ana.R" % os.path.basename(prefix),
             wdir=self.wdir)

        output = open(os.path.join(self.wdir, "CLUSTERS.txt")).readlines()
        for l in output[1:]:
            sp = l.split()
            clid, clheight, ids = sp[0], sp[2], sp[3:]
            self.clusters[int(clid)] = [float(clheight), map(int,ids)]
    # do_clustering()

    def cluster_completeness(self, clno, anomalous_flag, d_min, calc_redundancy=True):
        if clno not in self.clusters:
            print "Cluster No. %d not found" % clno
            return

        cls = self.clusters[clno][-1]
        msets = map(lambda x: self.miller_sets[self.arrays.keys()[x-1]], cls)
        #msets = map(lambda x: self.arrays.values()[x-1], cls)
        num_idx = sum(map(lambda x: x.size(), msets))
        all_idx = flex.miller_index()
        all_idx.reserve(num_idx)
        for mset in msets: all_idx.extend(mset.indices())

        # Calc median cell
        cells = numpy.array(map(lambda x: x.unit_cell().parameters(), msets))
        median_cell = map(lambda i: numpy.median(cells[:,i]), xrange(6))
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

    def show_cluster_summary(self, d_min, out=null_out()):
        tmp = []
        self.miller_sets = load_xds_data_only_indices(xac_files=self.arrays.keys(), d_min=d_min)

        for clno in self.clusters:
            cluster_height, IDs = self.clusters[clno]
            cmpl, redun = self.cluster_completeness(clno, anomalous_flag=False, d_min=d_min)
            acmpl, aredun = self.cluster_completeness(clno, anomalous_flag=True, d_min=d_min)
            tmp.append((clno, IDs, cluster_height, cmpl*100., redun, acmpl*100., aredun))

        self.miller_sets = None # clear memory

        tmp.sort(key=lambda x: (-x[4], -x[3])) # redundancy & completeness
        out.write("# d_min= %.3f\n" % (d_min))
        out.write("# Sorted by redundancy & completeness\n")
        out.write("Cluster Number   CLh   Cmpl Redun  ACmpl ARedun\n")
        for clno, IDs, clh, cmpl, redun, acmpl, aredun in tmp:
            out.write("%7d %6d %5.1f %6.2f %5.1f %6.2f %5.1f\n" % (clno, len(IDs), clh, cmpl, redun,
                                                                   acmpl, aredun))
        return tmp
    # show_cluster_summary()

# class CCClustering

