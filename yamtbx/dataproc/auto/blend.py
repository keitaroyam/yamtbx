"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os
import numpy
import shutil
from yamtbx.dataproc.xds import xds_ascii
from yamtbx.dataproc.xds import integrate_hkl_as_flex
from yamtbx import util
from cctbx.array_family import flex
from cctbx import miller
from libtbx.utils import null_out

blend_comm = "blend"

def run_blend(wdir, xds_ascii_files, logout="blend_a.log"):
    ofs_lst = open(os.path.join(wdir, "files.lst"), "w")
    map(lambda x: ofs_lst.write(os.path.relpath(x, wdir)+"\n"), xds_ascii_files)
    ofs_lst.close()
    util.call(blend_comm, "-aDO files.lst", stdin="tolerance 10\n",
              wdir=wdir,
              stdout=open(os.path.join(wdir, logout), "w"))
# run_blend()

def run_blend0R(wdir, xds_ascii_files, logout="blend0.log"):
    ofs_cell = open(os.path.join(wdir, "forR_macropar.dat"), "w")
    ofs_lst = open(os.path.join(wdir, "xds_lookup_table.txt"), "w")
    ofs_files = open(os.path.join(wdir, "NEW_list_of_files.dat"), "w") # Just to avoid error
    for i, xac in enumerate(xds_ascii_files):
        symm = xds_ascii.XDS_ASCII(xac, read_data=False).symm
        cell = " ".join(map(lambda x: "%7.3f"%x, symm.unit_cell().parameters()))
        ofs_lst.write("dataset_%.3d %s\n" % (i, xac))
        ofs_cell.write("%4d %s 0 0 0\n" % (i+1, cell))
        ofs_files.write("%s\n"%xac)

    ofs_lst.close()
    ofs_cell.close()
    ofs_files.close()

    shutil.copyfile(os.path.join(wdir, "forR_macropar.dat"),
                    os.path.join(wdir, "forR_macropar.dat.bak"))

    util.call("Rscript", os.path.join(os.environ["CCP4"], "share/blend/R/blend0.R"),
              wdir=wdir,
              stdout=open(os.path.join(wdir, logout), "w"))

    open(os.path.join(wdir, "hctojson.R"), "w").write("""\
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

load("BLEND0.RData",.GlobalEnv)
JSON<-HCtoJSON(npar.hc_ward)
cat(JSON, file="dendro.json")
""")
    util.call("Rscript", "hctojson.R",
              wdir=wdir,
              stdout=open(os.path.join(wdir, logout), "a"))

# run_blend0R()

def load_xds_data_only_indices(xac_files, d_min=None):
    miller_sets = {}
    for f in xac_files:
        if xds_ascii.is_xds_ascii(f):
            print "Loading", f
            ms = xds_ascii.XDS_ASCII(f, i_only=True).as_miller_set()
            miller_sets[f] = ms.resolution_filter(d_min=d_min)
        elif integrate_hkl_as_flex.is_integrate_hkl(f):
            print "Sorry, skipping", f
        else:
            print "Skipping unrecognized:", f
    return miller_sets
# load_xds_data_only_indices()

class BlendClusters:
    def __init__(self, workdir=None, d_min=None, load_results=True):
        self.workdir = workdir
        self.d_min = d_min

        self.clusters = {} # clno -> (cluster_height, LCV, aLCV, IDs)
        self.files = None

        if load_results:
            self.read_results()
            self.miller_sets = load_xds_data_only_indices(xac_files=self.files, d_min=self.d_min)
    # __init__()

    def read_results(self):
        self.clusters = {}

        lookup_table_txt = os.path.join(self.workdir, "xds_lookup_table.txt")
        self.files = map(lambda l: l.split()[1], open(lookup_table_txt))
        
        clusters_txt = os.path.join(self.workdir, "CLUSTERS.txt")
        ifs = open(clusters_txt)
        for i in xrange(4): ifs.readline()
        for l in ifs:
            sp = l.split()
            clno, num, clh, lcv, alcv = sp[:5]
            ids = map(int, sp[5:])
            assert int(num) == len(ids)
            assert max(ids) <= len(self.files)
            self.clusters[int(clno)] = (float(clh), float(lcv), float(alcv), ids)
    # read_results()

    def cluster_completeness(self, clno, anomalous_flag, calc_redundancy=True):
        if clno not in self.clusters:
            print "Cluster No. %d not found" % clno
            return

        cls = self.clusters[clno][3]
        msets = map(lambda x: self.miller_sets[self.files[x-1]], cls)
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
        all_set = all_set.resolution_filter(d_min=self.d_min)

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

    def show_cluster_summary(self, out=null_out()):
        tmp = []
        for clno in self.clusters:
            cluster_height, LCV, aLCV, IDs = self.clusters[clno]
            cmpl, redun = self.cluster_completeness(clno, anomalous_flag=False)
            acmpl, aredun = self.cluster_completeness(clno, anomalous_flag=True)
            tmp.append((clno, IDs, cluster_height, cmpl*100., redun, acmpl*100., aredun, LCV, aLCV))

        tmp.sort(key=lambda x: (-x[4], -x[3])) # redundancy & completeness
        out.write("# d_min= %.3f\n" % (self.d_min))
        out.write("# Sorted by redundancy and completeness\n")
        out.write("Cluster Number   CLh   Cmpl Redun  ACmpl ARedun  LCV  aLCV\n")
        for clno, IDs, clh, cmpl, redun, acmpl, aredun, LCV, aLCV in tmp:
            out.write("%7d %6d %5.1f %6.2f %5.1f %6.2f %5.1f %5.1f %5.1f\n" % (clno, len(IDs), clh, cmpl, redun,
                                                                               acmpl, aredun, LCV, aLCV))
        return tmp
    # show_cluster_summary()

# class BlendClusters


if __name__ == "__main__":
    import sys

    bc = BlendClusters(sys.argv[1], d_min=3.)
    bc.show_cluster_summary(sys.stdout)
