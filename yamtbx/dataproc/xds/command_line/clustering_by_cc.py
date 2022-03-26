"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import iotbx.phil
from cctbx.array_family import flex
from libtbx import easy_mp
from yamtbx.util import read_path_list
from yamtbx.dataproc.xds.xds_ascii import XDS_ASCII
import os
import numpy

master_params_str = """
d_min = 3
 .type = float
d_max = None
 .type = float
min_ios = 3
 .type = float
 .help = minimum I/sigma for input
nproc = 1
 .type = int
prefix = cctable
 .type = str
"""

def calc_cc(ari, arj):
  ari, arj = ari.common_sets(arj, assert_is_similar_symmetry=False)
  corr = flex.linear_correlation(ari.data(), arj.data())
  if corr.is_well_defined():
      return corr.coefficient(), ari.size()
  else:
      return float("nan"), ari.size()
# calc_cc()

def run(lstin, params):
    xac_files = read_path_list(lstin)

    common0 = len(os.path.commonprefix(xac_files))

    arrays = []

    for f in xac_files:
        xac = XDS_ASCII(f, i_only=True)
        xac.remove_rejected()
        a = xac.i_obs().resolution_filter(d_min=params.d_min, d_max=params.d_max)
        a = a.merge_equivalents(use_internal_variance=False).array()
        a = a.select(a.data()/a.sigmas()>=params.min_ios)
        arrays.append(a)

    # Prep 
    args = []
    for i in range(len(arrays)-1):
        for j in range(i+1, len(arrays)):
            args.append((i,j))
           
    # Calc all CC
    worker = lambda x: calc_cc(arrays[x[0]], arrays[x[1]])
    results = easy_mp.pool_map(fixed_func=worker,
                               args=args,
                               processes=params.nproc)

    # Make matrix
    mat = numpy.zeros(shape=(len(arrays), len(arrays)))
    for (i,j), (cc,nref) in zip(args, results):
        print(j,i,cc)
        mat[j,i] = cc
    
    open("%s.names"%params.prefix, "w").write("\n".join([os.path.dirname(x[common0:]) for x in xac_files]))
    open("%s.matrix"%params.prefix, "w").write(" ".join(["%.4f"%x for x in mat.flatten()]))

    ofs = open("%s.dat"%params.prefix, "w")
    ofs.write("i j cc nref\n")
    for (i,j), (cc,nref) in zip(args, results):
        ofs.write("%4d %4d %.4f %4d\n" % (i,j,cc,nref))

    open("%s_ana.R"%params.prefix, "w").write("""\
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

cc<-scan("%s.matrix")
md<-matrix(1-cc, ncol=%d, byrow=TRUE)
labs<-read.table("%s.names")
filenames<-read.table("%s")$V1
rownames(md)<-labs$V1
hc <- hclust(as.dist(md),method="ward")
pdf("tree.pdf")
plot(hc)
dev.off()

hc$labels <- 1:nrow(md)
groups <- treeToList2(hc)
cat("ClNumber             Nds         Clheight\\n",file="./CLUSTERS.txt")
for (i in 1:length(groups))
{
 sorted_groups <- sort(groups[[i]])
 linea <- paste(sprintf("     %%03d           %%3d         %%7.3f\\n",
                i,length(groups[[i]]),hc$height[i]),sep="")
 cat(linea, file="./CLUSTERS.txt", append=TRUE)
 write.table(filenames[sorted_groups], sprintf("cluster%%.3d.lst",i), quote=FALSE, row.names=FALSE, col.names=FALSE)
}

q(save="yes")
""" % (params.prefix, len(arrays),  params.prefix, lstin))
    print("R --vanilla < %s_ana.R" % params.prefix)
    

# run()
    

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    lstin = args[0]
    #run(lstin, params)
    from yamtbx.dataproc.auto.cc_clustering import CCClustering
    ccc = CCClustering(wdir=".",
                       xac_files=read_path_list(lstin),
                       d_min=params.d_min,
                       d_max=params.d_max,
                       min_ios=params.min_ios)

    ccc.do_clustering(nproc=params.nproc)
