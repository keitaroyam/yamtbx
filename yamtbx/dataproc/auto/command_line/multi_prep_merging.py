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
from libtbx import easy_mp
from libtbx.utils import multi_out
from cctbx.crystal import reindex
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex

import os
import time
import sys
import shutil
import numpy
import traceback
import tempfile
import glob
from cStringIO import StringIO

master_params_str = """
xdsdir = None
 .type = path
 .multiple = true
 .help = top directory containing xds results
workdir = None
 .type = path
 .help = ""

unit_cell = None
 .type = floats(size=6)
space_group = None
 .type = str
 .help = Choose the space group (name or number).
group_choice = None
 .type = int
 .help = Choose the group (run once and choose).

cell_method = *reindex refine
 .type = choice(multi=False)
 .help = ""
nproc = 1
 .type = int
prep_dials_files = True
 .type = bool
copy_into_workdir = True
 .type = bool

cell_grouping {
  tol_length = None
   .type = float
   .help = relative_length_tolerance
  tol_angle = None
   .type = float
   .help = absolute_angle_tolerance in degree
}
"""

def prepare_dials_files(wd, out, space_group=None, reindex_op=None, moveto=None):
    try:
        from yamtbx.dataproc.dials.command_line import import_xds_for_refine
        files = import_xds_for_refine.run(xds_inp=os.path.join(wd, "XDS.INP"),
                                          xparm=os.path.join(wd, "XPARM.XDS"),
                                          integrate_lp=os.path.join(wd, "INTEGRATE.LP"),
                                          integrate_hkl=os.path.join(wd, "INTEGRATE.HKL"),
                                          spot_xds=os.path.join(wd, "SPOT.XDS"),
                                          space_group=space_group, reindex_op=reindex_op,
                                          out_dir=wd)
        if moveto and wd!=moveto:
            for f in files: shutil.move(f, moveto)
    except:
        print >>out, "Error in generation of dials files in %s" % wd
        print >>out, traceback.format_exc()
# prepare_dials_files()

def rescale_with_specified_symm_worker(sym_wd_wdr, topdir, log_out, reference_symm, sgnum, sgnum_laue, prep_dials_files=False):
    # XXX Unsafe if multiple processes run this function for the same target directory at the same time

    sym, wd, wdr = sym_wd_wdr
    out = StringIO()
    print >>out,  os.path.relpath(wd, topdir),

    # Find appropriate data # XXX not works for DIALS data!!
    xac_file = util.return_first_found_file(("XDS_ASCII.HKL_noscale.org", "XDS_ASCII.HKL_noscale", 
                                             "XDS_ASCII_fullres.HKL.org", "XDS_ASCII_fullres.HKL",
                                             "XDS_ASCII.HKL.org", "XDS_ASCII.HKL"),
                                            wd=wd)
    if xac_file is None:
        print >>out, "Can't find XDS_ASCII file in %s" % wd
        log_out.write(out.getvalue())
        log_out.flush()
        return (wd, None)

    xac = XDS_ASCII(xac_file, read_data=False)
    print >>out, "%s %s (%s)" % (os.path.basename(xac_file), xac.symm.space_group_info(),
                                 ",".join(map(lambda x: "%.2f"%x, xac.symm.unit_cell().parameters())))

    if xac.symm.reflection_intensity_symmetry(False).space_group_info().type().number() == sgnum_laue:
        if xac.symm.unit_cell().is_similar_to(reference_symm.unit_cell(), 0.1, 10):
            print >>out,  "  Already scaled with specified symmetry"
            log_out.write(out.getvalue())
            log_out.flush()

            if wd != wdr: shutil.copy2(xac_file, wdr)

            if prep_dials_files: prepare_dials_files(wd, out, moveto=wdr)
            return (wdr, (numpy.array(xac.symm.unit_cell().parameters()),
                          os.path.join(wdr, os.path.basename(xac_file))))

    xdsinp = os.path.join(wd, "XDS.INP")
    cosets = reindex.reindexing_operators(reference_symm, xac.symm, 0.2, 20)

    if len(cosets.combined_cb_ops())==0:
        print >>out, "Can't find operator:"
        sym.show_summary(out, " ")
        reference_symm.show_summary(out, " ")
        log_out.write(out.getvalue())
        log_out.flush()
        return (wdr, None)

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
                                      ("MINIMUM_I/SIGMA", None), # use default
                                      ("REFINE(CORRECT)", None), # use default
                                      ])
    run_xds(wd)
    for f in ("XDS.INP", "CORRECT.LP", "XDS_ASCII.HKL", "GXPARM.XDS"):
        if os.path.exists(os.path.join(wd, f)):
            shutil.copyfile(os.path.join(wd, f), os.path.join(wdr, f+"_rescale"))

    revert_files(xds_files.generated_by_CORRECT, bk_prefix, wdir=wd, quiet=True)

    new_xac = os.path.join(wdr, "XDS_ASCII.HKL_rescale")

    if prep_dials_files:
        prepare_dials_files(wd, out,
                            space_group=reference_symm.space_group(),
                            reindex_op=cosets.combined_cb_ops()[0],
                            moveto=wdr)

    ret = None
    if os.path.isfile(new_xac):
        ret = (XDS_ASCII(new_xac, read_data=False).symm.unit_cell().parameters(), new_xac)
        print >>out, " OK:", ret[0]
    else:
        print >>out, "Error: rescaling failed (Can't find XDS_ASCII.HKL)"

    return (wd, ret)
# rescale_with_specified_symm_worker()

def rescale_with_specified_symm(topdir, dirs, symms, out, sgnum=None, reference_symm=None, nproc=1, prep_dials_files=False, copyto_root=None):
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
    out.flush()
    st_time = time.time()
    wd_ret = []

    if copyto_root:
        for wd in dirs:
            assert wd.startswith(os.path.join(topdir, ""))
            tmp = os.path.join(copyto_root, os.path.relpath(wd, topdir))
            if not os.path.exists(tmp): os.makedirs(tmp)
            wd_ret.append(tmp)
    else:
        wd_ret = dirs


    ret = easy_mp.pool_map(fixed_func=lambda x: rescale_with_specified_symm_worker(x, topdir, out, reference_symm, sgnum, sgnum_laue, prep_dials_files),
                           args=zip(symms, dirs, wd_ret), processes=nproc)
    cells = dict(filter(lambda x: x[1] is not None, ret)) # cell and file
    print >>out, "\nTotal wall-clock time for reindexing: %.2f sec (using %d cores)." % (time.time()-st_time, nproc)
    return cells, reference_symm
# rescale_with_specified_symm()

def reindex_with_specified_symm_worker(wd, wdr, topdir, log_out, reference_symm, sgnum_laue, prep_dials_files=False):
    """
    wd: directory where XDS file exists
    wdr: wd to return; a directory where transformed file should be saved.

    If wd!=wdr, files in wd/ are unchanged during procedure. Multiprocessing is unsafe when wd==wdr.
    """

    out = StringIO()
    print >>out, "%s:" % os.path.relpath(wd, topdir),

    # Find appropriate data
    xac_file = util.return_first_found_file(("XDS_ASCII.HKL_noscale.org", "XDS_ASCII.HKL_noscale", 
                                             "XDS_ASCII_fullres.HKL.org", "XDS_ASCII_fullres.HKL",
                                             "XDS_ASCII.HKL.org", "XDS_ASCII.HKL", "DIALS.HKL.org", "DIALS.HKL"),
                                            wd=wd)
    if xac_file is None:
        print >>out, "Can't find XDS_ASCII file in %s" % wd
        log_out.write(out.getvalue())
        log_out.flush()
        return (wdr, None)

    if xac_file.endswith(".org"): xac_file_out = xac_file[:-4]
    else: xac_file_out = xac_file

    xac = XDS_ASCII(xac_file, read_data=False)
    print >>out, "%s %s (%s)" % (os.path.basename(xac_file), xac.symm.space_group_info(),
                               ",".join(map(lambda x: "%.2f"%x, xac.symm.unit_cell().parameters())))

    if xac.symm.reflection_intensity_symmetry(False).space_group_info().type().number() == sgnum_laue:
        if xac.symm.unit_cell().is_similar_to(reference_symm.unit_cell(), 0.1, 10): # XXX Check unit cell consistency!!
            print >>out,  "  Already scaled with specified symmetry"
            log_out.write(out.getvalue())
            log_out.flush()

            if wd != wdr: shutil.copy2(xac_file, wdr)

            if prep_dials_files and "DIALS.HKL" not in xac_file:
                prepare_dials_files(wd, out, moveto=wdr)

            return (wdr, (numpy.array(xac.symm.unit_cell().parameters()), 
                          os.path.join(wdr, os.path.basename(xac_file))))
            

    cosets = reindex.reindexing_operators(reference_symm, xac.symm, 0.2, 20) # XXX ISN'T THIS TOO LARGE?

    if len(cosets.combined_cb_ops())==0:
        print >>out, "Can't find operator:"
        xac.symm.show_summary(out, " ")
        reference_symm.show_summary(out, " ")
        log_out.write(out.getvalue())
        log_out.flush()
        return (wdr, None)

    if wd == wdr:
        dest = tempfile.mkdtemp(prefix="multiprep", dir=wd)
    else:
        dest = wdr
        
    hklout = os.path.join(dest, os.path.basename(xac_file_out))

    newcell = xac.write_reindexed(op=cosets.combined_cb_ops()[0],
                                  space_group=reference_symm.space_group(),
                                  hklout=hklout)

    if "DIALS.HKL" in os.path.basename(xac_file):
        outstr = 'output.experiments="%sreindexed_experiments.json" ' % os.path.join(dest, "")
        outstr += 'output.reflections="%sreindexed_reflections.pickle" ' % os.path.join(dest, "")
        for f in ("experiments.json", "indexed.pickle"):
            if not os.path.isfile(os.path.join(os.path.dirname(xac_file), f)): continue
            util.call('dials.reindex %s change_of_basis_op=%s space_group="%s" %s'%(f, 
                                                                                    cosets.combined_cb_ops()[0].as_abc(), 
                                                                                    reference_symm.space_group_info(),
                                                                                    outstr),
                      wdir=os.path.dirname(xac_file))
    elif prep_dials_files:
        prepare_dials_files(wd, out,
                            space_group=reference_symm.space_group(),
                            reindex_op=cosets.combined_cb_ops()[0],
                            moveto=dest)

    newcell_str = " ".join(map(lambda x: "%.3f"%x, newcell.parameters()))
    print >>out,  "  Reindexed to transformed cell: %s with %s" % (newcell_str, cosets.combined_cb_ops()[0].as_hkl())
    log_out.write(out.getvalue())
    log_out.flush()

    if wd == wdr:
        for f in glob.glob(os.path.join(dest, "*")):
            f_in_wd = os.path.join(wd, os.path.basename(f))
            if os.path.exists(f_in_wd) and not os.path.exists(f_in_wd+".org"): os.rename(f_in_wd, f_in_wd+".org")
            os.rename(f, f_in_wd)

        shutil.rmtree(dest)
        ret = (numpy.array(newcell.parameters()), 
               os.path.join(wd, os.path.basename(xac_file_out)))
    else:
        ret = (numpy.array(newcell.parameters()), hklout)
           

    return (wdr, ret)
# reindex_with_specified_symm_worker()

def reindex_with_specified_symm(topdir, reference_symm, dirs, out, nproc=10, prep_dials_files=False, copyto_root=None):
    print >>out
    print >>out,  "Re-index to specified symmetry:"
    reference_symm.show_summary(out, "  ")
    print >>out
    print >>out
    out.flush()

    st_time = time.time()
    wd_ret = []

    if copyto_root:
        for wd in dirs:
            assert wd.startswith(os.path.join(topdir, ""))
            tmp = os.path.join(copyto_root, os.path.relpath(wd, topdir))
            if not os.path.exists(tmp): os.makedirs(tmp)
            wd_ret.append(tmp)
    else:
        wd_ret = dirs

    sgnum_laue = reference_symm.space_group().build_derived_reflection_intensity_group(False).type().number()

    ret = easy_mp.pool_map(fixed_func=lambda wd2: reindex_with_specified_symm_worker(wd2[0], wd2[1], topdir, out, reference_symm, sgnum_laue, prep_dials_files),
                           args=zip(dirs, wd_ret), processes=nproc)
    cells = dict(filter(lambda x: x[1] is not None, ret)) # cell and file

    print >>out, "\nTotal wall-clock time for reindexing: %.2f sec (using %d cores)." % (time.time()-st_time, nproc)

    return cells
# reindex_with_specified_symm()

class PrepMerging:
    def __init__(self, cell_graph):
        self.cell_graph = cell_graph
        self.cell_and_files = {}
        self.log_buffer = None
    # __init__()

    def find_groups(self):
        sio = StringIO()
        self.cell_graph.group_xds_results(sio)
        self.log_buffer = sio.getvalue()
        return self.log_buffer
    # find_groups()
    
    def prep_merging(self, workdir, group, symmidx=None, reference_symm=None, topdir=None, cell_method="reindex", nproc=1, prep_dials_files=True, into_workdir=True):
        assert (symmidx, reference_symm).count(None) == 1
        
        from yamtbx.util.xtal import format_unit_cell
        from cctbx.crystal import reindex

        cm = self.cell_graph
        
        prep_log_out = multi_out()
        prep_log_out.register("log", open(os.path.join(workdir, "prep_merge.log"), "w"), atexit_send_to=None)
        prep_log_out.register("stdout", sys.stdout)
        prep_log_out.write(self.log_buffer)
        prep_log_out.flush()

        if reference_symm is None: reference_symm = cm.get_reference_symm(group-1, symmidx)

        prep_log_out.write("\n\ngroup_choice= %d, symmetry= %s (%s)\n" % (group, reference_symm.space_group_info(),
                                                                          format_unit_cell(reference_symm.unit_cell())))
        prep_log_out.flush()

        # Scale with specified symmetry
        symms = map(lambda i: cm.symms[i], cm.groups[group-1])
        dirs = map(lambda i: cm.dirs[i], cm.groups[group-1])
        copyto_root = os.path.join(workdir, "input_files") if into_workdir else None

        if not topdir: topdir = os.path.dirname(os.path.commonprefix(dirs))
        
        if cell_method == "reindex":
            self.cell_and_files = reindex_with_specified_symm(topdir, reference_symm, dirs,
                                                         out=prep_log_out, nproc=nproc,
                                                         prep_dials_files=prep_dials_files, copyto_root=copyto_root)
        elif cell_method == "refine":
            self.cell_and_files, reference_symm = rescale_with_specified_symm(topdir, dirs, symms,
                                                                         reference_symm=reference_symm,
                                                                         out=prep_log_out, nproc=nproc,
                                                                         prep_dials_files=prep_dials_files,
                                                                         copyto_root=copyto_root)
        else:
            raise "Don't know this choice: %s" % cell_method

        prep_log_out.flush()

        cosets = reindex.reindexing_operators(reference_symm, reference_symm, max_delta=5)
        reidx_ops = cosets.combined_cb_ops()

        print >>prep_log_out, "\nReference symmetry:", reference_symm.unit_cell(), reference_symm.space_group_info().symbol_and_number()
        msg_reindex = "\n"
        if len(reidx_ops) > 1:
            msg_reindex += "!! ATTENTION !! Reindex operators found. You may need to reindex some files before merging.\n"
            for rop in reidx_ops:
                msg_reindex += " operator: %-16s Cell: (%s)\n" % (rop.as_hkl(),
                                                                format_unit_cell(reference_symm.unit_cell().change_basis(rop)))
            msg_reindex += "Try kamo.resolve_indexing_ambiguity command before merging!!"

        else:
            msg_reindex += "No reindex operators found. No need to run kamo.resolve_indexing_ambiguity."

        print >>prep_log_out, "%s\n\n" % msg_reindex
        prep_log_out.close()

        # Make list for merging
        ofs_lst = open(os.path.join(workdir, "formerge.lst"), "w")
        ofs_dat = open(os.path.join(workdir, "cells.dat"), "w")
        ofs_dat.write("file a b c al be ga\n")

        for wd in sorted(self.cell_and_files):
            cell, xas = self.cell_and_files[wd]
            ofs_lst.write(xas+"\n")
            ofs_dat.write(xas+" "+" ".join(map(lambda x:"%7.3f"%x, cell))+"\n")

        ofs_dat.close()
        ofs_lst.close()

        return msg_reindex, reidx_ops
    # prep_merging()

    def write_merging_scripts(self, workdir, sge_pe_name="par", prep_dials_files=True):
        open(os.path.join(workdir, "merge_blend.sh"), "w").write("""\
#!/bin/sh
# settings
dmin=2.8 # resolution
anomalous=false # true or false
lstin=formerge.lst # list of XDS_ASCII.HKL files
use_ramdisk=true # set false if there is few memory or few space in /tmp
# _______/setting

kamo.multi_merge \\
        workdir=blend_${dmin}A_framecc_b+B \\
        lstin=${lstin} d_min=${dmin} anomalous=${anomalous} \\
        space_group=None reference.data=None \\
        program=xscale xscale.reference=bmin xscale.degrees_per_batch=None \\
        reject_method=framecc+lpstats rejection.lpstats.stats=em.b+bfactor \\
        clustering=blend blend.min_cmpl=90 blend.min_redun=2 blend.max_LCV=None blend.max_aLCV=None \\
        max_clusters=None xscale.use_tmpdir_if_available=${use_ramdisk} \\
#        batch.engine=sge batch.par_run=merging batch.nproc_each=8 nproc=8 batch.sge_pe_name=%s
""" % sge_pe_name)
        os.chmod(os.path.join(workdir, "merge_blend.sh"), 0755)
        
        open(os.path.join(workdir, "merge_ccc.sh"), "w").write("""\
#!/bin/sh
# settings
dmin=2.8 # resolution
clustering_dmin=3.5  # resolution for CC calculation
anomalous=false # true or false
lstin=formerge.lst # list of XDS_ASCII.HKL files
use_ramdisk=true # set false if there is few memory or few space in /tmp
# _______/setting

kamo.multi_merge \\
        workdir=ccc_${dmin}A_framecc_b+B \\
        lstin=${lstin} d_min=${dmin} anomalous=${anomalous} \\
        space_group=None reference.data=None \\
        program=xscale xscale.reference=bmin xscale.degrees_per_batch=None \\
        reject_method=framecc+lpstats rejection.lpstats.stats=em.b+bfactor \\
        clustering=cc cc_clustering.d_min=${clustering_dmin} cc_clustering.b_scale=false cc_clustering.use_normalized=false \\
        cc_clustering.min_cmpl=90 cc_clustering.min_redun=2 \\
        max_clusters=None xscale.use_tmpdir_if_available=${use_ramdisk} \\
#        batch.engine=sge batch.par_run=merging batch.nproc_each=8 nproc=8 batch.sge_pe_name=%s
""" % sge_pe_name)
        os.chmod(os.path.join(workdir, "merge_ccc.sh"), 0755)

        open(os.path.join(workdir, "filter_cell.R"), "w").write(r"""iqrf <- 2.5
outliers <- function(x) {
 q1 <- quantile(x, 0.25)
 q3 <- quantile(x, 0.75)
 iqr <- q3 - q1
 return(x<q1-iqr*iqrf | x>q3+iqr*iqrf)
}

myhist <- function(v, title) {
 if(sd(v)==0) {
  plot.new()
  return()
 }
 hist(v, main=paste("Histogram of", title), xlab=title)
 q1 <- quantile(v, 0.25)
 q3 <- quantile(v, 0.75)
 iqr <- q3 - q1
 abline(v=c(q1-iqr*iqrf, q3+iqr*iqrf), col="blue")
}

cells <- read.table("cells.dat", h=T)
good <- subset(cells, ! (outliers(a) | outliers(b) | outliers(c) | outliers(al) | outliers(be) | outliers(ga)))
write.table(good$file, "formerge_goodcell.lst", quote=F, row.names=F, col.names=F)

pdf("hist_cells.pdf", width=14, height=7)
par(mfrow=c(2,3))
myhist(cells$a, "a")
myhist(cells$b, "b")
myhist(cells$c, "c")
myhist(cells$al,"alpha")
myhist(cells$be,"beta")
myhist(cells$ga,"gamma")
dev.off()
cat("See hist_cells.pdf\n\n")

cat(sprintf("%4d files given.\n", nrow(cells)))
cat(sprintf("mean: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", mean(cells$a), mean(cells$b), mean(cells$c), mean(cells$al), mean(cells$be), mean(cells$ga)))
cat(sprintf(" std: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", sd(cells$a), sd(cells$b), sd(cells$c), sd(cells$al), sd(cells$be), sd(cells$ga)))
cat(sprintf(" iqr: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", IQR(cells$a), IQR(cells$b), IQR(cells$c), IQR(cells$al), IQR(cells$be), IQR(cells$ga)))

cat(sprintf("\n%4d files removed.\n", nrow(cells)-nrow(good)))
cat(sprintf("mean: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", mean(good$a), mean(good$b), mean(good$c), mean(good$al), mean(good$be), mean(good$ga)))
cat(sprintf(" std: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", sd(good$a), sd(good$b), sd(good$c), sd(good$al), sd(good$be), sd(good$ga)))
cat(sprintf(" iqr: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", IQR(good$a), IQR(good$b), IQR(good$c), IQR(good$al), IQR(good$be), IQR(good$ga)))

cat("\nUse formerge_goodcell.lst instead!\n")
""")

        if prep_dials_files:
            wd_jref = os.path.join(workdir, "dials_joint_refine")
            os.mkdir(wd_jref)

            ofs_phil = open(os.path.join(wd_jref, "experiments_and_reflections.phil"), "w")
            ofs_phil.write("input {\n")
            for wd in sorted(self.cell_and_files):
                fe = os.path.join(wd, "experiments.json")
                fp = os.path.join(wd, "integrate_hkl.pickle")
                if os.path.isfile(fe) and os.path.isfile(fp):
                    ofs_phil.write(' experiments = "%s"\n' % fe)
                    ofs_phil.write(' reflections = "%s"\n' % fp)
            ofs_phil.write("}\n")
            ofs_phil.close()

            open(os.path.join(wd_jref, "joint_refine.sh"), "w").write("""\
#!/bin/sh

dials.combine_experiments experiments_and_reflections.phil reference_from_experiment.beam=0 reference_from_experiment.goniometer=0 average_detector=true compare_models=false
dials.refine combined_experiments.json combined_reflections.pickle auto_reduction.action=remove verbosity=9
""")

    # write_merging_scripts()
# class PrepMerging

def run(params):
    if not params.workdir:
        print "Give workdir="
        return
    if os.path.exists(params.workdir):
        print "workdir already exists:", params.workdir
        return

    params.workdir = os.path.abspath(params.workdir)
    
    if None not in (params.unit_cell, params.space_group):
        user_xs = crystal.symmetry(params.unit_cell, params.space_group)
    else:
        user_xs = None
        
    from yamtbx.dataproc.auto.command_line.multi_check_cell_consistency import CellGraph

    cm = CellGraph(tol_length=params.cell_grouping.tol_length,
                   tol_angle=params.cell_grouping.tol_angle)

    if len(params.xdsdir) == 1 and os.path.isfile(params.xdsdir[0]):
        params.xdsdir = util.read_path_list(params.xdsdir[0])
        
    xds_dirs = []
    for xd0 in params.xdsdir:
        for xd in glob.glob(xd0):
            xds_dirs.extend(map(lambda x: x[0], filter(lambda x: any(map(lambda y: y.startswith("XDS_ASCII.HKL"), x[2])) or "DIALS.HKL" in x[2],
                                                       os.walk(os.path.abspath(xd)))))
    
    for i, xd in enumerate(xds_dirs):
        cm.add_proc_result(i, xd)

    pm = PrepMerging(cm)
    print pm.find_groups()

    if len(cm.groups) == 0:
        print "Oh, no. No data."
        return

    if params.group_choice is None:
        while True:
            try:
                val = int(raw_input("Input group number [%d..%d]: " % (1, len(cm.groups))))
                if not 0 < val <= len(cm.groups): raise ValueError
                params.group_choice = val
                break
            except ValueError:
                continue

    symms = cm.get_selectable_symms(params.group_choice-1)
    symmidx = -1
    
    if user_xs:
        #for xs in cm.get_selectable_symms(params.group_choice):
        raise "Not supported now."

    while True:
        try:
            val = int(raw_input("Input symmetry number [%d..%d]: " % (0, len(symms)-1)))
            if not 0 <= val < len(symms): raise ValueError
            symmidx = val
            break
        except ValueError:
            continue

    os.mkdir(params.workdir)
                
    topdir = os.path.dirname(os.path.commonprefix(xds_dirs))
    
    pm.prep_merging(group=params.group_choice, symmidx=symmidx, workdir=params.workdir, topdir=topdir,
                    cell_method=params.cell_method, nproc=params.nproc, prep_dials_files=params.prep_dials_files, into_workdir=params.copy_into_workdir)
    pm.write_merging_scripts(params.workdir, "par", params.prep_dials_files)
# run()

def run_from_args(argv):
    if "-h" in argv or "--help" in argv:
        print """
kamo.multi_prep_merging is a helper program to prepare for merging multiple (small wedge) datasets.
Typical usage:
  kamo.multi_prep_merging workdir=merge_reidx xdsdir=../\*/

All parameters:
"""
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    
    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
