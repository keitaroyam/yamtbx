"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os
import shutil

from yamtbx.dataproc.xds import xscalelp
from yamtbx.dataproc import xds
from yamtbx import util
from yamtbx.util import batchjob

from libtbx import easy_mp

xscale_comm = "xscale_par"
"""
def run_xscale(d_min=None, d_max=100, nbins=9, anomalous_flag=None, min_ios=None, outfile="xscale.hkl", sg=None, cell=None, nproc=1):

    inp_out.write("MAXIMUM_NUMBER_OF_PROCESSORS= %d\n" % nproc)
    if min_ios is not None: inp_out.write("MINIMUM_I/SIGMA= %.2f\n")
    inp_out.write("OUTPUT_FILE= %s\n"%outfile)

    if anomalous_flag is not None:
        inp_out.write("FRIEDEL'S_LAW= %s\n" % ("FALSE" if anomalous_flag else "TRUE"))

    if d_min is not None:
        step = ( 1./(d_min**2) - 1./(d_max**2) ) / float(nbins)
        start = 1./(d_max**2)
        rshells = " ".join(map(lambda i: "%.2f" % (start + i * step)**(-1./2), xrange(1, nbins+1)))
        inp_out.write("RESOLUTION_SHELLS= %s\n\n" % rshells)

    if None not in (sg, cell):
        inp_out.write("SPACE_GROUP_NUMBER= %s\nUNIT_CELL_CONSTANTS= %s\n\n" % (sg, cell))
    """

def run_xscale(xscale_inp):
    ftable = {}
    outfile = None
    count = 0
    wdir = os.path.dirname(os.path.abspath(xscale_inp))

    os.rename(xscale_inp, xscale_inp+".org")
    ofs = open(xscale_inp, "w")

    # Check line length and make symlink if needed
    for l in open(xscale_inp+".org"):
        ll = l[:l.index("!")] if "!" in l else l
        if "OUTPUT_FILE=" in ll:
            outfile = ll[ll.index("=")+1:].strip()
        if "INPUT_FILE=" in ll and len(l) > 132: # one line is limited to 131 characters!
            filename = ll[ll.index("=")+1:].strip()
            if "*" in filename: filename = filename[filename.index("*")+1:].strip()
            assert " " not in filename
            lnkf = "lnk%.6d.hkl" % count
            assert not os.path.isfile(os.path.join(wdir, lnkf))
            filename_abs = os.path.normpath(os.path.join(wdir, filename)) if not os.path.isabs(filename) else filename
            os.symlink(filename_abs, os.path.join(wdir, lnkf))
            print "xscale: %s -> %s" % (lnkf, filename)
            count += 1
            ftable[lnkf] = filename
            l = l.replace(filename, lnkf)
        ofs.write(l)

    ofs.close()
    assert outfile is not None

    if len(ftable) == 0:
        os.rename(xscale_inp+".org", xscale_inp)

    # Run xscale
    util.call(xscale_comm, wdir=wdir)

    # Replace file names if needed
    if len(ftable) > 0:
        for f in "XSCALE.LP", outfile:
            f = os.path.join(wdir, f)
            if not os.path.isfile(f): continue

            os.rename(f, f+".org")
            ofs = open(f, "w")
            for l in open(f+".org"):
                for lfn in ftable: l = l.replace(lfn, ftable[lfn])
                ofs.write(l)
            ofs.close()

        for lfn in ftable:
            os.remove(os.path.join(wdir, lfn))

        os.rename(xscale_inp+".org", xscale_inp)
# run_xscale()

def _calc_cchalf_by_removing_worker_1(wdir, inp_head, inpfiles, iex, nproc=None):
    print "Doing", iex

    tmpdir = os.path.join(wdir, "work.%.3d"%iex)
    if not os.path.exists(tmpdir): os.mkdir(tmpdir)

    newinp = os.path.join(tmpdir, "XSCALE.INP")

    files = inpfiles[:iex] + inpfiles[iex+1:]
    ofs = open(newinp, "w")
    ofs.write(inp_head)
    ofs.write("\nSAVE_CORRECTION_IMAGES= FALSE\n")
    if nproc is not None: ofs.write("MAXIMUM_NUMBER_OF_PROCESSORS= %d\n\n"%nproc)
    else: ofs.write("\n")

    for f in files:
        tmp = min(os.path.relpath(f, tmpdir), f, key=lambda x:len(x))
        ofs.write(" INPUT_FILE= %s\n" % tmp)
    ofs.close()

    #util.call(xscale_comm, wdir=tmpdir)
    return tmpdir
# _calc_cchalf_by_removing_worker_1()

def _calc_cchalf_by_removing_worker_2(wdir, tmpdir, iex, stat_bin):
    assert stat_bin in ("total", "outer")
    newlp = os.path.join(tmpdir, "XSCALE.LP")
    newinp = os.path.join(tmpdir, "XSCALE.INP")
    table = xscalelp.read_stats_table(newlp)
    if table is None:
        shutil.rmtree(tmpdir)
        return iex, float("nan"), -1

    assert table["dmin"][-1] is None # None for total
    i_stat = -1 if stat_bin == "total" else -2
    cchalf_exi = table["cc_half"][i_stat]
    nuniq = table["nuniq"][i_stat]

    # backup .INP and .LP, and then remove directory.
    os.rename(newinp, os.path.join(wdir, "XSCALE.INP.ex%.3d"%iex))
    os.rename(newlp, os.path.join(wdir, "XSCALE.LP.ex%.3d"%iex))
    shutil.rmtree(tmpdir)

    return iex, cchalf_exi, nuniq
# _calc_cchalf_by_removing_worker_2()

def calc_cchalf_by_removing(wdir, inp_head, inpfiles, with_sigma=False, stat_bin="total", nproc=1, nproc_each=None, batchjobs=None):
    assert not with_sigma # Not supported now
    assert stat_bin in ("total", "outer")

    if not os.path.exists(wdir): os.makedirs(wdir)

    datout = open(os.path.join(wdir, "cchalf.dat"), "w")
    datout.write("idx exfile cc1/2(%s) Nuniq\n" % stat_bin)

    cchalf_list = [] # (i_ex, CC1/2, Nuniq)

    # Prep runs
    tmpdirs = map(lambda x: _calc_cchalf_by_removing_worker_1(wdir, inp_head, inpfiles, x, nproc_each),
                  xrange(len(inpfiles)))
    # Run XSCALE 
    if batchjobs is not None:
        jobs = []
        for tmpdir in tmpdirs:
            job = batchjob.Job(tmpdir, "xscale.sh", nproc=nproc_each)
            job.write_script(xscale_comm)
            batchjobs.submit(job)
            jobs.append(job)

        batchjobs.wait_all(jobs)
    else:
        easy_mp.pool_map(fixed_func=lambda x: util.call(xscale_comm, wdir=x),
                         args=tmpdirs,
                         processes=nproc)
    # Finish runs
    cchalf_list = map(lambda x: _calc_cchalf_by_removing_worker_2(wdir, x[1], x[0], stat_bin), enumerate(tmpdirs))

    for iex, cchalf_exi, nuniq in cchalf_list:
        datout.write("%3d %s %.4f %d\n" % (iex, inpfiles[iex], cchalf_exi, nuniq))

    cchalf_list.sort(key=lambda x: -x[1])
    print
    print "# Sorted table"
    for idx, cch, nuniq in cchalf_list:
        print "%3d %-.4f %4d %s" % (idx, cch, nuniq, inpfiles[idx])

    # Remove unuseful (failed) data
    cchalf_list = filter(lambda x: x[1]==x[1], cchalf_list)

    return cchalf_list
# calc_delta_cchalf()


def calc_delta_cchalf(prev_lp, tmpdir, with_sigma=False, precalc_cchalf_all=None):
    """
    Obsolete function. Maybe useful when starting with the last XSCALE.LP...
    """

    assert not with_sigma # Not supported now

    if not os.path.exists(tmpdir): os.makedirs(tmpdir)

    newinp = os.path.join(tmpdir, "XSCALE.INP")
    newlp = os.path.join(tmpdir, "XSCALE.LP")

    rel_org = os.path.relpath(os.path.dirname(prev_lp), tmpdir)

    if precalc_cchalf_all is None:
        # read CC1/2(all) from cwd
        orgtable = xscalelp.read_stats_table(prev_lp)
        assert orgtable["dmin"][-1] is None # None for total
        cchalf_all = orgtable["cc_half"][-1]
    else:
        cchalf_all = precalc_cchalf_all

    datout = open(os.path.join(tmpdir, "delta_cchalf.dat"), "w")
    datout.write("# CC1/2(all)= %.4f\n" % cchalf_all)
    datout.write("idx exfile cc1/2 delta_cc1/2\n")

    # Read inp and extract input files.
    # XXX What if reference file is included???
    orgkwds = xscalelp.read_control_cards(prev_lp)
    inpfiles = map(lambda x:x[1],
                   filter(lambda y: y[0]=="INPUT_FILE", orgkwds))
                          
    # XXX Need to take care of xscale specific inp manner - order matters!!

    delta_cchalf = []

    for iex in xrange(len(inpfiles)):
        print "Doing", iex
        files = inpfiles[:iex] + inpfiles[iex+1:]
        ofs = open(newinp, "w")
        for k, v in orgkwds:
            if k not in ("INPUT_FILE", "INCLUDE_RESOLUTION_RANGE"):
                ofs.write("%s= %s\n" % (k,v))

        for f in files:
            if not os.path.isabs(f): f = os.path.join(rel_org, f)
            ofs.write("INPUT_FILE= %s\n" % f)
        ofs.close()

        util.call(xscale_comm, wdir=tmpdir)
        table = xscalelp.read_stats_table(newlp)
        assert table["dmin"][-1] is None # None for total
        cchalf_exi = table["cc_half"][-1]
        delta_cchalf.append((iex, cchalf_exi - cchalf_all))

        os.rename(newinp, newinp+".ex%.3d"%iex)
        os.rename(newlp, newlp+".ex%.3d"%iex)

        datout.write("%3d %s %.4f %.4f\n" % (iex, inpfiles[iex], cchalf_exi, cchalf_exi-cchalf_all))

    delta_cchalf.sort(key=lambda x: -x[1])
    print
    print "# Sorted table"
    for idx, dch in delta_cchalf:
        print "%3d %-.4f %s" % (idx, dch, inpfiles[idx])

    return delta_cchalf, cchalf_all
# calc_delta_cchalf()

def decide_scaling_reference_based_on_bfactor(lpin, ref, return_as="index"):
    assert ref in ("bmed", "bmin", "bmax")
    assert return_as in ("index", "filename")

    KBs = map(lambda x: [x[0]]+x[1], enumerate(xscalelp.get_k_b(lpin))) # list of [i, K, B, filename]
    KBs.sort(key=lambda x: x[2])

    if len(KBs) == 0:
        return 0 # if only one data, scale and B values don't exist.

    if ref == "bmin":   idx = 0
    elif ref == "bmax": idx = -1
    else:               idx = len(KBs)//2

    if return_as == "index": return KBs[idx][0]
    else: return KBs[idx][-1]
# decide_scaling_reference_based_on_bfactor()
